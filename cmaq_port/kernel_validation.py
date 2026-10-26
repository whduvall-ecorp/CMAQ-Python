"""Compare Python plume helpers to unchanged, standalone CMAQ Fortran."""

import argparse
import ctypes
import json
import math
import os
import shutil
import subprocess
from dataclasses import astuple, dataclass
from pathlib import Path
from typing import TypedDict

from .inventory import PROJECT_ROOT, git, save_report, sha256
from .fire_rise import HCRIT as FIRE_HCRIT, GRAV as FIRE_GRAV, FirePlumeRise, fire_plmris
from .plume import delta_zs, plsprd, real32
from .rise import GRAV, HCRIT, plmris

DEFAULT_COMPILER = PROJECT_ROOT.parent / "AERMOD-Python" / "tools" / "mingw64" / "bin" / "gfortran.exe"

# Engineering limits fixed before execution; no automatic tolerance inflation.
TOLERANCES = {
    "PLSPRD": {"atol": 2e-5, "rtol": 2e-6, "units": "m"},
    "DELTA_ZS": {"atol": 1e-10, "rtol": 2e-6, "units": "ZSTK m; DDZF 1/m"},
    "PLMRIS": {"atol": 2e-4, "rtol": 2e-6, "units": "centerline m; WSTK m/s; exact top flag"},
    "CONNECTED": {"atol": 2e-4, "rtol": 2e-6,
                  "units": "rise; exact top flag; geometry uses DELTA_ZS tolerance; bounds use PLSPRD tolerance"},
    "FIRE_PLMRIS": {"atol": 2e-4, "rtol": 2e-6,
                    "units": "centerline m; reporting wind m/s; exact top flag and LOGDEV record"},
}


def compiler_path(path: Path) -> str:
    if os.name != "nt":
        return str(path)
    get_short = ctypes.WinDLL("kernel32", use_last_error=True).GetShortPathNameW
    get_short.argtypes = [ctypes.c_wchar_p, ctypes.c_wchar_p, ctypes.c_uint32]
    get_short.restype = ctypes.c_uint32
    buffer = ctypes.create_unicode_buffer(32768)
    length = get_short(str(path), buffer, len(buffer))
    if not length:
        raise ctypes.WinError(ctypes.get_last_error())
    if length >= len(buffer) or " " in buffer.value:
        raise ValueError("Compiler requires a usable Windows short path.")
    return buffer.value


def cases() -> list[tuple[tuple[float, float, float], float]]:
    return [
        (gradients, center)
        for gradients in ((0.0, 0.01, 0.02), (0.0098, 0.0098, 0.0098), (-0.01, 0.0, 0.03), (0.1, 0.1, 0.1))
        for center in (0.0, 1.0, 99.99, 100.0, 100.01, 499.99, 500.0, 500.01, 999.0, 1000.0, 1500.0)
    ]


class RiseDiagnostics(TypedDict):
    stability_path: tuple[int, ...]
    layer_transitions: int
    hot_layer_adjusted: bool


class FireDiagnostics(TypedDict):
    stability_path: tuple[int, ...]
    layer_transitions: int


@dataclass(frozen=True)
class KernelCase:
    label: str
    routine: str
    protocol: str
    expected: tuple[float, ...]
    diagnostics: RiseDiagnostics | FireDiagnostics | None = None
    top_log: str | None = None


def geometry_cases() -> list[KernelCase]:
    result: list[KernelCase] = []
    for grid in ((10.0,), (100.0, 500.0, 1000.0), (0.125, 0.5, 1.0), (99.99, 100.01, 1000.0)):
        for mapping in ((0,), (2, 0), (1, 1, 0), (2, 0, 1)):
            for stacks in ((0.0, 100.0, 400.0), (10.0, 100.0, 1500.0), (99.99, 500.01, 0.0), (1.0, 1.0, 1.0)):
                heights = tuple(tuple(value * (1 + 0.125 * s) for s in range(len(mapping))) for value in grid)
                geometry = delta_zs(mapping, stacks, heights)
                # Fortran list-directed arrays traverse layer index first.
                expected = tuple(
                    row[s] for matrix in (geometry.relative_heights, geometry.inverse_thicknesses)
                    for s in range(len(mapping)) for row in matrix
                )
                values = [f"{len(grid)} {len(mapping)} {len(stacks)}",
                          " ".join(str(index + 1) for index in mapping),
                          " ".join(format(value, ".17g") for value in stacks),
                          " ".join(format(row[s], ".17g") for s in range(len(mapping)) for row in heights)]
                result.append(KernelCase(
                    f"delta-zs-{len(result)}", "DELTA_ZS",
                    "DELTA_ZS\n" + "\n".join(values) + "\n", expected,
                ))
    return result


def rise_cases() -> list[KernelCase]:
    """Deterministic stratification/stack matrix; no random expected outputs."""
    result: list[KernelCase] = []
    grids = (
        (0.0, 20.0, 60.0, 140.0, 300.0, 600.0),
        (0.0, 100.0, 500.0, 1000.0),
    )
    profiles = (
        ("neutral", 0.0, (0.0,)),
        ("stable-flux", -0.03, (0.0001,)),
        ("stable-gradient", 0.0, (0.02,)),
        ("unstable", 0.03, (0.0,)),
        ("layer-changing", 0.03, (0.0, 0.02, -0.005, 0.0001, 0.01)),
        ("stable-to-neutral", -0.03, (0.02, 0.0, 0.0, 0.0, 0.0)),
    ) + tuple(
        (f"boundary-heat-{sign}-{factor}", sign * real32(
            real32(HCRIT * real32(289.7) / GRAV) * factor
        ), (0.0,))
        for sign in (-1, 1) for factor in (0.999999, 1.0, 1.000001)
    ) + (
        ("boundary-gradient-equal", 0.0, (0.001,)),
        ("boundary-gradient-above", 0.0, (0.001000001,)),
    )
    for grid in grids:
        n = len(grid) - 1
        centers = tuple((a + b) / 2 for a, b in zip(grid, grid[1:]))
        for profile, hfx, values in profiles:
            boundary = profile.startswith("boundary-")
            if boundary and n != 3:
                continue
            lapse = tuple(values[i % len(values)] for i in range(n))
            ta = tuple(290.0 - 0.006 * height for height in centers)
            wind = tuple((0.2, 1.0, 5.0, 8.0, 2.0)[i] for i in range(n))
            for ht in ((grid[1],) if boundary else (1.0, grid[1], grid[1] + 0.25, grid[-1] + 20.0)):
                layer = next((i for i in range(n) if ht <= grid[i + 1]), n - 1)
                for hmix in ((800.0,) if boundary else (0.0, 800.0)):
                    stack_cases = (
                        (280.0, 0.0, 1.0, 0.0, 0.2),
                        (290.0, 12.0, 2.0, 3.0, 0.3),
                        (291.0, 40.0, 4.0, 0.2, 0.3),
                        (450.0, 10.0, 2.0, 4.0, 0.2),
                        (1000.0, 20.0, 6.0, 1.0, 0.8),
                    )
                    for tk, ve, dm, ws, us in ((stack_cases[3],) if boundary else stack_cases):
                        arguments = (layer, hfx, hmix, dm, ht, tk, ve, 290.0, us,
                                     lapse, ta, wind, grid, centers, ws)
                        rise = plmris(*arguments)
                        geometry = delta_zs((0,), (ht,), tuple((v,) for v in grid[1:]))
                        zstk = tuple(row[0] for row in geometry.relative_heights)
                        scalars = (hfx, hmix, dm, ht, tk, ve, 290.0, us, ws)
                        lines = [f"{n} {layer + 1}"] + [
                            " ".join(format(v, ".17g") for v in row)
                            for row in (scalars, lapse, ta, wind, grid, centers)
                        ]
                        outputs = (rise.centerline, rise.stack_wind, float(rise.top_reached))
                        diagnostics: RiseDiagnostics = {
                            "stability_path": rise.stability_path,
                            "layer_transitions": rise.layer_transitions,
                            "hot_layer_adjusted": rise.hot_layer_adjusted,
                        }
                        label = f"{profile}-grid{n}-ht{ht}-hmix{hmix}-tk{tk}"
                        result.append(KernelCase(
                            label, "PLMRIS", "PLMRIS\n" + "\n".join(lines + [
                                " ".join(format(v, ".17g") for v in zstk),
                            ]) + "\n", outputs, diagnostics,
                        ))
                        bounds = plsprd(lapse, grid, rise.centerline)
                        result.append(KernelCase(
                            "connected-" + label, "CONNECTED",
                            "CONNECTED\n" + "\n".join(lines) + "\n",
                            outputs + zstk + tuple(row[0] for row in geometry.inverse_thicknesses)
                            + (bounds.top, bounds.bottom), diagnostics,
                        ))
    return result


def fire_cases() -> list[KernelCase]:
    """Synthetic FIRE_PLMRIS cases spanning stability, gradients, and layers."""
    grid = (0.0, 20.0, 60.0, 140.0, 300.0, 600.0)
    profiles = (
        ("neutral-gradient", (0.0,) * 5),
        ("stable-gradient", (0.02,) * 5),
        ("threshold-gradient", (0.001,) * 5),
        ("layer-changing", (0.0, 0.02, -0.005, 0.0001, 0.01)),
    )
    temperatures = tuple(290.0 - 0.006 * ((lower + upper) / 2)
                        for lower, upper in zip(grid, grid[1:]))
    winds = (0.0, 1.0, 5.0, 8.0, 2.0)
    result: list[KernelCase] = []

    def append_case(
        label: str, gradients: tuple[float, ...], heat: float, mixing: float,
        buoyancy: float, layer: int, surface_temperature: float | None = None,
    ) -> None:
        initial_wind = 0.2
        fire_temperature = 450.0
        friction = 0.3
        case_temperatures = temperatures
        if surface_temperature is not None:
            case_temperatures = (real32(surface_temperature), *temperatures[1:])
        rise: FirePlumeRise = fire_plmris(
            layer, heat, mixing, buoyancy, fire_temperature, friction,
            gradients, case_temperatures, winds, grid, initial_wind,
        )
        protocol = "\n".join((
            "FIRE_PLMRIS",
            f"{len(gradients)} {layer + 1}",
            " ".join(format(value, ".17g") for value in (
                heat, mixing, buoyancy, fire_temperature, friction, initial_wind,
            )),
            " ".join(format(value, ".17g") for value in gradients),
            " ".join(format(value, ".17g") for value in case_temperatures),
            " ".join(format(value, ".17g") for value in winds),
            " ".join(format(value, ".17g") for value in grid),
            "",
        ))
        diagnostics: FireDiagnostics = {
            "stability_path": rise.stability_path,
            "layer_transitions": rise.layer_transitions,
        }
        top_log = None
        if rise.top_log_height is not None:
            top_log = (
                "     Plume rise reached EMLAYS with ZPLM:"
                f"{len(gradients):3d}{rise.top_log_height:10.3f}"
            )
        result.append(KernelCase(
            label, "FIRE_PLMRIS", protocol,
            (rise.centerline, rise.fire_wind, float(rise.top_reached)),
            diagnostics, top_log,
        ))

    for profile, gradients in profiles:
        for heat in (-0.03, 0.0, 0.03):
            for mixing in (0.0, 800.0):
                for buoyancy in (0.0, 0.1, 100.0):
                    for layer in (0, 2, 4):
                        append_case(
                            f"{profile}-hfx{heat}-hmix{mixing}-bflx{buoyancy}-layer{layer}",
                            gradients, heat, mixing, buoyancy, layer,
                        )

    # Exercise strict HCRIT boundaries in the same default-REAL arithmetic as
    # the source, plus the independent strict DTHDZ > 0.001 classification.
    for sign in (-1, 1):
        boundary_heat = real32(real32(FIRE_HCRIT * real32(290.0)) / FIRE_GRAV)
        for factor in (0.999, 1.0, 1.001):
            append_case(
                f"boundary-hstar-{sign}-{factor}",
                (0.0,) * 5, real32(sign * real32(boundary_heat * real32(factor))),
                800.0, 1.0, 0, 290.0,
            )
    for label, gradient in (("equal", 0.001), ("above", 0.001001)):
        append_case(
            f"boundary-gradient-{label}", (gradient,) * 5, 0.0, 800.0, 1.0, 0,
        )
    return result


def validate(source: Path, compiler: Path) -> Path:
    source = source.resolve(strict=True)
    compiler = compiler.resolve(strict=True)
    originals = source / "CCTM" / "src" / "plrise" / "smoke"
    report: dict[str, object] = {
        "success": False,
        "schema_version": 2,
        "scope": "Isolated PLSPRD/DELTA_ZS/PLMRIS/FIRE_PLMRIS and connected calls; synthetic inputs, not a CMAQ run.",
        "tolerances": TOLERANCES,
    }
    destination = save_report(report, "kernels")
    run = destination.parent
    try:
        for filename in ("plsprd.f", "delta_zs.f", "plmris.F", "fire_plmris.F"):
            shutil.copy2(originals / filename, run / filename)
        shutil.copy2(source / "CCTM" / "src" / "ICL" / "fixed" / "const" / "CONST.EXT", run / "CONST.EXT")
        shutil.copy2(source / "license.md", run / "license.md")
        shutil.copy2(PROJECT_ROOT / "tools" / "plsprd_driver.f90", run / "driver.f90")
        shutil.copy2(PROJECT_ROOT / "tools" / "runtime_vars_test.f90", run / "runtime_vars_test.f90")
        report.update({
            "source_commit": git(source, "rev-parse", "HEAD").strip(),
            "source_status": git(source, "status", "--porcelain=v1"),
            "source_sha256": sha256(run / "plsprd.f"), "driver_sha256": sha256(run / "driver.f90"),
            "source_hashes": {
                name: sha256(run / name)
                for name in ("plsprd.f", "delta_zs.f", "plmris.F", "fire_plmris.F", "CONST.EXT")
            },
            "repository_license_sha256": sha256(run / "license.md"),
            "adapter_hashes": {name: sha256(run / name) for name in ("driver.f90", "runtime_vars_test.f90")},
            "runtime_adapter": "TEST-ONLY RUNTIME_VARS exports LOGDEV, directed to plmris.log; no physics stubs.",
            "compiler_sha256": sha256(compiler),
        })
        command = [compiler_path(compiler), "-O0", "-fcheck=all", "-fno-fast-math",
                   "-ffixed-line-length-none", "-static"]
        if os.name == "nt":
            command.append("-B" + compiler_path(compiler.parent.parent / "x86_64-w64-mingw32" / "lib") + os.sep)
        command.extend(["-cpp", '-DSUBST_CONST="CONST.EXT"', "runtime_vars_test.f90",
                        "plsprd.f", "delta_zs.f", "plmris.F", "fire_plmris.F", "driver.f90", "-o", "driver.exe"])
        report["compile_command"] = command
        report["compiler_version"] = subprocess.run(
            [compiler_path(compiler), "--version"], capture_output=True, text=True, check=True, timeout=30,
        ).stdout
        built = subprocess.run(command, cwd=run, capture_output=True, text=True, timeout=120)
        (run / "compile.log").write_text(built.stdout + built.stderr, encoding="utf-8")
        if built.returncode:
            raise RuntimeError(f"CMAQ kernel compile failed; see {run / 'compile.log'}")
        samples = [
            KernelCase(
                f"plsprd-{i}", "PLSPRD",
                "PLSPRD\n" + " ".join(format(v, ".17g") for v in (*g, h)) + "\n",
                astuple(plsprd(g, (0, 100, 500, 1000), h)),
            )
            for i, (g, h) in enumerate(cases())
        ] + geometry_cases() + rise_cases() + fire_cases()
        protocol = "".join(sample.protocol for sample in samples)
        (run / "inputs.txt").write_text(protocol, encoding="ascii")
        executed = subprocess.run(
            [str(run / "driver.exe")], cwd=run, input=protocol,
            text=True, capture_output=True, timeout=60,
        )
        (run / "outputs.txt").write_text(executed.stdout, encoding="ascii")
        (run / "stderr.txt").write_text(executed.stderr, encoding="utf-8")
        if executed.returncode:
            raise RuntimeError("CMAQ test driver failed; see stderr.txt in the run folder.")
        lines = executed.stdout.splitlines()
        if len(lines) != len(samples):
            raise ValueError("CMAQ driver returned an unexpected number of lines.")
        top_messages = (run / "plmris.log").read_text(encoding="ascii").splitlines()
        expected_messages = sum(
            sample.expected[2] == 1 for sample in samples if sample.routine in ("PLMRIS", "CONNECTED")
        ) + sum(sample.top_log is not None for sample in samples)
        if len(top_messages) != expected_messages or any(
            "Plume rise reached EMLAYS with ZPLM:" not in message for message in top_messages
        ):
            raise ValueError("PLMRIS LOGDEV diagnostics differ from Python top-reaching flags.")
        report["top_reaching_messages"] = len(top_messages)
        expected_fire_messages = [sample.top_log for sample in samples if sample.top_log is not None]
        if expected_fire_messages and top_messages[-len(expected_fire_messages):] != expected_fire_messages:
            raise ValueError("FIRE_PLMRIS LOGDEV records do not exactly match Python diagnostics.")
        report["fire_top_log_records_exact"] = True
        results = []
        for sample, line in zip(samples, lines, strict=True):
            expected = sample.expected
            observed = tuple(float(v) for v in line.split())
            if len(observed) != len(expected) or not all(math.isfinite(v) for v in observed):
                raise ValueError(f"Invalid {sample.routine} output: {line}")
            errors = tuple(abs(a - b) for a, b in zip(expected, observed, strict=True))
            absolute_tolerance = TOLERANCES[sample.routine]["atol"]
            limits = [float(absolute_tolerance)] * len(expected)
            relative_limits = [2e-6] * len(expected)
            if sample.routine in ("PLMRIS", "CONNECTED", "FIRE_PLMRIS"):
                limits[2], relative_limits[2] = 0.0, 0.0
            if sample.routine == "CONNECTED":
                limits[3:-2] = [1e-10] * (len(expected) - 5)
                limits[-2:] = [2e-5] * 2
            results.append({
                "label": sample.label, "routine": sample.routine, "protocol": sample.protocol,
                "python": expected, "fortran": observed,
                "python_diagnostics": sample.diagnostics, "expected_top_log": sample.top_log,
                "atol": absolute_tolerance, "rtol": 2e-6,
                "absolute_limits": limits, "relative_limits": relative_limits,
                "max_absolute_error": max(errors),
                "max_relative_error": max(e / abs(v) if v else 0.0 for e, v in zip(errors, observed, strict=True)),
                "passed": all(
                    e <= a + r * abs(v) for e, v, a, r in zip(errors, observed, limits, relative_limits, strict=True)
                ),
            })
        report.update({
            "cases": len(samples), "values": sum(len(sample.expected) for sample in samples), "results": results,
            "rtol": 2e-6, "precision": "default REAL (binary32)",
            "per_routine": {
                name: {
                    "cases": sum(row["routine"] == name for row in results),
                    "values": sum(len(row["python"]) for row in results if row["routine"] == name),
                    "max_absolute_error": max(row["max_absolute_error"] for row in results if row["routine"] == name),
                    "max_relative_error": max(row["max_relative_error"] for row in results if row["routine"] == name),
                } for name in ("PLSPRD", "DELTA_ZS", "PLMRIS", "CONNECTED", "FIRE_PLMRIS")
            },
            "max_absolute_error": max(row["max_absolute_error"] for row in results),
            "success": all(row["passed"] for row in results),
        })
        report["source_status_after"] = git(source, "status", "--porcelain=v1")
        unchanged = all(
            sha256(originals / name) == sha256(run / name)
            for name in ("plsprd.f", "delta_zs.f", "plmris.F", "fire_plmris.F")
        ) and sha256(source / "CCTM" / "src" / "ICL" / "fixed" / "const" / "CONST.EXT") == sha256(run / "CONST.EXT")
        unchanged = unchanged and sha256(source / "license.md") == sha256(run / "license.md")
        report["source_unchanged"] = unchanged and report["source_status"] == report["source_status_after"]
        if not report["source_unchanged"]:
            report["success"] = False
            raise ValueError("Original source changed during validation.")
        if not report["success"]:
            raise ValueError("CMAQ kernel comparison failed.")
    except (OSError, ValueError, RuntimeError, subprocess.SubprocessError) as exc:
        report["error"] = str(exc)
        raise
    finally:
        destination.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return destination


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=Path(r"C:\Dev\CMAQ"))
    parser.add_argument("--compiler", type=Path, default=DEFAULT_COMPILER)
    args = parser.parse_args()
    print(validate(args.source, args.compiler))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
