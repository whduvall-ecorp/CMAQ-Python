"""Compare the Python POLY port to the unchanged, standalone I/O API POLY.

The reference `reference/ioapi/poly.f` is the I/O API routine CMAQ's PREPLM
calls; it is self-contained, so it is compiled standalone with the project-local
GNU Fortran toolchain and compared against `cmaq_port.interpolation.poly` on a
deterministic case matrix. Synthetic inputs only; this is kernel-level evidence,
not a CMAQ run or a PREPLM port.
"""

import argparse
import json
import math
import os
import shutil
import subprocess
from pathlib import Path

from .interpolation import poly
from .inventory import PROJECT_ROOT, sha256, save_report
from .kernel_validation import DEFAULT_COMPILER, compiler_path
from .plume import real32

REFERENCE = PROJECT_ROOT / "reference" / "ioapi" / "poly.f"
DRIVER = PROJECT_ROOT / "tools" / "poly_driver.f90"
TOLERANCE = {"atol": 1e-3, "rtol": 2e-6, "units": "interpolated value (same units as YPTS)"}


def cases() -> list[tuple[str, float, tuple[float, ...], tuple[float, ...], int]]:
    """Deterministic (label, xpt, xpts, ypts, ndeg) tuples; no random outputs."""
    result: list[tuple[str, float, tuple[float, ...], tuple[float, ...], int]] = []

    # PREPLM-like profiles: monotonic heights with smooth temperature/wind.
    heights = (10.0, 30.0, 60.0, 120.0, 250.0, 500.0)
    temperature = tuple(290.0 - 0.0065 * h for h in heights)
    wind = tuple(0.5 + 0.02 * h for h in heights)
    for ndeg in (1, 2, 3, 4, 5):
        window = heights[: ndeg + 1]
        for field_name, field in (("temp", temperature), ("wind", wind)):
            values = field[: ndeg + 1]
            # Interpolate at the stack top inside and near the window edges.
            for xpt in (window[0], window[-1], window[0] + 1.0,
                        0.5 * (window[0] + window[-1]), window[-1] - 3.0, window[0] - 5.0):
                result.append((
                    f"preplm-{field_name}-deg{ndeg}-x{xpt}", xpt, window, values, ndeg,
                ))

    # Exact-polynomial stress cases: coefficients chosen so POLY reproduces them.
    abscissae = (-2.0, -0.5, 1.0, 2.5, 4.0, 7.0)
    polynomials = (
        ("linear", lambda x: 3.0 * x - 1.0),
        ("quadratic", lambda x: 2.0 * x * x - 3.0 * x + 5.0),
        ("cubic", lambda x: x * x * x - 2.0 * x + 1.0),
    )
    for name, fn in polynomials:
        for ndeg in (1, 2, 3):
            window = abscissae[: ndeg + 1]
            values = tuple(fn(x) for x in window)
            for xpt in (window[0], 0.5 * (window[0] + window[-1]), window[-1], 1.3, -1.7):
                result.append((f"poly-{name}-deg{ndeg}-x{xpt}", xpt, window, values, ndeg))

    # Non-monotonic abscissae and a descending profile exercise sign handling.
    result.append(("nonmono-deg2", 1.5, (0.0, 3.0, 1.0), (1.0, 4.0, 9.0), 2))
    result.append(("descending-deg3", 2.0, (5.0, 4.0, 2.0, 0.0), (0.0, 1.0, 3.0, 7.0), 3))
    return result


def validate(compiler: Path = DEFAULT_COMPILER) -> Path:
    compiler = compiler.resolve(strict=True)
    report: dict[str, object] = {
        "success": False,
        "schema_version": 1,
        "routine": "POLY",
        "scope": "Isolated I/O API POLY vs Python port; synthetic inputs, not a CMAQ or PREPLM run.",
        "reference": "I/O API poly.f (vendored unchanged); see reference/ioapi/PROVENANCE.md.",
        "tolerance": TOLERANCE,
        "precision": "default REAL (binary32)",
    }
    destination = save_report(report, "poly")
    run = destination.parent
    try:
        shutil.copy2(REFERENCE, run / "poly.f")
        shutil.copy2(DRIVER, run / "poly_driver.f90")
        shutil.copy2(PROJECT_ROOT / "reference" / "ioapi" / "PROVENANCE.md", run / "PROVENANCE.md")
        report.update({
            "reference_sha256": sha256(run / "poly.f"),
            "driver_sha256": sha256(run / "poly_driver.f90"),
            "provenance_sha256": sha256(run / "PROVENANCE.md"),
            "compiler_sha256": sha256(compiler),
        })
        command = [compiler_path(compiler), "-O0", "-fcheck=all", "-fno-fast-math",
                   "-ffixed-line-length-none", "-static"]
        if os.name == "nt":
            command.append("-B" + compiler_path(compiler.parent.parent / "x86_64-w64-mingw32" / "lib") + os.sep)
        command.extend(["poly.f", "poly_driver.f90", "-o", "poly_driver.exe"])
        report["compile_command"] = command
        report["compiler_version"] = subprocess.run(
            [compiler_path(compiler), "--version"], capture_output=True, text=True, check=True, timeout=30,
        ).stdout
        built = subprocess.run(command, cwd=run, capture_output=True, text=True, timeout=120)
        (run / "compile.log").write_text(built.stdout + built.stderr, encoding="utf-8")
        if built.returncode:
            raise RuntimeError(f"I/O API POLY compile failed; see {run / 'compile.log'}")

        samples = cases()
        protocol_lines: list[str] = []
        expected: list[float] = []
        for label, xpt, xpts, ypts, ndeg in samples:
            x = real32(xpt)
            xs = tuple(real32(v) for v in xpts)
            ys = tuple(real32(v) for v in ypts)
            expected.append(poly(x, xs, ys, ndeg))
            protocol_lines.extend([
                "POLY",
                f"{ndeg} {format(x, '.17g')}",
                " ".join(format(v, ".17g") for v in xs),
                " ".join(format(v, ".17g") for v in ys),
            ])
        protocol = "\n".join(protocol_lines) + "\n"
        (run / "inputs.txt").write_text(protocol, encoding="ascii")
        executed = subprocess.run(
            [str(run / "poly_driver.exe")], cwd=run, input=protocol,
            text=True, capture_output=True, timeout=60,
        )
        (run / "outputs.txt").write_text(executed.stdout, encoding="ascii")
        (run / "stderr.txt").write_text(executed.stderr, encoding="utf-8")
        if executed.returncode:
            raise RuntimeError("POLY test driver failed; see stderr.txt in the run folder.")
        lines = executed.stdout.splitlines()
        if len(lines) != len(samples):
            raise ValueError("POLY driver returned an unexpected number of lines.")

        atol = float(TOLERANCE["atol"])
        rtol = float(TOLERANCE["rtol"])
        results = []
        for (label, xpt, xpts, ypts, ndeg), value, line in zip(samples, expected, lines, strict=True):
            observed = float(line.strip())
            if not math.isfinite(observed):
                raise ValueError(f"Nonfinite POLY output for {label}: {line}")
            error = abs(value - observed)
            results.append({
                "label": label, "ndeg": ndeg, "xpt": xpt,
                "python": value, "fortran": observed,
                "max_absolute_error": error,
                "passed": error <= atol + rtol * abs(observed),
            })
        report.update({
            "cases": len(samples),
            "values": len(samples),
            "results": results,
            "max_absolute_error": max(row["max_absolute_error"] for row in results),
            "success": all(row["passed"] for row in results),
        })
        if not report["success"]:
            raise ValueError("POLY comparison failed.")
    except (OSError, ValueError, RuntimeError, subprocess.SubprocessError) as exc:
        report["error"] = str(exc)
        raise
    finally:
        destination.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return destination


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--compiler", type=Path, default=DEFAULT_COMPILER)
    args = parser.parse_args()
    print(validate(args.compiler))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
