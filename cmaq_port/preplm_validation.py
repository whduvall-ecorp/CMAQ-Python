"""Compare the Python PREPLM port to the unchanged CMAQ preplm.f.

The unchanged CCTM/src/plrise/smoke/preplm.f is compiled with a test-only
UTILIO_DEFN adapter (exposing the unchanged I/O API POLY and a trivial M3MSG2)
and compared against cmaq_port.preplm on a deterministic case matrix. Synthetic
inputs only; this is kernel-level evidence, not a CMAQ run.
"""

import argparse
import json
import math
import os
import shutil
import subprocess
from pathlib import Path

from .inventory import PROJECT_ROOT, git, sha256, save_report
from .kernel_validation import DEFAULT_COMPILER, compiler_path
from .plume import real32
from .preplm import preplm

REFERENCE_POLY = PROJECT_ROOT / "reference" / "ioapi" / "poly.f"
ADAPTER = PROJECT_ROOT / "tools" / "utilio_defn_test.f90"
DRIVER = PROJECT_ROOT / "tools" / "preplm_driver.f90"
TOLERANCE = {"atol": 2e-3, "rtol": 2e-5,
             "units": "TSTK K; WSTK m/s; DTHDZ K/m; WSPD m/s; LSTK/LPBL exact"}


def _profile(emlays: int):
    zf = tuple(float(sum(20 * (j + 1) for j in range(i + 1))) for i in range(emlays))
    zh = tuple(((0.0 if i == 0 else zf[i - 1]) + zf[i]) / 2.0 for i in range(emlays))
    ddzf = tuple(1.0 / (zf[i] - (0.0 if i == 0 else zf[i - 1])) for i in range(emlays))
    ta = tuple(293.0 - 0.0065 * zh[i] for i in range(emlays))
    qv = tuple(0.009 - 0.0005 * i for i in range(emlays))
    uw = tuple(1.0 + 0.7 * i for i in range(emlays))
    vw = tuple(0.5 + 0.2 * i for i in range(emlays))
    pres = tuple(1013.0 - 40.0 * i for i in range(emlays + 1))
    return ddzf, qv, ta, uw, vw, zh, zf, pres


def cases() -> list[dict]:
    result: list[dict] = []
    for emlays in (4, 6, 8):
        ddzf, qv, ta, uw, vw, zh, zf, pres = _profile(emlays)
        top = zf[-1]
        for fireflg in (False, True):
            for hts in (zf[0] * 0.5, zf[0], top * 0.3, top * 0.7):
                for hmix in (zf[0] * 0.5, top * 0.5, top * 0.95):
                    result.append({
                        "label": f"preplm-{emlays}-{int(fireflg)}-{hts:.1f}-{hmix:.1f}",
                        "fireflg": fireflg, "emlays": emlays, "hmix": hmix, "hts": hts,
                        "psfc": 1013.0, "ts": 291.5, "ddzf": ddzf, "qv": qv, "ta": ta,
                        "uw": uw, "vw": vw, "zh": zh, "zf": zf, "pres": pres,
                    })
    return result


def _protocol(case: dict) -> str:
    def row(values) -> str:
        return " ".join(format(real32(v), ".17g") for v in values)
    lines = [
        "PREPLM",
        f"{1 if case['fireflg'] else 0} {case['emlays']} "
        + " ".join(format(real32(case[k]), ".17g") for k in ("hmix", "hts", "psfc", "ts")),
        row(case["ddzf"]), row(case["qv"]), row(case["ta"]), row(case["uw"]),
        row(case["vw"]), row(case["zh"]), row(case["zf"]), row(case["pres"]),
    ]
    return "\n".join(lines) + "\n"


def validate(source: Path = Path(r"C:\Dev\CMAQ"), compiler: Path = DEFAULT_COMPILER) -> Path:
    source = source.resolve(strict=True)
    compiler = compiler.resolve(strict=True)
    original = source / "CCTM" / "src" / "plrise" / "smoke" / "preplm.f"
    report: dict[str, object] = {
        "success": False, "schema_version": 1, "routine": "PREPLM",
        "scope": "Isolated unchanged preplm.f vs Python port; synthetic inputs, not a CMAQ run.",
        "adapter": "TEST-ONLY UTILIO_DEFN exposes the unchanged I/O API POLY and a trivial M3MSG2.",
        "tolerance": TOLERANCE, "precision": "default REAL (binary32)",
    }
    destination = save_report(report, "preplm")
    run = destination.parent
    try:
        shutil.copy2(original, run / "preplm.f")
        shutil.copy2(REFERENCE_POLY, run / "poly.f")
        shutil.copy2(ADAPTER, run / "utilio_defn_test.f90")
        shutil.copy2(DRIVER, run / "preplm_driver.f90")
        report.update({
            "source_commit": git(source, "rev-parse", "HEAD").strip(),
            "source_status": git(source, "status", "--porcelain=v1"),
            "preplm_sha256": sha256(run / "preplm.f"),
            "poly_sha256": sha256(run / "poly.f"),
            "adapter_sha256": sha256(run / "utilio_defn_test.f90"),
            "driver_sha256": sha256(run / "preplm_driver.f90"),
            "compiler_sha256": sha256(compiler),
        })
        command = [compiler_path(compiler), "-O0", "-fcheck=all", "-fno-fast-math",
                   "-ffixed-line-length-none", "-static"]
        if os.name == "nt":
            command.append("-B" + compiler_path(compiler.parent.parent / "x86_64-w64-mingw32" / "lib") + os.sep)
        command.extend(["utilio_defn_test.f90", "poly.f", "preplm.f", "preplm_driver.f90",
                        "-o", "preplm_driver.exe"])
        report["compile_command"] = command
        report["compiler_version"] = subprocess.run(
            [compiler_path(compiler), "--version"], capture_output=True, text=True, check=True, timeout=30,
        ).stdout
        built = subprocess.run(command, cwd=run, capture_output=True, text=True, timeout=120)
        (run / "compile.log").write_text(built.stdout + built.stderr, encoding="utf-8")
        if built.returncode:
            raise RuntimeError(f"PREPLM compile failed; see {run / 'compile.log'}")

        samples = cases()
        protocol = "".join(_protocol(c) for c in samples)
        (run / "inputs.txt").write_text(protocol, encoding="ascii")
        executed = subprocess.run(
            [str(run / "preplm_driver.exe")], cwd=run, input=protocol,
            text=True, capture_output=True, timeout=60,
        )
        (run / "outputs.txt").write_text(executed.stdout, encoding="ascii")
        (run / "stderr.txt").write_text(executed.stderr, encoding="utf-8")
        if executed.returncode:
            raise RuntimeError("PREPLM driver failed; see stderr.txt in the run folder.")
        lines = executed.stdout.splitlines()
        if len(lines) != 4 * len(samples):
            raise ValueError("PREPLM driver returned an unexpected number of lines.")

        atol = float(TOLERANCE["atol"])
        rtol = float(TOLERANCE["rtol"])
        results = []
        max_abs = 0.0
        for i, case in enumerate(samples):
            block = lines[4 * i:4 * i + 4]
            f_lstk, f_lpbl = (int(v) for v in block[0].split())
            f_tstk, f_wstk = (float(v) for v in block[1].split())
            f_dthdz = tuple(float(v) for v in block[2].split())
            f_wspd = tuple(float(v) for v in block[3].split())
            py = preplm(
                case["fireflg"], case["emlays"], case["hmix"], case["hts"],
                case["psfc"], case["ts"], case["ddzf"], case["qv"], case["ta"],
                case["uw"], case["vw"], case["zh"], case["zf"], case["pres"],
            )
            scalar_pairs = [(py.stack_temperature, f_tstk), (py.stack_wind, f_wstk)]
            scalar_pairs += list(zip(py.theta_gradient, f_dthdz, strict=True))
            scalar_pairs += list(zip(py.wind_speed, f_wspd, strict=True))
            if not all(math.isfinite(f) for _, f in scalar_pairs):
                raise ValueError(f"Nonfinite PREPLM output for {case['label']}")
            errors = [abs(p - f) for p, f in scalar_pairs]
            case_max = max(errors)
            max_abs = max(max_abs, case_max)
            index_ok = py.lstk == f_lstk and py.lpbl == f_lpbl
            passed = index_ok and all(
                e <= atol + rtol * abs(f) for e, (_, f) in zip(errors, scalar_pairs)
            )
            results.append({
                "label": case["label"], "emlays": case["emlays"], "fireflg": case["fireflg"],
                "python_lstk_lpbl": [py.lstk, py.lpbl], "fortran_lstk_lpbl": [f_lstk, f_lpbl],
                "indices_match": index_ok, "max_absolute_error": case_max, "passed": passed,
            })
        report.update({
            "cases": len(samples),
            "values": sum(2 + 2 * c["emlays"] for c in samples),
            "results": results, "max_absolute_error": max_abs,
            "success": all(row["passed"] for row in results),
        })
        report["source_status_after"] = git(source, "status", "--porcelain=v1")
        report["source_unchanged"] = (
            sha256(original) == sha256(run / "preplm.f")
            and report["source_status"] == report["source_status_after"]
        )
        if not report["source_unchanged"]:
            report["success"] = False
            raise ValueError("Original preplm.f changed during validation.")
        if not report["success"]:
            raise ValueError("PREPLM comparison failed.")
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
