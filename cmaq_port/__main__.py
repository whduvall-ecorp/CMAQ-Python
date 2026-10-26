"""Command-line entry points for the staged CMAQ workflow."""

import argparse
import subprocess
import sys
from pathlib import Path

from .comparison import compare
from .fire_rise import fire_plmris
from .interpolation import poly
from .inventory import inventory, save_report
from .plume import delta_zs, plsprd
from .preplm import preplm
from .rise import plmris

LESSON = """CMAQ source walkthrough -- no model simulation is performed.

1. Meteorology: regional model -> PREP/mcip -> CCTM meteorological inputs.
2. Emissions: external processing plus configured inline streams/DESID.
3. Initial/boundary conditions: PREP/icon and PREP/bcon.
4. CCTM/src/driver/cmaq_main.F -> CMAQ_DRIVER in driver.F -> SCIPROC.
5. sciproc.F orchestrates transport, diffusion, photolysis, clouds, chemistry,
   and aerosols, subject to configuration and substeps. This is not an
   interchangeable list of independent plume factors.
6. CCTM/src/plrise/smoke contains point-source plume-rise implementations.
   Similar physical names in AERMOD do not establish identical algorithms.
7. CONC/ACONC and ELMO/AELMO output choices must be matched to time semantics,
   species definitions, units, and vertical layers before analysis.

Read C:\\Dev\\CMAQ\\DOCS\\Users_Guide\\Tutorials\\CMAQ_UG_tutorial_tracers.md.
Chemically inert does not automatically mean non-depositing: tracer surrogate
settings can enable dry deposition, scavenging, or aqueous processing.

AERMOD estimates source impacts at receptors; CMAQ evolves concentrations
on grid cells/layers over time. A receptor maximum is not a grid-cell mean.
Start comparison with matched passive-tracer source increments, not total
ozone/PM versus a single-source impact. Match horizontal AND vertical support.

Current tools inventory source and compare already aligned JSON datasets.
The Python plume routines PLSPRD, DELTA_ZS, PLMRIS, FIRE_PLMRIS, and PREPLM are
checked against standalone unchanged Fortran. The I/O API POLY interpolation
PREPLM uses is checked against the unchanged upstream poly.f.
No netCDF adapter, full CMAQ executable build, or full-model benchmark has been
validated in this project yet.
"""


def numerical_lesson() -> str:
    result = plsprd((0.0, 0.01, 0.02), (0, 100, 500, 1000), 150)
    clipped = plsprd((0.0, 0.01, 0.02), (0, 100, 500, 1000), 1500)
    geometry = delta_zs((2, 0), (10, 100, 300), ((100, 120), (500, 600)))
    grid = (0, 20, 60, 140, 300, 600)
    centers = (10, 40, 100, 220, 450)
    gradients = (0, 0.02, -0.005, 0.0001, 0.01)
    stack_geometry = delta_zs((0,), (20,), tuple((height,) for height in grid[1:]))
    rise = plmris(
        0, 0.03, 800, 2, 20, 450, 10, 290, 0.2, gradients,
        tuple(290 - 0.006 * height for height in centers), (0.2, 1, 5, 8, 2),
        grid, centers, 4, tuple(row[0] for row in stack_geometry.relative_heights),
    )
    spread = plsprd(gradients, grid, rise.centerline)
    fire = fire_plmris(
        0, 0.03, 800, 0.4, 450, 0.2, gradients,
        tuple(290 - 0.006 * height for height in centers), (0.2, 1, 5, 8, 2),
        grid, 0.2,
    )
    stack_top = 45.0
    heights_window = (10.0, 30.0, 60.0, 120.0)
    temperature_window = tuple(290.0 - 0.0065 * h for h in heights_window)
    stack_temperature = poly(stack_top, heights_window, temperature_window, 3)
    preplm_emlays = 6
    preplm_zf = (20.0, 60.0, 140.0, 300.0, 500.0, 800.0)
    preplm_zh = (10.0, 40.0, 100.0, 220.0, 400.0, 650.0)
    preplm_result = preplm(
        False, preplm_emlays, 500.0, 45.0, 1013.0, 291.5,
        tuple(1.0 / 40.0 for _ in range(preplm_emlays)),
        tuple(0.009 - 0.0005 * i for i in range(preplm_emlays)),
        tuple(293.0 - 0.0065 * h for h in preplm_zh),
        tuple(1.0 + 0.7 * i for i in range(preplm_emlays)),
        tuple(0.5 + 0.2 * i for i in range(preplm_emlays)),
        preplm_zh, preplm_zf,
        tuple(1013.0 - 40.0 * i for i in range(preplm_emlays + 1)),
    )
    return (
        "\n8. PLSPRD calculates initial plume bounds, not concentrations.\n"
        "   ZF=(0,100,500,1000) m; DTHDZ=(0,0.01,0.02) K/m; effective height=150 m.\n"
        f"   Calculated PLTOP={result.top:.9g} m; PLBOT={result.bottom:.9g} m.\n"
        f"   Height=1500 m: capped top={clipped.top:g} m, bottom={clipped.bottom:g} m.\n"
        "   This alternative is called when IPVERT is nonzero, not the default Turner branch.\n"
        "   PLSPRD retains default REAL rounding; its returned bounds are not AERMOD SIGZ.\n"
        "   Mixing height is not an active input; model-top clipping is not lid reflection.\n"
        "   Read docs/plume-spread.md for formulas, interfaces, and limitations.\n"
        "\n9. DELTA_ZS maps global stack heights onto local layer columns.\n"
        "   Zero-based map=(2,0); stack heights=(10,100,300) m.\n"
        "   Upper interfaces by layer=((100,120),(500,600)) m; no ground row.\n"
        f"   ZSTK by layer={geometry.relative_heights} m (negative values are retained).\n"
        f"   DDZF in first local column={geometry.inverse_thicknesses[0][0]:.9g},"
        f" {geometry.inverse_thicknesses[1][0]:.9g} 1/m.\n"
        "   First thickness is ZF(1), subsequent thicknesses are interface differences.\n"
        "   The source map indexes stacks, NOT the already-local meteorology columns.\n"
        "   Geometry only: no plume rise, layer fractions, or concentrations are calculated.\n"
        "   Read docs/stack-layer-geometry.md for indexing and original-code evidence.\n"
        "\n10. DELTA_ZS -> PLMRIS -> PLSPRD is a connected plume example.\n"
        "   HS=20 m, diameter=2 m, exhaust=450 K at 10 m/s; ambient stack temperature=290 K.\n"
        "   HFX=0.03, HMIX=800 m, USTAR=0.2 m/s, WSTK=4 m/s.\n"
        "   HFX is m K/s, not W/m2; PREPLM remains upstream of these explicit met inputs.\n"
        f"   PLMRIS centerline={rise.centerline:.9g} m; reporting WSTK={rise.stack_wind:g} m/s.\n"
        f"   IQ path={rise.stability_path}; layer transitions={rise.layer_transitions}.\n"
        f"   Alternative PLSPRD top={spread.top:.9g}, bottom={spread.bottom:.9g} m.\n"
        "   This is the nonzero-IPVERT spread branch, not the default Turner thickness.\n"
        "   IQ=1 unstable, 2 neutral, 3 stable, 4 momentum; path is a Python teaching diagnostic.\n"
        "   No emission fractions, fire plume rise, or concentrations are calculated.\n"
        "   Read docs/plume-rise.md for the loop, source quirks, and checked evidence.\n"
        "\n11. FIRE_PLMRIS is a distinct fire-source calculation.\n"
        "   It receives caller-supplied BFLX=0.4 m^4/s^3; it does not derive it from stack dimensions.\n"
        "   It receives HFX=0.03 m K/s, TSTK=450 K, USTAR=0.2 m/s, and fire wind=0.2 m/s.\n"
        f"   Calculated centerline={fire.centerline:.9g} m; reporting fire wind={fire.fire_wind:g} m/s.\n"
        f"   IQ path={fire.stability_path}; layer transitions={fire.layer_transitions}.\n"
        "   The fire-wind floor is 1 m/s; plume-centerline conversion is 2/3*ZPLM with no stack height.\n"
        "   HMIX<=0 selects the above-mixing-height branch; BFLX=0 is not given a positive floor.\n"
        "   The fire caller handles fire-specific BFLX; this routine does not invent preprocessing.\n"
        "   Read docs/fire-plume-rise.md for formulas, the distinct source path, and original-code evidence.\n"
        "\n12. POLY is the I/O API interpolation PREPLM uses, not a CMAQ-tree routine.\n"
        "   Newton divided-difference interpolation; here degree 3 on heights=(10,30,60,120) m.\n"
        f"   Ambient temperature interpolated to stack top {stack_top:g} m = {stack_temperature:.9g} K.\n"
        "   PREPLM calls POLY for TSTK and WSTK; this is the interpolation only, not PREPLM itself.\n"
        "   POLY needs NDEG+1 distinct points and NDEG>=1; high degree is numerically unstable.\n"
        "   The Python port matches the unchanged I/O API poly.f bitwise; see docs/interpolation.md.\n"
        "\n13. PREPLM prepares the stack-top meteorology PLMRIS consumes.\n"
        "   6-layer synthetic profile; stack height 45 m, mixing height 500 m, non-fire.\n"
        f"   LSTK={preplm_result.lstk}, LPBL={preplm_result.lpbl} (1-based layer indices).\n"
        f"   Interpolated stack-top TSTK={preplm_result.stack_temperature:.9g} K,"
        f" WSTK={preplm_result.stack_wind:.9g} m/s.\n"
        "   Wind speed, DTHDZ gradient, and POLY interpolation to stack top; single-precision throughout.\n"
        "   The commented-out surface-gradient term (ES/QSFC/TVSFC/THETG) is omitted as dead code.\n"
        "   The Python port matches the unchanged compiled preplm.f bitwise; see docs/preplm.md.\n"
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    source = commands.add_parser("inventory", help="Hash tracked local source; do not execute CMAQ scripts.")
    source.add_argument("--source", type=Path, default=Path(r"C:\Dev\CMAQ"))
    commands.add_parser("lesson", help="Print the source-guided learning walkthrough.")
    comparison = commands.add_parser("compare", help="Compare declared-aligned passive-tracer JSON samples.")
    comparison.add_argument("cmaq", type=Path)
    comparison.add_argument("aermod", type=Path)
    args = parser.parse_args()
    try:
        if args.command == "lesson":
            print(LESSON)
            print(numerical_lesson())
        elif args.command == "inventory":
            report = inventory(args.source)
            if report["warnings"]:
                print(f"WARNING: {report['warnings']}", file=sys.stderr)
            print(save_report(report, "inventory"))
        else:
            print(save_report(compare(args.cmaq, args.aermod), "comparison"))
    except (OSError, ValueError, subprocess.SubprocessError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
