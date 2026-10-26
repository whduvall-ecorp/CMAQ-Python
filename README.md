# CMAQ Python learning and validation workflow

An incremental companion to the Desktop AERMOD-Python project. This project
adds source provenance, a runnable learning guide, a staged Python workflow,
and an initial aligned passive-tracer comparison tool.

**PLSPRD, DELTA_ZS, PLMRIS, and FIRE_PLMRIS now match unchanged Fortran
within declared default-REAL tolerances in 1,308 synthetic cases (9,272
values); maximum absolute difference is 0.000030517578125. The I/O API `POLY`
interpolation and `PREPLM` (stack-top meteorology) match the unchanged Fortran
**bitwise** (0.0 difference) in 107 and 72 cases respectively — `PREPLM`
completes the science routines in `plrise/smoke`. No full CMAQ run has been
validated. The examples are component/tool tests, not model results.**

## Local source inspected

- Source checkout: `C:\Dev\CMAQ` (left unchanged).
- README identity: CMAQv5.5+, an evolving bugfix/development code line.
- Git branch: `cloud-cost-pilot`.
- Commit: `e10d0a04f319d6be6db086bc66ba13ce4044987b`.
- Git reports a clean working tree, but this is not a complete Linux checkout.
- Inventory: 1,731 tracked entries; 42 skip-worktree files absent on disk;
  104 Git symbolic links represented as ordinary Windows files rather than
  materialized symbolic links. Both limitations are explicit in the report.

A clean Git status does not prove all build inputs exist in usable form.
Do not assume this checkout is ready to compile CCTM. The source README's
version name is not a substitute for a pinned commit.

## Commands

From this folder, using Python 3.12 and Git (no extra Python dependencies):

```powershell
python -m cmaq_port inventory --source 'C:\Dev\CMAQ'
python -m cmaq_port lesson
python -m cmaq_port.kernel_validation
python -m unittest discover -s tests -v

# Synthetic metric demonstration only: neither model was run.
python -m cmaq_port compare '.\examples\synthetic-cmaq.json' '.\examples\synthetic-aermod.json'
```

Inventory and comparison reports go into new isolated folders under `runs`.
Each inventory hashes all tracked entries, distinguishes disk files from
stored Git blobs/link targets, and records branch/commit/status. It does not
hash untracked data. Unexpected missing tracked files fail explicitly.

The comparison CLI only accepts already harmonized JSON datasets. It does
not read CMAQ netCDF, convert units, perform spatial averaging, or run models.
See [the comparison protocol](docs/comparison-protocol.md) before using it.

## Learning and implementation documents

- [Understanding CMAQ](docs/understanding-cmaq.md): source map and exercises.
- [PLSPRD plume-spread lesson](docs/plume-spread.md): first Python science
  kernel, source attribution, formulas, and compiled evidence.
- [DELTA_ZS stack/layer geometry](docs/stack-layer-geometry.md): zero-based
  source maps, local columns, upper interfaces, and default-REAL arithmetic.
- [PLMRIS layered plume rise](docs/plume-rise.md): stability branches,
  residual-buoyancy iteration, source quirks, and the connected lesson.
- [FIRE_PLMRIS fire plume rise](docs/fire-plume-rise.md): caller-owned fire
  buoyancy flux, distinct stability/residual path, source diagnostics, and
  compiled validation.
- [POLY interpolation](docs/interpolation.md): the I/O API routine PREPLM uses,
  ported bitwise-faithfully and validated against the unchanged upstream `poly.f`.
- [PREPLM stack-top meteorology](docs/preplm.md): wind, potential-temperature
  gradient, and POLY interpolation to stack top, with side-by-side Fortran and
  a bitwise differential against the unchanged `preplm.f`.
- [Workflow plan](docs/workflow-plan.md): milestones, dependencies, verification gates.
- [CMAQ versus AERMOD comparison protocol](docs/comparison-protocol.md):
  what can be compared and what must be aligned first.

The [AERMOD learning project](../AERMOD-Python/README.md) remains separate.
Neither project's educational Python functions constitute a validated,
complete replacement for its original model.

## Baseline evidence

- **104 unit tests pass**, including 19 dedicated PLMRIS, 19 dedicated
  FIRE_PLMRIS, 10 POLY, and 12 PREPLM tests,
  PLSPRD boundaries/precision/guards,
  DELTA_ZS mapping/shapes/reciprocals,
  source exclusions/link stubs, provenance,
  exact pairing, invalid metadata/time/value rejection, and metric arithmetic.
  Python compilation checks pass; Pylance reports no errors in the new modules.
- [Source inventory](runs/20261024T172500-inventory-0b2743f9/inventory.json).
- [Five-operation Fortran comparison](runs/20261024T225921-kernels-e99396a9/kernels.json):
  1,308 cases, 9,272 values, maximum absolute difference 0.000030517578125,
  including 224 FIRE_PLMRIS cases and 488 connected geometry/rise/spread calls.
  Prior PLSPRD/DELTA_ZS/PLMRIS/connected results remain exact. Uses the existing
  AERMOD project-local compiler; no full CMAQ dependencies are required.
- POLY differential: `python -m cmaq_port.poly_validation` compiles the
  unchanged I/O API `poly.f` and matches the Python port **bitwise** (0.0
  maximum absolute difference) over 107 cases; the reference is vendored with
  provenance in [reference/ioapi/PROVENANCE.md](reference/ioapi/PROVENANCE.md).
- PREPLM differential: `python -m cmaq_port.preplm_validation` compiles the
  unchanged `preplm.f` with a test-only UTILIO_DEFN adapter (real POLY, trivial
  M3MSG2) and matches the Python port **bitwise** (0.0 difference) over 72
  cases / 1,008 values, with LSTK/LPBL exact and source verified unchanged.
- [Synthetic two-pair comparison](runs/20261024T172500-comparison-98e5b195/comparison.json):
  mean bias 0, MAE 2, RMSE 2, maximum absolute difference 2 ug/m3.
  These invented values test metric arithmetic only.

Build prerequisites in the local CMAQ documentation include a Linux-style
shell, compatible Fortran/C compilers, netCDF C/Fortran, I/O API, and MPI
for parallel runs. The user is handling WSL installation in another chat;
this project does not install or configure it. No WSL installation,
global dependency changes, benchmark download, or source edits were made.

Continue alternating bounded CMAQ and AERMOD routines, with executable
lessons, independent review, and original-code validation. This CMAQ update
adds FIRE_PLMRIS; recent AERMOD DISTF/platform work is separate and is not
included in these counts. Neither project is a complete,
validated replacement for its original model.
