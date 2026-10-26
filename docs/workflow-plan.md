# Staged CMAQ Python workflow and comparison roadmap

> Part of the program-level
> [AERMOD + CMAQ Python port and comparison roadmap](../../model-port-and-comparison-roadmap.md)
> (Track B). That master plan covers how this CMAQ work connects to the AERMOD
> port and the eventual cross-model comparison.

## Objective and boundaries

Follow the same evidence-driven pattern as AERMOD-Python: preserve reference
source, explain routines, run the original implementation, port bounded
numerical pieces, and verify each against the original before integration.
Do not attempt a whole-model automatic translation or infer regulatory
equivalence from matching a few kernels.

The foundation is implemented and the first bounded science port is checked.
The eventual first inter-model study is a passive-tracer source increment
with matched time, horizontal area, and vertical support.

## Stage 0: inspect and pin the local source -- implemented

- Source: `C:\Dev\CMAQ`; untouched by this project.
- README identifies CMAQv5.5+; branch is `cloud-cost-pilot`.
- Commit: `e10d0a04f319d6be6db086bc66ba13ce4044987b`.
- Baseline report hashes 1,731 tracked entries.
- 42 skip-worktree entries are absent from disk and hashed from Git index blobs.
- 104 tracked Git links are Windows regular files, not materialized symlinks.
- The manifest records `content_origin`, `available_in_worktree`, link targets,
  Git object IDs, actual content SHA-256, and explicit warnings.

These are facts about this checkout, not a supported-build certification.
Git-index link hashes describe stored target text, not target-file content.
Untracked data/builds are outside inventory scope. The source may change
later; reproduce the inventory and compare commit/hashes before validation.

## Stage 1: establish a supported reference environment -- in progress

A Linux environment is now available: `wsl --list --verbose` reports the
`Ubuntu` distribution in the `Running` state under WSL version 2, and
`wsl --status` confirms WSL2 as the default. This removes the earlier
"WSL not installed" blocker. The toolchain inside WSL still needs to be
inventoried and, where missing, installed before a build.

Interactive `wsl ... bash -lc` probes from this session have been slow and
sometimes hang; prefer running the environment checks directly in a WSL
terminal and recording the output, rather than driving long WSL commands
through this workflow.

Choose one of two reference-build hosts:

- **WSL2 (local, lowest cost):** suitable for building CMAQ, the serial
  benchmark, and capturing real caller states for kernel validation. Memory
  and I/O are bounded by the Windows host; large parallel runs may be slow.
- **Cloud spot instance (for larger/parallel runs):** e.g. a ~16-vCPU
  Google Cloud Spot VM. Indicative 2026 spot rates are roughly
  $0.32-0.33/hour for a 16-vCPU general-purpose VM, so a "<$50/month,
  ~$10/run" budget implies limited active hours, not continuous uptime
  (about 150 VM-hours per $50 at that rate). Spot VMs can be preempted at
  any time and carry no SLA; stopped VMs still incur persistent-disk charges,
  and benchmark input/output archives (tens of GB) add storage/egress cost.
  Do not create billable cloud resources without explicit approval.

Do not run the CMAQ configuration scripts in Windows PowerShell as if they
were native Windows build commands. The supplied C-shell scripts configure
and create directories/links; they are not read-only inspection tools.

In the selected Linux workspace (WSL2 or cloud):

1. Inventory the toolchain (git, gcc, gfortran, csh/tcsh, make, MPI,
   nc-config/nf-config) and install what is missing; record versions.
2. Obtain the pinned source with proper filenames and symbolic links.
   Keep this Windows checkout unchanged; do not reset it or flip branches.
   The Windows `C:\Dev\CMAQ` checkout is missing 42 skip-worktree files and
   stores 104 Git links as plain files, so a native Linux clone is required.
3. Build compatible C/Fortran compilers, netCDF C/Fortran, I/O API, and MPI
   as required by the selected serial/parallel configuration.
4. Record compiler/library versions, flags, mechanism, deposition scheme,
   namelists, link paths, CPU/platform, and MPI decomposition.
5. Build with the original scripts in a separate project/build location.
6. Preserve executable/library hashes and complete build/run logs.

Follow the local
[Linux configuration tutorial](C:/Dev/CMAQ/DOCS/Users_Guide/Tutorials/CMAQ_UG_tutorial_configure_linux_environment.md)
and [benchmark tutorial](C:/Dev/CMAQ/DOCS/Users_Guide/Tutorials/CMAQ_UG_tutorial_benchmark_cracmm2_stage.md).
The small Windows GNU Fortran setup used for isolated AERMOD tests does
not establish a compatible full CMAQ runtime.

## Stage 2: reproduce the official reference benchmark -- pending

Prerequisites: Stage 1 and agreed data download/storage budget.

The local documentation describes the 12NE3 benchmark for July 1-2, 2018:
100 columns, 105 rows, 35 layers, 12-km spacing. Its data catalog lists
approximately 6 GB of CRACMM2 input and 19 GB of output archives.
Those are documentation estimates, not a fresh check of remote sizes or
extracted storage. Alternative CB6r5 data/output combinations differ.

No benchmark data was downloaded here. Locate any existing local inputs
before acquiring more. Do not mix reference output mechanisms/deposition
schemes or assume a development branch is bitwise identical to release 5.5.

Required evidence:

- Version/configuration-matched inputs and expected outputs.
- Persistent per-run directories, timeouts, explicit run completion/errors.
- Input/executable/configuration hashes and captured logs.
- netCDF-aware checks of dimensions, grid, layers, time flags, species,
  units, masks/missing values, and selected numeric arrays.
- Per-variable engineering tolerances tied to precision/platform behavior.
  Report deviations honestly; no automatic tolerance inflation.
- Serial/parallel reproducibility check if using parallel CCTM.

Passing the benchmark validates that configuration on that case, not all
CMAQ science options and not equivalence with AERMOD.

## Stage 3: port and teach bounded CMAQ kernels -- started

Prerequisite: a compilable unchanged reference kernel and known caller state.

Initial investigation targets include plume rise/vertical allocation in
`CCTM\src\plrise\smoke`, followed by suitable transport helpers and species/
output conversions. PLSPRD is checked against unchanged standalone Fortran:
44 cases, 88 values, zero observed difference. See the
[plume-spread lesson](plume-spread.md). DELTA_ZS adds 64 cases and 720 values
with zero observed difference; see [stack/layer geometry](stack-layer-geometry.md).
PLMRIS adds 488 isolated and 488 connected geometry/rise/spread cases; see
[the plume-rise guide](plume-rise.md). FIRE_PLMRIS now adds 224 fire cases
and 672 checked values; see [the fire plume-rise guide](fire-plume-rise.md).
The current five-operation report covers 1,308 cases/9,272 values, with a
maximum absolute difference of 0.000030517578125 under declared default-REAL
tolerances. Prior routines remain exact. Other routines remain candidates.
Search prior utilities and library implementations before adding new logic.

The I/O API `POLY` interpolation and `PREPLM` (stack-top meteorology) are both
ported and validated. `PREPLM` calls `POLY` to interpolate ambient temperature
and wind to stack top. `POLY` is defined nowhere in the CMAQ source tree; it is
the Models-3/EDSS I/O API library function
`REAL FUNCTION POLY(X, XPTS, YPTS, NDEG)`, a documented Newton
divided-difference polynomial interpolation routine (degree `NDEG`, using
`NDEG+1` points; recommended `NDEG <= 5`). See the
[I/O API POLY documentation](https://www.cmascenter.org/ioapi/documentation/all_versions/html/POLY.html),
the [interpolation port guide](interpolation.md), and the
[PREPLM port guide](preplm.md).

`POLY` (`cmaq_port/interpolation.py`) matches the unchanged upstream `poly.f`
**bitwise** (0.0 maximum absolute difference) over 107 cases via
`python -m cmaq_port.poly_validation`; the reference is vendored with
provenance in [reference/ioapi/PROVENANCE.md](../reference/ioapi/PROVENANCE.md).

`PREPLM` (`cmaq_port/preplm.py`) matches the **unchanged compiled** `preplm.f`
bitwise (0.0 difference) over 72 cases / 1,008 values via
`python -m cmaq_port.preplm_validation`, which compiles the original with a
test-only `UTILIO_DEFN` adapter (exposing the unchanged POLY and a trivial
`M3MSG2`, no stubbed physics); `LSTK`/`LPBL` match exactly and the source is
verified unchanged by hashing. This completes the science routines in
`plrise/smoke` (PLSPRD, DELTA_ZS, PLMRIS, FIRE_PLMRIS, PREPLM); the remaining
files there are I/O/diagnostic. The `POLY3`/`POLY3B` routines in
`aero6/isocom.f` are cubic-equation solvers and are unrelated to this
interpolation `POLY`.

For each selected routine:

1. Read callers, types, units, module/IO dependencies, and error semantics.
2. Define explicit typed Python inputs and outputs; preserve unassigned state.
3. Add boundary/algebraic tests and explanatory source-variable notes.
4. Compile unchanged original code with a small, clearly test-only driver.
5. Compare deterministic synthetic cases and states captured from a real run.
6. Record compiler/source/driver identity and every comparison result.
7. Add Markdown explanations and a runnable lesson section.

Alternate one CMAQ routine then one AERMOD routine each round. FIRE_PLMRIS
and the concurrent AERMOD DISTF/platform work are validated in their separate
projects and reports. Use an independent
high-capability model review for new numerical work, then resolve findings
and rerun the original-code checks. The FIRE_PLMRIS review remains pending.
None of these isolated CMAQ routines needs the full CMAQ
runtime; Stage 1/2 still gate whole-model experiments. Captured real-run
states remain outstanding for all three routines. PLMRIS imports only LOGDEV
from a clearly test-only runtime adapter; unchanged CONST.EXT supplies physics
constants. This is not a production runtime implementation.

Only then connect kernels. Mass conservation, positivity/limiters, timestep
constraints, boundary exchange, and solver convergence need dedicated tests.
Do not quietly replace MPI/data/model modules with success-shaped stubs.

## Stage 4: output adapters and comparison contract -- partially implemented

Implemented: strict normalized JSON ingestion, pairing, metrics, hashed
input provenance, and explicit rejection of mismatched declarations.
Synthetic fixtures prove the tool mechanics only.
The full project now has 104 passing unit tests, including 11 PLSPRD, 11
DELTA_ZS, 19 PLMRIS, 19 FIRE_PLMRIS, 10 POLY, and 12 PREPLM tests;
the original 22 cover the foundation. Python compilation checks pass.

Pending:

- CMAQ netCDF extraction with verified CONC/ACONC or ELMO/AELMO time semantics.
- AERMOD output extraction tied to its existing reference harness.
- Species/mass conversions with molecular weights and applicable air state.
- Spatial/vertical aggregation and coordinate/projection checks.
- Machine-readable forcing/source configuration manifests linked to each
  contract ID; the current tool checks declarations, not the physics.

Reuse AERMOD's tested output parser where appropriate rather than making an
incompatible duplicate. Keep each adapter independently tested.

## Stage 5: matched passive-tracer study -- pending

Prerequisites: both original model baselines, adapters, and an agreed case.
Use original executables first; neither Python port is a full model today.

Define a simple matched source/forcing scenario, explicitly configure inert
transport with no removal, and compare source-on minus source-off increments
or an equivalently isolated tracer. Account for initial/boundary conditions,
spin-up, source plume rise and layer injection, boundary flux, and transit
time. AERMOD's hourly steady-state assumptions and CMAQ's transient transport
will remain different even with harmonized forcing.

Report paired concentration metrics, spatial patterns, integrated tracer
mass/flux where meaningful, and grid/timestep/sampling sensitivity. Use
observations or independent solutions to evaluate accuracy; similarity
between two models is not accuracy.

See [the detailed comparison protocol](comparison-protocol.md).

## Later studies

Only after passive transport is characterized: primary PM/removal,
NO/NO2 with documented chemistry, and regional background plus local
increments. Secondary ozone/aerosol chemistry, attribution, and sensitivities
need their own design; ISAM and DDM results are not automatically equivalent
to brute-force source-on/off changes.

## Resume point

Continue from this folder's README and inventory. Preserve the CMAQ checkout
and AERMOD validation work. WSL2 (Ubuntu) is now available as the Linux host;
the immediate next steps are (1) inventory/install the WSL build toolchain and
produce a native Linux clone of the pinned source, and (2) continue bounded CMAQ
kernel ports now that the `plrise/smoke` science routines (PLSPRD, DELTA_ZS,
PLMRIS, FIRE_PLMRIS, PREPLM) and the I/O API `POLY` are validated. Continue
paired bounded-kernel lessons in parallel. After the Linux build and data are in
place, establish the CMAQ reference benchmark before claiming any
model-result comparison.
