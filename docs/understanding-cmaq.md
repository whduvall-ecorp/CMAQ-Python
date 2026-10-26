# Understanding CMAQ alongside AERMOD

This walkthrough is grounded in the local source at the pinned inventory
commit. The supplied diagrams illustrate system structure, process budgets,
emissions, domain decomposition, and output choices. They are documentation
figures, not benchmark output produced by this computer. We have not modified
or redistributed the read-only attachments.

## The two models answer different questions

| Aspect | AERMOD | CMAQ CCTM |
| --- | --- | --- |
| Basic representation | Steady-state plume framework, driven hour by hour | Time-evolving Eulerian grid/layer concentrations |
| Typical interpretation | Source impacts at receptors | Regional concentration fields and process budgets |
| Meteorology | AERMET-processed surface/profile information | Gridded meteorology, commonly WRF processed through MCIP |
| Chemistry | Selected pollutant treatments/options | Mechanism-dependent gas, aerosol, and cloud chemistry |
| Outputs to reconcile | Receptor values/averages | Cell/layer values, instantaneous or time averaged |

This is a conceptual distinction, not a statement that every AERMOD/CMAQ
option has identical scope or the same required inputs. AERMOD point-source
impact cannot be directly ranked against a regional total ozone or PM map.

## Follow the original system

The local [CCTM overview](C:/Dev/CMAQ/CCTM/README.md)
describes the chemistry-transport component. Other source paths below
are relative to `C:\Dev\CMAQ`.

```text
Meteorology model -> PREP/mcip -------------------------+
Emissions processing / configured inline emissions ----+
PREP/icon + PREP/bcon (initial/boundary conditions) -----+
Mechanism, namelists, deposition, land/ocean inputs ----+
                                                       v
CCTM/src/driver/cmaq_main.F
  CMAQ_DRIVER in CCTM/src/driver/driver.F
    SCIPROC in CCTM/src/driver/sciproc.F
      transport / diffusion / photolysis / clouds / chemistry / aerosols
    concentration and configured diagnostic output
                                                       v
POST, PYTOOLS, and user analysis
```

This is a reading path, not a complete call graph or exact operator schedule.
SCIPROC contains configuration branches and substeps; transport and process
order must be read in context before porting anything.

## Source reading map

| Directory/file | What to inspect |
| --- | --- |
| `CCTM/src/driver/cmaq_main.F` | Program entry, initialization and time stepping |
| `CCTM/src/driver/driver.F` | Model-driver orchestration |
| `CCTM/src/driver/sciproc.F` | Science process calls and conditional paths |
| `CCTM/src/hadv`, `vadv`, `hdiff`, `vdiff` | Transport and mixing |
| `CCTM/src/plrise/smoke` | Point-source plume rise and vertical distribution |
| `CCTM/src/emis` | Emissions handling and DESID-related processing |
| `CCTM/src/gas`, `MECHS`, `phot` | Mechanisms, solvers, photolysis |
| `CCTM/src/aero`, `cloud`, `depv` | Aerosols, cloud processing, deposition |
| `CCTM/src/driver/ELMO_PROC.F` | Explicit/lumped output quantities |
| `CCTM/src/driver/AVG_CONC.F`, `wr_aconc.F` | Averaging and output handling |
| `CCTM/scripts/bldit_cctm.csh` | Build selections |
| `CCTM/scripts/run_cctm_Bench_2018_12NE3_CRACMM2.csh` | Example benchmark run configuration |

The `PYTOOLS` folder already contains Python utilities such as `shp2cmaq`
and `dmschlo`. These are useful preprocessing tools, not a Python CCTM
implementation. Reuse existing utilities where applicable rather than
mistaking a new wrapper for ported atmospheric science.

## Trace one process before translating it

For an initial bounded kernel, inspect the plume-rise directory and identify:

1. Its selected implementation in the build script.
2. Units, types/precision, constants, argument and module dependencies.
3. Callers and the state they establish.
4. Inactive/conditional outputs and error reporting.
5. A way to compile the unchanged original with an explicit test adapter.

A routine with the same physical label as an AERMOD routine need not use the
same formula, state, or resolution. Compare behavior only after documenting
those distinctions. The first port,
[PLSPRD](plume-spread.md), now provides a checked example: it derives initial
plume bounds in the caller's alternate nonzero-IPVERT branch, not the default
Turner branch or a full emissions-allocation calculation.
The second port, [DELTA_ZS](stack-layer-geometry.md), traces physical stack
heights into already-local meteorology columns. Its ZF omits the ground row
that PLSPRD requires; study that interface distinction before connecting them.
The full [PLMRIS port](plume-rise.md) now connects those helpers in executable
lesson section 10. Follow the initial rise branch, residual buoyancy as each
layer is crossed, and final plume-top-to-centerline conversion. The
[PREPLM port](preplm.md) now supplies that upstream meteorology (wind speed,
potential-temperature gradient, and stack-top temperature/wind); no layer
emission fractions are calculated.
The separate [FIRE_PLMRIS port](fire-plume-rise.md) has no stack geometry and
receives caller-computed `BFLX`; its stability tests, residual formula, wind
handling, and final centerline conversion are not copied from PLMRIS.

## Tracer does not necessarily mean passive in every process

Read `DOCS\Users_Guide\Tutorials\CMAQ_UG_tutorial_tracers.md`.
Its examples can use dry-deposition/scavenging surrogates and other mappings.
For the first paired transport study, explicitly configure:

- No gas/aqueous chemistry or aerosol conversion for the selected tracer.
- No dry/wet removal, unless deliberately matched in both models.
- Explicit source emissions, initial conditions, boundary conditions, and
  baseline subtraction; do not assume unspecified defaults equal zero.
- Transport on, with documented molecular weight and concentration units.

The repository's example tracers are instructional starting points, not a
ready-made matched AERMOD comparison.

## Exercises

1. Locate CMAQ_DRIVER's call from `cmaq_main.F`, then SCIPROC's call from
   `driver.F`. Separate declarations from executable calls.
2. Find HADV, ZADV, HDIFF, VDIFF, PHOT, CLDPROC, CHEM, and AERO in SCIPROC.
   Which are conditional or substepped?
3. Compare CONC versus ACONC timing and ELMO/AELMO species definitions.
4. Inspect the tracer tutorial's deposition/scavenging surrogate settings.
5. Compare one plume-rise routine's required state with AERMOD's PDF/HEFF/
   DELTAH reading path. Record differences before discussing equivalence.
6. Explain why averaging AERMOD's horizontal receptor grid at ground level
   still does not automatically create CMAQ's first-layer volume mean.

Run `python -m cmaq_port lesson` for the short executable walkthrough.
