# First CMAQ science lesson: PLSPRD

## Source and attribution

Python: [plume.py](../cmaq_port/plume.py). Original:
[plsprd.f](C:/Dev/CMAQ/CCTM/src/plrise/smoke/plsprd.f), from pinned commit
`e10d0a04f319d6be6db086bc66ba13ce4044987b`.
The routine describes a modification of Gillani's model, initially provided
by J. Godowitch (EPA, September 2003); D. Wong revised array declarations in
August 2015. Its SMOKE attribution is **COPYRIGHT (C) 2003, MCNC Environmental
Modeling Center. All Rights Reserved**, with a reference to a SMOKE COPYRIGHT
file. Retain this attribution and the source's EPA permission notice when
redistributing derived work. The upstream root license is preserved in
[LICENSE-CMAQ.md](../LICENSE-CMAQ.md). Each validation run retains the complete
unchanged source, including notices, and a copy of that license.

## What it does -- and does not do

PLSPRD calculates initial plume top and bottom, not a concentration,
Gaussian density, emission fraction, or complete plume-rise solution.
The actual caller is [PT3D_DEFN.F](C:/Dev/CMAQ/CCTM/src/emis/emis/PT3D_DEFN.F).
PLMRIS or FIRE_PLMRIS first establishes the effective height. The caller's
default `IPVERT == 0` branch uses a Turner thickness expression; **PLSPRD is
the alternate nonzero-IPVERT path**, not the universal CMAQ injection method.
The caller subsequently allocates emissions to layers under its
uniform-in-pressure assumption and performs additional error checks.

No mixing-height argument is active: HTMIX is commented out in the original.
Do not invent a mixing-lid reflection or equate these bounds to AERMOD SIGZ.

## Inputs and equations

| Python input | Fortran input | Meaning |
| --- | --- | --- |
| `gradients` | DTHDZ(1:KZ) | Potential-temperature gradients, K/m |
| `layer_heights` | ZF(0:KZ) | Full-layer interfaces in m, including ground |
| `effective_height` | CEFSTK | Upstream effective stack/plume height, m |

Select the first layer whose upper interface is at or above CEFSTK; use the
last layer when the effective height exceeds the model top. Equality belongs
to the lower layer.

```text
DTDZ = DTHDZ(selected layer) - 0.0098
SIGZ0 = max(10, 15 * exp(-117 * DTDZ))
DPTH = 3.545 * SIGZ0
PLTOP = min(ZF(KZ), CEFSTK + DPTH/2)
PLBOT = min(ZF(KZ)-1, max(0, CEFSTK - DPTH/2))
```

Clipping order follows the source. A center above the model can return top
1000 m and bottom 999 m on a 1000-m domain; that is not a physically validated
emissions allocation. Full caller checks are not replaced by this function.

The original uses default REAL, unlike AERMOD's mostly double-precision
kernels. Python explicitly rounds intermediate arithmetic/constants to
binary32 using the standard library. Its exponential uses Python's math
library followed by binary32 rounding; do not assume bitwise identity across
compilers and math libraries.

Python guards require finite, representable inputs, KZ gradients, KZ+1
interfaces increasing after rounding, ground interface zero, model top at
least 1 m, and nonnegative center. Invalid input/overflow raises ValueError;
these guards are not CMAQ's complete input-screening system.

## Try and explain

Run `python -m cmaq_port lesson`.

For gradients `(0, 0.01, 0.02)` K/m, interfaces `(0,100,500,1000)` m,
and center 150 m, the selected layer is the second. The calculated top is
175.97258 m and bottom 124.02742 m (rounded for display).

Exercises:

1. Set all gradients to 0.0098: why is SIGZ0 exactly 15 m?
2. Compare centers 100 and 100.01 m with different first/second gradients.
3. Use gradient 0.1 K/m: identify the 10-m minimum SIGZ0.
4. Move the center to the ground, then above the model. Track each clipping
   operation separately.
5. Explain why none of these outputs can be directly compared to an AERMOD
   receptor concentration.

## Fortran and Python side by side

The verbatim core computation of `SUBROUTINE PLSPRD`
([plsprd.f](C:/Dev/CMAQ/CCTM/src/plrise/smoke/plsprd.f), lines 86-109, with
the parameters `SZ0FAC = 3.545`, `SPRFAC = 15.`, and `GAMMA = -0.0098`
declared at lines 73-75):

```fortran
C Get ambient temperature above plume rise height (effective stack height)
      K = 0
      DO
         K = K + 1
         IF ( K .EQ. KZ .OR. CEFSTK .LE. ZF( K ) ) EXIT
      END DO
      DTDZ  = DTHDZ( K ) + GAMMA

C Compute initial vertical spread
      SIGZ0 = MAX( 10.0, SPRFAC * EXP( -117.0 * DTDZ ) )
      DPTH  = SZ0FAC * SIGZ0

C Compute plume top and bottom heights; plume is either completely
C within or outside mixing layer
      PLTOP = CEFSTK + DPTH / 2.0
      PLBOT = CEFSTK - DPTH / 2.0

C Make sure plume bottom is at least zero
      PLBOT = MAX( 0.0, PLBOT )

C Make sure that plume top and bottom height < ZF (KZ) .
      PLTOP = MIN( ZF( KZ ), PLTOP )
      PLBOT = MIN( ZF( KZ )-1.0, PLBOT )
```

The matching Python from `plsprd`
([plume.py](../cmaq_port/plume.py)):

```python
    index = 0
    while index < len(lapse) - 1 and center > heights[index + 1]:
        index += 1
    gradient = real32(lapse[index] + real32(-0.0098))
    exponent = real32(real32(-117.0) * gradient)
    try:
        exponential = real32(math.exp(exponent))
    except OverflowError as exc:
        raise ValueError("PLSPRD exponential overflow; check temperature gradients.") from exc
    sigma = max(real32(10.0), real32(real32(15.0) * exponential))
    depth = real32(real32(3.545) * sigma)
    half = real32(depth / real32(2.0))
    top = real32(center + half)
    bottom = max(real32(0.0), real32(center - half))
    top = min(heights[-1], top)
    bottom = min(real32(heights[-1] - real32(1.0)), bottom)
    return PlumeBounds(top, bottom)
```

Learning notes:

- CMAQ's original declares everything as default `REAL` (single precision),
  so *every* operation rounds to binary32. The Python port makes that explicit
  by wrapping each constant and intermediate in `real32(...)`. Contrast this
  with AERMOD's mostly double-precision kernels, where no such per-op rounding
  applies.
- `GAMMA = -0.0098` is the dry-adiabatic lapse rate (K/m). It is added to the
  selected potential-temperature gradient `DTHDZ(K)` to recover the actual
  ambient temperature gradient `DTDZ` before the exponential spread term.
- The vertical spread is exponential: `SIGZ0 = max(10, 15*exp(-117*DTDZ))`.
  The `10.0` is a hard 10 m floor on sigma, and the plume depth is
  `DPTH = 3.545 * SIGZ0` (`SZ0FAC = 3.545`, roughly ±half-depth about the
  center giving the familiar ~±1.77 sigma bounds).
- Top and bottom are clipped *separately* and in source order: bottom is first
  floored at `0` (ground), then top is capped at `ZF(KZ)` (model top) and
  bottom at `ZF(KZ)-1`. Order and the independent clips matter for reproducing
  edge cases.
- The layer-selection loop exits on `CEFSTK .LE. ZF(K)`, so a center sitting
  exactly on an interface selects the *lower* layer; the Python `while` uses
  the strict `center > heights[index + 1]` to preserve that interface-equality
  tie-break.

## Verified evidence and alternating workflow

`python -m cmaq_port.kernel_validation` copies and compiles the unchanged
original with an explicit assumed-shape interface, runtime checks, no fast
math, and the existing project-local GNU Fortran 16.2.0 compiler. This
standalone routine needs no MPI, netCDF, I/O API, or WSL.

[Validation report](../runs/20261024T173132-kernels-08d3afcf/kernels.json):
**44 synthetic cases, 88 values, zero observed absolute difference**.
The matrix covers four gradient profiles and eleven heights, including both
sides of interfaces and model-top clipping. Engineering thresholds are
`atol=2e-5 m`, `rtol=2e-6`, not complete-model scientific acceptance criteria.
Run directories retain source/compiler/driver hashes, source status/commit,
commands, logs, inputs, outputs, and per-case errors.

At this first-kernel milestone **33 CMAQ project unit tests passed**, including 11 plume tests.
These cases are synthetic, not states captured from a full CMAQ run.
The subsequent [DELTA_ZS lesson](stack-layer-geometry.md) adds geometry,
independent model review, and the current combined validation totals.

This round ports **CMAQ PLSPRD, then AERMOD FLUXES**, each with source notes,
lessons, and original-code checks. Continue one routine from each project
per round, choosing bounded callable science before integration. Full CMAQ
reference execution and real inter-model comparison remain separate gates.
The user is handling WSL installation in another chat; this workflow does
not install WSL or alter system configuration.
