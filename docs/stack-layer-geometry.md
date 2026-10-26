# DELTA_ZS: where are the layer interfaces relative to each stack?

This is the second CMAQ helper port. Python:
[delta_zs in plume.py](../cmaq_port/plume.py).
Original: [delta_zs.f](C:/Dev/CMAQ/CCTM/src/plrise/smoke/delta_zs.f),
pinned source commit `e10d0a04f319d6be6db086bc66ba13ce4044987b`.
Revision notes credit J. Young (November 2007) and D. Wong's assumed-shape
array update (August 2015). Retain the EPA permission notice in the unchanged
source and [upstream license](../LICENSE-CMAQ.md); see the
[plume-spread attribution](plume-spread.md#source-and-attribution).

## Why this is useful

Before interpreting plume height, learn the geometry surrounding a source.
The caller [PTMET.F](C:/Dev/CMAQ/CCTM/src/emis/emis/PTMET.F) fills local
meteorological columns and calls DELTA_ZS to derive:

- **ZSTK(L,S):** upper layer interface minus physical stack height, m.
- **DDZF(L,S):** reciprocal layer thickness, 1/m.

This routine does not calculate plume rise, emission fractions, or
concentrations. STKHT is a physical source height, not PLSPRD's effective
plume center. A negative ZSTK is meaningful: the interface is below the stack.
There is no clipping to ground or the model top.

## Indexing is the main lesson

The Python API is:

```python
from cmaq_port.plume import delta_zs

geometry = delta_zs(
    source_map=(2, 0),
    stack_heights=(10, 100, 300),
    layer_heights=((100, 120), (500, 600)),
)
assert geometry.relative_heights == ((-200, 110), (200, 590))
```

The Python map is **zero-based**, whereas Fortran SRC_MAP is one-based.
Each map entry selects from the stack list. It does **not** reorder the
meteorological columns: those columns are already in local-source order.
Maps may be permuted or repeated. In this example, local column zero uses
stack height 300 m but meteorological interfaces 100/500 m; local column one
uses stack height 10 m and its own interfaces 120/600 m.

Both inputs/outputs are **layer-major**: `[layer][local source]`. Unlike
PLSPRD, this ZF input contains upper interfaces **without a ground row**.
The number of layers/sources is inferred rather than passed separately.

## Source equations and precision

For each local source S, SRC=SRC_MAP(S):

```text
ZSTK(L,S) = ZF(L,S) - STKHT(SRC)
DDZF(1,S) = 1 / ZF(1,S)
DDZF(L,S) = 1 / (ZF(L,S) - ZF(L-1,S)), for L >= 2
```

The first thickness is measured from implicit ground zero. Later layers
use their own lower interface, not the stack height. Returned values retain
default REAL/binary32 rounding, including rounded input heights, each
subtraction, and each reciprocal.

Python rejects nonfinite/out-of-range values, negative stack heights,
noninteger/out-of-range maps, mismatched row sizes, nonpositive first
interfaces, and interfaces that cease increasing after binary32 rounding.
Overflow raises ValueError. These are API preconditions, not a port of all
CMAQ meteorological screening. At least one layer is required; zero local
sources returns an empty row per layer, matching the source's zero-iteration
loop. The production caller normally skips that call when no local sources
exist.

### Fortran and Python side by side

The original nested loop (`delta_zs.f`, lines 51-65):

```fortran
      DO S = 1, MY_NSRC
         SRC = SRC_MAP( S )

         ZF0 = ZF( 1,S )
         ZSTK( 1,S ) = ZF0 - STKHT( SRC )
         DDZF( 1,S ) = 1.0 / ZF0

         DO L = 2, EMLAYS
            ZF1 = ZF( L,S )
            ZSTK( L,S ) = ZF1 - STKHT( SRC )
            DDZF( L,S ) = 1.0 / ( ZF1 - ZF0 )
            ZF0 = ZF1
         END DO

      END DO
```

The matching Python loop (`plume.py`, lines 60-71):

```python
    relative: list[tuple[float, ...]] = []
    inverse: list[tuple[float, ...]] = []
    previous = (0.0,) * len(source_map)
    for row in heights:
        thicknesses = tuple(real32(upper - lower) for upper, lower in zip(row, previous, strict=True))
        if any(value <= 0 for value in thicknesses):
            raise ValueError("Upper interfaces must be positive and increase after REAL rounding.")
        relative.append(tuple(
            real32(height - stacks[index]) for height, index in zip(row, source_map, strict=True)
        ))
        inverse.append(tuple(real32(1.0 / value) for value in thicknesses))
        previous = row
```

**Learning notes**

- Fortran's `SRC_MAP( S )` is one-based and indexes into `STKHT`; the Python
  `source_map` entries are zero-based and index `stacks`, the stack list, not
  the local meteorological columns (`row`), which already arrive in
  local-source order.
- Fortran is source-major (outer `S`, inner `L`); Python is layer-major
  (outer `for row in heights`, i.e. per layer `L`). Both produce the same
  `[layer][local source]` outputs, so the loop nesting is transposed but the
  computed values match.
- The first thickness is special: Fortran sets `DDZF(1,S) = 1.0 / ZF( 1,S )`,
  measuring from implicit ground zero, while later layers use `ZF1 - ZF0`.
  Python reproduces this by seeding `previous` with zeros, so the first
  `upper - lower` subtraction is `ZF(1) - 0`, and subsequent rows subtract
  the actual lower interface.
- Every operation rounds to binary32: Fortran's default `REAL` arithmetic is
  mirrored by wrapping each subtraction and reciprocal in `real32(...)`,
  including the already-rounded input heights.

## Exercises

Run `python -m cmaq_port lesson`, section 9.

1. Explain the negative -200 m in the example.
2. Change the map to `(0,2)` without swapping meteorological columns.
3. Repeat a map entry: why can the resulting columns still differ?
4. Calculate the two thicknesses in the first local column: 100 and 400 m.
   Their reciprocals are approximately 0.01 and 0.0025 1/m.
5. Explain why including a ground row, as in PLSPRD, is invalid here.
6. Move the stack above every interface; note that ZSTK stays negative.

## Original-code evidence and review

The shared validator compiles unchanged PLSPRD and DELTA_ZS with an explicit
assumed-shape test interface, bounds/runtime checks, and no fast math.
DELTA_ZS's declaration extends beyond fixed-form column 72; the compiler
therefore uses `-ffixed-line-length-none`, recorded in the run report.
The original file was not reformatted. The initial failed compile/log is
retained rather than hidden.

[Two-helper milestone report](../runs/20261024T174254-kernels-9ed417a8/kernels.json):

| Routine | Synthetic cases | Scalar values | Maximum observed difference |
| --- | --- | --- | --- |
| PLSPRD | 44 | 88 | 0 |
| DELTA_ZS | 64 | 720 | 0 |
| Total | 108 | 808 | 0 |

The 64 geometry cases combine four grids, four maps, and four stack lists,
including nonidentity/repeated maps and source-dependent columns.
The adapter serializes arrays in Fortran's layer-index-first storage order;
Python still exposes layer-major tuples. DELTA_ZS tolerances are
`atol=1e-10` in each output's units and `rtol=2e-6`; PLSPRD keeps its existing
`atol=2e-5 m`, `rtol=2e-6`. Reports retain per-case criteria, inputs, hashes,
commands, results, and errors. Errors across unlike units are not a single
physical metric.

At this milestone **44 CMAQ unit tests passed**, including 11 new geometry tests. An independent
model review checked arithmetic, mapping, array serialization, and the
original-code comparison method; no blocking correctness issues were found.
These are synthetic tests, not captured full-model states. Full CMAQ execution,
realistic-state coverage, and concentration-level validation remain pending.
The subsequent [PLMRIS lesson](plume-rise.md) completes a connected
geometry/rise/spread check and records the current totals.

This round pairs DELTA_ZS with
[AERMOD HSPRIM](../../AERMOD-Python/docs/stack-tip-downwash.md).
Do not equate layer-relative geometry to an adjusted plume release height.
