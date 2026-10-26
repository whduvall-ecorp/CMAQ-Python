# POLY: the I/O API interpolation behind PREPLM

Python: `poly` in `cmaq_port/interpolation.py`.
Original: `REAL FUNCTION POLY (XPT, XPTS, YPTS, NDEG)` in the EDSS/Models-3
I/O API (`ioapi/poly.f`, version `$Id: poly.f 219 2015-08-17 ...$`).

## Why this routine

CMAQ's `CCTM/src/plrise/smoke/preplm.f` interpolates ambient temperature and
wind to stack top with `POLY`:

```fortran
TSTK =      POLY( HTS, ZH( M:EMLAYS ), TA( M:EMLAYS ), DEG )
WSTK = MAX( POLY( HTS, ZH( M:EMLAYS ), WSPD( M:EMLAYS ), DEG ), 0.1 )
```

`POLY` is defined nowhere in the CMAQ source tree; it is an I/O API library
function. Porting it is the prerequisite for a faithful `PREPLM` port, which is
the next meteorological-preparation kernel in the
[workflow plan](workflow-plan.md). The `POLY3`/`POLY3B` routines in
`CCTM/src/aero/aero6` are unrelated cubic-equation solvers.

## Algorithm

`POLY` evaluates the unique polynomial of degree `NDEG` through the `NDEG+1`
points `(XPTS(k), YPTS(k))` using Newton divided differences in the form

```text
P(x) = f[x1] + f[x1,x2](x-x1) + f[x1,x2,x3](x-x1)(x-x2) + ...
       + f[x1..x_{NDEG+1}] (x-x1)...(x-x_{NDEG})
```

Each divided difference is computed by the explicit sum

```text
f[x1,...,x_{i+1}] = sum_{k=1}^{i+1} YPTS(k) / prod_{j != k} (XPTS(k) - XPTS(j))
```

exactly as the source does: it unrolls the linear (degree-1) term, then the
`K = 1` denominator, and preserves the original operation order. `NDEG >= 1`
is required (the source indexes `XPTS(2)`), the `NDEG+1` abscissae must be
distinct (the denominators are differences of abscissae), and `NDEG <= 5` is
recommended upstream because high-degree polynomial interpolation is unstable.

## Port fidelity

Every elementary REAL operation rounds to binary32 via `real32`, matching the
single-precision Fortran result. Inputs are rounded to REAL on entry. The port
adds explicit guards the library omits: `NDEG >= 1`, matching point counts,
and distinct abscissae raise `ValueError` instead of indexing out of bounds or
dividing by zero. It does not read meteorology, run `PREPLM`, or apply the
`MAX(..., 0.1)` wind floor — that floor belongs to the `PREPLM` caller.

## Validation

`python -m cmaq_port.poly_validation` compiles the vendored, unchanged
`reference/ioapi/poly.f` standalone with the project-local GNU Fortran
toolchain and compares it to the Python port on a deterministic matrix of 107
cases: PREPLM-like temperature/wind profiles at degrees 1-5, exact
linear/quadratic/cubic polynomials, and non-monotonic/descending abscissae.

- **Maximum absolute difference: 0.0 (bitwise agreement) across all 107 cases.**
- Unit tests in `tests/test_interpolation.py` independently prove data-point
  reproduction, polynomial exactness within float32, extrapolation, and the
  degree/count/distinctness guards.

Provenance and license for the vendored reference are in
[reference/ioapi/PROVENANCE.md](../reference/ioapi/PROVENANCE.md). This is
kernel-level evidence on synthetic inputs; it is not a `PREPLM` port, a CMAQ
run, or any model-equivalence claim. The runnable example is section 12 of
`python -m cmaq_port lesson`.
