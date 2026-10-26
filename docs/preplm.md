# PREPLM: stack-top meteorology for plume rise (Fortran and Python)

Python: `preplm` in [`cmaq_port/preplm.py`](../cmaq_port/preplm.py).
Original: `SUBROUTINE PREPLM` in `CCTM/src/plrise/smoke/preplm.f`.

`PREPLM` turns per-source meteorology into the exact inputs `PLMRIS` needs:
wind speed per layer, the potential-temperature gradient `DTHDZ`, the stack and
mixing-layer indices `LSTK`/`LPBL`, and the temperature and wind interpolated to
stack top (`TSTK`, `WSTK`). It is the meteorological preparation step — not a
plume rise, concentration, or CMAQ run.

CMAQ computes in single precision (`REAL`), so the port rounds **every**
elementary operation to binary32 with `real32`, exactly as the other CMAQ
kernels do. It also reuses the already-ported I/O API [`POLY`](interpolation.md)
for the stack-top interpolation.

---

## Wind speed and virtual temperature

**Fortran:**
```fortran
      DO L = 1, EMLAYS
         P = UW( L )
         Q = VW( L )
         WSPD( L ) = SQRT( P * P + Q * Q )
         TV( L ) = TA( L ) * ( 1.0 + 0.622 * ( QV( L ) / ( 1.0 + QV( L ) ) ) )
      END DO
```

**Python:**
```python
for layer in range(emlays):
    p, q = uw_r[layer], vw_r[layer]
    wspd.append(real32(math.sqrt(_add(_mul(p, p), _mul(q, q)))))
    ratio = _div(qv_r[layer], _add(real32(1.0), qv_r[layer]))
    tv.append(_mul(ta_r[layer], _add(real32(1.0), _mul(real32(0.622), ratio))))
```

Learning note: `_add`, `_mul`, `_div` each round to binary32. A plain
`math.sqrt(p*p + q*q)` in float64 would drift from the Fortran `REAL` result by
a few ulps across many layers; the per-operation rounding keeps the port bit for
bit identical.

---

## Potential-temperature gradient

**Fortran:**
```fortran
      THV1  = TF( 1 ) * ( 1000.0 / PRES( 1 ) ) ** 0.286
      DO L = 2, EMLAYS
         IF ( HMIX .GT. ZF( L-1 ) ) LPBL = L
         IF ( HTS  .GT. ZF( L-1 ) ) LSTK = L
         THVK = TF( L ) * ( 1000.0 / PRES( L ) ) ** 0.286
         DTHDZ( L ) = DDZF( L ) * ( THVK - THV1 )
         THV1 = THVK
      END DO
C Set the 1st level gradient to the 2nd layer value
      DTHDZ( 1 ) = DTHDZ( 2 )
```

**Python:**
```python
thv1 = _mul(tf[0], _pow(_div(real32(1000.0), pres_r[1]), real32(0.286)))
for fortran_l in range(2, emlays + 1):
    if hmix_r > zf_r[fortran_l - 2]:
        lpbl = fortran_l
    if hts_r > zf_r[fortran_l - 2]:
        lstk = fortran_l
    thvk = _mul(tf[fortran_l - 1], _pow(_div(real32(1000.0), pres_r[fortran_l]), real32(0.286)))
    dthdz[fortran_l - 1] = _mul(ddzf_r[fortran_l - 1], _sub(thvk, thv1))
    thv1 = thvk
dthdz[0] = dthdz[1]
```

Learning notes:
- `LSTK`/`LPBL` are **1-based** layer indices; the port keeps them 1-based.
- `PRES` is declared `PRES(0:)` in Fortran. The active code uses `PRES(1)` and
  `PRES(L)`, so the Python `pres` sequence mirrors that with `pres[0]` present
  but unused by this path.
- The first-level gradient is deliberately overwritten with the second layer's
  value (`DTHDZ(1) = DTHDZ(2)`).

---

## Stack-top interpolation with POLY

**Fortran:**
```fortran
      IF ( .NOT. FIREFLG ) THEN
         M    = MAX( 1, LSTK - DEG - 1 )
         TSTK =      POLY( HTS, ZH( M:EMLAYS ), TA( M:EMLAYS ), DEG )
         WSTK = MAX( POLY( HTS, ZH( M:EMLAYS ), WSPD( M:EMLAYS ), DEG ), 0.1 )
      ELSE
         TSTK = TS
         WSTK = WSPD( 1 )
      END IF
```

**Python:**
```python
if not fireflg:
    m = max(1, lstk - DEG - 1)
    start = m - 1
    window_zh = zh_r[start:start + DEG + 1]
    tstk = poly(hts_r, window_zh, ta_r[start:start + DEG + 1], DEG)
    wstk = real32(max(poly(hts_r, window_zh, tuple(wspd[start:start + DEG + 1]), DEG), real32(0.1)))
else:
    tstk = ts_r
    wstk = wspd[0]
```

Learning notes:
- `DEG = 3`, so `POLY` uses four points. The Fortran passes the whole slice
  `ZH(M:EMLAYS)` but `POLY` reads only the first `DEG+1` entries; the port
  slices exactly those four points.
- The `MAX( ..., 0.1 )` wind floor lives in `PREPLM`, not in `POLY`.
- The fire branch bypasses interpolation entirely, so `PREPLM` needs at least
  `DEG+1 = 4` layers only when `FIREFLG` is false.

---

## A dead-code teaching point

The source computes `ES`, `QSFC`, `TVSFC`, and `THETG` (a surface virtual
potential temperature) for this line:

```fortran
!     DTHDZ( 1 ) = ( THV1 - THETG ) / ZF( 1 )
```

That line is **commented out**, and `DTHDZ(1)` is overwritten with `DTHDZ(2)`
anyway, so the whole `ES -> QSFC -> TVSFC -> THETG` chain has no effect on any
output. The port omits it and documents why, rather than copying dead
arithmetic. `PSFC` therefore survives in the signature (interface fidelity) but
does not influence results.

---

## Validation

`python -m cmaq_port.preplm_validation` compiles the **unchanged** `preplm.f`
from the local CMAQ checkout with a test-only `UTILIO_DEFN` adapter — which
exposes the unchanged I/O API `POLY` and a trivial `M3MSG2`, with no stubbed
physics — and compares it against the Python port:

- **72 cases / 1,008 values, maximum absolute difference 0.0 (bitwise).**
- `LSTK`/`LPBL` match exactly; `source_unchanged` is verified by hashing.

Unit tests are in [`tests/test_preplm.py`](../tests/test_preplm.py). PREPLM
completes the science routines in `plrise/smoke` (PLSPRD, DELTA_ZS, PLMRIS,
FIRE_PLMRIS, PREPLM). This is kernel-level evidence on synthetic inputs, not a
CMAQ run or model-equivalence claim.
