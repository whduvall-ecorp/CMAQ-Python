# FIRE_PLMRIS: trace a fire plume without borrowing stack physics

Python implementation: [`cmaq_port/fire_rise.py`](../cmaq_port/fire_rise.py).
Unchanged reference:
[`fire_plmris.F`](C:/Dev/CMAQ/CCTM/src/plrise/smoke/fire_plmris.F), from
source commit `e10d0a04f319d6be6db086bc66ba13ce4044987b`.

The source credits a December 1995 prototype by CJC, based on RADM 2.6
`PLUMER`, and a March 1999 copy from `plmris.F` by M. Houyoux. Its header
retains the U.S. EPA permission notice and the SMOKE 2002 MCNC Environmental
Modeling Center attribution, all rights reserved, with conditions in the
SMOKE `COPYRIGHT` file. That referenced file is not present in this local
checkout. The original source header and repository-root license are retained
in the validation run, and the missing terms are not inferred from the root
MIT license. See the companion [CMAQ attribution](../LICENSE-CMAQ.md), and
obtain and review the applicable full source terms before redistributing
derivatives.

## Actual caller and distinct inputs

The fire/stack split is in
[`PT3D_DEFN.F`](C:/Dev/CMAQ/CCTM/src/emis/emis/PT3D_DEFN.F). Before calling
this routine, the caller runs PREPLM, applies its `USTARMIN`, and converts
surface sensible heat from W/m² to `HFX` in m K/s. For a fire it passes the
already computed `BFLUX` array value as `BFLX`; FIRE_PLMRIS itself does not
derive buoyancy flux from stack diameter, exit velocity, or temperature.
The routine has no stack height or stack diameter argument.

| Python input | Original argument | Meaning |
| --- | --- | --- |
| `fire_layer` | `LSTK` | Zero-based Python layer containing the fire top; Fortran is one-based |
| `heat_flux` | `HFX` | Caller-converted sensible heat, m K/s |
| `mixing_height` | `HMIX` | Mixing height in m; source tests `HMIX <= 0`, not `HMIX <= fire height` |
| `buoyancy_flux` | `BFLX` | Caller-supplied m⁴/s³; never calculated or given a positive floor here |
| `fire_temperature` | `TSTK` | Temperature at top of fire, K |
| `friction_velocity` | `USTAR` | Caller-supplied value after its own minimum policy |
| `gradients` | `DTHDZ` | Layer potential-temperature gradients, K/m |
| `temperatures` | `TA` | Layer temperatures, K |
| `winds` | `WSPD` | Layer wind speeds, m/s |
| `interfaces` | `ZF(0:)` | Ground-inclusive layer interfaces, m |
| `fire_wind` | `WSTK` | Fire-top wind input and reporting output, m/s |

Every value is rounded to default Fortran REAL/binary32 precision. The Python
API requires consistent array lengths, a valid zero-based fire layer,
`ZF[0] == 0`, strictly increasing interfaces after REAL rounding, positive
`TSTK`, layer `TA`, and `USTAR`, and nonnegative input winds. It does not
impose a caller-specific `USTARMIN`. These explicit physical preconditions
are not extra CMAQ preprocessing.

`BFLX <= 0` is deliberately not replaced by a stack-style positive minimum.
With zero flux the original equations produce zero rise. Negative flux can
produce negative neutral rise if no fractional power is evaluated; a
negative base passed to a fractional power has no real-valued result, so the
Python port raises an explicit `ValueError` at that operation instead of
silently clipping it. The caller owns interpretation and screening of the
fire-specific flux.

## Follow the source

The first calculations use the surface layer temperature and caller inputs:

```text
HSTAR = GRAV * HFX / TA(1)
LPLM  = LSTK
WSTK  = max(WSTK, 1 m/s)
DHM   = 0
```

`GRAV=9.80622`; `HCRIT=REAL(1e-4)*REAL(0.03)` and `SMALL=3e-5`.
The fire routine's statement functions are not the stack PLMRIS functions:

```text
NEUTRL(B,U,US) = 1.3327804 * B / (U * US * US)
STABLE(B,U,S)  = 2.6 * (B / (U * S)) ** (1/3)
UNSTBL(B,U)    = 30 * (B / U) ** 0.6
```

There is no height-dependent cap in `NEUTRL`. All operations, constants,
and exponents are rounded as default REAL; Python `math.pow` is rounded back
to REAL and may differ by a few ulps from the compiler's single-precision
power routine.

### Initial stability choice

- If `HSTAR > HCRIT`, use the unstable branch. For `HMIX > 0`, compare
  unstable and neutral rise. For `HMIX <= 0`, compare stable-gradient and
  neutral rise. `HMIX` is not reduced by a fire/stack height.
- Otherwise, if `HSTAR < -HCRIT` **or** `DTHDZ(LSTK) > 0.001`, use stable
  treatment.
- Otherwise, use neutral treatment. Both HCRIT comparisons and the gradient
  test are strict.

The source initializes `DHM=0`; its later momentum override tests
`DHM > DH AND WSTK > 1`. For positive buoyancy rise this cannot select a
momentum path. The source has a wind reset guarded by `LPLM != LSTK`, but
`LPLM` was just assigned `LSTK`, so that initial reset is unreachable.
Neither behavior is “fixed” by the port.

### Layer-by-layer residual buoyancy

`ZPLM` initially contains plume-top rise. Compare it against `ZF(LPLM)`,
the upper interface of the current fire layer. If it crosses that interface,
calculate residual `RBFLX` using the active stability code:

```text
IQ=1: RX = RX/45
     RBFLX = WPLM * RX**(5/3)
IQ=2: PX = (2/3)*ZPLM
     RBFLX = (1/2.664) * RX**(5/3) * WPLM * USTAR**2 / PX**(2/3)
IQ=3: RBFLX = (1/59.319) * WPLM * SX * RX**3
```

The first residual calculation uses the fire wind and `TSTK`; later passes
use the wind/temperature for the current layer. Advance one layer, replace a
zero `WSPD` with 1 m/s, recalculate `SX`, and select stable or neutral/unstable
rise from that layer's gradient. The loop only advances upward.

At the final layer, if the estimated plume top crosses `ZF(EMLAYS)`, the
routine clips **rise** to that interface and writes this source diagnostic:

```text
Plume rise reached EMLAYS with ZPLM: <layer><height to 0.001 m>
```

`ZPLM` is then always converted to the fire centerline with `2/3*ZPLM`.
Unlike stack PLMRIS, no stack height is added, and there is no hot-first-layer
centerline adjustment. The returned `top_log_height` records the pre-conversion
clipped value; the validation checks the original formatted `LOGDEV` line
exactly.

## Runnable numerical lesson

Run `python -m cmaq_port lesson`, section 11. It supplies a synthetic
`BFLX=0.4 m⁴/s³`, `HFX=0.03 m K/s`, `TSTK=450 K`, `USTAR=0.2 m/s`, and
`WSTK=0.2 m/s` on a five-layer grid. The example computes an approximately
**13.3278 m** centerline and reports the fire-wind floor of **1 m/s**.
These are illustrative inputs, not a CMAQ run, plume emission fractions, or
concentration results.

## Fortran and Python side by side

The two excerpts below trace the same logic: FIRE_PLMRIS's distinct stability
selection, the residual-buoyancy loop over layers, and the final
`2/3*ZPLM` centerline conversion that adds **no** stack height.

### Stability selection and residual loop (Fortran)

The following is condensed pseudocode drawn from
`fire_plmris.F` (`SUBROUTINE FIRE_PLMRIS`, stability selection at source
lines 155–239 and the residual-buoyancy loop at source lines 245–301).

```text
Pseudocode (from fire_plmris.F)

HSTAR = GRAV * HFX / TA(1)          ! surface temperature is correct
LPLM  = LSTK
WSTK  = MAX( WSTK, 1.0 )            ! fire-wind floor 1 m/s
DHM   = 0.0

IF ( HSTAR > HCRIT ) THEN                       ! buoyant (hot) case
   IF ( HMIX <= 0.0 ) THEN                       ! fire above mixing height
      SX  = MAX( GRAV * DTHDZ(LPLM) / TSTK, SMALL )
      IF ( LPLM /= LSTK ) WSTK = WSPD(LPLM)      ! unreachable on first pass
      DHSM = 0.0  (if DTHDZ(LPLM) > 0.001)  else  DHM
      DHM  = MIN( DHSM, DHM )
      DHN  = NEUTRL( BFLX, WSTK, USTAR )
      DH   = STABLE( BFLX, WSTK, SX )
      IQ   = 2 (if DHN < DH) else 3;  DH = MIN( DHN, DH )
      IF ( DHM > DH .AND. WSTK > 1.0 ) THEN DH = DHM; IQ = 4
      DHT  = 1.5 * DH
   ELSE                                           ! unstable (within mixed layer)
      DH   = UNSTBL( BFLX, WSTK )
      DHN  = NEUTRL( BFLX, WSTK, USTAR )
      IQ   = 2 (if DHN < DH) else 1;  DH = MIN( DHN, DH )
      IF ( DHM > DH .AND. WSTK > 1.0 ) THEN DH = DHM; IQ = 4
      DHT  = 1.5 * DH
   END IF
ELSE IF ( HSTAR < -HCRIT .OR. DTHDZ(LSTK) > 0.001 ) THEN   ! stable case
   SX  = MAX( GRAV * DTHDZ(LSTK) / TSTK, SMALL )
   DHT = 1.5 * STABLE( BFLX, WSTK, SX )
   DHN = 1.5 * NEUTRL( BFLX, WSTK, USTAR )
   IF ( DHN < DHT ) THEN DHT = DHN; IQ = 2  ELSE  IQ = 3
ELSE                                             ! neutral case
   DHT = 1.5 * NEUTRL( BFLX, WSTK, USTAR )
   IQ  = 2
END IF

NN = 0;  RBFLX = BFLX;  ZPLM = DHT
IF ( IQ == 4 ) GO TO 199        ! momentum rise: skip the buoyancy loop

DO
   RX = ZPLM - ZF(LPLM)
   IF ( RX <= 0.0 ) EXIT
   IF ( LPLM == EMLAYS ) THEN   ! clip at model top and log, then exit
      ZPLM = MIN( ZPLM, ZF(EMLAYS) );  WRITE(LOGDEV,...);  EXIT
   END IF
   NN = NN + 1
   WPLM = (NN > 1) ? WSPD(LPLM) : WSTK           ! first pass uses fire wind
   IF ( IQ == 1 ) RBFLX = WPLM * (D45*RX)**FIVETHD
   IF ( IQ == 2 ) RBFLX = D2664 * RX**FIVETHD * WPLM * USTAR**2 / (TWOTHD*ZPLM)**TWOTHD
   IF ( IQ == 3 ) RBFLX = D59319 * WPLM * SX * RX**3
   LPLM = LPLM + 1
   WPLM = WSPD(LPLM);  TPLM = TA(LPLM);  IF ( WPLM == 0.0 ) WPLM = 1.0
   SX   = GRAV * DTHDZ(LPLM) / TPLM
   IF ( SX > SMALL ) THEN                          ! stable above
      DHT = 1.5*STABLE(RBFLX,WPLM,SX);  DHN = 1.5*NEUTRL(RBFLX,WPLM,USTAR)
      IF ( DHN < DHT ) THEN DHT = DHN; IQ = 2  ELSE  IQ = 3
   ELSE                                            ! neutral/unstable above
      DHN = NEUTRL(RBFLX,WPLM,USTAR);  DH = UNSTBL(RBFLX,WPLM);  IQ = 1
      IF ( DHN < DH ) THEN DH = DHN; IQ = 2
      DHT = 1.5 * DH
   END IF
   ZPLM = ZF(LPLM-1) + DHT
END DO
```

The final conversion is short enough to quote verbatim from the routine
(source lines 308–312):

```fortran
199   CONTINUE

C Compute plume rise amount (DH) and actual final plume centerline height (ZPLM):
      DH  = TWOTHD * ZPLM
      ZPLM = DH
```

### Stability selection and residual loop (Python)

The port in `cmaq_port/fire_rise.py` (function `fire_plmris`) mirrors each
branch. The stability selection:

```python
hstar = _div(_mul(GRAV, hfx), ta[0])
layer = fire_layer
wstk = max(wstk, real32(1.0))
dhm = real32(0.0)
sx = real32(0.0)

if hstar > HCRIT:
    if hmix <= 0:
        sx = max(_div(_mul(GRAV, lapse[layer]), tstk), SMALL)
        dhsm = real32(0.0) if lapse[layer] > real32(0.001) else dhm
        dhm = min(dhsm, dhm)
        dhn = _neutral(bflx, wstk, ustar)
        dh = _stable(bflx, wstk, sx)
        iq = 2 if dhn < dh else 3
        dh = min(dhn, dh)
        if dhm > dh and wstk > 1.0:
            dh, iq = dhm, 4
        dht = _mul(real32(1.5), dh)
    else:
        dh = _unstable(bflx, wstk)
        dhn = _neutral(bflx, wstk, ustar)
        iq = 2 if dhn < dh else 1
        dh = min(dhn, dh)
        if dhm > dh and wstk > 1.0:
            dh, iq = dhm, 4
        dht = _mul(real32(1.5), dh)
elif hstar < -HCRIT or lapse[layer] > real32(0.001):
    sx = max(_div(_mul(GRAV, lapse[layer]), tstk), SMALL)
    dht = _mul(real32(1.5), _stable(bflx, wstk, sx))
    dhn = _mul(real32(1.5), _neutral(bflx, wstk, ustar))
    if dhn < dht:
        dht, iq = dhn, 2
    else:
        iq = 3
else:
    dht = _mul(real32(1.5), _neutral(bflx, wstk, ustar))
    iq = 2
```

The residual-buoyancy loop over layers:

```python
buoyancy = bflx
zplm = dht
transitions = 0
top_log_height: float | None = None
path = [iq]
if iq != 4:
    while True:
        residual_height = _sub(zplm, zf[layer + 1])
        if residual_height <= 0:
            break
        if layer == layers - 1:
            zplm = min(zplm, zf[layers])
            top_log_height = zplm
            break

        transitions += 1
        if transitions > 1:
            plume_wind = wspd[layer]
        else:
            plume_wind = wstk

        if iq == 1:
            residual_height = _mul(D45, residual_height)
            buoyancy = _mul(plume_wind, _pow(residual_height, FIVETHD))
        elif iq == 2:
            px = _mul(TWOTHD, zplm)
            buoyancy = _div(
                _mul(_mul(_mul(D2664, _pow(residual_height, FIVETHD)), plume_wind),
                     _mul(ustar, ustar)),
                _pow(px, TWOTHD),
            )
        else:
            buoyancy = _mul(
                _mul(_mul(D59319, plume_wind), sx),
                _mul(_mul(residual_height, residual_height), residual_height),
            )

        layer += 1
        plume_wind = wspd[layer]
        plume_temperature = ta[layer]
        if plume_wind == 0:
            plume_wind = real32(1.0)
        sx = _div(_mul(GRAV, lapse[layer]), plume_temperature)
        if sx > SMALL:
            dht = _mul(real32(1.5), _stable(buoyancy, plume_wind, sx))
            dhn = _mul(real32(1.5), _neutral(buoyancy, plume_wind, ustar))
            if dhn < dht:
                dht, iq = dhn, 2
            else:
                iq = 3
            dh = _div(dht, real32(1.5))
        else:
            dhn = _neutral(buoyancy, plume_wind, ustar)
            dh = _unstable(buoyancy, plume_wind)
            iq = 1
            if dhn < dh:
                dh, iq = dhn, 2
            dht = _mul(real32(1.5), dh)

        path.append(iq)
        zplm = _add(zf[layer], dht)
```

The final two-thirds conversion, with no stack height added:

```python
zplm = _mul(TWOTHD, zplm)
return FirePlumeRise(
    zplm, wstk, top_log_height is not None, top_log_height, tuple(path), transitions,
)
```

### Learning notes

- The caller supplies `BFLX` directly (in m⁴/s³); FIRE_PLMRIS never derives
  buoyancy from stack diameter, exit velocity, or temperature as stack PLMRIS
  does — there are no stack dimensions in the signature.
- The fire wind is floored at 1 m/s (`WSTK = MAX(WSTK, 1.0)` /
  `wstk = max(wstk, real32(1.0))`) before any rise is computed.
- The final centerline is `2/3 * ZPLM` with **no** stack height added and no
  hot-first-layer adjustment, unlike PLMRIS which adds back the stack height.
- `HMIX <= 0` selects the above-mixing-height branch (stable/neutral momentum
  logic), while `HMIX > 0` under a buoyant `HSTAR` selects the unstable branch.
- All arithmetic is single-precision binary32: the Fortran `REAL` values are
  reproduced in Python via `real32`, which is why the fire comparison tolerates
  differences at the binary32 rounding scale.

## Differential evidence and limits

Latest compiled report:
[`kernels.json`](../runs/20261024T225921-kernels-e99396a9/kernels.json).

| Operation | Cases | Values | Max absolute difference |
| --- | ---: | ---: | ---: |
| PLSPRD | 44 | 88 | 0 |
| DELTA_ZS | 64 | 720 | 0 |
| PLMRIS | 488 | 1,464 | 0 |
| Connected DELTA_ZS → PLMRIS → PLSPRD | 488 | 6,328 | 0 |
| FIRE_PLMRIS | 224 | 672 | 0.000030517578125 |
| **Total** | **1,308** | **9,272** | **0.000030517578125** |

All 136 original top-reaching `LOGDEV` lines across PLMRIS and FIRE_PLMRIS
were present; the eight FIRE_PLMRIS records matched the exact Fortran-formatted
line. The fire matrix spans stable/neutral/unstable gradients and heat fluxes,
positive and zero BFLX, fire-layer indices, zero-wind floors, mixing-height
branches, top clipping, and exact/adjacent HCRIT and 0.001-gradient boundaries.
The fire comparison uses `atol=2e-4` (m or m/s) and `rtol=2e-6`, with the
top flag and diagnostic line checked exactly. The maximum observed difference
was 3.05×10⁻⁵ in the tested values; this is below the predeclared tolerance,
not evidence of cross-compiler bitwise identity.

Original source, `CONST.EXT`, and test-only logging adapter are copied into
the persistent run directory with hashes and compile metadata. The standalone
driver supplies only the `RUNTIME_VARS.LOGDEV` needed by the unchanged
routine. No production runtime, PREPLM, fire-BFLX calculation, layer
emissions allocation, full CMAQ build/run, or actual CMAQ/AERMOD comparison
is validated here. The local source checkout remains unchanged. Independent
review of this new port is still required before treating the implementation
as reviewed.
