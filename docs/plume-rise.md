# PLMRIS: follow an elevated plume through meteorological layers

This completes the interrupted third CMAQ port. Python:
[rise.py](../cmaq_port/rise.py). Original:
[plmris.F](C:/Dev/CMAQ/CCTM/src/plrise/smoke/plmris.F), from pinned commit
`e10d0a04f319d6be6db086bc66ba13ce4044987b`.

The source describes a Briggs algorithm adapted from RADM 2.6, with a
December 1995 prototype by CJC and a March 1999 SMOKE copy by M. Houyoux.
Its notice states **COPYRIGHT (C) 2002, MCNC Environmental Modeling Center.
All Rights Reserved**, referring to SMOKE COPYRIGHT conditions. Retain
that attribution, the EPA permission notice in the original, and
[LICENSE-CMAQ.md](../LICENSE-CMAQ.md). Validation runs preserve the complete
unchanged source and license.

## Where it fits

The actual non-fire caller is
[PT3D_DEFN.F](C:/Dev/CMAQ/CCTM/src/emis/emis/PT3D_DEFN.F):

```text
Meteorology + physical stack parameters
  PREPLM -> layer index, temperature/wind at stack, temperature gradients
  DELTA_ZS -> interface-minus-stack heights
  PLMRIS -> final effective plume centerline height
  default IPVERT=0: Turner thickness
  alternate IPVERT!=0: PLSPRD bounds
  pressure-based layer allocation and caller error checks
```

PREPLM, fire plume rise, caller emission fractions, and concentrations are
not ported here. The test chain uses the alternate PLSPRD branch explicitly;
it does not silently replace the default Turner method.

## Inputs: distinguish geometry, temperatures, and winds

| Input | Meaning |
| --- | --- |
| `stack_layer` | Zero-based LSTK; equality with an upper interface selects the lower layer |
| `heat_flux` | HFX already converted to m K/s, **not W/m2** |
| `mixing_height` | HMIX, m |
| `stack_diameter`, `stack_height` | Physical stack dimensions, m |
| `exhaust_temperature`, `exhaust_velocity` | STKTK in K and STKVE in m/s |
| `stack_temperature`, `stack_wind` | Ambient TSTK/WSTK at stack height, normally from PREPLM |
| `friction_velocity` | USTAR, m/s; caller supplies its own minimum |
| `gradients` | DTHDZ potential-temperature gradients, K/m |
| `temperatures`, `winds` | Layer temperatures in K and wind speeds in m/s |
| `interfaces` | Ground-inclusive ZF, n+1 entries |
| `centers` | ZH, n entries strictly inside their layers |
| `relative_heights` | Optional ZSTK, n entries; otherwise derived by DELTA_ZS |

Every numerical input is rounded to the original default REAL/binary32.
The optional ZSTK must exactly equal ZF minus stack height after that rounding.
The API requires at least two layers, consistent shapes/index/geometry,
positive temperatures/friction velocity, and nonnegative physical stack
dimensions/exit velocity. Inputs must be finite and representable; invalid
values fail explicitly. The source's WSTK floor is preserved, not replaced
by a generic positive-wind guard. These are library preconditions, not a
replacement for complete production meteorological screening.

## First calculate buoyancy and momentum

The original [CONST.EXT](C:/Dev/CMAQ/CCTM/src/ICL/fixed/const/CONST.EXT)
defines GRAV=9.80622 as default REAL. Do not substitute AERMOD's gravity.

```text
HSTAR = GRAV * HFX / TA(1)
BFLX = 0.25 * GRAV * (STKTK-TSTK) * STKVE * STKDM * STKDM / STKTK
WSTK = max(WSTK, 1)
DHM = 3 * STKDM * STKVE / WSTK
```

**BFLX<=0 returns immediately:** centerline=STKHT+max(DHM,2). This is not
the later plume-top conversion, and it skips the hot-first-layer adjustment.
Do not reuse AERMOD FLUXES's positive buoyancy floor here; CMAQ retains
nonpositive buoyancy to select this branch.

## Initial rise and the four IQ codes

The source's statement functions are:

```text
NEUTRL(H,B,U,US) =
  min(10*H, 1.2*(B/(U*US*US))^0.6 * (H+1.3*B/(U*US*US))^0.4)
STABLE(B,U,S) = 2.6*(B/(U*S))^(1/3)
UNSTBL(B,U) = 30*(B/U)^0.6
```

Exponent constants and arithmetic order follow default REAL rounding.
HCRIT is the rounded product `1e-4*0.03`, and SMALL=3e-5.

- HSTAR>HCRIT selects unstable treatment. Below HMIX, take the lower
  neutral/unstable buoyant rise, then allow momentum to dominate. At or
  above HMIX, use the stable-gradient treatment and its stable momentum cap.
- HSTAR<-HCRIT or DTHDZ(LSTK)>0.001 selects the stable branch.
- Otherwise use neutral treatment. Equal threshold values do not satisfy
  the corresponding strict inequality.

IQ identifies the selected calculation: **1 unstable, 2 neutral, 3 stable,
4 momentum**. It can change as the plume enters new layers. Momentum IQ=4
skips the iterative buoyancy loop; it is not a fourth atmospheric stability.

## The layer loop

For each current layer, compare plume-top rise ZPLM to its upper-interface
height above the stack, ZSTK. If the plume does not cross the interface,
stop. If it crosses the last interface, limit **rise** to ZSTK(last), write
the original top-reaching diagnostic, and stop.

Otherwise compute residual buoyancy from the remaining rise using the
current IQ expression, advance one layer, and use that layer's wind,
temperature, and gradient for the next rise. Upper-layer transport wind
is floored at 1 m/s. Reconstruct total plume-top rise from the preceding
interface plus the new rise increment. The layer index only moves upward;
the loop is bounded by the number of input layers.

WSTK returned for reporting retains its stack-wind floor; it is not
overwritten on every upper-layer wind update. The source has an initial
LPLM!=LSTK check that cannot be true immediately after LPLM=LSTK.

### Fortran and Python side by side

The original is a bare `DO ... END DO` with `EXIT`s; the port is a `while True`
loop with `break`s. The structure, the IQ-dependent residual-buoyancy
expressions, the 1 m/s wind floor, and the `ZPLM = ZSTK(LPLM-1) + DHT`
reconstruction line up one-to-one.

**Fortran (`plmris.F`):**
```fortran
      DO       ! infinite loop computing further plume rise
         RX = ZPLM - ZSTK( LPLM )
         IF ( RX .LE. 0.0 ) THEN
            EXIT
         END IF
         IF ( LPLM .EQ. EMLAYS ) THEN   ! we're finished
            ZPLM = MIN( ZPLM, ZSTK( EMLAYS ) )
            WRITE( LOGDEV,'(5X, A, I3, F10.3)' )
     &             'Plume rise reached EMLAYS with ZPLM:', EMLAYS, ZPLM
            EXIT
         END IF
         NN = NN + 1
         IF ( NN .GT. 1 ) THEN
            WPLM = WSPD( LPLM )
         ELSE
            WPLM = WSTK
         END IF
         IF ( IQ .EQ. 1 ) THEN
            RX = D45 * RX
            RBFLX = WPLM * ( RX ** FIVETHD )
         ELSE IF ( IQ .EQ. 2 ) THEN
            PX = STKHT + TWOTHD * ZPLM
            RBFLX = D2664 * ( RX ** FIVETHD ) * WPLM * ( USTAR * USTAR ) / PX ** TWOTHD
         ELSE
            RBFLX = D59319 * WPLM * SX * RX ** 3
         END IF
         IF ( LPLM .LT. EMLAYS ) LPLM = LPLM + 1
         WPLM = MAX( WSPD( LPLM ), 1.0 )
         SX = GRAV * DTHDZ( LPLM ) / TPLM
         ... stable/neutral choice of DHT and IQ ...
         ZPLM = ZSTK( LPLM-1 ) + DHT
      END DO
```

**Python (`cmaq_port/rise.py`):**
```python
while True:
    rx = _sub(zplm, zstk[layer])
    if rx <= 0:
        break
    if layer == n - 1:
        zplm, top = min(zplm, zstk[-1]), True
        break
    nn += 1
    wplm = wind[layer] if nn > 1 else ws
    if iq == 1:
        rx = _mul(D45, rx)
        rbflx = _mul(wplm, _pow(rx, FIVETHD))
    elif iq == 2:
        px = _add(ht, _mul(TWOTHD, zplm))
        rbflx = _div(_mul(_mul(_mul(D2664, _pow(rx, FIVETHD)), wplm),
                          _mul(us, us)), _pow(px, TWOTHD))
    else:
        rbflx = _mul(_mul(_mul(D59319, wplm), sx), _mul(_mul(rx, rx), rx))
    layer += 1
    wplm = max(wind[layer], 1.0)
    sx = _div(_mul(GRAV, lapse[layer]), ta[layer])
    # ... stable/neutral choice of dht and iq ...
    zplm = _add(zstk[layer - 1], dht)
```

Learning notes:
- `D45`, `D2664`, `D59319`, `FIVETHD`, `TWOTHD` are the pre-rounded constants
  `1/45`, `1/2.664`, `1/59.319`, `5/3`, `2/3` from the source; each Python
  `_mul`/`_div`/`_pow` rounds to binary32 so the residual-buoyancy flux matches
  bit for bit.
- The Fortran `GO TO 199` that skips the whole loop for momentum rise (`IQ==4`)
  becomes the Python `if iq != 4:` guard around the `while` loop.
- The top-reaching `WRITE(LOGDEV, ...)` diagnostic is reproduced exactly by the
  harness and checked against the original formatted record.

## Preserve two surprising source behaviors

1. If the estimated centerline is at/below the first interface and the exit
   is hotter than TA(1), the source assigns **ZPLM=ZH(2)**. It then still
   computes final centerline **STKHT+(2/3)*ZPLM**. The result is not simply
   ZH(2); the Python port deliberately does not "correct" this formula.
2. Model-top clipping applies to plume-top **rise**, not directly to the
   final absolute centerline. For a stack above the domain, a negative
   ZSTK(last) can yield a final centerline below the stack. This reproduces
   the routine, not a certified emissions allocation. The production caller
   has additional checks; no success-shaped layer fractions are invented.

`PlumeRise` returns the final centerline, updated reporting WSTK, and a
`top_reached` flag. It also provides `stability_path`, `layer_transitions`,
and `hot_layer_adjusted` for teaching. Those three are Python internal
diagnostics, **not independently observed Fortran module outputs**.
The top-reaching flag is independently checked against the actual original
routine's log output for each call.

## Runnable lesson

Run `python -m cmaq_port lesson`, section 10:

```text
DELTA_ZS -> PLMRIS -> PLSPRD
```

The synthetic example uses ZF=(0,20,60,140,300,600) m, layer centers
(10,40,100,220,450) m, gradients (0,0.02,-0.005,0.0001,0.01) K/m,
TA=290-0.006*ZH K, and winds (0.2,1,5,8,2) m/s.
The stack is 20 m tall and 2 m wide, emitting at 450 K and 10 m/s.
Ambient stack temperature is 290 K, WSTK=4 m/s, USTAR=0.2 m/s,
HFX=0.03 m K/s, and HMIX=800 m.

- PLMRIS centerline: approximately **100.348648 m**.
- Reporting WSTK: **4 m/s**.
- Python IQ path: **(1,3,1,2)** through three layer transitions.
- Alternate PLSPRD bounds: approximately **250.559677 m** and **0 m**.

This example describes plume geometry only, not concentration or layer
emissions. A larger height is not proof of a smaller ground-level impact.

### Exercises

1. Set STKTK=TSTK. Which computations does the early return bypass?
2. Change WSTK below 1 m/s; distinguish reporting wind from upper-layer wind.
3. Track IQ as the plume crosses the stable second layer into the third.
4. Explain why PLSPRD can return bounds in a different layer from the stack.
5. Compare the hot-first-layer result with ZH(2); show the final conversion.
6. Contrast BFLX<=0 here with AERMOD's positive FLUXES floor.

## Verification and independent review

[Latest combined report](../runs/20261024T225921-kernels-e99396a9/kernels.json):

| Operation | Cases | Checked values | Maximum absolute difference |
| --- | --- | --- | --- |
| Original PLSPRD | 44 | 88 | 0 |
| Original DELTA_ZS | 64 | 720 | 0 |
| Original PLMRIS | 488 | 1,464 | 0 |
| Connected DELTA_ZS -> PLMRIS -> PLSPRD adapter | 488 | 6,328 | 0 |
| Original FIRE_PLMRIS | 224 | 672 | 0.000030517578125 |
| **Total** | **1,308** | **9,272** | **0.000030517578125** |

The connected adapter is not a fourth source routine. The matrix spans
neutral, stable-heat-flux, stable-gradient, unstable, and changing-layer
profiles; interface equality, above-domain stacks, cold/equal/hot exhaust,
wind floors, hot-first-layer adjustments, momentum shortcuts, and model-top
diagnostics. Eight targeted inputs also check exact/adjacent +/-HCRIT and
0.001-gradient thresholds. There are 136 actual top-reaching log messages
across PLMRIS, connected, and FIRE_PLMRIS calls; the eight fire logs match
the exact original formatted records. All **104 unit tests pass**, including
**19 dedicated PLMRIS tests**, **19 FIRE_PLMRIS tests**, and **12 PREPLM tests**.

The earlier independent high-capability model review traced PLMRIS equations,
loop, caller, adapters, and tests; it did not review the new FIRE_PLMRIS port.
The independent review of FIRE_PLMRIS remains pending. Deliberately
constructed exact ties between competing rise formulas remain a coverage
limitation.

Original routines and CONST.EXT are copied unchanged into each run; source,
include, compiler, driver, and logging-adapter hashes, commands, inputs,
outputs, logs, and per-case criteria persist alongside the report.
The test-only RUNTIME_VARS module exports **only LOGDEV**, directed to a
persistent log. It does not implement the production runtime, replace
physics, or establish a full CMAQ installation.

Compilation uses runtime checks, no fast math, extended fixed lines, and a
preprocessor definition mapping SUBST_CONST to the unchanged CONST.EXT.
The existing local GNU Fortran toolchain runs natively on Windows.
PLMRIS comparison thresholds are `atol=2e-4`, `rtol=2e-6`; the top flag
must match exactly, and connected geometry/bounds retain their tighter
DELTA_ZS/PLSPRD criteria. These thresholds are engineering checks in each
output's units, not scientific acceptance criteria. PLMRIS and all earlier
operations remain exact in this run; FIRE_PLMRIS observed 3.05e-5 maximum
absolute difference against `atol=2e-4`, `rtol=2e-6`. Cross-compiler
bitwise identity is not claimed.

No full CMAQ run, real-state capture, emissions allocation, or actual
CMAQ/AERMOD concentration comparison has been completed. (PREPLM itself is now
ported and validated bitwise against the unchanged `preplm.f`; see
[PREPLM](preplm.md).) Original source remains unchanged; WSL is untouched.
