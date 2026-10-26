"""CMAQ fire-source plume rise; source attribution is documented separately."""

import math
from collections.abc import Sequence
from dataclasses import dataclass

from .plume import real32

GRAV = real32(9.80622)
HCRIT = real32(real32(1.0e-4) * real32(0.03))
SMALL = real32(3.0e-5)
D3 = real32(1.0 / 3.0)
D6 = real32(1.0 / 6.0)
D45 = real32(1.0 / 45.0)
D2664 = real32(1.0 / 2.664)
D59319 = real32(1.0 / 59.319)
TWOTHD = real32(2.0 / 3.0)
FIVETHD = real32(5.0 / 3.0)


def _mul(left: float, right: float) -> float:
    return real32(real32(left) * real32(right))


def _div(left: float, right: float) -> float:
    return real32(real32(left) / real32(right))


def _add(left: float, right: float) -> float:
    return real32(real32(left) + real32(right))


def _sub(left: float, right: float) -> float:
    return real32(real32(left) - real32(right))


def _pow(base: float, exponent: float) -> float:
    try:
        return real32(math.pow(base, exponent))
    except (OverflowError, ValueError) as exc:
        raise ValueError("FIRE_PLMRIS requires nonnegative buoyancy for fractional powers.") from exc


def _neutral(buoyancy: float, wind: float, friction: float) -> float:
    denominator = _mul(_mul(wind, friction), friction)
    return _div(_mul(real32(1.3327804), buoyancy), denominator)


def _stable(buoyancy: float, wind: float, stability: float) -> float:
    return _mul(real32(2.6), _pow(_div(buoyancy, _mul(wind, stability)), D3))


def _unstable(buoyancy: float, wind: float) -> float:
    return _mul(real32(30.0), _pow(_div(buoyancy, wind), real32(0.6)))


@dataclass(frozen=True)
class FirePlumeRise:
    """FIRE_PLMRIS centerline and reporting-wind outputs, in m and m/s."""

    centerline: float
    fire_wind: float
    top_reached: bool
    top_log_height: float | None
    stability_path: tuple[int, ...]
    layer_transitions: int


def fire_plmris(
    fire_layer: int,
    heat_flux: float,
    mixing_height: float,
    buoyancy_flux: float,
    fire_temperature: float,
    friction_velocity: float,
    gradients: Sequence[float],
    temperatures: Sequence[float],
    winds: Sequence[float],
    interfaces: Sequence[float],
    fire_wind: float,
) -> FirePlumeRise:
    """Port FIRE_PLMRIS, using a zero-based fire layer and ground-inclusive ZF.

    ``heat_flux`` is already HFX in m K/s and ``buoyancy_flux`` is caller-
    supplied BFLX in m**4/s**3. Unlike stack PLMRIS this routine receives no
    stack dimensions, derives no BFLX, and returns its final two-thirds
    conversion without adding a stack height.
    """
    layers = len(gradients)
    if layers < 1 or len(temperatures) != layers or len(winds) != layers or len(interfaces) != layers + 1:
        raise ValueError("FIRE_PLMRIS requires matching meteorology arrays and n+1 ZF interfaces.")
    if isinstance(fire_layer, bool) or not isinstance(fire_layer, int) or not 0 <= fire_layer < layers:
        raise ValueError("Fire layer must be a zero-based integer in the model layer range.")

    hfx, hmix, bflx, tstk, ustar, wstk = map(real32, (
        heat_flux, mixing_height, buoyancy_flux, fire_temperature, friction_velocity, fire_wind,
    ))
    lapse, ta, wspd, zf = (
        tuple(map(real32, values)) for values in (gradients, temperatures, winds, interfaces)
    )
    if tstk <= 0 or ustar <= 0 or any(value <= 0 for value in ta):
        raise ValueError("Fire and ambient temperatures and USTAR must be positive.")
    if wstk < 0 or any(value < 0 for value in wspd):
        raise ValueError("Fire and layer wind speeds must be nonnegative.")
    if zf[0] != 0 or any(lower >= upper for lower, upper in zip(zf, zf[1:])):
        raise ValueError("ZF must start at zero and increase strictly after REAL rounding.")

    hstar = _div(_mul(GRAV, hfx), ta[0])
    layer = fire_layer
    wstk = max(wstk, real32(1.0))
    dhm = real32(0.0)
    sx = real32(0.0)

    if hstar > HCRIT:
        if hmix <= 0:
            sx = max(_div(_mul(GRAV, lapse[layer]), tstk), SMALL)
            # LPLM is initialized to LSTK above, so the original wind-reset
            # condition (LPLM != LSTK) is unreachable on this first pass.
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

    buoyancy = bflx
    zplm = dht
    transitions = 0
    top_log_height: float | None = None
    path = [iq]
    if iq != 4:
        while True:
            # Fortran LPLM is one-based and indexes ZF(LPLM), the upper
            # interface of that layer; Python layer is zero-based.
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
            # ZF( LPLM-1 ) in the zero-based-grid source corresponds to
            # Python ZF[layer] after the layer index has advanced.
            zplm = _add(zf[layer], dht)

    zplm = _mul(TWOTHD, zplm)
    return FirePlumeRise(
        zplm, wstk, top_log_height is not None, top_log_height, tuple(path), transitions,
    )
