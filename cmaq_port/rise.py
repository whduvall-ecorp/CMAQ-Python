"""Full educational PLMRIS port; CMAQ/SMOKE attribution: docs/plume-rise.md.

Each elementary REAL operation rounds to binary32, including statement
function intermediates and the source's noninteger exponent constants.
"""

import math
from collections.abc import Sequence
from dataclasses import dataclass

from .plume import delta_zs, real32


def _add(a: float, b: float) -> float:
    return real32(a + b)


def _sub(a: float, b: float) -> float:
    return real32(a - b)


def _mul(a: float, b: float) -> float:
    return real32(a * b)


def _div(a: float, b: float) -> float:
    return real32(a / b)


def _pow(a: float, b: float) -> float:
    return real32(math.pow(a, b))


GRAV = real32(9.80622)
HCRIT = _mul(real32(1e-4), real32(0.03))
SMALL = real32(3e-5)
D3, D6 = _div(1.0, 3.0), _div(1.0, 6.0)
D45, D2664, D59319 = _div(1.0, 45.0), _div(1.0, real32(2.664)), _div(1.0, real32(59.319))
TWOTHD, FIVETHD = _div(2.0, 3.0), _div(5.0, 3.0)


def _neutral(h: float, b: float, u: float, us: float) -> float:
    ratio = _div(b, _mul(_mul(u, us), us))
    return min(_mul(10.0, h), _mul(real32(1.2), _mul(
        _pow(ratio, real32(0.6)),
        _pow(_add(h, _div(_mul(real32(1.3), b), _mul(_mul(u, us), us))), real32(0.4)),
    )))


def _stable(b: float, u: float, s: float) -> float:
    return _mul(real32(2.6), _pow(_div(b, _mul(u, s)), D3))


def _unstable(b: float, u: float) -> float:
    return _mul(30.0, _pow(_div(b, u), real32(0.6)))


@dataclass(frozen=True)
class PlumeRise:
    """Centerline (m), updated reporting wind (m/s), and source top diagnostic."""

    centerline: float
    stack_wind: float
    top_reached: bool
    stability_path: tuple[int, ...]
    layer_transitions: int
    hot_layer_adjusted: bool


def plmris(
    stack_layer: int, heat_flux: float, mixing_height: float,
    stack_diameter: float, stack_height: float, exhaust_temperature: float,
    exhaust_velocity: float, stack_temperature: float, friction_velocity: float,
    gradients: Sequence[float], temperatures: Sequence[float], winds: Sequence[float],
    interfaces: Sequence[float], centers: Sequence[float], stack_wind: float,
    relative_heights: Sequence[float] | None = None,
) -> PlumeRise:
    """PLMRIS: zero-based LSTK; ZF includes ground, other arrays have EMLAYS.

    HFX is already m K/s (not W/m2). TSTK/WSTK and LSTK normally come from
    PREPLM, which is not ported. Optional ZSTK must exactly match DELTA_ZS
    after REAL rounding; otherwise derive it. No caller-side USTARMIN is
    invented here. Diagnostic path uses the original IQ codes (1..4).
    """
    n = len(gradients)
    if n < 2 or any(len(v) != n for v in (temperatures, winds, centers)) or len(interfaces) != n + 1:
        raise ValueError("PLMRIS requires >=2 layers, matching met/center arrays and n+1 interfaces.")
    if isinstance(stack_layer, bool) or not isinstance(stack_layer, int) or not 0 <= stack_layer < n:
        raise ValueError("LSTK must be a zero-based integer layer index.")
    hfx, hmix, dm, ht, tk, ve, ts, us, ws = map(real32, (
        heat_flux, mixing_height, stack_diameter, stack_height, exhaust_temperature,
        exhaust_velocity, stack_temperature, friction_velocity, stack_wind,
    ))
    lapse, ta, wind, zf, zh = (tuple(map(real32, v)) for v in (
        gradients, temperatures, winds, interfaces, centers,
    ))
    if min(tk, ts, us, *ta) <= 0 or min(dm, ht, ve) < 0:
        raise ValueError("Temperatures and USTAR must be positive; diameter, height and velocity nonnegative.")
    if zf[0] != 0 or any(a >= b for a, b in zip(zf, zf[1:])):
        raise ValueError("ZF must start at zero and strictly increase after REAL rounding.")
    if any(not lower < center < upper for lower, center, upper in zip(zf, zh, zf[1:])):
        raise ValueError("Each ZH must lie strictly inside its layer.")
    expected_layer = next((i for i in range(n) if ht <= zf[i + 1]), n - 1)
    if stack_layer != expected_layer:
        raise ValueError("LSTK is inconsistent with stack height and interfaces.")
    geometry = delta_zs((0,), (ht,), tuple((z,) for z in zf[1:]))
    zstk = tuple(row[0] for row in geometry.relative_heights)
    if relative_heights is not None and tuple(map(real32, relative_heights)) != zstk:
        raise ValueError("ZSTK must equal DELTA_ZS(ZF, STKHT), including its shape.")

    hstar = _div(_mul(GRAV, hfx), ta[0])
    bflx = _div(_mul(_mul(_mul(_mul(_mul(0.25, GRAV), _sub(tk, ts)), ve), dm), dm), tk)
    layer = stack_layer
    ws = max(ws, 1.0)
    dhm = _div(_mul(_mul(3.0, dm), ve), ws)
    if bflx <= 0:
        return PlumeRise(_add(ht, max(dhm, 2.0)), ws, False, (), 0, False)
    sx = 0.0
    if hstar > HCRIT:
        if _sub(hmix, ht) <= 0:
            sx = max(_div(_mul(GRAV, lapse[layer]), ts), SMALL)
            if lapse[layer] > real32(0.001):
                product = _mul(_mul(_mul(ve, ve), dm), dm)
                dhsm = _div(_mul(_mul(real32(0.646), _pow(
                    _div(product, _mul(tk, ws)), D3)), real32(math.sqrt(ts))),
                    _pow(lapse[layer], D6))
            else:
                dhsm = dhm
            dhm = min(dhsm, dhm)
            dhn, dh = _neutral(ht, bflx, ws, us), _stable(bflx, ws, sx)
            iq = 2 if dhn < dh else 3
        else:
            dhn, dh = _neutral(ht, bflx, ws, us), _unstable(bflx, ws)
            iq = 2 if dhn < dh else 1
        dh = min(dhn, dh)
        if dh < dhm:
            dh, iq = dhm, 4
        dht = _mul(1.5, dh)
    elif hstar < -HCRIT or lapse[stack_layer] > real32(0.001):
        sx = max(_div(_mul(GRAV, lapse[stack_layer]), ts), SMALL)
        dhn = _mul(1.5, _neutral(ht, bflx, ws, us))
        dht = _mul(1.5, _stable(bflx, ws, sx))
        iq = 2 if dhn < dht else 3
        dht = min(dhn, dht)
    else:
        dht, iq = _mul(1.5, _neutral(ht, bflx, ws, us)), 2

    zplm, nn, top = dht, 0, False
    path = [iq]
    if iq != 4:
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
            if sx > SMALL:
                dhn = _mul(1.5, _neutral(ht, rbflx, wplm, us))
                dht = _mul(1.5, _stable(rbflx, wplm, sx))
                iq = 2 if dhn < dht else 3
                dht = min(dhn, dht)
            else:
                dhn, dh = _neutral(ht, rbflx, wplm, us), _unstable(rbflx, wplm)
                iq = 2 if dhn < dh else 1
                dht = _mul(1.5, min(dhn, dh))
            path.append(iq)
            zplm = _add(zstk[layer - 1], dht)
    hot = _add(ht, _mul(TWOTHD, zplm)) <= zf[1] and tk > ta[0]
    if hot:
        zplm = zh[1]
    return PlumeRise(_add(ht, _mul(TWOTHD, zplm)), ws, top, tuple(path), nn, hot)
