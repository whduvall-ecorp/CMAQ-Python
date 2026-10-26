"""CMAQ/SMOKE PREPLM: stack-top meteorology for PLMRIS. See docs/preplm.md.

Faithful port of CCTM/src/plrise/smoke/preplm.f. It derives wind speed, the
potential-temperature gradient, and stack-top temperature/wind from per-source
meteorology, interpolating to stack top with the I/O API POLY routine. Every
elementary REAL operation rounds to binary32 to match the single-precision
Fortran result. This is the meteorological preparation PLMRIS consumes; it is
not PLMRIS, a concentration, or a CMAQ run.
"""

import math
from collections.abc import Sequence
from dataclasses import dataclass

from .interpolation import poly
from .plume import real32

DEG = 3
CTOK = real32(273.15)


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


@dataclass(frozen=True)
class StackMeteorology:
    """PREPLM outputs; LSTK/LPBL are 1-based layer indices like the Fortran."""

    lstk: int
    lpbl: int
    stack_temperature: float
    stack_wind: float
    theta_gradient: tuple[float, ...]
    wind_speed: tuple[float, ...]


def preplm(
    fireflg: bool,
    emlays: int,
    hmix: float,
    hts: float,
    psfc: float,
    ts: float,
    ddzf: Sequence[float],
    qv: Sequence[float],
    ta: Sequence[float],
    uw: Sequence[float],
    vw: Sequence[float],
    zh: Sequence[float],
    zf: Sequence[float],
    pres: Sequence[float],
) -> StackMeteorology:
    """Compute PLMRIS stack-top inputs from per-source meteorology.

    ``ddzf``/``qv``/``ta``/``uw``/``vw``/``zh``/``zf`` hold ``emlays`` values;
    ``pres`` mirrors the Fortran ``PRES(0:)`` with ``emlays + 1`` entries, of
    which ``pres[0]`` (``PRES(0)``) is unused by this code path. ``emlays >= 2``
    is required, and ``emlays >= DEG + 1`` (4) for the non-fire polynomial
    interpolation to stack top.
    """
    if isinstance(fireflg, int) and not isinstance(fireflg, bool):
        raise ValueError("fireflg must be a bool.")
    if not isinstance(fireflg, bool):
        raise ValueError("fireflg must be a bool.")
    if isinstance(emlays, bool) or not isinstance(emlays, int) or emlays < 2:
        raise ValueError("emlays must be an integer >= 2.")
    if not fireflg and emlays < DEG + 1:
        raise ValueError("Non-fire PREPLM needs emlays >= 4 for POLY interpolation.")
    profiles = {"ddzf": ddzf, "qv": qv, "ta": ta, "uw": uw, "vw": vw, "zh": zh, "zf": zf}
    for name, values in profiles.items():
        if len(values) != emlays:
            raise ValueError(f"{name} must contain emlays values.")
    if len(pres) != emlays + 1:
        raise ValueError("pres must contain emlays + 1 values (PRES(0:emlays)).")

    hmix_r = real32(hmix)
    hts_r = real32(hts)
    # psfc is validated for interface fidelity; in the source it only fed the
    # commented-out surface-gradient term (see the note below), so it does not
    # affect any output here.
    real32(psfc)
    ts_r = real32(ts)
    ddzf_r = tuple(real32(v) for v in ddzf)
    qv_r = tuple(real32(v) for v in qv)
    ta_r = tuple(real32(v) for v in ta)
    uw_r = tuple(real32(v) for v in uw)
    vw_r = tuple(real32(v) for v in vw)
    zh_r = tuple(real32(v) for v in zh)
    zf_r = tuple(real32(v) for v in zf)
    pres_r = tuple(real32(v) for v in pres)

    # Wind speed and virtual temperature per layer.
    wspd: list[float] = []
    tv: list[float] = []
    for layer in range(emlays):
        p, q = uw_r[layer], vw_r[layer]
        wspd.append(real32(math.sqrt(_add(_mul(p, p), _mul(q, q)))))
        ratio = _div(qv_r[layer], _add(real32(1.0), qv_r[layer]))
        tv.append(_mul(ta_r[layer], _add(real32(1.0), _mul(real32(0.622), ratio))))

    # The source computes ES, QSFC, TVSFC, and THETG for a surface-gradient
    # term (DTHDZ(1) = (THV1 - THETG)/ZF(1)) that is commented out in preplm.f;
    # DTHDZ(1) is overridden with DTHDZ(2) below, so that chain has no effect on
    # any PREPLM output and is intentionally omitted here.
    lpbl = 1 if hmix_r <= zf_r[0] else 0
    lstk = 1 if hts_r <= zf_r[0] else 0

    # Virtual temperature interpolated to the full-layer face heights.
    tf: list[float] = [0.0] * emlays
    for layer in range(emlays - 1):
        delz = _sub(zh_r[layer + 1], zh_r[layer])
        tf[layer] = _add(
            tv[layer],
            _div(_mul(_sub(tv[layer + 1], tv[layer]), _sub(zf_r[layer], zh_r[layer])), delz),
        )
    last = emlays - 1
    delz = _sub(zh_r[last], zh_r[last - 1])
    tf[last] = _add(
        tv[last],
        _div(_mul(_sub(tv[last], tv[last - 1]), _sub(zf_r[last], zh_r[last])), delz),
    )

    # Potential-temperature gradient; PRES is 1-based PRES(1..emlays).
    thv1 = _mul(tf[0], _pow(_div(real32(1000.0), pres_r[1]), real32(0.286)))
    dthdz: list[float] = [0.0] * emlays
    for fortran_l in range(2, emlays + 1):
        if hmix_r > zf_r[fortran_l - 2]:
            lpbl = fortran_l
        if hts_r > zf_r[fortran_l - 2]:
            lstk = fortran_l
        thvk = _mul(tf[fortran_l - 1], _pow(_div(real32(1000.0), pres_r[fortran_l]), real32(0.286)))
        dthdz[fortran_l - 1] = _mul(ddzf_r[fortran_l - 1], _sub(thvk, thv1))
        thv1 = thvk
    dthdz[0] = dthdz[1]

    if not fireflg:
        # POLY uses DEG+1 points starting at M (1-based); it reads only the
        # first DEG+1 entries of the ZH(M:EMLAYS) slice.
        m = max(1, lstk - DEG - 1)
        start = m - 1
        window_zh = zh_r[start:start + DEG + 1]
        tstk = poly(hts_r, window_zh, ta_r[start:start + DEG + 1], DEG)
        wstk = real32(max(poly(hts_r, window_zh, tuple(wspd[start:start + DEG + 1]), DEG), real32(0.1)))
    else:
        tstk = ts_r
        wstk = wspd[0]

    return StackMeteorology(lstk, lpbl, tstk, wstk, tuple(dthdz), tuple(wspd))
