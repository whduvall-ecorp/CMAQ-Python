"""I/O API POLY interpolation used by CMAQ PREPLM; see docs/interpolation.md.

Faithful port of the EDSS/Models-3 I/O API ``REAL FUNCTION POLY`` (Newton
divided-difference interpolation). Every elementary REAL operation rounds to
binary32 to match the single-precision Fortran result, including the source's
unrolled degree-1 term and the explicit divided-difference sums.
"""

from collections.abc import Sequence

from .plume import real32


def poly(xpt: float, xpts: Sequence[float], ypts: Sequence[float], ndeg: int) -> float:
    """Interpolate YPTS(XPTS) at XPT with a degree-NDEG polynomial.

    ``xpts`` and ``ypts`` must each hold ``ndeg + 1`` points and ``ndeg >= 1``,
    exactly as the I/O API routine requires (it unrolls the linear term and
    indexes ``XPTS(2)``). The abscissae must be distinct after REAL rounding;
    the source divides by their differences and does not guard against ties.
    """
    if isinstance(ndeg, bool) or not isinstance(ndeg, int):
        raise ValueError("NDEG must be an integer.")
    if ndeg < 1:
        raise ValueError("POLY requires NDEG >= 1 (needs at least two points).")
    x = real32(xpt)
    xs = tuple(real32(value) for value in xpts)
    ys = tuple(real32(value) for value in ypts)
    if len(xs) != ndeg + 1 or len(ys) != ndeg + 1:
        raise ValueError("XPTS and YPTS must each contain NDEG + 1 values.")
    if len(set(xs)) != len(xs):
        raise ValueError("XPTS must be distinct after REAL rounding.")

    # Degree-1 (linear) term, unrolled exactly as the Fortran source.
    xscr = xs[0]
    coeff = real32(x - xscr)
    dscr = real32(coeff / real32(xscr - xs[1]))
    yscr = ys[0]
    result = real32(yscr + real32(dscr * real32(yscr - ys[1])))

    # Higher-order terms: result += COEFF_i * f[x_1, ..., x_{i+1}].
    for i in range(2, ndeg + 1):
        # Divided difference f[x_1, ..., x_{i+1}] via the explicit sum; the
        # K = 1 denominator is unrolled first, then each remaining point K.
        xscr = xs[0]
        dscr = real32(xscr - xs[1])
        for j in range(3, i + 2):
            dscr = real32(dscr * real32(xscr - xs[j - 1]))
        tdiff = real32(ys[0] / dscr)
        for k in range(2, i + 2):
            xscr = xs[k - 1]
            dscr = real32(xscr - xs[0])
            for j in range(2, k):
                dscr = real32(dscr * real32(xscr - xs[j - 1]))
            for j in range(k + 1, i + 2):
                dscr = real32(dscr * real32(xscr - xs[j - 1]))
            tdiff = real32(tdiff + real32(ys[k - 1] / dscr))
        coeff = real32(coeff * real32(x - xs[i - 1]))
        result = real32(result + real32(coeff * tdiff))

    return result
