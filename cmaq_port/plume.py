"""CMAQ/SMOKE plume helpers; see docs/plume-spread.md for attribution."""

import math
import struct
from collections.abc import Sequence
from dataclasses import dataclass


def real32(value: float) -> float:
    """Round a value to source default REAL precision without a NumPy dependency."""
    if not math.isfinite(value):
        raise ValueError("Default REAL values must be finite.")
    try:
        rounded = struct.unpack("=f", struct.pack("=f", value))[0]
    except OverflowError as exc:
        raise ValueError("Value exceeds default REAL range.") from exc
    if not math.isfinite(rounded):
        raise ValueError("Value exceeds default REAL range.")
    return rounded


@dataclass(frozen=True)
class PlumeBounds:
    """PLTOP/PLBOT (m), not a Gaussian density or layer emission fraction."""

    top: float
    bottom: float


@dataclass(frozen=True)
class StackLayerGeometry:
    """Layer-major ZSTK (m) and DDZF (1/m), columns in local-source order."""

    relative_heights: tuple[tuple[float, ...], ...]
    inverse_thicknesses: tuple[tuple[float, ...], ...]


def delta_zs(
    source_map: Sequence[int], stack_heights: Sequence[float],
    layer_heights: Sequence[Sequence[float]],
) -> StackLayerGeometry:
    """DELTA_ZS with zero-based SRC_MAP and ZF[layer][local source].

    ZF contains upper interfaces only, unlike PLSPRD's ground-inclusive grid.
    Source-map entries index the stack list, not the local meteorology columns.
    """
    if not layer_heights:
        raise ValueError("DELTA_ZS requires at least one layer.")
    stacks = tuple(real32(value) for value in stack_heights)
    if any(value < 0 for value in stacks):
        raise ValueError("Stack heights must be nonnegative.")
    if any(
        isinstance(index, bool) or not isinstance(index, int) or not 0 <= index < len(stacks)
        for index in source_map
    ):
        raise ValueError("Source map must contain zero-based integer stack indices.")
    heights = tuple(tuple(real32(value) for value in row) for row in layer_heights)
    if any(len(row) != len(source_map) for row in heights):
        raise ValueError("Each layer must contain one height per local source.")
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
    return StackLayerGeometry(tuple(relative), tuple(inverse))


def plsprd(
    gradients: Sequence[float], layer_heights: Sequence[float], effective_height: float,
) -> PlumeBounds:
    """Initial plume bounds using DTHDZ (K/m) and ZF (m, including index zero).

    KZ is inferred from gradients. Equality with an interface selects the
    lower layer. Preserve source default REAL rounding and clipping order.
    """
    if not gradients or len(layer_heights) != len(gradients) + 1:
        raise ValueError("Provide KZ gradients and KZ+1 interfaces including ZF(0).")
    lapse = tuple(real32(value) for value in gradients)
    heights = tuple(real32(value) for value in layer_heights)
    center = real32(effective_height)
    if center < 0 or heights[0] != 0 or heights[-1] < 1:
        raise ValueError("Require nonnegative effective height, ZF(0)=0, and model top >=1 m.")
    if any(upper <= lower for lower, upper in zip(heights, heights[1:])):
        raise ValueError("Layer interfaces must increase strictly after REAL rounding.")
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
