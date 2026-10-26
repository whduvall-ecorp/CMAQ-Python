import math
import unittest

from cmaq_port.plume import plsprd, real32
from cmaq_port.kernel_validation import cases
from cmaq_port.__main__ import numerical_lesson


class PlumeSpreadTests(unittest.TestCase):
    def test_zero_exponent_depth_and_default_real_constants(self) -> None:
        result = plsprd((0.0098,), (0, 1000), 100)
        depth = real32(real32(3.545) * 15)
        self.assertEqual(result.top, real32(100 + real32(depth / 2)))
        self.assertEqual(result.bottom, real32(100 - real32(depth / 2)))

    def test_spread_floor_is_ten_meters_before_depth_scaling(self) -> None:
        result = plsprd((0.1,), (0, 1000), 100)
        self.assertEqual(result.top, real32(100 + real32(real32(real32(3.545) * 10) / 2)))

    def test_ground_clipping(self) -> None:
        self.assertEqual(plsprd((0.0098,), (0, 1000), 0).bottom, 0)

    def test_model_top_clips_top_and_bottom_separately(self) -> None:
        result = plsprd((0.0098,), (0, 1000), 1500)
        self.assertEqual((result.top, result.bottom), (1000, 999))

    def test_interface_equality_selects_lower_layer(self) -> None:
        equal = plsprd((0.0098, 0.1), (0, 100, 1000), 100)
        above = plsprd((0.0098, 0.1), (0, 100, 1000), 100.01)
        self.assertGreater(equal.top, above.top)

    def test_last_layer_used_for_height_above_grid(self) -> None:
        self.assertEqual(
            plsprd((0, 0.1), (0, 100, 1000), 1200),
            plsprd((0.1,), (0, 1000), 1200),
        )

    def test_no_mixing_height_argument_or_lid_reflection_is_invented(self) -> None:
        result = plsprd((0.0098,), (0, 1000), 100)
        self.assertLess(result.bottom, 100)
        self.assertGreater(result.top, 100)

    def test_invalid_grids_raise(self) -> None:
        for gradients, heights in (((), (0,)), ((0,), (0,)), ((0, 0), (0, 100, 100)), ((0,), (1, 100)), ((0,), (0, 0.5))):
            with self.subTest(heights=heights), self.assertRaises(ValueError):
                plsprd(gradients, heights, 10)

    def test_nonfinite_or_invalid_height_raises(self) -> None:
        for height in (-1, math.nan, math.inf, 1e100):
            with self.subTest(height=height), self.assertRaises(ValueError):
                plsprd((0,), (0, 1000), height)
        with self.assertRaises(ValueError):
            plsprd((math.nan,), (0, 1000), 10)

    def test_gradient_overflow_is_explicit(self) -> None:
        with self.assertRaises(ValueError):
            plsprd((-100,), (0, 1000), 10)

    def test_comparison_matrix_is_deterministic(self) -> None:
        self.assertEqual(cases(), cases())
        self.assertEqual(len(cases()), 44)
        text = numerical_lesson()
        self.assertIn("bottom=999 m", text)
        self.assertIn("not the default Turner branch", text)


if __name__ == "__main__":
    unittest.main()
