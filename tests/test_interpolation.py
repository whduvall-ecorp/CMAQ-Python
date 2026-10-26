import unittest

from cmaq_port.interpolation import poly
from cmaq_port.plume import real32


class PolyInterpolationTests(unittest.TestCase):
    def test_reproduces_data_points_exactly(self) -> None:
        xs = (0.0, 1.0, 2.0, 3.0)
        ys = (1.0, 3.0, 9.0, 31.0)
        for ndeg in (1, 2, 3):
            window, values = xs[: ndeg + 1], ys[: ndeg + 1]
            for x, y in zip(window, values, strict=True):
                self.assertAlmostEqual(poly(x, window, values, ndeg), y, places=3)

    def test_linear_interpolation_matches_closed_form(self) -> None:
        # P(x) = y1 + (x - x1)/(x1 - x2) * (y1 - y2) in default REAL arithmetic.
        x1, x2, y1, y2, xp = 0.0, 1.0, 10.0, 20.0, 0.25
        dscr = real32(real32(xp - x1) / real32(x1 - x2))
        expected = real32(y1 + real32(dscr * real32(y1 - y2)))
        self.assertEqual(poly(xp, (x1, x2), (y1, y2), 1), expected)

    def test_quadratic_is_reproduced_within_float32(self) -> None:
        fn = lambda x: 2.0 * x * x - 3.0 * x + 5.0
        window = (-1.0, 0.5, 4.0)
        values = tuple(fn(x) for x in window)
        for x in (-1.0, 0.0, 1.7, 4.0, 2.3):
            self.assertAlmostEqual(poly(x, window, values, 2), fn(x), places=2)

    def test_cubic_is_reproduced_within_float32(self) -> None:
        fn = lambda x: x ** 3 - 2.0 * x + 1.0
        window = (-2.0, -0.5, 1.0, 3.0)
        values = tuple(fn(x) for x in window)
        worst = max(abs(poly(x, window, values, 3) - fn(x))
                    for x in (-2.0, -1.3, 0.0, 0.7, 1.9, 3.0))
        self.assertLess(worst, 5e-3)

    def test_extrapolation_follows_the_polynomial(self) -> None:
        fn = lambda x: 3.0 * x - 1.0
        self.assertAlmostEqual(poly(5.0, (0.0, 1.0), (fn(0.0), fn(1.0)), 1), fn(5.0), places=3)

    def test_requires_degree_at_least_one(self) -> None:
        with self.assertRaises(ValueError):
            poly(0.0, (1.0,), (2.0,), 0)
        with self.assertRaises(ValueError):
            poly(0.0, (1.0, 2.0), (3.0, 4.0), -1)

    def test_rejects_non_integer_degree(self) -> None:
        with self.assertRaises(ValueError):
            poly(0.0, (0.0, 1.0), (0.0, 1.0), True)  # type: ignore[arg-type]
        with self.assertRaises(ValueError):
            poly(0.0, (0.0, 1.0), (0.0, 1.0), 1.0)  # type: ignore[arg-type]

    def test_rejects_point_count_mismatch(self) -> None:
        with self.assertRaises(ValueError):
            poly(0.0, (0.0, 1.0, 2.0), (0.0, 1.0), 2)
        with self.assertRaises(ValueError):
            poly(0.0, (0.0, 1.0), (0.0, 1.0, 2.0), 2)

    def test_rejects_duplicate_abscissae(self) -> None:
        with self.assertRaises(ValueError):
            poly(0.5, (1.0, 1.0), (2.0, 3.0), 1)
        with self.assertRaises(ValueError):
            poly(0.5, (0.0, 2.0, 2.0), (1.0, 2.0, 3.0), 2)

    def test_inputs_are_rounded_to_default_real(self) -> None:
        # A value distinct in float64 but equal in float32 must be accepted and
        # must collapse to the same result as its rounded form.
        near = 1.0 + 1e-9
        self.assertEqual(real32(near), real32(1.0))
        self.assertEqual(
            poly(near, (0.0, 2.0), (0.0, 4.0), 1),
            poly(1.0, (0.0, 2.0), (0.0, 4.0), 1),
        )


if __name__ == "__main__":
    unittest.main()
