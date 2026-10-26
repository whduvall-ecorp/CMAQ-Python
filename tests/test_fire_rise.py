import math
import unittest

from cmaq_port.__main__ import numerical_lesson
from cmaq_port.fire_rise import FirePlumeRise, GRAV, HCRIT, fire_plmris
from cmaq_port.kernel_validation import fire_cases
from cmaq_port.plume import real32


GRID = (0.0, 20.0, 60.0, 140.0, 300.0, 600.0)
CENTERS = tuple((lower + upper) / 2 for lower, upper in zip(GRID, GRID[1:]))
TEMPERATURES = tuple(290.0 - 0.006 * height for height in CENTERS)


class FirePlumeRiseTests(unittest.TestCase):
    def calculate(
        self, *, layer: int = 0, heat: float = 0.0, mixing: float = 800.0,
        buoyancy: float = 1.0, fire_temperature: float = 450.0,
        friction: float = 0.3, gradients: tuple[float, ...] = (0.0,) * 5,
        temperatures: tuple[float, ...] = TEMPERATURES,
        winds: tuple[float, ...] = (0.0, 1.0, 5.0, 8.0, 2.0),
        interfaces: tuple[float, ...] = GRID, fire_wind: float = 0.2,
    ) -> FirePlumeRise:
        return fire_plmris(
            layer, heat, mixing, buoyancy, fire_temperature, friction,
            gradients, temperatures, winds, interfaces, fire_wind,
        )

    def test_zero_buoyancy_flux_is_not_given_a_stack_style_floor(self) -> None:
        result = self.calculate(buoyancy=0.0)
        self.assertEqual(result.centerline, 0.0)
        self.assertEqual(result.stability_path, (2,))
        self.assertFalse(result.top_reached)

    def test_nonpositive_buoyancy_is_not_given_a_stack_style_floor(self) -> None:
        negative = self.calculate(buoyancy=-0.1)
        self.assertLess(negative.centerline, 0.0)
        self.assertEqual(self.calculate(buoyancy=0.0).centerline, 0.0)

    def test_negative_buoyancy_fractional_power_is_reported_explicitly(self) -> None:
        with self.assertRaisesRegex(ValueError, "nonnegative buoyancy"):
            self.calculate(heat=0.03, buoyancy=-0.1)

    def test_neutral_branch_uses_fire_statement_function_and_two_thirds_conversion(self) -> None:
        result = self.calculate(layer=2)
        expected = real32(real32(1.3327804) / real32(real32(0.3) * real32(0.3)))
        self.assertAlmostEqual(result.centerline, expected, places=5)
        self.assertEqual(result.stability_path, (2,))

    def test_unstable_and_stable_initial_branches_are_distinct(self) -> None:
        unstable = self.calculate(heat=0.03, friction=0.01)
        stable = self.calculate(heat=-0.03, friction=0.01)
        self.assertNotEqual(unstable.centerline, stable.centerline)
        self.assertIn(unstable.stability_path[0], (1, 2))
        self.assertIn(stable.stability_path[0], (2, 3))

    def test_above_mixing_height_is_selected_by_nonpositive_hmix_not_fire_altitude(self) -> None:
        at_positive_mix = self.calculate(
            layer=4, heat=0.03, mixing=20.0, friction=0.01,
        )
        at_nonpositive_mix = self.calculate(
            layer=4, heat=0.03, mixing=0.0, friction=0.01,
        )
        self.assertNotEqual(at_positive_mix.centerline, at_nonpositive_mix.centerline)

    def test_hstar_threshold_is_strict_and_uses_surface_temperature(self) -> None:
        surface = 290.0
        temperatures = (surface, *TEMPERATURES[1:])
        boundary_hfx = real32(real32(HCRIT * real32(surface)) / GRAV)
        equal_positive = self.calculate(
            heat=boundary_hfx, friction=0.01, temperatures=temperatures,
        )
        above_positive = self.calculate(
            heat=real32(boundary_hfx * real32(1.001)), friction=0.01,
            temperatures=temperatures,
        )
        equal_negative = self.calculate(
            heat=-boundary_hfx, friction=0.01, temperatures=temperatures,
        )
        below_negative = self.calculate(
            heat=-real32(boundary_hfx * real32(1.001)), friction=0.01,
            temperatures=temperatures,
        )
        self.assertEqual(equal_positive.stability_path[0], 2)
        self.assertEqual(equal_negative.stability_path[0], 2)
        self.assertEqual(above_positive.stability_path[0], 1)
        self.assertEqual(below_negative.stability_path[0], 3)

    def test_gradient_classification_uses_strict_point_zero_zero_one(self) -> None:
        equal = self.calculate(gradients=(0.001,) * 5, buoyancy=100.0)
        above = self.calculate(gradients=(0.001001,) * 5, buoyancy=100.0)
        self.assertEqual(equal.stability_path[0], 2)
        self.assertEqual(above.stability_path[0], 3)

    def test_fire_wind_is_floored_for_reporting_without_stack_geometry(self) -> None:
        result = self.calculate(fire_wind=0.2, winds=(25, 1, 5, 8, 2))
        self.assertEqual(result.fire_wind, 1.0)
        self.assertNotEqual(result.centerline, 20.0)
        self.assertEqual(GRID, (0, 20, 60, 140, 300, 600))

    def test_returned_fire_wind_does_not_follow_layer_winds(self) -> None:
        result = self.calculate(fire_wind=4.0, winds=(0.0, 0.0, 0.0, 0.0, 0.0))
        self.assertEqual(result.fire_wind, 4.0)

    def test_zero_layer_wind_is_floored_when_a_layer_is_entered(self) -> None:
        result = self.calculate(
            heat=0.03, buoyancy=100.0, winds=(0.0,) * 5,
        )
        self.assertGreaterEqual(result.layer_transitions, 1)
        self.assertIn(result.stability_path[-1], (1, 2, 3))

    def test_top_reaching_clips_rise_and_reports_the_preconversion_height(self) -> None:
        result = self.calculate(heat=0.03, buoyancy=1e6)
        self.assertTrue(result.top_reached)
        self.assertEqual(result.top_log_height, GRID[-1])
        self.assertEqual(result.centerline, real32(real32(2.0 / 3.0) * GRID[-1]))

    def test_fire_centerline_does_not_add_a_stack_height(self) -> None:
        result = self.calculate(buoyancy=0.01)
        self.assertLess(result.centerline, 20.0)

    def test_one_layer_profile_is_supported(self) -> None:
        result = fire_plmris(
            0, 0.0, 800.0, 1.0, 450.0, 0.3, (0.0,), (290.0,), (0.0,),
            (0.0, 100.0), 0.2,
        )
        self.assertTrue(math.isfinite(result.centerline))

    def test_bad_shapes_and_layer_indices_are_rejected(self) -> None:
        with self.assertRaises(ValueError):
            self.calculate(temperatures=(290.0,))
        for index in (-1, 5, True):
            with self.subTest(index=index), self.assertRaises(ValueError):
                self.calculate(layer=index)

    def test_nonfinite_unrepresentable_and_nonphysical_inputs_are_rejected(self) -> None:
        for value in (math.nan, math.inf, 1e100):
            with self.subTest(value=value), self.assertRaises(ValueError):
                self.calculate(fire_temperature=value)
        for temperature in (0.0, -1.0):
            with self.subTest(temperature=temperature), self.assertRaises(ValueError):
                self.calculate(fire_temperature=temperature)
        with self.assertRaises(ValueError):
            self.calculate(friction=0)
        with self.assertRaises(ValueError):
            self.calculate(fire_wind=-1)
        with self.assertRaises(ValueError):
            self.calculate(winds=(0, -1, 0, 0, 0))

    def test_invalid_interfaces_and_layer_temperatures_are_rejected(self) -> None:
        for interfaces in (
            (1, 20, 60, 140, 300, 600),
            (0, 20, 20, 140, 300, 600),
            (0, 20, 20.0000001, 140, 300, 600),
        ):
            with self.subTest(interfaces=interfaces), self.assertRaises(ValueError):
                self.calculate(interfaces=interfaces)
        with self.assertRaises(ValueError):
            self.calculate(temperatures=(0, 280, 280, 280, 280))

    def test_deterministic_fortran_matrix_covers_key_branches_and_top_log(self) -> None:
        samples = fire_cases()
        self.assertEqual(samples, fire_cases())
        self.assertEqual(len(samples), 224)
        self.assertEqual(len({sample.label for sample in samples}), 224)
        paths = {
            code for sample in samples if sample.diagnostics is not None
            for code in sample.diagnostics["stability_path"]
        }
        self.assertEqual(paths, {1, 2, 3})
        self.assertTrue(any(sample.top_log is not None for sample in samples))
        self.assertTrue(any(sample.expected[2] == 0 for sample in samples))
        self.assertEqual(sum(len(sample.expected) for sample in samples), 672)
        for sign in (-1, 1):
            exact = next(
                sample for sample in samples
                if sample.label == f"boundary-hstar-{sign}-1.0"
            )
            lines = exact.protocol.splitlines()
            heat = real32(float(lines[2].split()[0]))
            surface_temperature = real32(float(lines[4].split()[0]))
            self.assertEqual(
                real32(real32(GRAV * heat) / surface_temperature),
                sign * HCRIT,
            )
        gradients = {
            sample.label: real32(float(sample.protocol.splitlines()[3].split()[0]))
            for sample in samples if sample.label.startswith("boundary-gradient-")
        }
        self.assertEqual(gradients["boundary-gradient-equal"], real32(0.001))
        self.assertGreater(gradients["boundary-gradient-above"], real32(0.001))

    def test_lesson_explains_distinct_fire_inputs_and_caller_owned_buoyancy(self) -> None:
        lesson = numerical_lesson()
        self.assertIn("FIRE_PLMRIS", lesson)
        self.assertIn("caller-supplied BFLX=", lesson)
        self.assertIn("with no stack height", lesson)


if __name__ == "__main__":
    unittest.main()
