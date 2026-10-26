import math
import unittest

from cmaq_port.__main__ import numerical_lesson
from cmaq_port.kernel_validation import rise_cases
from cmaq_port.plume import delta_zs, real32
from cmaq_port.rise import GRAV, HCRIT, PlumeRise, plmris


GRID = (0.0, 20.0, 60.0, 140.0, 300.0, 600.0)
CENTERS = (10.0, 40.0, 100.0, 220.0, 450.0)
TEMPERATURES = tuple(290 - 0.006 * height for height in CENTERS)


class PlumeRiseTests(unittest.TestCase):
    def calculate(
        self, *, layer: int = 0, heat: float = 0, mixing: float = 800,
        diameter: float = 2, height: float = 20, exhaust: float = 450,
        velocity: float = 10, ambient: float = 290, friction: float = 0.2,
        gradients: tuple[float, ...] = (0, 0, 0, 0, 0),
        temperatures: tuple[float, ...] = TEMPERATURES,
        winds: tuple[float, ...] = (0.2, 1, 5, 8, 2),
        interfaces: tuple[float, ...] = GRID, centers: tuple[float, ...] = CENTERS,
        wind: float = 4, relative: tuple[float, ...] | None = None,
    ) -> PlumeRise:
        return plmris(layer, heat, mixing, diameter, height, exhaust, velocity, ambient,
                      friction, gradients, temperatures, winds, interfaces, centers, wind, relative)

    def test_nonpositive_buoyancy_returns_momentum_height_and_updated_wind(self) -> None:
        result = self.calculate(exhaust=280, wind=0.2)
        self.assertEqual(result.centerline, 80)
        self.assertEqual(result.stack_wind, 1)
        self.assertEqual(result.stability_path, ())
        self.assertEqual(result.layer_transitions, 0)
        self.assertFalse(result.top_reached)
        self.assertFalse(result.hot_layer_adjusted)

    def test_equal_temperatures_take_early_return(self) -> None:
        self.assertEqual(self.calculate(exhaust=290).centerline, 35)
        self.assertEqual(self.calculate(exhaust=290).stability_path, ())

    def test_minimum_momentum_rise_is_two_meters(self) -> None:
        self.assertEqual(self.calculate(exhaust=280, velocity=0).centerline, 22)

    def test_zero_buoyancy_early_return_does_not_apply_hot_layer_adjustment(self) -> None:
        result = self.calculate(height=1, velocity=0)
        self.assertEqual(result.centerline, 3)
        self.assertFalse(result.hot_layer_adjusted)

    def test_wind_floor_changes_returned_wind_not_input_arrays(self) -> None:
        self.assertEqual(self.calculate(exhaust=280, wind=-1).stack_wind, 1)
        self.assertEqual(self.calculate(exhaust=280, wind=0).stack_wind, 1)
        self.assertEqual(self.calculate(exhaust=280, wind=8).stack_wind, 8)
        self.assertEqual(GRID, (0, 20, 60, 140, 300, 600))

    def test_hot_first_layer_preserves_source_assignment_then_center_conversion(self) -> None:
        result = self.calculate(height=1, exhaust=291, velocity=40, diameter=4, wind=0.2, friction=0.3)
        expected = real32(real32(1) + real32(real32(2 / 3) * real32(40)))
        self.assertEqual(result.centerline, expected)
        self.assertNotEqual(result.centerline, CENTERS[1])
        self.assertTrue(result.hot_layer_adjusted)

    def test_momentum_dominated_unstable_branch_skips_layer_iteration(self) -> None:
        result = self.calculate(height=1, exhaust=291, velocity=40, diameter=4,
                                wind=0.2, friction=0.3, heat=0.03)
        self.assertEqual(result.stability_path, (4,))
        self.assertEqual(result.layer_transitions, 0)
        self.assertEqual(result.centerline, 481)
        self.assertFalse(result.top_reached)

    def test_stable_and_neutral_starting_branches(self) -> None:
        self.assertEqual(self.calculate().stability_path[0], 2)
        self.assertEqual(self.calculate(heat=-0.03, gradients=(0.0001,) * 5).stability_path[0], 3)
        self.assertEqual(self.calculate(gradients=(0.02,) * 5).stability_path[0], 3)

    def test_unstable_stack_below_and_above_mixing_height(self) -> None:
        self.assertEqual(self.calculate(heat=0.03).stability_path[0], 1)
        self.assertEqual(self.calculate(heat=0.03, mixing=0).stability_path[0], 3)

    def test_layer_changing_profile_uses_multiple_branches(self) -> None:
        result = self.calculate(heat=0.03, gradients=(0, 0.02, -0.005, 0.0001, 0.01))
        self.assertEqual(result.stability_path, (1, 3, 1, 2))
        self.assertEqual(result.layer_transitions, 3)
        self.assertFalse(result.top_reached)
        self.assertEqual(result.stack_wind, 4)

    def test_above_domain_stack_retains_source_top_and_center_conversion(self) -> None:
        result = self.calculate(layer=4, height=620)
        self.assertTrue(result.top_reached)
        self.assertEqual(result.centerline, real32(620 + real32(real32(2 / 3) * -20)))
        self.assertLess(result.centerline, 620)

    def test_explicit_geometry_matches_derived_geometry(self) -> None:
        geometry = delta_zs((0,), (20,), tuple((height,) for height in GRID[1:]))
        relative = tuple(row[0] for row in geometry.relative_heights)
        self.assertEqual(self.calculate(relative=relative), self.calculate())
        with self.assertRaises(ValueError):
            self.calculate(relative=(1, 2, 3, 4, 5))

    def test_interface_equality_selects_lower_stack_layer(self) -> None:
        self.assertTrue(math.isfinite(self.calculate(height=20, layer=0).centerline))
        self.assertTrue(math.isfinite(self.calculate(height=20.25, layer=1).centerline))
        with self.assertRaises(ValueError):
            self.calculate(height=20, layer=1)

    def test_invalid_shapes_and_indices_raise(self) -> None:
        with self.assertRaises(ValueError):
            self.calculate(gradients=(0,))
        with self.assertRaises(ValueError):
            self.calculate(winds=(1, 2))
        for layer in (-1, 5, True):
            with self.subTest(layer=layer), self.assertRaises(ValueError):
                self.calculate(layer=layer)

    def test_invalid_grid_and_centers_raise(self) -> None:
        for interfaces in ((1, 20, 60, 140, 300, 600), (0, 20, 20, 140, 300, 600),
                           (0, 20, 20.0000001, 140, 300, 600)):
            with self.subTest(interfaces=interfaces), self.assertRaises(ValueError):
                self.calculate(interfaces=interfaces)
        with self.assertRaises(ValueError):
            self.calculate(centers=(20, 40, 100, 220, 450))

    def test_nonfinite_and_nonphysical_parameters_raise(self) -> None:
        for value in (math.nan, math.inf, 1e100):
            with self.subTest(value=value), self.assertRaises(ValueError):
                self.calculate(heat=value)
        with self.assertRaises(ValueError):
            self.calculate(friction=0)
        with self.assertRaises(ValueError):
            self.calculate(ambient=0)
        with self.assertRaises(ValueError):
            self.calculate(velocity=-1)
        with self.assertRaises(ValueError):
            self.calculate(diameter=-1)

    def test_matrix_covers_all_science_paths_and_connected_outputs(self) -> None:
        samples = rise_cases()
        self.assertEqual(samples, rise_cases())
        self.assertEqual(len(samples), 976)
        self.assertEqual(len({sample.label for sample in samples}), 976)
        isolated = [sample for sample in samples if sample.routine == "PLMRIS"]
        paths = {code for sample in isolated if sample.diagnostics is not None
                 for code in sample.diagnostics["stability_path"]}
        self.assertEqual(paths, {1, 2, 3, 4})
        self.assertTrue(any(sample.expected[2] == 1 for sample in isolated))
        self.assertTrue(any(sample.diagnostics["hot_layer_adjusted"] for sample in isolated
                            if sample.diagnostics is not None))
        self.assertTrue(any(sample.diagnostics["layer_transitions"] >= 3 for sample in isolated
                            if sample.diagnostics is not None))
        self.assertEqual(sum(len(sample.expected) for sample in samples), 7792)

    def test_original_matrix_hits_exact_and_adjacent_heat_and_gradient_thresholds(self) -> None:
        samples = [sample for sample in rise_cases()
                   if sample.routine == "PLMRIS" and sample.label.startswith("boundary-")]
        self.assertEqual(len(samples), 8)
        for sign in (-1, 1):
            sample = next(sample for sample in samples if sample.label.startswith(f"boundary-heat-{sign}-1.0-"))
            lines = sample.protocol.splitlines()
            heat = real32(float(lines[2].split()[0]))
            temperature = real32(float(lines[4].split()[0]))
            self.assertEqual(real32(real32(GRAV * heat) / temperature), sign * HCRIT)
        sample = next(sample for sample in samples if sample.label.startswith("boundary-gradient-equal"))
        self.assertEqual(real32(float(sample.protocol.splitlines()[3].split()[0])), real32(0.001))

    def test_lesson_calculates_plume_rise_and_distinguishes_alternate_spread(self) -> None:
        text = numerical_lesson()
        self.assertIn("DELTA_ZS -> PLMRIS -> PLSPRD", text)
        self.assertIn("HFX is m K/s, not W/m2", text)
        self.assertIn("PREPLM remains upstream", text)
        self.assertIn("nonzero-IPVERT", text)


if __name__ == "__main__":
    unittest.main()
