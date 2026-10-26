import math
import unittest

from cmaq_port.plume import delta_zs, real32
from cmaq_port.kernel_validation import geometry_cases
from cmaq_port.__main__ import numerical_lesson


class StackGeometryTests(unittest.TestCase):
    def test_first_layer_measures_thickness_from_ground(self) -> None:
        result = delta_zs((0,), (30,), ((100,), (500,)))
        self.assertEqual(result.relative_heights, ((70,), (470,)))
        self.assertEqual(result.inverse_thicknesses, ((real32(0.01),), (real32(0.0025),)))

    def test_map_selects_global_stack_not_meteorology_column(self) -> None:
        result = delta_zs((2, 0), (10, 100, 300), ((100, 120), (500, 600)))
        self.assertEqual(result.relative_heights, ((-200, 110), (200, 590)))
        self.assertEqual(result.inverse_thicknesses[1], (real32(1 / 400), real32(1 / 480)))

    def test_repeated_source_map_uses_distinct_local_columns(self) -> None:
        self.assertEqual(delta_zs((1, 1), (0, 150), ((100, 200),)).relative_heights, ((-50, 50),))

    def test_stack_above_domain_is_not_clipped(self) -> None:
        self.assertEqual(delta_zs((0,), (1500,), ((100,), (1000,))).relative_heights, ((-1400,), (-500,)))

    def test_no_local_sources_preserves_layer_shape(self) -> None:
        result = delta_zs((), (), ((), ()))
        self.assertEqual(result.relative_heights, ((), ()))
        self.assertEqual(result.inverse_thicknesses, ((), ()))

    def test_subtraction_and_reciprocal_are_rounded_to_real32(self) -> None:
        result = delta_zs((0,), (0.1,), ((0.3,), (0.7,)))
        self.assertEqual(result.relative_heights[0][0], real32(real32(0.3) - real32(0.1)))
        self.assertEqual(result.inverse_thicknesses[1][0], real32(1 / real32(real32(0.7) - real32(0.3))))

    def test_invalid_maps_raise(self) -> None:
        for mapping in ((-1,), (1,), (True,), (0.5,)):
            with self.subTest(mapping=mapping), self.assertRaises(ValueError):
                delta_zs(mapping, (10,), ((100,),))  # type: ignore[arg-type]

    def test_invalid_grids_raise(self) -> None:
        for grid in ((), ((0,),), ((100,), (100,)), ((100,), (99,)), ((100, 200),), ((math.nan,),)):
            with self.subTest(grid=grid), self.assertRaises(ValueError):
                delta_zs((0,), (10,), grid)

    def test_real32_collapsed_interfaces_raise(self) -> None:
        with self.assertRaises(ValueError):
            delta_zs((0,), (10,), ((100,), (100.0000001,)))

    def test_invalid_stack_and_overflow_raise(self) -> None:
        for stack in (-1, math.inf, 1e100):
            with self.subTest(stack=stack), self.assertRaises(ValueError):
                delta_zs((0,), (stack,), ((100,),))
        with self.assertRaises(ValueError):
            delta_zs((0,), (0,), ((1e-40,),))

    def test_reference_matrix_shape_mapping_and_counts(self) -> None:
        samples = geometry_cases()
        self.assertEqual(samples, geometry_cases())
        self.assertEqual(len(samples), 64)
        self.assertEqual(sum(len(sample.expected) for sample in samples), 720)
        self.assertEqual(len({sample.label for sample in samples}), 64)
        self.assertTrue(all(sample.protocol.startswith("DELTA_ZS\n") for sample in samples))
        lesson = numerical_lesson()
        self.assertIn("NOT the already-local meteorology columns", lesson)
        self.assertIn("no ground row", lesson)
