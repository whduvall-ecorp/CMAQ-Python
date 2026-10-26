import copy
import json
import math
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from cmaq_port.__main__ import LESSON, main
from cmaq_port.comparison import compare, load_dataset
from cmaq_port.inventory import inventory, save_report, sha256


def dataset(model: str) -> dict:
    return {
        "schema_version": 1, "model": model,
        "provenance": "SYNTHETIC TEST ONLY: no model was run",
        "contract": {
            "tracer_id": "synthetic-inert-gas", "source_id": "test-stack",
            "emissions_id": "test-emissions-v1", "meteorology_case_id": "test-met-v1",
            "spatial_support_id": "test-area-average-v1",
            "vertical_support_id": "test-layer-average-v1",
            "processes_id": "test-no-chemistry-or-deposition",
            "units": "ug/m3", "quantity": "passive_tracer_source_increment",
        },
        "records": [
            {"start_utc": "2018-07-01T00:00:00Z", "end_utc": "2018-07-01T01:00:00Z",
             "support_id": "cell-1-layer-1", "value": 1.0},
            {"start_utc": "2018-07-01T01:00:00Z", "end_utc": "2018-07-01T02:00:00Z",
             "support_id": "cell-1-layer-1", "value": 3.0},
        ],
    }


class ComparisonTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.left = self.root / "cmaq.json"
        self.right = self.root / "aermod.json"
        self.cmaq = dataset("CMAQ")
        self.aermod = dataset("AERMOD")
        self.write()

    def write(self) -> None:
        self.left.write_text(json.dumps(self.cmaq), encoding="utf-8")
        self.right.write_text(json.dumps(self.aermod), encoding="utf-8")

    def test_identical_values_have_zero_errors(self) -> None:
        report = compare(self.left, self.right)
        self.assertEqual(report["pairs"], 2)
        self.assertEqual(report["root_mean_square_error"], 0)
        self.assertEqual(report["mean_bias"], 0)
        inputs = report["inputs"]
        self.assertIsInstance(inputs, list)
        assert isinstance(inputs, list)
        self.assertEqual(inputs[0]["sha256"], sha256(self.left))
        self.assertNotIn("passed", report)

    def test_known_metrics_and_direction(self) -> None:
        self.cmaq["records"][0]["value"] = 3
        self.cmaq["records"][1]["value"] = 1
        self.write()
        report = compare(self.left, self.right)
        self.assertEqual(report["mean_bias"], 0)
        self.assertEqual(report["mean_absolute_error"], 2)
        self.assertEqual(report["root_mean_square_error"], 2)
        self.assertEqual(report["maximum_absolute_difference"], 2)
        self.assertEqual(report["difference_direction"], "CMAQ minus AERMOD")

    def test_row_order_and_equivalent_utc_notation_do_not_change_pairing(self) -> None:
        self.cmaq["records"].reverse()
        for row in self.cmaq["records"]:
            row["start_utc"] = row["start_utc"].replace("Z", "+00:00")
        self.write()
        self.assertEqual(compare(self.left, self.right)["root_mean_square_error"], 0)

    def test_every_contract_mismatch_is_rejected(self) -> None:
        for key in self.cmaq["contract"]:
            with self.subTest(key=key):
                original = self.cmaq["contract"][key]
                self.cmaq["contract"][key] = "different"
                self.write()
                with self.assertRaises(ValueError):
                    compare(self.left, self.right)
                self.cmaq["contract"][key] = original

    def test_missing_pairs_are_not_silently_intersected(self) -> None:
        self.cmaq["records"].pop()
        self.write()
        with self.assertRaisesRegex(ValueError, "differ"):
            compare(self.left, self.right)

    def test_spatial_and_time_keys_must_match_even_if_contract_ids_match(self) -> None:
        for key, value in (("support_id", "other-cell"), ("end_utc", "2018-07-01T00:30:00Z")):
            with self.subTest(key=key):
                self.cmaq = dataset("CMAQ")
                self.cmaq["records"][0][key] = value
                self.write()
                with self.assertRaises(ValueError):
                    compare(self.left, self.right)

    def test_model_order_is_enforced(self) -> None:
        with self.assertRaisesRegex(ValueError, "CMAQ first"):
            compare(self.right, self.left)

    def test_duplicate_keys_are_rejected(self) -> None:
        self.cmaq["records"].append(copy.deepcopy(self.cmaq["records"][0]))
        self.write()
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            compare(self.left, self.right)

    def test_overlapping_intervals_are_rejected(self) -> None:
        self.cmaq["records"][1]["start_utc"] = "2018-07-01T00:30:00Z"
        self.write()
        with self.assertRaisesRegex(ValueError, "Overlapping"):
            compare(self.left, self.right)

    def test_invalid_or_non_utc_time_is_rejected(self) -> None:
        for timestamp in ("not-a-time", "2018-07-01T00:00:00", "2018-07-01T00:00:00-08:00"):
            with self.subTest(timestamp=timestamp):
                self.cmaq["records"][0]["start_utc"] = timestamp
                self.write()
                with self.assertRaises(ValueError):
                    load_dataset(self.left)

    def test_zero_or_negative_duration_is_rejected(self) -> None:
        self.cmaq["records"][0]["end_utc"] = self.cmaq["records"][0]["start_utc"]
        self.write()
        with self.assertRaisesRegex(ValueError, "positive duration"):
            load_dataset(self.left)

    def test_invalid_values_are_rejected(self) -> None:
        for value in (math.nan, math.inf, True, "1", None):
            with self.subTest(value=value):
                self.cmaq["records"][0]["value"] = value
                self.write()
                with self.assertRaises(ValueError):
                    load_dataset(self.left)

    def test_signed_increments_are_not_clipped(self) -> None:
        self.cmaq["records"][0]["value"] = -1
        self.write()
        self.assertEqual(compare(self.left, self.right)["mean_bias"], -1)

    def test_large_finite_values_do_not_overflow_squared_errors(self) -> None:
        for row in self.cmaq["records"]:
            row["value"] = 1e200
        for row in self.aermod["records"]:
            row["value"] = 0
        self.write()
        rmse = compare(self.left, self.right)["root_mean_square_error"]
        assert isinstance(rmse, float)
        self.assertAlmostEqual(rmse / 1e200, 1)

    def test_unknown_schema_and_empty_records_are_rejected(self) -> None:
        for key, value in (("schema_version", 2), ("schema_version", True), ("records", [])):
            with self.subTest(key=key, value=value):
                self.cmaq = dataset("CMAQ")
                self.cmaq[key] = value
                self.write()
                with self.assertRaises(ValueError):
                    load_dataset(self.left)

    def test_empty_provenance_or_contract_field_is_rejected(self) -> None:
        self.cmaq["provenance"] = ""
        self.write()
        with self.assertRaises(ValueError):
            load_dataset(self.left)
        self.cmaq = dataset("CMAQ")
        self.cmaq["contract"]["source_id"] = ""
        self.write()
        with self.assertRaises(ValueError):
            load_dataset(self.left)

    def test_cli_reports_errors_and_nonzero_exit(self) -> None:
        with patch("sys.argv", ["cmaq_port", "compare", str(self.right), str(self.left)]):
            with patch("sys.stderr") as error:
                self.assertEqual(main(), 1)
                self.assertTrue(error.write.called)


class InventoryTests(unittest.TestCase):
    def test_skipped_files_and_windows_link_stubs_have_explicit_origins(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            driver = root / "CCTM" / "src" / "driver" / "cmaq_main.F"
            driver.parent.mkdir(parents=True)
            driver.write_text("driver", encoding="ascii")
            (root / "README.md").write_text("version", encoding="ascii")
            (root / "link.f").write_text("target.f", encoding="ascii")
            index = (
                "H 100644 " + "0" * 40 + " 0\tCCTM/src/driver/cmaq_main.F\0"
                "H 100644 " + "1" * 40 + " 0\tREADME.md\0"
                "S 100644 " + "2" * 40 + " 0\tskipped.md\0"
                "H 120000 " + "3" * 40 + " 0\tlink.f\0"
            )

            def fake_git(_path: Path, *args: str) -> str:
                if args[0] == "rev-parse":
                    return "commit\n"
                if args[0] == "status":
                    return ""
                if args[0] == "branch":
                    return "branch\n"
                return index

            with patch("cmaq_port.inventory.git", side_effect=fake_git):
                with patch("cmaq_port.inventory.git_blob", return_value=b"target.f"):
                    report = inventory(root)
            self.assertFalse(report["worktree_complete"])
            self.assertEqual(report["excluded_worktree_paths"], ["skipped.md"])
            self.assertEqual(report["tracked_link_paths"], ["link.f"])
            files = report["files"]
            assert isinstance(files, list)
            link = next(row for row in files if row["path"] == "link.f")
            self.assertEqual(link["content_origin"], "git_index_link_target")
            self.assertFalse(link["materialized_as_symlink"])
            with patch("cmaq_port.inventory.git", side_effect=fake_git):
                (root / "README.md").unlink()
                with self.assertRaisesRegex(ValueError, "Tracked file missing"):
                    inventory(root)

    def test_inventory_hashes_tracked_working_files_without_changing_source(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            driver = root / "CCTM" / "src" / "driver" / "cmaq_main.F"
            driver.parent.mkdir(parents=True)
            driver.write_text("test driver", encoding="ascii")
            (root / "README.md").write_text("test version\n", encoding="ascii")

            def fake_git(path: Path, *args: str) -> str:
                self.assertEqual(path, root)
                if args[0] == "rev-parse":
                    return "commit-test\n"
                if args[0] == "status":
                    return ""
                if args[0] == "branch":
                    return "branch-test\n"
                return (
                    "H 100644 " + "0" * 40 + " 0\tREADME.md\0"
                    "H 100644 " + "1" * 40 + " 0\tCCTM/src/driver/cmaq_main.F\0"
                )

            with patch("cmaq_port.inventory.git", side_effect=fake_git):
                report = inventory(root)
            self.assertEqual(report["tracked_file_count"], 2)
            self.assertEqual(report["commit"], "commit-test")
            self.assertEqual(report["branch"], "branch-test")
            self.assertTrue(report["clean"])
            files = report["files"]
            assert isinstance(files, list)
            self.assertEqual(files[0]["sha256"], sha256(driver))
            self.assertEqual(driver.read_text(encoding="ascii"), "test driver")

    def test_unrecognized_source_errors(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            with self.assertRaises(ValueError):
                inventory(Path(temporary))

    def test_unique_persistent_report_paths(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            first = save_report({"test": 1}, "inventory", project=Path(temporary))
            second = save_report({"test": 2}, "inventory", project=Path(temporary))
            self.assertNotEqual(first, second)
            self.assertEqual(json.loads(first.read_text(encoding="utf-8")), {"test": 1})

    def test_lesson_distinguishes_models_and_current_limits(self) -> None:
        self.assertIn("A receptor maximum is not a grid-cell mean", LESSON)
        self.assertIn("no model simulation", LESSON)
        self.assertIn("Chemically inert does not automatically mean non-depositing", LESSON)
        self.assertIn("No netCDF adapter, full CMAQ executable build", LESSON)


if __name__ == "__main__":
    unittest.main()
