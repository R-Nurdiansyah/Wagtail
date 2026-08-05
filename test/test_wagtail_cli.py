#!/usr/bin/env python3
"""
Unit tests for wagtail.py argument parsing (v1.2 refactor).

Covers parse_options() / _pop_archive_intermediate() / Options.config_overrides:
  - --amplicon canonicalisation + marker restriction
  - the hyphenated --archive-result / --archive-intermediate / --archive-all flags
  - leftover args forwarded to Snakemake
  - --config override assembly

wagtail.py is loaded directly from its file path so this never picks up an
unrelated installed `wagtail` package.
"""

import importlib.util
import os
import sys
import unittest

_ROOT = os.path.join(os.path.dirname(os.path.realpath(__file__)), "..")
_spec = importlib.util.spec_from_file_location(
    "wagtail_cli", os.path.join(_ROOT, "wagtail.py"))
wagtail = importlib.util.module_from_spec(_spec)
# Register before exec_module: @dataclass resolves cls.__module__ via
# sys.modules, which would be None for an unregistered importlib module.
sys.modules[_spec.name] = wagtail
_spec.loader.exec_module(wagtail)


class TestAmplicon(unittest.TestCase):

    def test_amplicon_restricts_markers_and_forwards_rest(self):
        opts = wagtail.parse_options(["--amplicon", "ITS", "-c", "8"])
        self.assertEqual(opts.amplicon, "ITS")
        self.assertEqual(opts.active_markers, ["ITS"])
        self.assertEqual(opts.snakemake_args, ["-c", "8"])

    def test_amplicon_is_canonicalised_to_upper(self):
        self.assertEqual(wagtail.parse_options(["--amplicon", "its"]).amplicon, "ITS")

    def test_no_amplicon_means_all_markers(self):
        opts = wagtail.parse_options(["-c", "4"])
        self.assertIsNone(opts.amplicon)
        self.assertEqual(opts.active_markers, wagtail.MARKERS)

    def test_invalid_amplicon_exits(self):
        with self.assertRaises(SystemExit):
            wagtail.parse_options(["--amplicon", "XYZ"])


class TestArchiveFlags(unittest.TestCase):

    def test_result_and_all_flags(self):
        opts = wagtail.parse_options(["--archive-result", "--archive-all"])
        self.assertTrue(opts.archive_result)
        self.assertTrue(opts.archive_all)
        self.assertEqual(opts.snakemake_args, [])

    def test_intermediate_all(self):
        opts = wagtail.parse_options(["--archive-intermediate", "all"])
        self.assertEqual(opts.archive_intermediate, "all")

    def test_intermediate_space_separated_list_is_sorted_unique(self):
        opts = wagtail.parse_options(["--archive-intermediate", "3", "1", "5"])
        self.assertEqual(opts.archive_intermediate, "1,3,5")

    def test_intermediate_comma_separated_list(self):
        opts = wagtail.parse_options(["--archive-intermediate", "2,4"])
        self.assertEqual(opts.archive_intermediate, "2,4")

    def test_intermediate_stops_at_next_flag(self):
        opts = wagtail.parse_options(["--archive-intermediate", "1", "2", "-c", "8"])
        self.assertEqual(opts.archive_intermediate, "1,2")
        self.assertEqual(opts.snakemake_args, ["-c", "8"])

    def test_intermediate_out_of_range_exits(self):
        with self.assertRaises(SystemExit):
            wagtail.parse_options(["--archive-intermediate", "9"])

    def test_intermediate_missing_value_exits(self):
        with self.assertRaises(SystemExit):
            wagtail.parse_options(["--archive-intermediate"])


class TestConfigOverrides(unittest.TestCase):

    def test_single_config_block_with_all_overrides(self):
        opts = wagtail.parse_options(
            ["--amplicon", "16S", "--archive-result", "--archive-intermediate", "1,2"])
        overrides = opts.config_overrides()
        self.assertEqual(overrides[0], "--config")
        self.assertIn("amplicon=16S", overrides)
        self.assertIn("archive_result=true", overrides)          # config key stays underscore
        self.assertIn("archive_intermediate=1,2", overrides)

    def test_no_overrides_is_empty(self):
        opts = wagtail.parse_options(["-c", "4"])
        self.assertEqual(opts.config_overrides(), [])


class TestMiscFlags(unittest.TestCase):

    def test_gg2_flag_and_legacy_alias(self):
        self.assertTrue(wagtail.parse_options(["--gg2"]).use_gg2)
        self.assertTrue(wagtail.parse_options(["--greengenes2"]).use_gg2)

    def test_db_dir_override_consumed(self):
        opts = wagtail.parse_options(["--db-dir", "/custom/db", "-c", "2"])
        self.assertEqual(opts.db_dir, "/custom/db")
        self.assertEqual(opts.snakemake_args, ["-c", "2"])

    def test_check_and_list_db_flags(self):
        self.assertTrue(wagtail.parse_options(["--check-db"]).check_db)
        self.assertTrue(wagtail.parse_options(["--list-db"]).list_db)


if __name__ == "__main__":
    unittest.main()
