#!/usr/bin/env python3
"""Unit tests for bin/check_sqlite_status.py"""

import os
import sqlite3
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = PROJECT_ROOT / "bin" / "check_sqlite_status.py"


class TestCheckSqliteStatus(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.db_path = os.path.join(self.temp_dir.name, "logs.db")
        conn = sqlite3.connect(self.db_path)
        conn.execute(
            """
            CREATE TABLE logs (
                id TEXT,
                rule TEXT,
                outcome TEXT,
                log TEXT,
                note TEXT
            )
            """
        )
        conn.execute(
            "INSERT INTO logs (id, rule, outcome, log, note) VALUES (?, ?, ?, ?, ?)",
            ("sample1", "ruleA", "SUCCESS", "log", "note"),
        )
        conn.commit()
        conn.close()

    def tearDown(self):
        self.temp_dir.cleanup()

    def _run_script(self, filename, rule):
        result = subprocess.run(
            [sys.executable, str(SCRIPT_PATH), self.db_path, filename, rule],
            check=True,
            capture_output=True,
            text=True,
        )
        return result.stdout.strip()

    def test_returns_logged_outcome(self):
        self.assertEqual(self._run_script("sample1", "ruleA"), "SUCCESS")

    def test_returns_none_when_entry_missing(self):
        self.assertEqual(self._run_script("missing", "ruleA"), "NONE")


if __name__ == "__main__":
    unittest.main()
