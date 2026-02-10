#!/usr/bin/env python3
"""Unit tests for sqlite-related helper scripts."""

import csv
import os
import sqlite3
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)), '..')] + sys.path

from bin.sql_merge import merge_sqlite_logs

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = PROJECT_ROOT / 'test' / 'test_data'
SQLITE_TO_CSV = PROJECT_ROOT / 'bin' / 'sqlite_to_csv.py'
SQLITE_LOGGER = PROJECT_ROOT / 'bin' / 'wagtail_sqlite_logger.py'


class TestSqliteTools(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.temp_dir.cleanup()

    def _create_db_with_logs(self, path, rows):
        conn = sqlite3.connect(path)
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
        conn.executemany("INSERT INTO logs VALUES (?, ?, ?, ?, ?)", rows)
        conn.commit()
        conn.close()

    def test_merge_sqlite_logs_combines_unique_entries(self):
        db1 = os.path.join(self.temp_dir.name, 'one.db')
        db2 = os.path.join(self.temp_dir.name, 'two.db')
        merged = os.path.join(self.temp_dir.name, 'merged.db')
        self._create_db_with_logs(db1, [("s1", "rule", "OK", "log1", "note1")])
        self._create_db_with_logs(db2, [("s2", "rule", "FAIL", "log2", "note2")])

        merge_sqlite_logs(merged, [db1, db2])

        conn = sqlite3.connect(merged)
        rows = conn.execute('SELECT id, outcome FROM logs ORDER BY id').fetchall()
        conn.close()
        self.assertEqual(rows, [("s1", "OK"), ("s2", "FAIL")])

    def test_sqlite_to_csv_exports_rows(self):
        db_path = os.path.join(self.temp_dir.name, 'export.db')
        self._create_db_with_logs(db_path, [("s3", "rule", "OK", "log", "note")])
        csv_path = os.path.join(self.temp_dir.name, 'logs.csv')

        subprocess.run(
            [sys.executable, str(SQLITE_TO_CSV), db_path, csv_path],
            check=True,
            capture_output=True,
            text=True,
        )

        with open(csv_path, newline='') as handle:
            reader = csv.reader(handle)
            header = next(reader)
            row = next(reader)

        self.assertIn('id', header)
        self.assertEqual(row[0], 's3')

    def test_wagtail_sqlite_logger_appends_rows(self):
        db_path = os.path.join(self.temp_dir.name, 'logger.db')
        log_path = DATA_DIR / 'log_generic.txt'

        subprocess.run(
            [
                sys.executable,
                str(SQLITE_LOGGER),
                db_path,
                'sampleX',
                'ruleY',
                'PASS',
                str(log_path),
                'note',
            ],
            check=True,
            capture_output=True,
            text=True,
        )

        conn = sqlite3.connect(db_path)
        stored = conn.execute('SELECT id, rule, outcome, note FROM logs').fetchone()
        conn.close()
        self.assertEqual(stored, ('sampleX', 'ruleY', 'PASS', 'note'))


if __name__ == '__main__':
    unittest.main()
