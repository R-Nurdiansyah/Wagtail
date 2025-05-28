import sys
import sqlite3
import os
import time

def connect_with_retry(db_path, retries=300, delay=1):
    for i in range(retries):
        try:
            conn = sqlite3.connect(db_path)
            conn.execute("PRAGMA journal_mode=WAL;")
            return conn
        except sqlite3.OperationalError as e:
            if "database is locked" in str(e):
                time.sleep(delay)
            else:
                raise
    raise RuntimeError("Could not connect to database after retries")

def execute_with_retry(cursor, query, params=(), retries=300, delay=1):
    for i in range(retries):
        try:
            cursor.execute(query, params)
            return
        except sqlite3.OperationalError as e:
            if "database is locked" in str(e):
                time.sleep(delay)
            else:
                raise
    raise RuntimeError("Could not execute query after retries")

db_path = sys.argv[1]
id_ = sys.argv[2]
rule = sys.argv[3]
outcome = sys.argv[4]
log_path = sys.argv[5]  # now only one log file

log = ""
if os.path.exists(log_path):
    with open(log_path, "r") as f:
        log = f.read()

conn = connect_with_retry(db_path)
c = conn.cursor()
execute_with_retry(
    c,
    """
    CREATE TABLE IF NOT EXISTS logs (
        id TEXT,
        rule TEXT,
        outcome TEXT,
        log TEXT
    )
    """
)
execute_with_retry(
    c,
    "INSERT INTO logs (id, rule, outcome, log) VALUES (?, ?, ?, ?)",
    (id_, rule, outcome, log)
)
conn.commit()
conn.close()
