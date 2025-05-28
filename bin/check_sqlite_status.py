import sys
import sqlite3
import time

def connect_with_retry(db_path, retries=300, delay=2):
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

def execute_with_retry(cursor, query, params=(), retries=300, delay=2):
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

if len(sys.argv) != 4:
    print("NONE")
    sys.exit(0)

db_path, filename, rule = sys.argv[1:4]
try:
    conn = connect_with_retry(db_path)
    c = conn.cursor()
    execute_with_retry(
        c,
        'SELECT status FROM logs WHERE filename=? AND rule=? ORDER BY id DESC LIMIT 1',
        (filename, rule)
    )
    row = c.fetchone()
    if row:
        print(row[0])
    else:
        print("NONE")
except Exception:
    print("NONE")
finally:
    try:
        conn.close()
    except Exception:
        pass
