import sys
import sqlite3
import os

def merge_sqlite_logs(merged_db, db_files):
    # Remove merged_db if exists
    if os.path.exists(merged_db):
        os.remove(merged_db)

    # Create merged_db and logs table with UNIQUE constraint
    conn = sqlite3.connect(merged_db)
    conn.execute("""
        CREATE TABLE IF NOT EXISTS logs (
            id TEXT,
            rule TEXT,
            outcome TEXT,
            log TEXT,
            note TEXT,
            UNIQUE(id, rule, outcome, log)
        );
    """)
    conn.commit()

    for db in db_files:
        if not os.path.exists(db) or os.path.getsize(db) == 0 or db == merged_db:
            continue
        print(f"Merging from {db} ...")
        try:
            # Check if logs table exists
            c = sqlite3.connect(db)
            count = c.execute("SELECT count(*) FROM sqlite_master WHERE type='table' AND name='logs';").fetchone()[0]
            c.close()
            if count == 1:
                conn.execute("ATTACH DATABASE ? AS to_merge;", (db,))
                conn.execute("BEGIN;")
                conn.execute("INSERT OR REPLACE INTO logs SELECT * FROM to_merge.logs;")
                conn.execute("COMMIT;")
                conn.execute("DETACH DATABASE to_merge;")
        except Exception as e:
            print(f"Error merging {db}: {e}")
    conn.close()
    print(f"Merge complete. Combined log stored in {merged_db}.")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        #print("Usage: python sql_merge.py merged_db.db db1.db db2.db ...")
        sys.exit(1)
    merged_db = sys.argv[1]
    db_files = sys.argv[2:]
    merge_sqlite_logs(merged_db, db_files)
