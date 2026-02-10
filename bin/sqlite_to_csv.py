import sys
import sqlite3
import csv

if len(sys.argv) != 3:
    print("Usage: sqlite_to_csv.py <input.db> <output.csv>")
    sys.exit(1)

db_path, csv_path = sys.argv[1:3]

conn = sqlite3.connect(db_path)
cur = conn.cursor()
cur.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='logs'")
if not cur.fetchone():
    print("No 'logs' table found in database.")
    sys.exit(1)

cur.execute("SELECT * FROM logs")
rows = cur.fetchall()
colnames = [desc[0] for desc in cur.description]

with open(csv_path, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(colnames)
    writer.writerows(rows)

conn.close()
