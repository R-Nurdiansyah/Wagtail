import sys
import sqlite3
import os

db_path = sys.argv[1]
id_ = sys.argv[2]
rule = sys.argv[3]
outcome = sys.argv[4]
stdout_path = sys.argv[5]
stderr_path = sys.argv[6]

stdout = ""
stderr = ""

if os.path.exists(stdout_path):
    with open(stdout_path, "r") as f:
        stdout = f.read()
if os.path.exists(stderr_path):
    with open(stderr_path, "r") as f:
        stderr = f.read()

conn = sqlite3.connect(db_path)
c = conn.cursor()
c.execute("""
CREATE TABLE IF NOT EXISTS logs (
    id TEXT,
    rule TEXT,
    outcome TEXT,
    stdout TEXT,
    stderr TEXT
)
""")
c.execute("INSERT INTO logs (id, rule, outcome, stdout, stderr) VALUES (?, ?, ?, ?, ?)",
          (id_, rule, outcome, stdout, stderr))
conn.commit()
conn.close()
