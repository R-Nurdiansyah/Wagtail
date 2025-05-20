import sys
import sqlite3

if len(sys.argv) != 4:
    print("NONE")
    sys.exit(0)

db_path, filename, rule = sys.argv[1:4]
try:
    conn = sqlite3.connect(db_path)
    c = conn.cursor()
    c.execute('SELECT status FROM logs WHERE filename=? AND rule=? ORDER BY id DESC LIMIT 1', (filename, rule))
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
