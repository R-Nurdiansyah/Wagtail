#!/usr/bin/env python3
"""
wagtail.py  –  Wagtail pipeline v0.15 wrapper

Usage
-----
Run the pipeline (all unrecognised flags are forwarded to Snakemake):
    python wagtail.py --use-conda -c 8 --configfile pipeline/config.yaml

Use GreenGenes2 database variant:
    python wagtail.py --gg2 --use-conda -c 8 --configfile pipeline/config.yaml

Scan database directory for required files and validate taxonomy TSVs:
    python wagtail.py --check-db
    python wagtail.py --check-db --db-dir /custom/path/to/database

List all database files detected in the database directory:
    python wagtail.py --list-db
    python wagtail.py --list-db --db-dir /custom/path/to/database

Print Wagtail version:
    python wagtail.py --version
"""

import glob
import os
import sys

# ── paths ─────────────────────────────────────────────────────────────────────
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DB_DIR   = os.path.join(BASE_DIR, "database")
SMK_DIR  = os.path.join(BASE_DIR, "pipeline")
SMK_FILE = os.path.join(SMK_DIR, "wagtail.smk")
SMK_GG2  = os.path.join(SMK_DIR, "wagtail_gg2.smk")
BIN_DIR  = os.path.join(BASE_DIR, "bin")

VERSION  = "0.15"

# ── database specification ────────────────────────────────────────────────────
MARKERS = ["16S", "18S", "ITS", "CO1"]

# (db_type, glob_suffix, skip_for_16S, description)
DB_SPECS = [
    ("identification", "_database*",  False, "Marker identification (identify_amplicon.py)"),
    ("deblur_ref",     "_deblur*",    True,  "Deblur reference (denoise-other; built-in for 16S)"),
    ("taxonomy",       "_taxonomy*",  False, "Taxonomy reference (mappy + extract_taxonomy)"),
]


# ── helpers ───────────────────────────────────────────────────────────────────

def _find(db_dir: str, marker: str, suffix: str) -> list[str]:
    return sorted(glob.glob(os.path.join(db_dir, f"{marker}{suffix}")))


def _pop_flag(args: list[str], flag: str, has_value: bool = False):
    """Remove *flag* (and optionally its following value) from args in-place.
    Returns the value if has_value=True, else True if flag was present.
    """
    if flag not in args:
        return None
    idx = args.index(flag)
    args.pop(idx)
    if has_value:
        if idx >= len(args):
            print(f"{flag} requires an argument", file=sys.stderr)
            sys.exit(1)
        return args.pop(idx)
    return True


# ── --list-db ─────────────────────────────────────────────────────────────────

def list_databases(db_dir: str):
    """Print a formatted table of all detected database files."""
    col = 16
    print(f"\nDatabase directory: {db_dir}\n")
    print(f"  {'Marker':<8} {'Type':<{col}} {'File'}")
    print(f"  {'-'*6}   {'-'*col}   {'-'*50}")

    for marker in MARKERS:
        for db_type, suffix, skip_16s, _ in DB_SPECS:
            if marker == "16S" and skip_16s:
                print(f"  {marker:<8} {db_type:<{col}} (built-in QIIME2 reference)")
                continue
            matches = _find(db_dir, marker, suffix)
            if matches:
                for m in matches:
                    print(f"  {marker:<8} {db_type:<{col}} {os.path.basename(m)}")
            else:
                print(f"  {marker:<8} {db_type:<{col}} [NOT FOUND]")

    print()
    sys.exit(0)


# ── --check-db ────────────────────────────────────────────────────────────────

def check_databases(db_dir: str):
    """
    1. Verify that every required database file is present (exactly one match
       per pattern).
    2. Validate the taxonomy TSV files using check_db_taxonomy.py.
    Exits 0 on pass, 1 on any issue.
    """
    print(f"\nChecking databases in: {db_dir}\n")

    problems    = []
    taxonomy_ok = []   # paths of taxonomy TSVs that passed presence check

    # ── presence / uniqueness check ───────────────────────────────────────────
    for marker in MARKERS:
        marker_ok = True
        for db_type, suffix, skip_16s, description in DB_SPECS:
            if marker == "16S" and skip_16s:
                print(f"  [{marker}] {db_type:<16} OK  (built-in QIIME2 reference)")
                continue

            matches = _find(db_dir, marker, suffix)

            if len(matches) == 0:
                msg = f"[{marker}] {db_type}: no file matching '{marker}{suffix}'"
                print(f"  ✗ {msg}")
                problems.append(msg)
                marker_ok = False

            elif len(matches) > 1:
                names = [os.path.basename(m) for m in matches]
                msg = f"[{marker}] {db_type}: multiple matches (expected 1): {names}"
                print(f"  ✗ {msg}")
                problems.append(msg)
                marker_ok = False

            else:
                path = matches[0]
                size = os.path.getsize(path)
                if size == 0:
                    msg = f"[{marker}] {db_type}: file is empty: {os.path.basename(path)}"
                    print(f"  ✗ {msg}")
                    problems.append(msg)
                    marker_ok = False
                else:
                    print(f"  [{marker}] {db_type:<16} OK  {os.path.basename(path)}")
                    if db_type == "taxonomy":
                        taxonomy_ok.append(path)

    # ── taxonomy format check ─────────────────────────────────────────────────
    # Import early so any problem with check_db_taxonomy.py is caught immediately,
    # not only when TSV files happen to be present.
    sys.path.insert(0, BIN_DIR)
    from check_db_taxonomy import check_file, summarise   # noqa: E402

    # Only run on plain TSV files (not .gz / .fasta / .qza) because
    # check_db_taxonomy.py uses plain open().
    tsv_files = [p for p in taxonomy_ok
                 if not any(p.endswith(ext) for ext in (".gz", ".fasta", ".fa", ".qza"))]

    if tsv_files:
        print(f"\nValidating taxonomy format ({len(tsv_files)} file(s))...")
        for path in tsv_files:
            results = check_file(path)
            summarise(results, os.path.basename(path))
            bad = [r for r in results if not r.ok]
            if bad:
                problems.append(
                    f"Taxonomy format errors in {os.path.basename(path)} "
                    f"({len(bad)} row(s) with issues)"
                )
    elif taxonomy_ok:
        skipped = [os.path.basename(p) for p in taxonomy_ok]
        print(f"\n  (Taxonomy format check skipped for compressed/non-TSV files: "
              f"{', '.join(skipped)})")
        print("   Run bin/check_db_taxonomy.py directly to check these files.")

    # ── summary ───────────────────────────────────────────────────────────────
    print(f"\n{'='*70}")
    if problems:
        print(f"  ISSUES FOUND — {len(problems)} problem(s):")
        for p in problems:
            print(f"    • {p}")
    else:
        print("  PASS — all databases present and taxonomy formats valid")
    print(f"{'='*70}\n")

    sys.exit(1 if problems else 0)


# ── pipeline runner ───────────────────────────────────────────────────────────

def run_pipeline(args: list[str], use_gg2: bool = False):
    if not any(a in ("-s", "--snakefile") for a in args):
        smk = SMK_GG2 if use_gg2 else SMK_FILE
        if not os.path.isfile(smk):
            print(f"Snakefile not found: {smk}", file=sys.stderr)
            sys.exit(1)
        args = ["-s", smk] + args

    sys.argv = ["snakemake"] + args
    import runpy
    try:
        runpy.run_module("snakemake", run_name="__main__")
    except Exception as e:
        print("Could not run Snakemake CLI. Check your Snakemake installation.",
              file=sys.stderr)
        print(str(e), file=sys.stderr)
        sys.exit(1)


# ── entry point ───────────────────────────────────────────────────────────────

def main():
    args = list(sys.argv[1:])

    # --version
    if "--version" in args:
        print(f"Wagtail {VERSION}")
        sys.exit(0)

    # --db-dir (shared by --check-db and --list-db; consumed before dispatch)
    db_dir = _pop_flag(args, "--db-dir", has_value=True) or DB_DIR

    # --check-db
    if _pop_flag(args, "--check-db"):
        check_databases(db_dir)
        return   # check_databases calls sys.exit

    # --list-db
    if _pop_flag(args, "--list-db"):
        list_databases(db_dir)
        return   # list_databases calls sys.exit

    # --gg2 (pass-through to pipeline; remove the flag itself before snakemake)
    use_gg2 = bool(_pop_flag(args, "--gg2"))
    # Also strip any legacy aliases that were accepted in v0.14
    for alias in ("--greengenes", "--greengenes2", "--green_genes"):
        if _pop_flag(args, alias):
            use_gg2 = True

    run_pipeline(args, use_gg2=use_gg2)


if __name__ == "__main__":
    main()
