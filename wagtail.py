#!/usr/bin/env python3
"""
wagtail.py  –  Wagtail pipeline v1.2 wrapper

A thin launcher around the Snakemake workflow. It consumes a small set of
Wagtail-specific flags and forwards everything else, unchanged, to Snakemake.

Usage
-----
Run the pipeline (all unrecognised flags are forwarded to Snakemake):
    python wagtail.py --use-conda -c 8 --configfile pipeline/config.yaml

Run for a single amplicon only (16S, 18S, ITS, or CO1):
    python wagtail.py --amplicon 16S --use-conda -c 8 --configfile pipeline/config.yaml
    python wagtail.py --amplicon ITS --use-conda -c 8 --configfile pipeline/config.yaml

Use GreenGenes2 database variant:
    python wagtail.py --gg2 --use-conda -c 8 --configfile pipeline/config.yaml

Archive run outputs during cleanup (any combination):
    # final results only (0_/6_/7_) → <run>_final_results.zip
    python wagtail.py --archive-result --use-conda -c 8 --configfile pipeline/config.yaml
    # selected intermediates (1-5) → <run>_intermediate.zip
    python wagtail.py --archive-intermediate all --use-conda -c 8 ...
    python wagtail.py --archive-intermediate 1,3,5 --use-conda -c 8 ...
    python wagtail.py --archive-intermediate 2 4   --use-conda -c 8 ...
    # everything (0_..7_) → <run>_all_results.zip
    python wagtail.py --archive-all --use-conda -c 8 --configfile pipeline/config.yaml

Scan database directory for required files and validate taxonomy TSVs:
    python wagtail.py --check-db
    python wagtail.py --check-db --amplicon ITS
    python wagtail.py --check-db --db-dir /custom/path/to/database

List all database files detected in the database directory:
    python wagtail.py --list-db
    python wagtail.py --list-db --amplicon CO1
    python wagtail.py --list-db --db-dir /custom/path/to/database

Print Wagtail version:
    python wagtail.py --version
"""

from __future__ import annotations

import glob
import os
import sys
from dataclasses import dataclass
from typing import NoReturn

# ── paths ─────────────────────────────────────────────────────────────────────
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DB_DIR   = os.path.join(BASE_DIR, "database")
SMK_DIR  = os.path.join(BASE_DIR, "pipeline")
SMK_FILE = os.path.join(SMK_DIR, "wagtail.smk")
SMK_GG2  = os.path.join(SMK_DIR, "wagtail_gg2.smk")
BIN_DIR  = os.path.join(BASE_DIR, "bin")

VERSION = "1.2"

# ── database specification ────────────────────────────────────────────────────
MARKERS = ["16S", "18S", "ITS", "CO1"]

# (db_type, glob_suffix, skip_for_16S, description)
DB_SPECS = [
    ("identification", "_database*", False, "Marker identification (identify_amplicon.py)"),
    ("deblur_ref",     "_deblur*",   True,  "Deblur reference (denoise-other; built-in for 16S)"),
    ("taxonomy",       "_taxonomy*",  False, "Taxonomy reference (mappy + extract_taxonomy)"),
]

# Legacy aliases that have always meant --gg2.
_GG2_ALIASES = ("--gg2", "--greengenes", "--greengenes2", "--green_genes")


def _die(msg: str) -> NoReturn:
    """Print a usage error to stderr and exit 1 (argparse-style)."""
    print(msg, file=sys.stderr)
    sys.exit(1)


# ── low-level flag consumption ────────────────────────────────────────────────
# These mutate the args list in place, popping Wagtail-specific flags so that
# whatever remains can be forwarded verbatim to Snakemake.

def _pop_flag(args: list[str], flag: str, has_value: bool = False):
    """Remove *flag* (and optionally its following value) from args in-place.
    Returns the value if has_value=True, else True if the flag was present.
    """
    if flag not in args:
        return None
    idx = args.index(flag)
    args.pop(idx)
    if has_value:
        if idx >= len(args):
            _die(f"{flag} requires an argument")
        return args.pop(idx)
    return True


def _pop_any_flag(args: list[str], flags) -> bool:
    """Pop every occurrence of any flag in *flags*; return True if any matched."""
    found = False
    for flag in flags:
        if _pop_flag(args, flag):
            found = True
    return found


def _is_archive_token(tok: str) -> bool:
    """True if *tok* is a valid --archive-intermediate value token:
    'all', or a (comma-separated) group of digits."""
    if tok == "all":
        return True
    parts = [p for p in tok.split(",") if p]
    return bool(parts) and all(p.isdigit() for p in parts)


def _pop_archive_intermediate(args: list[str]) -> str | None:
    """Remove --archive-intermediate and its following value token(s) in-place.

    Accepts 'all' or integers 1-5, given space- and/or comma-separated, e.g.
        --archive-intermediate all
        --archive-intermediate 1,3,5
        --archive-intermediate 1 3 5
    Returns None if absent, 'all', or a comma-joined sorted-unique digit string.
    """
    flag = "--archive-intermediate"
    if flag not in args:
        return None
    idx = args.index(flag)
    args.pop(idx)

    raw = []
    while idx < len(args) and _is_archive_token(args[idx]):
        raw.append(args.pop(idx))
    if not raw:
        _die(f"{flag} requires 'all' or integers 1-5")

    parts = [p for tok in raw for p in tok.split(",") if p]
    if "all" in parts:
        return "all"

    nums = sorted({int(p) for p in parts})
    for n in nums:
        if not 1 <= n <= 5:
            _die(f"{flag}: {n} out of range (choose 1-5 or 'all')")
    return ",".join(str(n) for n in nums)


# ── parsed options ────────────────────────────────────────────────────────────

@dataclass
class Options:
    """Wagtail-specific options parsed from the command line, plus the
    leftover arguments to forward to Snakemake."""
    db_dir: str
    amplicon: str | None              # canonical marker, or None for auto-detect
    active_markers: list[str]         # markers the db commands operate on
    check_db: bool
    list_db: bool
    use_gg2: bool
    archive_result: bool
    archive_all: bool
    archive_intermediate: str | None  # None | 'all' | '1,3,5'
    snakemake_args: list[str]

    def config_overrides(self) -> list[str]:
        """`--config KEY=VALUE ...` entries to inject for the Snakemake run.

        Emitted as a single --config block so later entries can't clobber
        earlier ones. The .smk reads these via config.get(...).
        """
        overrides = []
        if self.amplicon:
            overrides.append(f"amplicon={self.amplicon}")
        if self.archive_result:
            overrides.append("archive_result=true")
        if self.archive_all:
            overrides.append("archive_all=true")
        if self.archive_intermediate:
            overrides.append(f"archive_intermediate={self.archive_intermediate}")
        return ["--config", *overrides] if overrides else []


def parse_options(argv: list[str]) -> Options:
    """Consume every Wagtail-specific flag from *argv* and return an Options.
    Whatever is left in the working list becomes Options.snakemake_args."""
    args = list(argv)

    db_dir = _pop_flag(args, "--db-dir", has_value=True) or DB_DIR

    amplicon = _parse_amplicon(args)
    active_markers = [amplicon] if amplicon else MARKERS

    return Options(
        db_dir               = db_dir,
        amplicon             = amplicon,
        active_markers       = active_markers,
        check_db             = bool(_pop_flag(args, "--check-db")),
        list_db              = bool(_pop_flag(args, "--list-db")),
        use_gg2              = _pop_any_flag(args, _GG2_ALIASES),
        archive_result       = bool(_pop_flag(args, "--archive-result")),
        archive_all          = bool(_pop_flag(args, "--archive-all")),
        archive_intermediate = _pop_archive_intermediate(args),
        snakemake_args       = args,
    )


def _parse_amplicon(args: list[str]) -> str | None:
    """Pop --amplicon and return its canonical (upper-case) marker, or None."""
    raw = _pop_flag(args, "--amplicon", has_value=True)
    if raw is None:
        return None
    amplicon = raw.upper()
    if amplicon not in MARKERS:
        _die(f"Unknown amplicon {raw!r}. Choose from: {', '.join(MARKERS)}")
    return amplicon


# ── database commands ─────────────────────────────────────────────────────────

def _find(db_dir: str, marker: str, suffix: str) -> list[str]:
    return sorted(glob.glob(os.path.join(db_dir, f"{marker}{suffix}")))


def _amplicon_filter_note(markers: list[str]) -> None:
    if markers != MARKERS:
        print(f"  Amplicon filter: {', '.join(markers)}\n")


def list_databases(db_dir: str, markers: list[str]) -> int:
    """Print a formatted table of all detected database files. Returns 0."""
    col = 16
    print(f"\nDatabase directory: {db_dir}\n")
    _amplicon_filter_note(markers)
    print(f"  {'Marker':<8} {'Type':<{col}} {'File'}")
    print(f"  {'-'*6}   {'-'*col}   {'-'*50}")

    for marker in markers:
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
    return 0


def _check_db_file(db_dir: str, marker: str, db_type: str, suffix: str):
    """Presence/uniqueness/size check for one (marker, db_type).
    Returns (ok, message, taxonomy_path_or_None)."""
    matches = _find(db_dir, marker, suffix)
    if not matches:
        return False, f"[{marker}] {db_type}: no file matching '{marker}{suffix}'", None
    if len(matches) > 1:
        names = [os.path.basename(m) for m in matches]
        return False, f"[{marker}] {db_type}: multiple matches (expected 1): {names}", None
    path = matches[0]
    if os.path.getsize(path) == 0:
        return False, f"[{marker}] {db_type}: file is empty: {os.path.basename(path)}", None
    tax = path if db_type == "taxonomy" else None
    return True, f"[{marker}] {db_type:<16} OK  {os.path.basename(path)}", tax


def check_databases(db_dir: str, markers: list[str]) -> int:
    """
    1. Verify that every required database file is present (exactly one match
       per pattern) and non-empty.
    2. Validate the taxonomy TSV files using check_db_taxonomy.py.
    Returns 0 on pass, 1 on any issue.
    """
    print(f"\nChecking databases in: {db_dir}\n")
    _amplicon_filter_note(markers)

    problems: list[str] = []
    taxonomy_ok: list[str] = []   # taxonomy TSVs that passed the presence check

    # ── presence / uniqueness check ───────────────────────────────────────────
    for marker in markers:
        for db_type, suffix, skip_16s, _ in DB_SPECS:
            if marker == "16S" and skip_16s:
                print(f"  [{marker}] {db_type:<16} OK  (built-in QIIME2 reference)")
                continue
            ok, message, tax = _check_db_file(db_dir, marker, db_type, suffix)
            print(f"  {message}" if ok else f"  ✗ {message}")
            if not ok:
                problems.append(message)
            elif tax:
                taxonomy_ok.append(tax)

    # ── taxonomy format check ─────────────────────────────────────────────────
    problems.extend(_check_taxonomy_formats(taxonomy_ok))

    # ── summary ───────────────────────────────────────────────────────────────
    print(f"\n{'='*70}")
    if problems:
        print(f"  ISSUES FOUND — {len(problems)} problem(s):")
        for p in problems:
            print(f"    • {p}")
    else:
        print("  PASS — all databases present and taxonomy formats valid")
    print(f"{'='*70}\n")

    return 1 if problems else 0


def _check_taxonomy_formats(taxonomy_ok: list[str]) -> list[str]:
    """Validate the plain-TSV taxonomy files; return a list of problem strings."""
    # Import here (but unconditionally) so a problem with check_db_taxonomy.py is
    # caught even when no TSV files happen to be present.
    sys.path.insert(0, BIN_DIR)
    from check_db_taxonomy import check_file, summarise   # noqa: E402

    # check_db_taxonomy.py uses plain open(), so only validate uncompressed TSVs.
    tsv_files = [p for p in taxonomy_ok
                 if not any(p.endswith(ext) for ext in (".gz", ".fasta", ".fa", ".qza"))]

    if not tsv_files:
        if taxonomy_ok:
            skipped = ", ".join(os.path.basename(p) for p in taxonomy_ok)
            print(f"\n  (Taxonomy format check skipped for compressed/non-TSV files: {skipped})")
            print("   Run bin/check_db_taxonomy.py directly to check these files.")
        return []

    print(f"\nValidating taxonomy format ({len(tsv_files)} file(s))...")
    problems = []
    for path in tsv_files:
        results = check_file(path)
        summarise(results, os.path.basename(path))
        bad = [r for r in results if not r.ok]
        if bad:
            problems.append(
                f"Taxonomy format errors in {os.path.basename(path)} "
                f"({len(bad)} row(s) with issues)"
            )
    return problems


# ── pipeline runner ───────────────────────────────────────────────────────────

def run_pipeline(snakemake_args: list[str], use_gg2: bool = False) -> None:
    """Invoke the Snakemake CLI in-process with the chosen Snakefile."""
    args = list(snakemake_args)
    if not any(a in ("-s", "--snakefile") for a in args):
        smk = SMK_GG2 if use_gg2 else SMK_FILE
        if not os.path.isfile(smk):
            _die(f"Snakefile not found: {smk}")
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

def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)

    if "--version" in argv:
        print(f"Wagtail {VERSION}")
        return 0

    opts = parse_options(argv)

    # Database-inspection commands short-circuit the pipeline.
    if opts.check_db:
        return check_databases(opts.db_dir, opts.active_markers)
    if opts.list_db:
        return list_databases(opts.db_dir, opts.active_markers)

    # Forward to Snakemake, injecting the config overrides built from our flags.
    run_pipeline(opts.snakemake_args + opts.config_overrides(), use_gg2=opts.use_gg2)
    return 0


if __name__ == "__main__":
    sys.exit(main())
