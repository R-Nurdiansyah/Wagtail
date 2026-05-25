#!/usr/bin/env python3
"""
Taxonomy QC checker for multiple database formats.

Checks:
  1. Number of taxonomic levels
  2. Proper prefix notation (d__/k__, p__, c__, o__, f__, g__, s__, optionally sh__)
  3. No whitespace except immediately after semicolons
  4. Formatting consistency: feature_ID <TAB> taxonomy, 
     starts with d__/k__, no trailing semicolon
"""

import re
import sys
from pathlib import Path
from dataclasses import dataclass, field
from typing import Optional

# ── constants ────────────────────────────────────────────────────────────────

# Accepted rank prefixes in canonical order (stop at s__; sh__ is optional extra)
RANK_ORDER = ["d", "k", "p", "c", "o", "f", "g", "s"]
OPTIONAL_EXTRA_RANKS = {"sh"}          # UNITE's species-hypothesis rank
ROOT_PREFIXES = {"d", "k"}            # both are valid root ranks

# Regex: a valid level looks like  prefix__anything  (no whitespace in name)
LEVEL_RE = re.compile(r'^([a-z]+)__(.*)$')

# ── data classes ─────────────────────────────────────────────────────────────

@dataclass
class LevelIssue:
    level_index: int
    raw: str
    problems: list[str]

@dataclass
class RowResult:
    line_number: int
    feature_id: str
    taxonomy_raw: str
    levels: list[str]
    n_levels: int
    issues: list[str] = field(default_factory=list)
    level_issues: list[LevelIssue] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return not self.issues and not self.level_issues

# ── parsing helpers ───────────────────────────────────────────────────────────

def split_taxonomy(tax: str) -> list[str]:
    """Split on semicolons, stripping only trailing/leading whitespace from each token."""
    return [t.strip() for t in tax.split(";") if t.strip()]


def detect_format(raw_tax: str) -> str:
    """Heuristic: which DB format is this line's taxonomy field?"""
    if re.search(r'\bd__', raw_tax) or re.search(r'\bk__', raw_tax):
        return "prefix_semicolon"          # QIIME-style with __ prefixes
    # plain semicolon-separated names (9-level PR²-style)
    return "plain_semicolon"


def check_whitespace(tax: str) -> list[str]:
    """
    Whitespace rule: the only allowed whitespace is a single space
    immediately after a semicolon (between levels), not inside a level name.
    """
    problems = []
    # Check for whitespace INSIDE level tokens (after splitting on ';')
    for tok in split_taxonomy(tax):
        if re.search(r'\s', tok):
            problems.append(f"Whitespace inside level token: {tok!r}")
    # Check for whitespace BEFORE a semicolon
    if re.search(r'\s;', tax):
        problems.append("Whitespace before semicolon")
    return problems


def check_prefix_row(levels: list[str]) -> tuple[list[str], list[LevelIssue]]:
    """
    Validate prefix-style taxonomy (d__/k__ … s__).
    Returns (row-level issues, per-level issues).
    """
    row_issues = []
    lv_issues  = []

    # --- root rank ---
    if not levels:
        return ["Empty taxonomy"], []

    first_match = LEVEL_RE.match(levels[0])
    if not first_match:
        row_issues.append(f"First level has no __ prefix: {levels[0]!r}")
        return row_issues, lv_issues
    root_prefix = first_match.group(1)
    if root_prefix not in ROOT_PREFIXES:
        row_issues.append(
            f"Root prefix must be d__ or k__, got {root_prefix!r}: {levels[0]!r}"
        )

    # Build expected prefix sequence starting from root
    if root_prefix == "d":
        expected = ["d", "p", "c", "o", "f", "g", "s"]
    else:  # k
        expected = ["k", "p", "c", "o", "f", "g", "s"]

    # Separate core levels from optional extras (sh__)
    core_levels  = [lv for lv in levels if not (LEVEL_RE.match(lv) and
                    LEVEL_RE.match(lv).group(1) in OPTIONAL_EXTRA_RANKS)]
    extra_levels = [lv for lv in levels if (LEVEL_RE.match(lv) and
                    LEVEL_RE.match(lv).group(1) in OPTIONAL_EXTRA_RANKS)]

    # --- level count (core only) ---
    n_core = len(core_levels)
    if n_core != len(expected):
        row_issues.append(
            f"Expected {len(expected)} core levels ({'/'.join(expected)}__), "
            f"found {n_core}"
        )

    # --- per-level prefix order ---
    for i, lv in enumerate(core_levels):
        probs = []
        m = LEVEL_RE.match(lv)
        if not m:
            probs.append(f"No __ prefix found")
        else:
            prefix = m.group(1)
            if i < len(expected) and prefix != expected[i]:
                probs.append(
                    f"Expected prefix {expected[i]}__, got {prefix}__"
                )
        if probs:
            lv_issues.append(LevelIssue(i, lv, probs))

    # --- optional sh__ must come after s__ ---
    if extra_levels:
        last_core_prefix = LEVEL_RE.match(core_levels[-1]).group(1) if core_levels and LEVEL_RE.match(core_levels[-1]) else None
        if last_core_prefix != "s":
            row_issues.append("sh__ level present but s__ is not the last core level")

    return row_issues, lv_issues


def check_formatting(feature_id: str, tax_raw: str, levels: list[str]) -> list[str]:
    """Row-level formatting checks (point 4 in requirements)."""
    issues = []

    # Must be two-column TSV (checked at parse time, but flag if ID empty)
    if not feature_id.strip():
        issues.append("Empty feature ID")

    # Starts with d__ or k__
    if levels and LEVEL_RE.match(levels[0]):
        root_p = LEVEL_RE.match(levels[0]).group(1)
        if root_p not in ROOT_PREFIXES:
            issues.append(f"Taxonomy must start with d__ or k__, starts with {root_p}__")
    elif levels:
        # plain format — no prefix at all
        issues.append("Taxonomy does not use __ prefix notation")

    # Trailing semicolon
    stripped = tax_raw.rstrip()
    if stripped.endswith(";"):
        issues.append("Trailing semicolon present")

    return issues

# ── main row checker ──────────────────────────────────────────────────────────

def check_row(line_number: int, line: str) -> Optional[RowResult]:
    """Parse and check one TSV line. Returns None for blank/header lines."""
    line = line.rstrip("\n")
    if not line.strip() or line.startswith("#"):
        return None

    parts = line.split("\t")
    if len(parts) < 2:
        # single-column — can't check
        return RowResult(
            line_number=line_number,
            feature_id=parts[0] if parts else "",
            taxonomy_raw="",
            levels=[],
            n_levels=0,
            issues=["Row has only one column (no TAB-separated taxonomy field)"],
        )

    feature_id = parts[0]
    # Some files have a header row "Feature ID\tTaxon"
    if feature_id.strip().lower() in {"feature id", "feature_id", "#otuid", "otu id"}:
        return None

    taxonomy_raw = parts[1]
    fmt = detect_format(taxonomy_raw)

    # Split into levels
    levels = split_taxonomy(taxonomy_raw)
    n_levels = len(levels)

    # Whitespace check (applies to all formats)
    ws_issues = check_whitespace(taxonomy_raw)

    # Prefix / order checks
    if fmt == "prefix_semicolon":
        row_issues, lv_issues = check_prefix_row(levels)
    else:
        # plain format: we can only count levels and note missing prefixes
        row_issues = ["Plain (no __ prefix) format — cannot validate rank order"]
        lv_issues  = []

    # Formatting check
    fmt_issues = check_formatting(feature_id, taxonomy_raw, levels)

    all_row_issues = ws_issues + row_issues + fmt_issues

    return RowResult(
        line_number=line_number,
        feature_id=feature_id,
        taxonomy_raw=taxonomy_raw,
        levels=levels,
        n_levels=n_levels,
        issues=all_row_issues,
        level_issues=lv_issues,
    )

# ── file-level processing ─────────────────────────────────────────────────────

def check_file(path: str) -> list[RowResult]:
    results = []
    with open(path, encoding="utf-8") as fh:
        for ln, line in enumerate(fh, start=1):
            r = check_row(ln, line)
            if r is not None:
                results.append(r)
    return results


def summarise(results: list[RowResult], label: str):
    print(f"\n{'='*70}")
    print(f"  {label}")
    print(f"{'='*70}")
    if not results:
        print("  (no rows found)")
        return

    from collections import Counter, defaultdict

    ok_rows  = [r for r in results if r.ok]
    bad_rows = [r for r in results if not r.ok]
    level_counts = Counter(r.n_levels for r in results)
    levels_str = ", ".join(f"{n} levels" for n in sorted(level_counts))

    print(f"  Rows: {len(results)}  |  Levels: {levels_str}")

    if not bad_rows:
        print(f"  ✓ PASS — all {len(ok_rows)} rows OK")
        return

    print(f"  ✗ FAIL — {len(bad_rows)}/{len(results)} rows have issues\n")

    # Group problems by category (strip variable parts) so e.g. whitespace per
    # level doesn't print 6 times — keep one representative example per category.
    def problem_category(prob: str) -> str:
        if prob.startswith("Whitespace inside level token"):
            return "Whitespace inside level token(s)"
        if prob.startswith("Expected prefix"):
            return "Wrong rank prefix order"
        return prob

    seen_categories: dict[str, tuple[str, RowResult]] = {}
    for r in bad_rows:
        all_problems = r.issues + [p for li in r.level_issues for p in li.problems]
        for prob in all_problems:
            cat = problem_category(prob)
            if cat not in seen_categories:
                seen_categories[cat] = (prob, r)

    for cat, (prob, r) in seen_categories.items():
        print(f"  • {cat}")
        print(f"    e.g. line {r.line_number}: {r.feature_id}")
        print(f"         {r.taxonomy_raw[:100]}{'...' if len(r.taxonomy_raw) > 100 else ''}")

# ── inline data for the four databases shown ─────────────────────────────────

INLINE_DATABASES = {
    "control (7-level, d__ prefix)": """\
FLASV1.1346;tax=d:Bacteria,p:Actinobacteriota,c:Actinobacteria,o:Propionibacteriales,f:Propionibacteriaceae,g:Cutibacterium,s:Cutibacterium_acnes;\td__Bacteria; p__Actinobacteriota; c__Actinobacteria; o__Propionibacteriales; f__Propionibacteriaceae; g__Cutibacterium; s__Cutibacterium_acnes
FLASV2.1310;tax=d:Bacteria,p:Proteobacteria,c:Alphaproteobacteria,o:Rhizobiales,f:Xanthobacteraceae,g:Bradyrhizobium,s:MFD_s_2;\td__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rhizobiales; f__Xanthobacteraceae; g__Bradyrhizobium; s__MFD_s_2
FLASV3.1372;tax=d:Bacteria,p:Firmicutes,c:Bacilli,o:Bacillales,f:Bacillaceae,g:Bacillus,s:MFD_s_3;\td__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Bacillaceae; g__Bacillus; s__MFD_s_3
""",

    "9-level plain semicolon (PR2-style)": """\
AB353770.1.1740_U\tEukaryota;TSAR;Alveolata;Dinoflagellata;Dinophyceae;Peridiniales;Kryptoperidiniaceae;Unruhdinium;Unruhdinium_kevei;
AB284159.1.1765_U\tEukaryota;TSAR;Alveolata;Dinoflagellata;Dinophyceae;Peridiniales;Protoperidiniaceae;Protoperidinium;Protoperidinium_bipes;
FJ355953.1.1907_U\tEukaryota;Obazoa;Opisthokonta;Fungi;Ascomycota;Pezizomycotina;Eurotiomycetes;Knufia;Knufia_epidermidis;
""",

    "8-level UNITE (k__ + sh__)": """\
Feature ID\tTaxon
SH1227328.10FU_MT153946_refs\tk__Fungi;p__Ascomycota;c__Dothideomycetes;o__Abrothallales;f__Abrothallaceae;g__Abrothallus;s__Abrothallus_subhalei;sh__SH1227328.10FU
SH1227742.10FU_JN206177_refs\tk__Fungi;p__Mucoromycota;c__Mucoromycetes;o__Mucorales;f__Mucoraceae;g__Mucor;s__Mucor_inaequisporus;sh__SH1227742.10FU
SH1232203.10FU_KY102517_refs\tk__Fungi;p__Ascomycota;c__Saccharomycetes;o__Saccharomycetales;f__Saccharomycetales_fam_Incertae_sedis;g__Candida;s__Candida_vrieseae;sh__SH1232203.10FU
""",

    "7-level eukaryote (k__ prefix)": """\
MH535936.1.<1.>890\tk__Eukaryota_2759;p__phylum_class_order_family_genus_Amoebozoa sp._1892891;c__class_order_family_genus_Amoebozoa sp._1892891;o__order_family_genus_Amoebozoa sp._1892891;f__family_genus_Amoebozoa sp._1892891;g__genus_Amoebozoa sp._1892891;s__Amoebozoa sp._1892891
MH535935.1.<1.>889\tk__Eukaryota_2759;p__phylum_class_order_family_genus_Amoebozoa sp._1892891;c__class_order_family_genus_Amoebozoa sp._1892891;o__order_family_genus_Amoebozoa sp._1892891;f__family_genus_Amoebozoa sp._1892891;g__genus_Amoebozoa sp._1892891;s__Amoebozoa sp._1892891
MG559732.1.<1.>690\tk__Eukaryota_2759;p__Discosea_555280;c__Flabellinia_1485085;o__order_Vannellidae_95227;f__Vannellidae_95227;g__Clydonella_218657;s__Clydonella sawyeri_2201168
""",
}

# ── entry point ───────────────────────────────────────────────────────────────

import io

def main():
    if len(sys.argv) > 1:
        # File mode: check_taxonomy.py file1.tsv [file2.tsv ...]
        for path in sys.argv[1:]:
            results = check_file(path)
            summarise(results, Path(path).name)
    else:
        # Inline demo mode
        print("Taxonomy QC Checker — inline demo mode")
        print("(Pass TSV file paths as arguments to check real files)\n")
        for label, data in INLINE_DATABASES.items():
            results = []
            for ln, line in enumerate(io.StringIO(data), start=1):
                r = check_row(ln, line)
                if r is not None:
                    results.append(r)
            summarise(results, label)

    print(f"\n{'='*70}\nDone.\n")


if __name__ == "__main__":
    main()