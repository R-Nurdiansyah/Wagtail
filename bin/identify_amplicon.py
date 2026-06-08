#!/usr/bin/env python3
"""
identify_marker.py
Subsample reads and competitively align against 16S, 18S, ITS, and CO1
databases using mappy. Writes the predicted marker name to a .marker TSV.

Two modes:
  Single sample:  --reads FILE  (or --file-map + --sample)  --output FILE
  Batch:          --sample-list FILE  --file-map FILE  --output-dir DIR
                  Loads each database ONCE and loops over every sample, so the
                  multi-GB minimap2 index is built a single time per batch
                  instead of once per sample. Each sample's result is written
                  to <output-dir>/<sample>.marker.

Exit codes:
    0  ran to completion (per-sample status is in each .marker file)
    1  fatal error (bad arguments, unreadable sample list, etc.)
       In single-sample mode, also returns 1 when the call is not "OK"
       (AMBIGUOUS / UNKNOWN), so the Snakemake shell can catch it via $status.
"""

from __future__ import annotations

import argparse
import gc
import gzip
import logging
import sys
from pathlib import Path
from dataclasses import dataclass
import random
import csv

import mappy as mp


# ── CLI ───────────────────────────────────────────────────────────────────────
def parse_args():
    p = argparse.ArgumentParser(description="Identify amplicon marker gene from raw reads.")
    # Single-sample inputs
    p.add_argument("--reads",           help="Input FASTQ(.gz). Mutually exclusive with --file-map/--sample.")
    p.add_argument("--sample",          help="Sample accession to look up in --file-map (single-sample mode)")
    p.add_argument("--output",          help="Path to write the predicted-marker TSV (single-sample mode)")
    # Batch inputs
    p.add_argument("--sample-list",     help="File with one accession per line. Triggers batch mode "
                                             "(requires --file-map and --output-dir).")
    p.add_argument("--output-dir",      help="Directory to write <sample>.marker files (batch mode)")
    # Shared
    p.add_argument("--file-map",        help="Headerless 2-column TSV: accession<TAB>file_path")
    p.add_argument("--db-16S",          required=True, help="16S FASTA(.gz)/.mmi for 16S")
    p.add_argument("--db-ITS",          required=True, help="ITS FASTA(.gz)/.mmi for ITS")
    p.add_argument("--db-CO1",          required=True, help="CO1 FASTA(.gz)/.mmi for CO1")
    p.add_argument("--db-18S",          required=True, help="18S FASTA(.gz)/.mmi for 18S")
    p.add_argument("--n-subsample",     type=int, default=1000,
                   help="Number of reads to subsample (default: 1000)")
    p.add_argument("--min-confidence",  type=float, default=0.5,
                   help="Min combined score (hit_fraction × mean_identity) for a confident call "
                        "(default: 0.5)")
    p.add_argument("--preset",          default="sr",
                   help="Minimap2 preset: 'sr' for short reads")
    p.add_argument("--minimizer-w",    type=int, default=19,
                   help="Minimap2 minimizer window. Larger = smaller index / less RAM, "
                        "slightly lower sensitivity. Ignored when loading a prebuilt .mmi "
                        "(its baked-in w wins). Default: 19")
    p.add_argument("--seed",            type=int, default=42,
                   help="Random seed for subsampling (default: 42)")
    p.add_argument("--min-identity",   type=float, default=0.50,
                   help="Minimum sequence identity for a read to count as a hit (default: 0.50)")
    p.add_argument("--min-aln-length", type=int, default=25,
                   help="Minimum alignment length for a read to count as a hit (default: 25)")
    p.add_argument("--min-margin",     type=float, default=0.05,
                   help="Minimum combined-score margin between best and second-best marker "
                        "for an unambiguous call (default: 0.05)")
    p.add_argument("--min-16s-score", type=float, default=0.80,
                   help="If the 16S combined score (hit_fraction × mean_identity) meets "
                        "this threshold, always call 16S (overrides ranking; guards "
                        "against 16S/18S cross-homology). Default: 0.80")
    return p.parse_args()

# ── Data container ────────────────────────────────────────────────────────────
@dataclass
class AlignmentStats:
    marker:          str
    n_reads:         int          # total reads attempted
    n_hits:          int          # reads with at least one valid hit
    hit_fraction:    float        # n_hits / n_reads
    mean_identity:   float        # mean % identity across all valid hits
    median_identity: float        # median % identity across all valid hits

def make_empty(marker: str, n_reads: int) -> AlignmentStats:
    """A zero-hit stats record (database missing, unloadable, or no alignments)."""
    return AlignmentStats(
        marker=marker, n_reads=n_reads, n_hits=0,
        hit_fraction=0.0, mean_identity=0.0, median_identity=0.0
    )

# ── Helpers ───────────────────────────────────────────────────────────────────
def reads_from_filemap(filemap_path: str, sample_name: str) -> str:
    """Return the file path for sample_name from a headerless TSV filemap."""
    with open(filemap_path, newline="") as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if len(row) >= 2 and row[0].strip() == sample_name:
                return row[1].strip()
    raise ValueError(f"Sample '{sample_name}' not found in filemap '{filemap_path}'")

def load_filemap(filemap_path: str) -> dict[str, str]:
    """Load the whole filemap into a dict once — avoids re-scanning per sample."""
    mapping: dict[str, str] = {}
    with open(filemap_path, newline="") as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if len(row) >= 2 and row[0].strip():
                mapping[row[0].strip()] = row[1].strip()
    return mapping

def subsample_reads(fastq_path: str, n: int, seed: int = 42) -> list[str]:
    """
    Reservoir-sample n reads in a single streaming pass — O(n) peak memory.
    Handles arbitrarily large FASTQ files without loading them into RAM first.
    """
    rng = random.Random(seed)
    reservoir: list[str] = []
    opener = gzip.open if fastq_path.endswith(".gz") else open
    try:
        with opener(fastq_path, "rt") as fh:
            i = 0
            while True:
                header = fh.readline()
                if not header:
                    break
                seq = fh.readline().strip()
                fh.readline()   # +
                fh.readline()   # qual
                if not seq:
                    continue
                i += 1
                if len(reservoir) < n:
                    reservoir.append(seq)
                else:
                    j = rng.randint(0, i - 1)
                    if j < n:
                        reservoir[j] = seq
    except Exception as exc:
        logging.error(f"Failed to read {fastq_path}: {exc}")
        return []
    return reservoir

# ── Alignment ─────────────────────────────────────────────────────────────────
def build_aligner(db_path: str, marker: str, preset: str, w: int):
    """Load a minimap2 index. Returns an mp.Aligner or None on failure.

    A prebuilt .mmi already encodes its own k/w, so minimap2 ignores the w
    passed here for those files — build the .mmi with the desired w instead.
    """
    if not Path(db_path).exists():
        logging.warning(f"  [{marker}] Database not found, skipping: {db_path}")
        return None
    try:
        return mp.Aligner(db_path, preset=preset, best_n=1, w=w)
    except Exception as exc:
        logging.warning(f"  [{marker}] Could not load {db_path}: {exc}")
        return None

def align_reads(
    reads:          list[str],
    aligner,
    marker:         str,
    min_identity:   float,
    min_aln_length: int,
) -> AlignmentStats:
    """Align one sample's reads against an already-loaded aligner.

    Records only the best primary hit per read.
    Identity = mlen / blen (matching bases / alignment block length).
    Does NOT load or free the aligner — the caller owns its lifecycle, which
    is what lets batch mode reuse one index across many samples.
    """
    identities = []
    for seq in reads:
        for hit in aligner.map(seq):
            if not hit.is_primary:
                continue
            identity = hit.mlen / hit.blen
            aln_len  = hit.q_en - hit.q_st
            if identity >= min_identity and aln_len >= min_aln_length:
                identities.append(identity)
            break   # one best hit per read only

    n_hits = len(identities)
    if n_hits == 0:
        return make_empty(marker, len(reads))

    sorted_ids = sorted(identities)
    mid = len(sorted_ids) // 2
    median = (
        sorted_ids[mid]
        if len(sorted_ids) % 2 == 1
        else (sorted_ids[mid - 1] + sorted_ids[mid]) / 2
    )
    return AlignmentStats(
        marker          = marker,
        n_reads         = len(reads),
        n_hits          = n_hits,
        hit_fraction    = n_hits / len(reads),
        mean_identity   = sum(identities) / n_hits,
        median_identity = median,
    )

def align_to_db(
    reads:          list[str],
    db_path:        str,
    marker:         str,
    preset:         str,
    min_identity:   float,
    min_aln_length: int,
    w:              int = 19,
) -> AlignmentStats:
    """Single-sample convenience wrapper: load index, align, free index."""
    aligner = build_aligner(db_path, marker, preset, w)
    if aligner is None:
        return make_empty(marker, len(reads))
    try:
        return align_reads(reads, aligner, marker, min_identity, min_aln_length)
    finally:
        # Free the minimap2 index from C-level memory before the next DB loads.
        del aligner
        gc.collect()

# ── Decision ──────────────────────────────────────────────────────────────────
def identify_marker(
    stats:            dict[str, AlignmentStats],
    min_confidence:   float,
    min_margin:       float,
    min_16s_score:    float,
) -> tuple[str, str, AlignmentStats, AlignmentStats | None]:
    """
    Primary metric:   combined score = hit_fraction × mean_identity
    Tiebreaker:       hit_fraction

    Rationale: mean_identity alone is misleading when two markers have
    different numbers of hits — you'd be comparing means over different read
    populations.  The combined score is analogous to BLAST (query_coverage ×
    identity): a marker that aligns 60 % of reads at 94.8 % identity scores
    0.569, clearly beating one that aligns only 13 % at 96.5 % (0.126).

    Priority rule: if the 16S combined score (hit_fraction × mean_identity)
    >= min_16s_score, call 16S unconditionally.  16S and 18S share enough
    SSU homology that 16S reads routinely align to 18S at nearly identical
    identity; a strong 16S combined score is a more reliable discriminator.

    Returns: (status, predicted_marker, best_stats, second_stats)
      status: "OK" | "AMBIGUOUS" | "UNKNOWN"
    """
    # Sort by combined score descending; break ties by hit fraction
    ranked = sorted(
        stats.values(),
        key=lambda s: (s.hit_fraction * s.mean_identity, s.hit_fraction),
        reverse=True,
    )
    best   = ranked[0]
    second = ranked[1] if len(ranked) > 1 else None

    # ── 16S priority rule ────────────────────────────────────────────────────
    sixteen_s = stats.get("16S")
    if sixteen_s and (sixteen_s.hit_fraction * sixteen_s.mean_identity) >= min_16s_score:
        sixteen_s_score = sixteen_s.hit_fraction * sixteen_s.mean_identity
        non_16s = next((s for s in ranked if s.marker != "16S"), None)
        logging.info(
            f"  16S priority rule triggered: "
            f"16S combined_score={sixteen_s_score:.4f} >= {min_16s_score}"
        )
        return "OK", "16S", sixteen_s, non_16s

    # ── Standard ranking gates ────────────────────────────────────────────────
    best_score   = best.hit_fraction * best.mean_identity
    second_score = (second.hit_fraction * second.mean_identity) if second else 0.0
    margin       = round(best_score - second_score, 4)

    if best_score < min_confidence:
        return "UNKNOWN", best.marker, best, second
    if margin < min_margin:
        return "AMBIGUOUS", best.marker, best, second
    return "OK", best.marker, best, second

# ── Output ────────────────────────────────────────────────────────────────────
MARKER_HEADER = [
    "status", "predicted_marker",
    "best_combined_score", "best_mean_identity", "best_hit_fraction",
    "second_marker", "second_combined_score", "second_mean_identity",
    "second_hit_fraction", "score_margin",
]

def write_marker_file(
    output_path: str | Path,
    status:      str,
    predicted:   str,
    best:        AlignmentStats,
    second:      AlignmentStats | None,
) -> None:
    """Write the 10-column marker TSV for one sample."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    best_score    = best.hit_fraction * best.mean_identity
    second_score  = (second.hit_fraction * second.mean_identity) if second else 0.0
    score_margin  = best_score - second_score

    with open(output_path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(MARKER_HEADER)
        writer.writerow([
            status,
            predicted,
            f"{best_score:.4f}",
            f"{best.mean_identity:.4f}",
            f"{best.hit_fraction:.3f}",
            second.marker                  if second else "None",
            f"{second_score:.4f}"          if second else "NA",
            f"{second.mean_identity:.4f}"  if second else "NA",
            f"{second.hit_fraction:.3f}"   if second else "NA",
            f"{score_margin:.4f}",
        ])

# ── Single-sample mode ────────────────────────────────────────────────────────
def run_single(args) -> int:
    # Resolve reads path: either --reads directly or --file-map + --sample lookup
    if args.reads and args.file_map:
        logging.error("--reads and --file-map/--sample are mutually exclusive")
        return 1
    if args.reads:
        reads_path = args.reads
    elif args.file_map and args.sample:
        try:
            reads_path = reads_from_filemap(args.file_map, args.sample)
        except (ValueError, OSError) as exc:
            logging.error(f"Filemap lookup failed: {exc}")
            return 1
    else:
        logging.error("Provide either --reads or both --file-map and --sample")
        return 1

    if not args.output:
        logging.error("--output is required in single-sample mode")
        return 1

    # 1. Subsample
    logging.info(f"Subsampling {args.n_subsample} reads from {reads_path}")
    reads = subsample_reads(reads_path, args.n_subsample, args.seed)
    if not reads:
        logging.error("No reads could be loaded from input file.")
        return 1
    logging.info(f"  → {len(reads)} reads loaded")

    # 2. Align against each database (16S first so the early-exit can fire
    #    before the larger CO1/18S indices are loaded).
    databases = [
        ("16S", args.db_16S),
        ("ITS", args.db_ITS),
        ("18S", args.db_18S),
        ("CO1", args.db_CO1),
    ]
    stats: dict[str, AlignmentStats] = {}
    for marker, db_path in databases:
        logging.info(f"Aligning against {marker}: {db_path}")
        result = align_to_db(
            reads, db_path, marker, args.preset,
            args.min_identity, args.min_aln_length, args.minimizer_w,
        )
        stats[marker] = result
        logging.info(
            f"  → {marker}: hits={result.n_hits}/{result.n_reads} "
            f"({result.hit_fraction:.3f})  mean_id={result.mean_identity:.4f}  "
            f"median_id={result.median_identity:.4f}"
        )
        if marker == "16S" and (result.hit_fraction * result.mean_identity) >= args.min_16s_score:
            score = result.hit_fraction * result.mean_identity
            logging.info(f"  16S combined score {score:.4f} >= {args.min_16s_score} "
                         f"— skipping remaining databases")
            break

    # 3. Identify + log + write
    status, predicted, best, second = identify_marker(
        stats, args.min_confidence, args.min_margin, args.min_16s_score
    )
    _log_decision(status, predicted, best, second, args)
    write_marker_file(args.output, status, predicted, best, second)
    logging.info(f"Written to {args.output}")

    return 0 if status == "OK" else 1

def _log_decision(status, predicted, best, second, args) -> None:
    logging.info(f"Status:    {status}")
    logging.info(f"Predicted: {predicted}")
    best_score   = best.hit_fraction * best.mean_identity
    second_score = (second.hit_fraction * second.mean_identity) if second else 0.0
    logging.info(f"  → Best:        {best.marker}  score={best_score:.4f}  "
                 f"mean_id={best.mean_identity:.4f}  hit_frac={best.hit_fraction:.3f}")
    if second:
        logging.info(f"  → Second-best: {second.marker}  score={second_score:.4f}  "
                     f"mean_id={second.mean_identity:.4f}  hit_frac={second.hit_fraction:.3f}")
    if status == "AMBIGUOUS":
        logging.warning(
            f"Marker ambiguous between best ({best.marker} score={best_score:.4f}) "
            f"and second-best ({second.marker} score={second_score:.4f}); "
            f"score margin {best_score - second_score:.4f} below threshold {args.min_margin}"
        )
    elif status == "UNKNOWN":
        logging.warning(
            f"Amplicon is not 16S, 18S, ITS, or CO1 — best combined score "
            f"{best_score:.4f} (hit_frac={best.hit_fraction:.3f} × "
            f"mean_id={best.mean_identity:.4f}) below threshold {args.min_confidence}"
        )

# ── Batch mode ────────────────────────────────────────────────────────────────
def run_batch(args) -> int:
    if not args.file_map:
        logging.error("--file-map is required in batch mode")
        return 1
    if not args.output_dir:
        logging.error("--output-dir is required in batch mode")
        return 1

    try:
        samples = [ln.strip() for ln in open(args.sample_list) if ln.strip()]
    except OSError as exc:
        logging.error(f"Could not read --sample-list {args.sample_list}: {exc}")
        return 1
    if not samples:
        logging.error("Sample list is empty.")
        return 1

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    filemap = load_filemap(args.file_map)
    logging.info(f"Batch mode: {len(samples)} samples → {out_dir}")

    # 1. Subsample reads for every sample up front. Failed lookups/empty files
    #    are recorded and given an empty .marker (mirrors the old touch-on-fail).
    sample_reads: dict[str, list[str]] = {}
    failed: dict[str, str] = {}
    for s in samples:
        path = filemap.get(s)
        if path is None:
            failed[s] = "not found in filemap"
            continue
        reads = subsample_reads(path, args.n_subsample, args.seed)
        if not reads:
            failed[s] = "no reads loaded"
            continue
        sample_reads[s] = reads
    logging.info(f"  Subsampled OK: {len(sample_reads)}  |  failed: {len(failed)}")

    # 2. Align against each database, loading each index exactly once.
    #    Samples that trigger the 16S priority rule are 'resolved' and skip the
    #    remaining (larger) databases — and if every sample resolves on 16S,
    #    ITS/18S/CO1 are never even loaded.
    databases = [
        ("16S", args.db_16S),
        ("ITS", args.db_ITS),
        ("18S", args.db_18S),
        ("CO1", args.db_CO1),
    ]
    sample_stats: dict[str, dict[str, AlignmentStats]] = {s: {} for s in sample_reads}
    resolved: set[str] = set()

    for marker, db_path in databases:
        todo = (list(sample_reads.keys()) if marker == "16S"
                else [s for s in sample_reads if s not in resolved])
        if not todo:
            logging.info(f"All samples resolved before {marker}; skipping its index load.")
            continue

        logging.info(f"Loading {marker} index ({db_path}) — aligning {len(todo)} sample(s)")
        aligner = build_aligner(db_path, marker, args.preset, args.minimizer_w)
        if aligner is None:
            for s in todo:
                sample_stats[s][marker] = make_empty(marker, len(sample_reads[s]))
            continue

        for n, s in enumerate(todo, 1):
            st = align_reads(sample_reads[s], aligner, marker,
                             args.min_identity, args.min_aln_length)
            sample_stats[s][marker] = st
            if marker == "16S" and (st.hit_fraction * st.mean_identity) >= args.min_16s_score:
                resolved.add(s)
            if n % 500 == 0:
                logging.info(f"    {marker}: {n}/{len(todo)} aligned")

        # Free this index before the next (possibly larger) one loads.
        del aligner
        gc.collect()
        logging.info(f"  {marker} done — resolved so far: {len(resolved)}/{len(sample_reads)}")

    # 3. Decide and write one .marker per sample.
    counts: dict[str, int] = {}
    for s in sample_reads:
        status, predicted, best, second = identify_marker(
            sample_stats[s], args.min_confidence, args.min_margin, args.min_16s_score
        )
        write_marker_file(out_dir / f"{s}.marker", status, predicted, best, second)
        counts[status] = counts.get(status, 0) + 1

    # 4. Failed samples get an empty marker so the DAG can still build (the
    #    downstream SQLite check soft-fails them).
    for s, reason in failed.items():
        (out_dir / f"{s}.marker").touch()
        logging.warning(f"  {s}: {reason} — wrote empty marker")

    logging.info(f"Batch complete. Calls: {counts} | failed: {len(failed)}")
    return 0

# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [identify_marker] %(levelname)s %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    args = parse_args()
    if args.sample_list:
        sys.exit(run_batch(args))
    sys.exit(run_single(args))


if __name__ == "__main__":
    main()
