#!/usr/bin/env python3
"""
identify_marker.py
Subsample reads and competitively align against 16S, ITS, and CO1 databases
using mappy. Writes the predicted marker name to --output.

Exit codes:
    0  marker identified (check --output for the name)
    1  error during alignment or confidence below threshold
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
    p.add_argument("--reads",           help="Input FASTQ(.gz). Mutually exclusive with --file-map/--sample.")
    p.add_argument("--file-map",        help="Headerless 2-column TSV: accession<TAB>file_path")
    p.add_argument("--sample",          help="Sample accession to look up in --file-map")
    p.add_argument("--db-16S",          required=True, help="16S FASTA(.gz) for 16S")
    p.add_argument("--db-ITS",          required=True, help="ITS FASTA(.gz) for ITS")
    p.add_argument("--db-CO1",          required=True, help="CO1 FASTA(.gz) for CO1")
    p.add_argument("--db-18S",          required=True, help="18S FASTA(.gz) for 18S")
    p.add_argument("--output",          required=True, help="Path to write predicted marker name")
    p.add_argument("--n-subsample",     type=int, default=1000,
                   help="Number of reads to subsample (default: 1000)")
    p.add_argument("--min-confidence",  type=float, default=0.5,
                   help="Min combined score (hit_fraction × mean_identity) for a confident call "
                        "(default: 0.5)")
    p.add_argument("--preset",          default="sr",
                   help="Minimap2 preset: 'sr' for short reads")
    p.add_argument("--seed",            type=int, default=42,
                   help="Random seed for subsampling (default: 42)")
    p.add_argument("--min-identity",   type=float, default=0.50,
                   help="Minimum sequence identity for a read to count as a hit (default: 0.50)")
    p.add_argument("--min-aln-length", type=int, default=25,
                   help="Minimum alignment length for a read to count as a hit (default: 25)")
    p.add_argument("--min-margin",     type=float, default=0.05,
                   help="Minimum combined-score margin between best and second-best marker "
                        "for an unambiguous call (default: 0.05)")
    p.add_argument("--min-16s-fraction", type=float, default=0.90,
                   help="If 16S hit fraction meets this threshold, always call 16S "
                        "(overrides ranking; guards against 16S/18S cross-homology). "
                        "Default: 0.90")
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
    identities:      list[float]  # per-hit identities (one per read, best hit only)

# ── Helpers ───────────────────────────────────────────────────────────────────
def reads_from_filemap(filemap_path: str, sample_name: str) -> str:
    """Return the file path for sample_name from a headerless TSV filemap."""
    with open(filemap_path, newline="") as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if len(row) >= 2 and row[0].strip() == sample_name:
                return row[1].strip()
    raise ValueError(f"Sample '{sample_name}' not found in filemap '{filemap_path}'")

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

def align_to_db(
    reads:          list[str],
    db_path:        str,
    marker:         str,
    preset:         str,
    min_identity:   float,
    min_aln_length: int,
) -> AlignmentStats:
    """
    Align reads against one database and return per-marker alignment stats.
    For each read, only the best primary hit is recorded.
    Identity = mlen / blen (matching bases / alignment block length).
    """
    empty = AlignmentStats(
        marker=marker, n_reads=len(reads), n_hits=0,
        hit_fraction=0.0, mean_identity=0.0, median_identity=0.0, identities=[]
    )
 
    if not Path(db_path).exists():
        logging.warning(f"  [{marker}] Database not found, skipping: {db_path}")
        return empty
 
    try:
        aligner = mp.Aligner(db_path, preset=preset, best_n=2)
    except Exception as exc:
        logging.warning(f"  [{marker}] Could not load {db_path}: {exc}")
        return empty
 
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
 
    # Free the minimap2 index from C-level memory before returning.
    # Python's GC will not release mmap'd index buffers promptly without this.
    del aligner

    n_hits = len(identities)
    if n_hits == 0:
        return empty

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
        identities      = identities,
    )

def identify_marker(
    stats:            dict[str, AlignmentStats],
    min_confidence:   float,
    min_margin:       float,
    min_16s_fraction: float,
) -> tuple[str, str, AlignmentStats, AlignmentStats | None]:
    """
    Primary metric:   combined score = hit_fraction × mean_identity
    Tiebreaker:       hit_fraction

    Rationale: mean_identity alone is misleading when two markers have
    different numbers of hits — you'd be comparing means over different read
    populations.  The combined score is analogous to BLAST (query_coverage ×
    identity): a marker that aligns 60 % of reads at 94.8 % identity scores
    0.569, clearly beating one that aligns only 13 % at 96.5 % (0.126).

    Priority rule: if 16S hit fraction >= min_16s_fraction, call 16S
    unconditionally.  16S and 18S share enough SSU homology that 16S reads
    routinely align to 18S at nearly identical identity; a strong 16S hit
    fraction is a more reliable discriminator.

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
    # If 16S hit fraction is high, short-circuit to 16S regardless of ranking.
    # This guards against the 16S/18S cross-homology scenario where 18S edges
    # out 16S by a fraction of a percent.
    sixteen_s = stats.get("16S")
    if sixteen_s and sixteen_s.hit_fraction >= min_16s_fraction:
        non_16s = next((s for s in ranked if s.marker != "16S"), None)
        logging.info(
            f"  16S priority rule triggered: "
            f"16S hit_frac={sixteen_s.hit_fraction:.3f} >= {min_16s_fraction}"
        )
        return "OK", "16S", sixteen_s, non_16s

    # ── Standard ranking gates ────────────────────────────────────────────────
    best_score   = best.hit_fraction * best.mean_identity
    second_score = (second.hit_fraction * second.mean_identity) if second else 0.0
    margin       = round(best_score - second_score, 4)

    # Gate 1: best combined score must exceed min_confidence
    if best_score < min_confidence:
        return "UNKNOWN", best.marker, best, second

    # Gate 2: best combined score must be clearly better than second
    if margin < min_margin:
        return "AMBIGUOUS", best.marker, best, second

    return "OK", best.marker, best, second

# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [identify_marker] %(levelname)s %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S"
    )
    args = parse_args()

    # Resolve reads path: either --reads directly or --file-map + --sample lookup
    if args.reads and args.file_map:
        logging.error("--reads and --file-map/--sample are mutually exclusive")
        sys.exit(1)
    if args.reads:
        reads_path = args.reads
    elif args.file_map and args.sample:
        try:
            reads_path = reads_from_filemap(args.file_map, args.sample)
        except (ValueError, OSError) as exc:
            logging.error(f"Filemap lookup failed: {exc}")
            sys.exit(1)
    else:
        logging.error("Provide either --reads or both --file-map and --sample")
        sys.exit(1)

    # 1. Subsample
    logging.info(f"Subsampling {args.n_subsample} reads from {reads_path}")
    reads = subsample_reads(reads_path, args.n_subsample)
    if not reads:
        logging.error("No reads could be loaded from input file.")
        sys.exit(1)
    logging.info(f"  → {len(reads)} reads loaded")

    # 2. Align against each database.
    # 16S is listed first so the early-exit check below can fire before the
    # larger CO1/18S databases are loaded into memory.
    databases = [
        ("16S", args.db_16S),
        ("ITS", args.db_ITS),
        ("18S", args.db_18S),
        ("CO1", args.db_CO1),   # largest DB last — skipped for most 16S samples
    ]
    stats: dict[str, AlignmentStats] = {}
    for marker, db_path in databases:
        logging.info(f"Aligning against {marker}: {db_path}")
        result = align_to_db(
            reads, db_path, marker, args.preset,
            args.min_identity, args.min_aln_length,
        )
        gc.collect()   # reclaim C-level memory freed by del aligner inside align_to_db
        stats[marker] = result
        logging.info(
            f"  → {marker}: hits={result.n_hits}/{result.n_reads} "
            f"({result.hit_fraction:.3f})  "
            f"mean_id={result.mean_identity:.4f}  "
            f"median_id={result.median_identity:.4f}"
        )

        # Early-exit: if 16S priority rule already fires, skip remaining databases.
        # For a typical 16S sample this avoids loading the CO1 MIDORI2 index
        # (~1.5 GB uncompressed → several GB RAM for the minimap2 index).
        if marker == "16S" and result.hit_fraction >= args.min_16s_fraction:
            logging.info(
                f"  16S fraction {result.hit_fraction:.3f} >= {args.min_16s_fraction} "
                f"— skipping remaining databases"
            )
            break

    # 3. Identify marker
    status, predicted, best, second = identify_marker(
        stats, args.min_confidence, args.min_margin, args.min_16s_fraction
    )
 
    logging.info(f"Status:    {status}")
    logging.info(f"Predicted: {predicted}")

    best_score   = best.hit_fraction * best.mean_identity
    second_score = (second.hit_fraction * second.mean_identity) if second else 0.0

    logging.info(
        f"  → Best:        {best.marker}  "
        f"score={best_score:.4f}  "
        f"mean_id={best.mean_identity:.4f}  "
        f"hit_frac={best.hit_fraction:.3f}"
    )
    if second:
        logging.info(
            f"  → Second-best: {second.marker}  "
            f"score={second_score:.4f}  "
            f"mean_id={second.mean_identity:.4f}  "
            f"hit_frac={second.hit_fraction:.3f}"
        )

    if status == "AMBIGUOUS":
        logging.warning(
            f"Marker ambiguous between best ({best.marker} "
            f"score={best_score:.4f}) and second-best "
            f"({second.marker} score={second_score:.4f}); "
            f"score margin {best_score - second_score:.4f} "
            f"is below threshold {args.min_margin}"
        )
    elif status == "UNKNOWN":
        logging.warning(
            f"Amplicon is not 16S, 18S, ITS, or CO1 — "
            f"best combined score {best_score:.4f} "
            f"(hit_frac={best.hit_fraction:.3f} × mean_id={best.mean_identity:.4f}) "
            f"is below confidence threshold {args.min_confidence}"
        )
 
    # 4. Write TSV — 8 columns
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
 
    second_marker    = second.marker                        if second else "None"
    second_mean_id   = f"{second.mean_identity:.4f}"       if second else "NA"
    second_hit_frac  = f"{second.hit_fraction:.3f}"        if second else "NA"
    second_score_str = f"{second_score:.4f}"               if second else "NA"
    score_margin     = best_score - second_score

    with open(args.output, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow([
            "status",
            "predicted_marker",
            "best_combined_score",
            "best_mean_identity",
            "best_hit_fraction",
            "second_marker",
            "second_combined_score",
            "second_mean_identity",
            "second_hit_fraction",
            "score_margin",
        ])
        writer.writerow([
            status,
            predicted,
            f"{best_score:.4f}",
            f"{best.mean_identity:.4f}",
            f"{best.hit_fraction:.3f}",
            second_marker,
            second_score_str,
            second_mean_id,
            second_hit_frac,
            f"{score_margin:.4f}",
        ])
 
    logging.info(f"Written to {args.output}")
 
    # Exit 1 on non-OK so Snakemake shell block catches it via $status
    if status != "OK":
        sys.exit(1)


if __name__ == "__main__":
    main()
