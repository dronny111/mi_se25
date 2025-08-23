#!/usr/bin/env python3
"""
Thin wrapper around the Bowtie2 command‑line tool.

*   The user supplies a Bowtie2 index that points to the tomato genome
    (e.g. the ITAG4.0 FASTA you downloaded from SGN or Ensembl).
*   The wrapper creates a temporary FASTA that contains every guide
    (protospacer + PAM) and runs Bowtie2 allowing up to N mismatches.
*   The SAM output is parsed to count how many alignments each guide
    obtained.  The *intended* on‑target hit is always counted, so the
    number of **off‑target** hits = total – 1.
"""

import subprocess
import tempfile
from pathlib import Path
from collections import defaultdict
import sys


def _run_bowtie2(fasta_path: Path,
                 index_prefix: Path,
                 max_mismatches: int = 3,
                 seed_len: int = 20) -> str:
    """
    Executes Bowtie2 and returns the raw SAM text.
    """
    cmd = [
        "bowtie2",
        "-x", str(index_prefix),          # Bowtie2 index prefix (no suffix)
        "-f",                             # input is FASTA
        "-U", str(fasta_path),            # single‑ended reads (our guides)
        "-N", str(max_mismatches),        # max mismatches in the seed
        "-L", str(seed_len),              # seed length (default 20)
        "-a",                             # report *all* alignments
        "--quiet",                        # silence the alignment summary
        "-S", "-"                          # write SAM to STDOUT
    ]

    # Run Bowtie2 and capture STDOUT (the SAM)
    try:
        result = subprocess.run(
            cmd, check=True, capture_output=True, text=True
        )
    except FileNotFoundError:
        sys.exit("❌  Bowtie2 not found – please install it and ensure the binary is on your $PATH.")
    except subprocess.CalledProcessError as e:
        sys.exit(f"❌  Bowtie2 failed:\n{e.stderr}")

    return result.stdout


def count_off_targets(guides_df,
                     genome_index_prefix: Path,
                     max_mismatches: int = 3,
                     max_off_targets: int = 0) -> dict:
    """
    Takes a DataFrame that contains a column ``protospacer`` (20 nt) and ``pam`` (3 nt).
    Returns a dictionary ``guide_id → off_target_count``.
    Guides that exceed *max_off_targets* are **removed** from the result.
    """
    # -------------------------------------------------
    # 1.  Write a temporary FASTA with all guides (protospacer+PAM)
    # -------------------------------------------------
    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".fa") as fa:
        for i, row in guides_df.iterrows():
            guide_id = f"g{i}"
            seq = row["protospacer"] + row["pam"]  # 23 nt (20 + 3)
            fa.write(f">{guide_id}\n{seq}\n")
        tmp_fa_path = Path(fa.name)

    # -------------------------------------------------
    # 2.  Run Bowtie2
    # -------------------------------------------------
    sam_text = _run_bowtie2(tmp_fa_path,
                            genome_index_prefix,
                            max_mismatches=max_mismatches,
                            seed_len=20)

    # -------------------------------------------------
    # 3.  Parse SAM – count alignments per guide
    # -------------------------------------------------
    # SAM fields: QNAME (guide id), FLAG, RNAME, POS, MAPQ, CIGAR, …
    off_target_counts = defaultdict(int)
    for line in sam_text.splitlines():
        if line.startswith("@"):
            continue  # header line
        fields = line.split("\t")
        qname = fields[0]          # guide identifier we gave (g0, g1, …)
        off_target_counts[qname] += 1

    # -------------------------------------------------
    # 4.  Convert to “off‑target = total‑1” (subtract the on‑target hit)
    # -------------------------------------------------
    filtered_guides = {}
    for i, row in guides_df.iterrows():
        gid = f"g{i}"
        total_hits = off_target_counts.get(gid, 0)
        off_targets = max(total_hits - 1, 0)   # cannot be negative
        if off_targets <= max_off_targets:
            filtered_guides[i] = off_targets   # keep the index of the original row

    # Clean up the temporary FASTA file
    tmp_fa_path.unlink(missing_ok=True)

    return filtered_guides
