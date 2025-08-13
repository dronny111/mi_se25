#!/usr/bin/env python3
"""
Generate every 20‑nt protospacer that satisfies the NGG PAM (SpCas9) and score it
with the real Doench/Azimuth model that is available from the Microsoft
Research package “azimuth”.
"""

import re
import numpy as np
from pathlib import Path
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord  # ← keep the earlier import for reverse complements


from azimuth_ruleset2 import rule_set2_score


# ------------------------------------------------------------
# Primer utilities
# ------------------------------------------------------------
PAM_RE = re.compile(r"[ATGC]GG")   # 5'NGG

def reverse_complement(seq: str) -> str:
    return str(SeqRecord(seq).reverse_complement())

def find_guides(seq: str, strand: int = 1):
    """Return a list of dicts describing each protospacer on a single strand."""
    guides = []
    for m in PAM_RE.finditer(seq):
        pam_start, pam_end = m.span()

        if strand == 1:                       # + strand → spacer upstream of PAM
            start = pam_start - 20
            if start < 0:                    # not enough nt upstream
                continue
            spacer = seq[start:pam_start]
        else:                                 # – strand → spacer downstream of PAM
            end = pam_end + 20
            if end > len(seq):
                continue
            spacer = seq[pam_start:pam_end]
            spacer = reverse_complement(spacer)

        guides.append({
            "spacer": spacer,
            "pam": seq[pam_start:pam_end],
            "start": start if strand == 1 else pam_end,
            "end": pam_start if strand == 1 else end,
            "strand": "+" if strand == 1 else "-",
        })
    return guides

import math
def logistic(score):
    return 1/(1+math.exp(-score))

# ------------------------------------------------------------
# Real Azimuth scoring
# ------------------------------------------------------------
def azimuth_score(spacer: str) -> float:
    score = rule_set2_score(spacer)

    score = logistic(score)
    
    return float(score) if score is not None else -float("inf")

# ------------------------------------------------------------
# Generate guides + score
# ------------------------------------------------------------
def generate_guides_from_fasta(in_path: Path, out_path: Path):
    """Creates a CSV with every valid guide and its Azimuth score."""
    with open(in_path) as f:
        record = SeqIO.read(f, "fasta")
    seq = str(record.seq)

    all_guides = []
    for strand in (1, -1):
        all_guides.extend(find_guides(seq, strand=strand))

    df = pd.DataFrame(all_guides)
    # Add a proper one‑column score column
    df["score"] = df["spacer"].apply(azimuth_score)
    df.to_csv(out_path, index=False)
    print(f"Saved {len(df)} guides to {out_path}")

# ------------------------------------------------------------
# CLI driver (kept the same as before)
# ------------------------------------------------------------
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Generate SpCas9 guides for SlAREB1.")
    parser.add_argument("in_fasta", type=Path, help="CDS FASTA file")
    parser.add_argument("out_csv", type=Path, help="Output CSV")
    args = parser.parse_args()
    generate_guides_from_fasta(args.in_fasta, args.out_csv)
