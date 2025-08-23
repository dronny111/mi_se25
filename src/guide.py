#!/usr/bin/env python3
"""
Generate every NGG guide from a CDS, compute a Doench Rule‑Set 2 score,
(optionally) apply a sigmoid transform, and write a CSV that now also
contains a micro‑homology repair‑outcome prediction.
"""

import re
from pathlib import Path
import pandas as pd

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# -------------------------------------------------
# 1.  Scoring (Rule‑Set 2) and repair outcome
# -------------------------------------------------
from azimuth_ruleset2 import rule_set2_score
from repair_pred import predict_repair_outcome, serialize_outcome

# -------------------------------------------------
# 2.  Core helpers (unchanged)
# -------------------------------------------------
PAM_REGEX = re.compile(r"[ATGC]GG")   # NGG

def _rev_comp(seq: str) -> str:
    return str(Seq(seq).reverse_complement())

def _find_guides_on_strand(seq: str, strand: str) -> list[dict]:
    guides = []
    for m in PAM_REGEX.finditer(seq):
        pam_start = m.start()
        if pam_start < 20:
            continue
        protospacer = seq[pam_start - 20:pam_start]
        guides.append({
            "protospacer": protospacer,
            "pam": seq[pam_start:pam_start + 3],
            "genomic_start": pam_start - 20,
            "genomic_end":   pam_start + 3,
            "strand": strand
        })
    return guides


# -------------------------------------------------
# 3.  Main driver – generate, score, (optional) sigmoid,
#     predict repair outcome, write CSV
# -------------------------------------------------
def generate_guides_from_cds_fasta(cds_fasta: Path,
                                   out_csv: Path
                                   ) -> None:
    """
    Parameters
    ----------
    cds_fasta : Path
        FASTA containing the *coding* DNA of SlAREB1.
    out_csv : Path
        Where the guide table (including scores & repair prediction) will be written.
    apply_sigmoid : bool
        If True the raw Rule‑Set 2 score is passed through a sigmoid so the
        column ``score`` lies in (0, 1).  Set False to keep the raw value.
    """
    # -----------------------------------------------------------------
    # 1) Load the CDS (robust version that also works if the FASTA has >1 record)
    # -----------------------------------------------------------------
    cds_path_str = str(cds_fasta)
    try:
        record = SeqIO.read(cds_path_str, "fasta")
    except ValueError:                     # >1 record → just take the first
        record = next(SeqIO.parse(cds_path_str, "fasta"))
    seq = str(record.seq).upper()

    # -----------------------------------------------------------------
    # 2) Enumerate guides on both strands
    # -----------------------------------------------------------------
    guides = _find_guides_on_strand(seq, "+")
    rev_seq = _rev_comp(seq)
    guides += _find_guides_on_strand(rev_seq, "-")

    # -----------------------------------------------------------------
    # 3) Score each guide (Rule‑Set 2) and predict the repair outcome
    # -----------------------------------------------------------------
    for g in guides:
        raw = rule_set2_score(g["protospacer"])


        g["score"] = raw


        # ---- repair outcome (micro‑homology model) -----------------
        outcome = predict_repair_outcome(g["protospacer"])
        # Store the JSON string – it will appear as a single column in the CSV.
        g["predicted_outcome"] = serialize_outcome(outcome)

    df = pd.DataFrame(guides)
    
    num_scored = df["score"].notna().sum()

    df[df["score"].notna()].to_csv(out_csv, index=False)

    print(f"✅  Saved {len(df)} guides (scored {num_scored}) to {out_csv}")
