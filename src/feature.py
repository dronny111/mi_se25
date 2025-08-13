#!/usr/bin/env python3
"""
Feature extraction for the guide set.  We create a DataFrame that
will serve as the training data for a deep‑learning model.
"""

import pandas as pd
import numpy as np
from collections import Counter
import re
from typing import List

# ------------------------------------
# Helper functions
# ------------------------------------
def gc_content(seq: str) -> float:
    return (seq.count("G") + seq.count("C")) / len(seq)

def kmer_counts(seq: str, k: int=4) -> dict:
    """Return dict of kmer frequencies (normalized)."""
    counts = Counter([seq[i:i+k] for i in range(len(seq)-k+1)])
    total = sum(counts.values())
    return {k: v/total for k,v in counts.items()}

def one_hot(seq: str) -> np.ndarray:
    """Return 20‑mer as one‑hot array (4x20)."""
    mapping = {"A":0,"C":1,"G":2,"T":3}
    vec = np.zeros((4, len(seq)), dtype=int)
    for i,n in enumerate(seq):
        vec[mapping[n], i] = 1
    return vec

# ------------------------------------
# Main
# ------------------------------------
def features_from_guides(df: pd.DataFrame) -> pd.DataFrame:
    """Take guide df from guide.py and produce a numeric feature matrix."""
    feats = []
    for _, row in df.iterrows():
        sp = row["spacer"]
        row_feats = {
            "gc": gc_content(sp),
            "len": len(sp),
            **kmer_counts(sp, k=4),   # 4-mers
        }
        feats.append(row_feats)
    feat_df = pd.DataFrame(feats)
    # add one‑hot for each guide (may be used later in a CNN)
    return feat_df, df["score"].values

if __name__ == "__main__":
    guides = pd.read_csv("data/guides/SlAREB1_guides.csv")
    feat_df, y = features_from_guides(guides)
    feat_df.to_csv("data/features/SlAREB1_features.csv", index=False)
