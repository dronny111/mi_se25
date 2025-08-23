
#!/usr/bin/env python3
"""
Feature extraction for the guide set.  We create a DataFrame that
will serve as the training data for a machine/deep‑learning model.
"""

import pandas as pd
import numpy as np
from collections import Counter
import re
from typing import List
import time
import itertools
# ------------------------------------
# Helper functions
# ------------------------------------


def gc_content(s):
    '''
    GC content for only the 20mer, as per the Doench paper/code
    '''

    return len(s[4:24].replace('A', '').replace('T', ''))



def normalize_features(data,axis):
    '''
    input: Pandas.DataFrame of dtype=np.float64 array, of dimensions
    mean-center, and unit variance each feature
    '''
    data -= data.mean(axis)
    data /= data.std(axis)
    return data



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
        sp = row["protospacer"]

        gc = gc_content(sp)

        gc_above_10 = (gc > 10)*1
        gc_below_10 = (gc < 10) *1
        row_feats = {
            "gc": gc,
            "GC > 10": gc_above_10,
            "GC < 10": gc_below_10,
            **kmer_counts(sp, k=4),   # 4-mers
        }
        feats.append(row_feats)
    feat_df = pd.DataFrame(feats)

    feat_df[['gc', 'GC > 10', 'GC < 10']] = normalize_features(feat_df[['gc', 'GC > 10', 'GC < 10']], axis=0)
    # add one‑hot for each guide (may be used later in a CNN)
    return feat_df, df["score"].values

if __name__ == "__main__":

    guides = pd.read_csv("data/guides/SlAREB1_guides.csv")

    # guides = guides[guides.score.isna() == False].reset_index(drop=True)

    feat_df, y = features_from_guides(guides)
    
    feat_df.to_csv("data/features/SlAREB1_features.csv", index=False)
