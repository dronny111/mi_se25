#!/usr/bin/env python3
"""
Driver that pulls everything together:

1. Fetch data bundle
2. Generate guides
3. Extract features
4. Train & evaluate
"""

import sys
sys.path.append("src")   # allow imports

from fetcher import run_fetcher
from guide import generate_guides_from_fasta
from feature import features_from_guides
from model import train_tabular, train_cnn
from pathlib import Path
import pandas as pd

# ------------------------------------
def main():
    # 1. get the data
    bundle = run_fetcher()

    # 2. pick the first CDS – change if you prefer a different transcript
    cds_path = next(bundle.glob("data/cDNA/*.fasta"))
    guides_csv = Path("data/guides/SlAREB1_guides.csv")
    guides_csv.parent.mkdir(parents=True, exist_ok=True)
    generate_guides_from_fasta(cds_path, guides_csv)

    # 3. extract features & target
    guides_df = pd.read_csv(guides_csv)
    feat_df, y = features_from_guides(guides_df)
    feat_df.to_csv("data/features/SlAREB1_features.csv", index=False)

    # # 4. train models
    # train_tabular()
    # train_cnn()

if __name__ == "__main__":
    main()
