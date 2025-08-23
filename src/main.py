#!/usr/bin/env python3
import sys
sys.path.append("src")  # add src to path for imports

from pathlib import Path
from guide import generate_guides_from_cds_fasta
from feature import features_from_guides
from model import train_tabular, train_cnn, load_data
from rl_guide_selection import run_thompson_sampling
import pandas as pd

import random

random.seed(42)


# 2)  Pick the SlAREB1 CDS
cds = "data/bundle/SlAREB1_bundle_SGN_first/CDS/Solyc04g078840.2.1.fasta"
# cds = "data/Slycopersicum_691_ITAG4.0.cds_primaryTranscriptOnly.fa"

# 3)  Generate guides + scores + repair prediction
raw_guides = Path("data/guides/SlAREB1_guides_rule2.csv")
generate_guides_from_cds_fasta(cds, raw_guides)

# 4)  Optional ML demo (kept from the original project)
guides = pd.read_csv(raw_guides)
print(f"📊  {len(guides)} guides generated from {cds}")



feat, _ = features_from_guides(guides)
feat_path = Path("data/features/SlAREB1_features.csv")
feat.to_csv(feat_path, index=False)

X_tab, X_onehot, y = load_data()
train_tabular(X_tab, y)
train_cnn(X_onehot, y)

# 5)  **RL bandit** – propose top‑N guides for the wet‑lab
print("\n🎯  Running Thompson‑Sampling bandit …")
rl_top = run_thompson_sampling(pd.read_csv(raw_guides),
                               n_pulls=200, top_n=5)
rl_path = raw_guides.with_name(raw_guides.stem + "_rl_selected.csv")
rl_top.to_csv(rl_path, index=False)
print(f"✅  RL‑selected guides saved to {rl_path}")


# ------------------------------------------------------------
# 1) Plot 1: “Guide landscape” (score + position)
# ------------------------------------------------------------
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


all_guides= pd.read_csv(raw_guides)
top5 = pd.read_csv(rl_path)
# ------------------------------------------------------------
# 2)  Determine the length of the gene (coding‑strand coordinate system)
# ------------------------------------------------------------
# The longest possible start coordinate is the max genomic_start value;
# the gene length is roughly the sum of start+20 (protospacer length)
gene_len = all_guides['genomic_end'].max()        # end includes PAM (20+3 bp)
# Convert the *0‑based* start on the reverse strand to a *1‑based* coding‑strand coordinate:
#   coding_position = gene_len - genomic_start
all_guides['pos_coding'] = gene_len - all_guides['genomic_start']
top5['pos_coding']      = gene_len - top5['genomic_start']

# ------------------------------------------------------------
# 3)  Make the scatter plot (all guides in gray, top‑5 in red)
# ------------------------------------------------------------
plt.figure(figsize=(10,4))

# background – all guides (light gray)
sns.scatterplot(
    data=all_guides,
    x='pos_coding',
    y='score',
    hue='strand',                     # optional: colour by strand
    palette={'+' : '#1f77b4', '-':'#ff7f0e'},
    alpha=0.3,
    edgecolor=None,
    legend=False
)

# foreground – the 5 RL‑selected guides (large red X)
sns.scatterplot(
    data=top5,
    x='pos_coding',
    y='score',
    color='red',
    s=120,
    marker='X',
    edgecolor='black',
    linewidth=1.2,
    label='RL‑selected (top‑5)'
)

# annotate each of the five points with its protospacer and posterior mean
for _, row in top5.iterrows():
    plt.text(
        row['pos_coding'],
        row['score'] + 0.02,
        f"{row['protospacer'][:6]}…\nμ={row['posterior_mean']:.2f}",
        fontsize=9,
        ha='center',
        va='bottom',
        color='darkred'
    )

plt.title("All NGG guides on SlAREB1 – Rule‑Set 2 efficiency")
plt.xlabel("Position on coding strand (bp from start of CDS)")
plt.ylabel("Sigmoid‑scaled Rule‑Set 2 score (0 → 1)")
plt.ylim(0, 1.05)
plt.xlim(0, gene_len)
plt.tight_layout()
plt.show()


## Plot 2
import json
import numpy as np

# ------------------------------------------------------------
# 1) Load the top‑5 table (same as above)
# ------------------------------------------------------------
top5 = pd.read_csv("data/guides/SlAREB1_guides_rule2_rl_selected.csv")

# ------------------------------------------------------------
# 2) Parse JSON and build a table: guide → indel size → prob
# ------------------------------------------------------------
records = []
for idx, row in top5.iterrows():
    outcomes = json.loads(row['predicted_outcome'])
    for o in outcomes:
        records.append({
            "guide_id": row['protospacer'],
            "size": o['size'],           # -N for deletions, +N for insertions
            "prob": o['prob']
        })

out_df = pd.DataFrame(records)

# ------------------------------------------------------------
# 3) Pivot to a matrix (guides as columns, sizes as rows)
# ------------------------------------------------------------
pivot = out_df.pivot_table(index='size',
                           columns='guide_id',
                           values='prob',
                           aggfunc='sum',
                           fill_value=0)

# Sort sizes from most negative (large deletions) to most positive
pivot = pivot.sort_index()

# ------------------------------------------------------------
# 4) Stacked bar plot
# ------------------------------------------------------------
fig, ax = plt.subplots(figsize=(9,5))

bottom = np.zeros(len(pivot.columns))
labels = pivot.columns.tolist()

for size in pivot.index:
    probs = pivot.loc[size].values
    ax.bar(labels, probs, bottom=bottom,
           label=f"{size} bp", width=0.6)
    bottom += probs

ax.set_ylabel("Predicted probability")
ax.set_title("Micro‑homology repair spectrum for the 5 RL‑selected guides")
ax.legend(title="Indel size", bbox_to_anchor=(1.05, 1), loc='upper left')
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.show()


# Plot 3

from sklearn.linear_model import LinearRegression
import numpy as np

# ------------------------------------------------------------
# 1) Load the full guide table and the RL table
# ------------------------------------------------------------
all_guides = pd.read_csv("data/guides/SlAREB1_guides_rule2.csv")
rl_guides  = pd.read_csv("data/guides/SlAREB1_guides_rule2_rl_selected.csv")

# ------------------------------------------------------------
# 2) Merge to have both "score" and "posterior_mean" in the same frame
# ------------------------------------------------------------
merged = pd.merge(all_guides,
                  rl_guides[['protospacer','posterior_mean']],
                  on='protospacer',
                  how='left')               # only the top‑5 will have a value

# ------------------------------------------------------------
# 3) Simple regression (all points)
# ------------------------------------------------------------
X = merged['score'].values.reshape(-1,1)   # predictor = Rule‑Set 2 score
y = merged['posterior_mean'].fillna(0).values  # missing → 0 (non‑selected)

model = LinearRegression()
model.fit(X, y)
pred = model.predict(X)

# ------------------------------------------------------------
# 4) Plot
# ------------------------------------------------------------
plt.figure(figsize=(7,5))
sns.scatterplot(data=merged,
                x='score',
                y='posterior_mean',
                hue=merged['posterior_mean'].notna(),
                palette={True:'red', False:'gray'},
                alpha=0.6,
                legend=False)

# regression line
xs = np.linspace(0,1,200).reshape(-1,1)
plt.plot(xs, model.predict(xs), '--', color='black',
         label=f"y = {model.coef_[0]:.2f}·x + {model.intercept_: .2f}")

# highlight the top‑5 points
top = merged[merged['posterior_mean'].notna()]
for _, r in top.iterrows():
    plt.text(r['score']+0.02, r['posterior_mean']+0.02,
             r['protospacer'][:6]+"…", fontsize=8, color='darkred')

plt.xlabel("Sigmoid‑scaled Rule‑Set 2 score")
plt.ylabel("Posterior mean (Thompson bandit)")
plt.title("How the RL agent incorporates efficiency into its belief")
plt.xlim(0,1.05)
plt.ylim(-0.05,1.05)
plt.legend()
plt.tight_layout()
plt.show()


print("\n🚀  Full pipeline finished – you now have a RL selected shortlist for wet-lab experiments.")
