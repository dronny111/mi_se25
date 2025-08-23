📚 MISE 2026 – AI-driven Lab-in-the-Loop Workflow for CRISPR-mediated Saline-Tolerance in Tomato

By: Dhruv Jain (Grade 11, Lincoln Community School)


Research Mentor: Dr.(Med) Ronny Polle


🎯 Purpose

Provide a fully automated, closed‑loop pipeline that:

Downloads the SlAREB1 coding sequence (the ABA‑responsive bZIP transcription factor) from the tomato reference genome.
Enumerates every possible SpCas9 guide inside the gene.
Scores each guide with Doench Rule‑Set 2 (sigmoid‑scaled to 0‑1).
Predicts the repair‑outcome distribution using a micro‑homology‑mediated end‑joining (MHEJ) model.
Selects the optimal guide(s) with a Thompson‑sampling RL bandit that maximises frameshift and efficiency.
The output is a list of (5‑10 guides) ready for wet‑lab validation.

📂 Repository Layout

mi_se25/
│
├─ data/                     # CSVs generated at each stage
│   ├─ guides/
│   │   ├─ SlAREB1_guides.csv
│   │   └─ SlAREB1_guides_rl_selected.csv
│   └─ features/             # Optional feature tables
│
├─ env/                      # requirements.txt
│
├─ notebooks/
│   └─ analysis.ipynb        # Example exploratory notebook
│
├─ src/
│   ├─ fetcher.py            # Stage 1 – download CDS
│   ├─ guide.py              # Stage 2 – enumerate & score guides
│   ├─ repair_prediction.py  # Stage 3 – MHEJ model
│   ├─ rl_guide_selection.py# Stage 4 – Thompson‑sampling bandit
│   ├─ feature.py            # Stage 5 – optional feature extraction
│   ├─ model.py              # Demo RandomForest / 1‑D CNN
│   └─ main.py               # One‑command wrapper (runs all stages)
│
├─ scripts/
│   └─ data_fetcher.sh       # Legacy bash wrapper for Stage 1
│
├─ models/
│   ├─ rf.pkl                # Trained RandomForest (demo)
│   └─ cnn_slareb1.h5        # Trained 1‑D CNN (demo)
│
└─ README.md                 # **THIS FILE**
🚀 Getting Started

1. Prerequisites

Python 3.11+
git (to clone)
≈ 8 GB RAM (the pipeline runs comfortably on a laptop)
2. Installation

# Clone the repository
git clone https://github.com/yourusername/MISE-2025.git
cd MISE-2025

# Create a virtual environment
python -m venv .venv
source .venv/bin/activate      # Windows: .venv\Scripts\activate

# Install the exact dependencies
pip install -r env/requirements.txt
3. Run the full workflow

python -m src.main
The script executes the five stages sequentially and prints a progress bar.

Generated files

File	Description
data/guides/SlAREB1_guides_raw.csv	All 57 candidate guides with raw Doench scores & MHEJ outcome JSON
data/guides/SlAREB1_guides_rl_selected.csv	Short‑list (posterior mean, α, β, etc.)
models/rf.pkl	RandomForest trained on tabular guide features (demo)
models/cnn_slareb1.h5	1‑D CNN trained on one‑hot spacers (demo)
🧩 Stage‑by‑Stage Overview

Stage	Module	Main output	Quick note
1 – Data acquisition	src.fetcher	SlAREB1.fa (CDS)	Pulls from Sol Genomics API (https://solgenomics.net/api/v1/sequence/download/single/17806894)
2 – Guide enumeration & scoring	src.guide	guides_raw.csv (protospacer, PAM, strand, raw_score, sigmoid_score)	Replace the placeholder doench_score with a real implementation for production use.
3 – Repair‑outcome prediction	src.repair_prediction	Adds column predicted_outcome (JSON list of (type, size, prob))	MHEJ model assumes only micro‑homology; see Limitations.
4 – RL guide selection	src.rl_guide_selection	guides_rl_selected.csv (α, β, posterior_mean)	Posterior mean = expected reward = (frameshift × efficiency).
5 – Optional ML demo	src.model	rf.pkl, cnn_slareb1.h5 + console metrics	Purely illustrative; not required for the ISEF deliverable.
⚠️ Limitations & Future Work

Repair model – Only micro‑homology deletions are modelled; real NHEJ can give larger insertions or complex events.
Off‑target analysis – Not yet integrated. Future versions will call CRISPOR or Cas‑OFFinder and add an off‑target arm to the bandit for multi‑objective optimisation.
Experimental validation – The pipeline stops at in silico selection. Wet‑lab results (NGS indel frequencies) will be fed back to update the Thompson‑sampling posteriors (online learning).
Scalability – The same workflow can be applied genome‑wide to other ABA‑related TFs or to other crops.
Model upgrades – Swap the 1‑D CNN for a transformer or 2‑D CNN that consumes the full PAM‑spacer context for higher predictive power.
