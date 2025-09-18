Designing Halophytic S. Lycopersicum using CRISPR-Cas9 and AI Optimized gRNA Selection

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

🚀 Getting Started

1. Prerequisites

Python 3.11+

Git installed

2. Installation

# Clone the repository

git clone https://github.com/yourusername/mi-se25.git

cd mi_se25

# Create a virtual environment

python -m venv .venv

source .venv/bin/activate      # Windows: .venv\Scripts\activate

# Install the exact dependencies
pip install -r env/requirements.txt


3. Run the full workflow

python -m src.main

The script executes the five stages sequentially and plots all graphs.


🧩 Stage‑by‑Stage Overview
1 – Data acquisition	src.fetcher	SlAREB1.fasta (CDS)	Pulls from Sol Genomics API (https://solgenomics.net/api/v1/sequence/download/single/17806894)
2 – Guide enumeration & scoring	src.guide	guides_raw.csv (protospacer, PAM, strand, raw_score, sigmoid_score)	Replace the placeholder doench_score with a real implementation for production use.
3 – Repair‑outcome prediction	src.repair_prediction	Adds column predicted_outcome (JSON list of (type, size, prob))	MHEJ model assumes only micro‑homology; see Limitations.
4 – RL guide selection	src.rl_guide_selection	guides_rl_selected.csv (α, β, posterior_mean)	Posterior mean = expected reward = (frameshift × efficiency).

