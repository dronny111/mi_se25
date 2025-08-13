# MISE 2025 CRISPR/ML Project for Tomato Saline Tolerance

## Purpose
Predict and optimise SpCas9‑guide RNA sequences for editing the *SlAREB1* gene in *Solanum lycopersicum*.

## Workflow
1. **Download reference data**  
   `python -m src.fetcher`

2. **Generate Guides**  
   `python -m src.guide`

3. **Construct Features**  
   `python -m src.feature`

4. **Train Models**  
   `python -m src.model --train both`

5. **Results**  
   – Random Forest: `models/rf.pkl`  
   – 1‑D CNN: `models/cnn_slareb1.h5`  
   – Evaluation metrics printed to console.

## What is included
- Bash script `data_fetcher.sh` (already in the repo)  
- Python modules in `src/` for each step  
- Sample Jupyter notebook in `/notebooks/analysis.ipynb`

## Extending the pipeline
- Replace the dummy `doench_score` in `src/guide.py` with a real Doench or Azimuth score.  
- Load experimentally labelled edit‑efficiency data from a public CRISPR library and re‑train the models.  
- Expand the CNN architecture to a transformer or 2‑D CNN if you include the full PAM+spacer context.  

## Dependencies
Create a virtual environment and install with:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r env/requirements.txt
