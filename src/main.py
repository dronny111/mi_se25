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

# Plot 4

import json
import pathlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import colors
from BCBio import GFF
from typing import Tuple, List, Dict

# ----------------------------------------------------------------------
# Helper: parse a GFF3 and return a dict with exon/UTR/CDS blocks
# ----------------------------------------------------------------------
def parse_gene(gff_path: pathlib.Path,
               gene_id_prefix: str) -> Dict[str, List[Tuple[int, int]]]:
    """
    Reads a GFF3 with BCBio.GFF and extracts all sub‑features of the
    gene whose ``ID`` starts with ``gene_id_prefix``.  Returns a dict
    with the keys ``exon``, ``CDS``, ``five_prime_UTR`` and
    ``three_prime_UTR``.  Each value is a list of tuples ``(start,
    length)`` that are **0‑based, half‑open** (as used by matplotlib).

    Parameters
    ----------
    gff_path : pathlib.Path
        Path to the .gff3 file.
    gene_id_prefix : str
        e.g. ``'Solyc04g078840'`` – the prefix that uniquely identifies
        the gene of interest.

    Returns
    -------
    dict
        ``{'exon': [(s,l),…], 'CDS': [(s,l),…], ...}``
    """
    exon_features = []
    try:
        with open(gff_path, 'r') as file:
            for rec in GFF.parse(file):
                if rec.id.startswith('SL4.0ch04'):  # Check chromosome 4
                    for feature in rec.features:


                        # Check if this is our gene of interest
                        if (hasattr(feature, 'qualifiers') and 
                            'ID' in feature.qualifiers and 
                            any(gene_id_prefix in id_val for id_val in feature.qualifiers['ID'])):
                            gene_found = True
                            print(f"Found gene: {feature.qualifiers['ID'][0]}")

                            
                        
                        if (feature.type == 'exon' and 
                            hasattr(feature, 'qualifiers') and 
                            'Parent' in feature.qualifiers and
                            any(gene_id_prefix in parent for parent in feature.qualifiers['Parent'])):
                            
                            start_pos = int(feature.location.start)
                            end_pos = int(feature.location.end)
                            length = end_pos - start_pos
                            
                            exon_features.append((start_pos, length))

                            print(f"Found Exon: {start_pos}-{end_pos} (length: {length})")

    except FileNotFoundError:
        print(f"Error: GFF3 file '{gff_path}' not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading GFF3 file: {e}")
        sys.exit(1)
    
    blocks = {"exon": exon_features,}

    return blocks

    # # Initialise containers
    # blocks = {"exon": [], "CDS": [], "five_prime_UTR": [], "three_prime_UTR": []}
    # for feat in gene_record.features:
    #     ftype = feat.type
    #     if ftype not in blocks:
    #         continue
    #     # BCBio returns location as FeatureLocation (0‑based, half‑open)
    #     start = int(feat.location.start.position)
    #     length = int(feat.location.end.position - feat.location.start.position)
    #     blocks[ftype].append((start, length))

    # # Sort every list left‑to‑right (genomic order)
    # for k in blocks:
    #     blocks[k].sort(key=lambda x: x[0])
    # return blocks


def plot_slareb1(
    gff_path: pathlib.Path,
    gene_id_prefix: str,
    domain_aa: Tuple[int, int],
    guides_path: pathlib.Path,
    top5_path: pathlib.Path,
    figsize: Tuple[int, int] = (14, 3),
    cmap: str = "viridis",
    interactive: bool = False,
    outfile: pathlib.Path | None = None,
):
    """
    Create a linear map of *SlAREB1* (or any gene you point at).

    Parameters
    ----------
    gff_path           : Path to the exon‑CDS‑UTR GFF3.
    gene_id_prefix    : e.g. "Solyc04g078840".
    domain_aa          : (start_aa, end_aa) inclusive – e.g. (122, 180).
    guides_path        : CSV with columns ``genomic_start`` (1‑based)
                         and any guide‑score columns you like.
    top5_path          : CSV with the *RL‑selected* guides – same format.
    figsize            : Figure size.
    cmap               : Matplotlib colormap used for colouring all guides.
    interactive        : If True, produces a Plotly figure instead of Matplotlib.
    outfile            : If given, saves the figure (png / pdf for mpl,
                         html for Plotly) and also returns the figure object.
    """
    # ------------------------------------------------------------------
    # 1️⃣ Load gene structure
    # ------------------------------------------------------------------
    blocks = parse_gene(gff_path, gene_id_prefix)

    print(blocks)

    # Gene orientation (assume all sub‑features share the same strand)
    strand = "-"  # for SlAREB1 – you could infer from the GFF if you wish

    # Coding‑strand coordinates (0 = first base of CDS)
    # cds_blocks = blocks["exon"]
    # The first CDS block *in genomic order* is the **most 3'** one on the
    # minus strand, so we have to reverse the coordinates.
    gene_start = min(s for s, _ in blocks["exon"])
    gene_end   = max(s + l for s, l in blocks["exon"])
    gene_len   = gene_end - gene_start

    # Helper: convert a genomic coordinate (1‑based) to a position on the
    # coding strand (0‑based, start of CDS = 0)
    def genomic_to_coding(pos_1based: int) -> int:
        """Return distance from CDS start on the coding strand."""
        # Position on the coding strand = (gene_end - pos)  because we are on the - strand
        return gene_end - pos_1based

    # ------------------------------------------------------------------
    # 2️⃣ Load guide data
    # ------------------------------------------------------------------
    all_guides = pd.read_csv(guides_path)
    top5 = pd.read_csv(top5_path)

    # Compute the coding‑strand coordinate of the **center** of each guide
    # (most tools report the 5′‑most base of the protospacer, i.e. the start)
    all_guides["pos_coding"] = all_guides["genomic_start"].apply(genomic_to_coding)
    top5["pos_coding"]      = top5["genomic_start"].apply(genomic_to_coding)

    # ------------------------------------------------------------------
    # 3️⃣ Define the bZIP‑1 domain on the coding strand
    # ------------------------------------------------------------------
    aa_start, aa_end = domain_aa
    # Convert AA → nt (0‑based, half‑open)
    nt_start = (aa_start - 1) * 3
    nt_end   = aa_end * 3               # exclusive
    domain_rel = (nt_start, nt_end - nt_start)

    # ------------------------------------------------------------------
    # 4️⃣ Plot (Matplotlib)
    # ------------------------------------------------------------------
    if interactive:
        # ---------------------------------------------------------------
        #   Plotly version – a few extra lines but very handy for
        #   exploring guide‑specific data with tooltips.
        # ---------------------------------------------------------------
        import plotly.graph_objects as go
        fig = go.Figure()

        # Gene backbone (light gray)
        fig.add_shape(
            type="rect",
            x0=0, x1=gene_len,
            y0=0.35, y1=0.65,
            fillcolor="lightgray",
            line=dict(width=0),
        )
        # Strand arrow (points left)
        fig.add_annotation(
            x=gene_len * 0.97,
            y=0.5,
            axref="x",
            ayref="y",
            ax=gene_len * 0.03,
            ay=0.5,
            showarrow=True,
            arrowhead=3,
            arrowsize=2,
            arrowwidth=2,
            arrowcolor="black",
        )

        # Exons / UTR / CDS coloured separately
        colour_map = {
            "exon": "black",
            "five_prime_UTR": "orange",
            "three_prime_UTR": "green",
            "CDS": "steelblue",
        }
        ybottom, yheight = 0.35, 0.30
        for ftype, colour in colour_map.items():
            for start, length in blocks[ftype]:
                # Convert genomic start → coding‑strand start
                if strand == "-":
                    rel_start = gene_end - (start + length)
                else:
                    rel_start = start - gene_start
                fig.add_shape(
                    type="rect",
                    x0=rel_start,
                    x1=rel_start + length,
                    y0=ybottom,
                    y1=ybottom + yheight,
                    fillcolor=colour,
                    line=dict(width=0),
                )

        # bZIP‑1 domain (semi‑transparent)
        fig.add_shape(
            type="rect",
            x0=domain_rel[0],
            x1=domain_rel[0] + domain_rel[1],
            y0=ybottom,
            y1=ybottom + yheight,
            fillcolor="royalblue",
            opacity=0.4,
            line=dict(width=1, dash="dash"),
        )
        fig.add_annotation(
            x=domain_rel[0] + domain_rel[1] / 2,
            y=ybottom + yheight + 0.05,
            text=f"bZIP‑1 (AA {aa_start}-{aa_end})",
            showarrow=False,
            font=dict(color="royalblue"),
        )

        # All guides – colour by “on_target_score” if present, otherwise gray
        guide_score_col = None
        for col in ["on_target_score", "score", "cutting_efficiency"]:
            if col in all_guides.columns:
                guide_score_col = col
                break
        if guide_score_col:
            norm = colors.Normalize(vmin=all_guides[guide_score_col].min(),
                                    vmax=all_guides[guide_score_col].max())
            cmap_obj = plt.get_cmap(cmap)
            guide_colors = [colors.to_hex(cmap_obj(norm(v))) for v in all_guides[guide_score_col]]
        else:
            guide_colors = ["gray"] * len(all_guides)

        fig.add_trace(
            go.Scatter(
                x=all_guides["pos_coding"],
                y=[0.6] * len(all_guides),
                mode="markers",
                marker=dict(symbol="line-ns-open", size=8, color=guide_colors),
                hovertemplate=(
                    "Genomic start: %{customdata[0]}<br>"
                    "Coding pos: %{x}<br>"
                    "%{customdata[1]}"
                ),
                customdata=np.stack(
                    [all_guides["genomic_start"],
                     all_guides[guide_score_col] if guide_score_col else [""] * len(all_guides)],
                    axis=-1,
                ),
                name="All guides",
            )
        )

        # Top‑5 RL guides – red arrows + hover info
        for _, row in top5.iterrows():
            fig.add_annotation(
                x=row["pos_coding"],
                y=0.8,
                ax=row["pos_coding"],
                ay=0.6,
                axref="x",
                ayref="y",
                showarrow=True,
                arrowhead=2,
                arrowsize=1.2,
                arrowwidth=2,
                arrowcolor="red",
                text=row["protospacer"][:6] + "…",
                font=dict(color="darkred", size=9),
                align="center",
                yanchor="bottom",
            )
            # Tooltip for the arrow
            fig.add_trace(
                go.Scatter(
                    x=[row["pos_coding"]],
                    y=[0.8],
                    mode="markers",
                    marker=dict(opacity=0),
                    hovertemplate=(
                        "Top‑5 guide<br>"
                        f"Protospacer: {row['protospacer']}<br>"
                        f"PAM: {row['pam']}<br>"
                        f"Genomic start: {row['genomic_start']}"
                    ),
                    showlegend=False,
                )
            )

        fig.update_yaxes(visible=False, range=[0, 1])
        fig.update_xaxes(title_text="Position on coding strand (bp from CDS start)")
        fig.update_layout(
            title=dict(
                text=(
                    f"Linear map of {gene_id_prefix} (bZIP‑1 domain highlighted)<br>"
                    "Grey lines = all NGG guides, red arrows = RL‑selected top‑5"
                ),
                x=0.5,
                xanchor="center",
            ),
            height=250,
            margin=dict(l=40, r=40, t=80, b=40),
            legend=dict(yanchor="top", y=1.02, xanchor="right", x=1),
            showlegend=False,
        )

        if outfile:
            fig.write_html(str(outfile))
        return fig

    # ------------------------------------------------------------------
    #   Matplotlib (static) – this is what you asked for originally.
    # ------------------------------------------------------------------
    import matplotlib.pyplot as plt  # already imported, kept for clarity

    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)

    # ---- Gene backbone -------------------------------------------------
    ax.broken_barh([(0, gene_len)], (0.4, 0.2), facecolors="lightgray")

    # ---- Strand indicator arrow ----------------------------------------
    arrow_x = gene_len * 0.97
    ax.arrow(
        arrow_x,
        0.5,
        -gene_len * 0.94,
        0,
        width=0.01,
        head_width=0.08,
        head_length=gene_len * 0.02,
        length_includes_head=True,
        fc="black",
        ec="black",
    )

    # ---- Feature blocks (exon / UTR / CDS) -----------------------------
    colour_map = {
        "exon": "black",
        "five_prime_UTR": "orange",
        "three_prime_UTR": "green",
        "CDS": "steelblue",
    }
    for ftype, colour in colour_map.items():
        for s, l in blocks[ftype]:
            # Convert to coding‑strand coordinate system
            if strand == "-":
                rel_start = gene_end - (s + l)
            else:
                rel_start = s - gene_start
            ax.broken_barh([(rel_start, l)], (0.4, 0.2), facecolors=colour)

    # ---- bZIP‑1 domain (semi‑transparent) -------------------------------
    ax.broken_barh([domain_rel], (0.4, 0.2), facecolors="royalblue", alpha=0.4)
    ax.text(
        domain_rel[0] + domain_rel[1] / 2,
        0.65,
        f"bZIP‑1 (AA {aa_start}-{aa_end})",
        ha="center",
        va="bottom",
        fontsize=9,
        color="royalblue",
    )

    # ---- All guides (tiny tick marks) ----------------------------------
    # Colour by a score column if present, otherwise uniform gray.
    score_col = None
    for col in ["on_target_score", "score", "cutting_efficiency"]:
        if col in all_guides.columns:
            score_col = col
            break

    if score_col:
        norm = colors.Normalize(vmin=all_guides[score_col].min(),
                                vmax=all_guides[score_col].max())
        cmap_obj = plt.get_cmap(cmap)
        guide_colors = [cmap_obj(norm(v)) for v in all_guides[score_col]]
    else:
        guide_colors = ["gray"] * len(all_guides)

    ax.scatter(
        all_guides["pos_coding"],
        np.full_like(all_guides["pos_coding"], 0.6),
        marker="|",
        s=80,
        c=guide_colors,
        edgecolors="none",
        alpha=0.7,
        linewidth=0,
        zorder=2,
    )

    # ---- Top‑5 RL guides (red arrows + label) -------------------------
    for _, row in top5.iterrows():
        pos = row["pos_coding"]
        ax.arrow(
            pos,
            0.75,
            0,
            -0.15,
            width=0.02,
            head_width=5,
            head_length=0.04,
            length_includes_head=True,
            fc="red",
            ec="red",
            zorder=5,
        )
        ax.text(
            pos,
            0.80,
            row["protospacer"][:6] + "…",
            ha="center",
            va="bottom",
            fontsize=8,
            color="darkred",
        )

    # ---- Aesthetic tweaks -----------------------------------------------
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.set_xlabel("Position on coding strand (bp from start of CDS)")
    ax.set_title(
        f"Linear map of {gene_id_prefix} (bZIP‑1 domain highlighted)\n"
        "Grey ticks = all NGG guides; red arrows = RL‑selected top‑5"
    )

    # ---- Legend (feature colours) ---------------------------------------
    legend_patches = [
        mpatches.Patch(color=col, label=label)
        for label, col in [
            ("CDS", colour_map["CDS"]),
            ("5′ UTR", colour_map["five_prime_UTR"]),
            ("3′ UTR", colour_map["three_prime_UTR"]),
            ("Exon (non‑CDS)", colour_map["exon"]),
            ("bZIP‑1 domain", "royalblue"),
        ]
    ]
    ax.legend(handles=legend_patches, loc="upper right", fontsize=8)

    # ------------------------------------------------------------------
    #   Save / return
    # ------------------------------------------------------------------
    if outfile:
        fig.savefig(outfile, dpi=300, bbox_inches="tight")
        print(f"Saved figure → {outfile}")

    plt.show()
    return fig



# ------------------------------------------------------------
# 1) Load GFF3 to get exon coordinates of SlAREB1
# ------------------------------------------------------------
gff_path = Path("data/Slycopersicum_691_ITAG4.0.gene_exons.gff3")
guides = Path("data/guides/SlAREB1_guides_rule2.csv")
top5 = Path("data/guides/SlAREB1_guides_rule2_rl_selected.csv")

plot_slareb1(
    gff_path=gff_path,
    gene_id_prefix="Solyc04g078840",
    domain_aa=(122, 180),                # bZIP‑1 domain (AA 122‑180)
    guides_path=guides,
    top5_path=top5,
    figsize=(15, 2.5),
    cmap="viridis",                     # colour map for guide scores
    interactive=False,                  # ← set True for Plotly
    outfile=Path("data/figs/SlAREB1_linear_map.png"),
)



print("\n🚀  Full pipeline finished – you now have a data‑driven shortlist for experiments.")
