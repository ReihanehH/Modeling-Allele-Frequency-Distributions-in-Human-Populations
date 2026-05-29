import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import os

# ── Setup ──────────────────────────────────────────────────────────────────────
DATA_PATH = "data/processed/all_genes_tidy.csv"
PLOT_DIR  = "data/plots/exploratory"
os.makedirs(PLOT_DIR, exist_ok=True)

GENES = ["APOE", "CYP2C19", "HLA-B"]
GENE_COLORS = {"APOE": "#E63946", "CYP2C19": "#457B9D", "HLA-B": "#2A9D8F"}
POP_SHORT = {
    "African/African American": "AFR", "Admixed American": "AMR",
    "Ashkenazi Jewish": "ASJ",         "East Asian": "EAS",
    "European (Finnish)": "FIN",       "Middle Eastern": "MID",
    "European (non-Finnish)": "NFE",   "Amish": "AMI",
    "South Asian": "SAS",              "Remaining": "REM",
}

df      = pd.read_csv(DATA_PATH)
overall = df[df["Population"] == "Overall"].dropna(subset=["AF"])
by_pop  = df[df["Population"] != "Overall"].dropna(subset=["AF"])

# ── Plot 4: ECDF — cumulative proportion of variants by AF ─────────────────────
# Shows clearly what fraction of variants fall below any AF threshold
fig, ax = plt.subplots(figsize=(10, 6))

for gene in GENES:
    data = np.sort(overall[overall["Gene"] == gene]["AF"].values)
    ecdf = np.arange(1, len(data) + 1) / len(data)
    ax.plot(data, ecdf, color=GENE_COLORS[gene], linewidth=2, label=gene)

# Mark key thresholds
for thresh, label in [(0.001, "0.1%"), (0.01, "1%"), (0.05, "5%")]:
    ax.axvline(thresh, color="grey", linestyle="--", linewidth=0.8, alpha=0.7)
    ax.text(thresh * 1.1, 0.05, label, color="grey", fontsize=9)

ax.set_xscale("log")
ax.set_xlabel("Allele Frequency (log scale)", fontsize=12)
ax.set_ylabel("Cumulative Proportion of Variants", fontsize=12)
ax.set_title("ECDF of Allele Frequency per Gene\n(shows what fraction of variants fall below each AF threshold)",
             fontsize=13, fontweight="bold")
ax.legend(fontsize=11)
ax.set_ylim(0, 1.02)
ax.yaxis.set_major_formatter(ticker.PercentFormatter(xmax=1))
plt.tight_layout()
plt.savefig(f"{PLOT_DIR}/04_ecdf_AF.png", dpi=150)
plt.close()
print(f"Saved → {PLOT_DIR}/04_ecdf_AF.png")

# ── Plot 5: Heatmap — median log₁₀(AF) per gene × population ─────────────────
pop_labels = list(POP_SHORT.keys())

median_matrix = np.zeros((len(GENES), len(pop_labels)))
for i, gene in enumerate(GENES):
    for j, pop in enumerate(pop_labels):
        vals = by_pop[(by_pop["Gene"] == gene) &
                      (by_pop["Population"] == pop) &
                      (by_pop["AF"] > 0)]["AF"]
        median_matrix[i, j] = np.log10(vals.median()) if len(vals) > 0 else np.nan

fig, ax = plt.subplots(figsize=(13, 4))
im = ax.imshow(median_matrix, cmap="RdYlGn", aspect="auto")

ax.set_xticks(range(len(pop_labels)))
ax.set_xticklabels([POP_SHORT[p] for p in pop_labels], fontsize=11)
ax.set_yticks(range(len(GENES)))
ax.set_yticklabels(GENES, fontsize=12, fontweight="bold")
ax.set_title("Median log₁₀(AF) per Gene and Population\n(greener = higher AF, redder = lower AF)",
             fontsize=13, fontweight="bold")

# Annotate cells with values
for i in range(len(GENES)):
    for j in range(len(pop_labels)):
        val = median_matrix[i, j]
        if not np.isnan(val):
            ax.text(j, i, f"{val:.1f}", ha="center", va="center",
                    fontsize=9, color="black", fontweight="bold")

plt.colorbar(im, ax=ax, label="Median log₁₀(AF)", shrink=0.8)
plt.tight_layout()
plt.savefig(f"{PLOT_DIR}/05_heatmap_median_AF.png", dpi=150)
plt.close()
print(f"Saved → {PLOT_DIR}/05_heatmap_median_AF.png")

# ── Plot 6: Rare variant proportion by population (grouped bars) ───────────────
RARE_THRESH = 0.01
fig, ax = plt.subplots(figsize=(13, 6))

n_pops  = len(pop_labels)
n_genes = len(GENES)
width   = 0.25
x       = np.arange(n_pops)

for i, gene in enumerate(GENES):
    proportions = []
    for pop in pop_labels:
        vals  = by_pop[(by_pop["Gene"] == gene) & (by_pop["Population"] == pop)]["AF"].dropna()
        if len(vals) > 0:
            proportions.append((vals < RARE_THRESH).sum() / len(vals) * 100)
        else:
            proportions.append(np.nan)
    offset = (i - 1) * width
    bars = ax.bar(x + offset, proportions, width, label=gene,
                  color=GENE_COLORS[gene], alpha=0.85, edgecolor="white")

ax.set_xticks(x)
ax.set_xticklabels([POP_SHORT[p] for p in pop_labels], fontsize=11)
ax.set_ylabel("% of Variants with AF < 1%", fontsize=12)
ax.set_title("Proportion of Rare Variants (AF < 1%) by Population and Gene",
             fontsize=13, fontweight="bold")
ax.set_ylim(0, 105)
ax.axhline(100, color="grey", linestyle="--", linewidth=0.7, alpha=0.5)
ax.legend(fontsize=11)
ax.yaxis.set_major_formatter(ticker.PercentFormatter())
plt.tight_layout()
plt.savefig(f"{PLOT_DIR}/06_rare_variant_proportion_by_population.png", dpi=150)
plt.close()
print(f"Saved → {PLOT_DIR}/06_rare_variant_proportion_by_population.png")

print("\nDone.")