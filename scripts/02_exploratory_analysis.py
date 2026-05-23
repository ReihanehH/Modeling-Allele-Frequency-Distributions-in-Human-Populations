import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

# ── Setup ──────────────────────────────────────────────────────────────────────
DATA_PATH  = "data/processed/all_genes_tidy.csv"
PLOT_DIR   = "plots/exploratory"
OUTPUT_DIR = "output"
os.makedirs(PLOT_DIR, exist_ok=True)
os.makedirs(OUTPUT_DIR, exist_ok=True)

GENES = ["APOE", "CYP2C19", "HLA-B"]
GENE_COLORS = {"APOE": "#E63946", "CYP2C19": "#457B9D", "HLA-B": "#2A9D8F"}
POP_SHORT = {
    "African/African American": "AFR", "Admixed American": "AMR",
    "Ashkenazi Jewish": "ASJ", "East Asian": "EAS",
    "European (Finnish)": "FIN", "Middle Eastern": "MID",
    "European (non-Finnish)": "NFE", "Amish": "AMI",
    "South Asian": "SAS", "Remaining": "REM",
}

# ── Load data ──────────────────────────────────────────────────────────────────
df      = pd.read_csv(DATA_PATH)
overall = df[df["Population"] == "Overall"].dropna(subset=["AF"])
by_pop  = df[df["Population"] != "Overall"].dropna(subset=["AF"])

# ── 1. Summary statistics → output/ ───────────────────────────────────────────
stats = (
    overall.groupby("Gene")["AF"]
    .describe(percentiles=[0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99])
    .round(8)
)
stats.to_csv(f"{OUTPUT_DIR}/summary_statistics_overall.csv")
print(f"Saved → {OUTPUT_DIR}/summary_statistics_overall.csv")

pop_stats = (
    by_pop.groupby(["Gene", "Population"])["AF"]
    .describe(percentiles=[0.25, 0.5, 0.75])
    .round(8)
)
pop_stats.to_csv(f"{OUTPUT_DIR}/summary_statistics_by_population.csv")
print(f"Saved → {OUTPUT_DIR}/summary_statistics_by_population.csv")

print("\n" + "=" * 60)
print("RARE VARIANT PREVALENCE (AF < 0.01)")
print("=" * 60)
for gene in GENES:
    g     = overall[overall["Gene"] == gene]
    rare  = (g["AF"] < 0.01).sum()
    total = len(g)
    print(f"  {gene:8s}: {rare:>4d} / {total} variants are rare  ({100*rare/total:.1f}%)")

# ── 2. Plot 1: linear overview + per-gene log-scale histograms (FIXED) ─────────
fig = plt.figure(figsize=(16, 10))
fig.suptitle("Overall AF Distribution per Gene", fontsize=14, fontweight="bold")

# Top: linear scale combined
ax_lin = fig.add_subplot(2, 1, 1)
for gene in GENES:
    data = overall[overall["Gene"] == gene]["AF"]
    ax_lin.hist(data, bins=80, alpha=0.6, color=GENE_COLORS[gene], label=gene)
ax_lin.set_title("Linear Scale — all three genes")
ax_lin.set_xlabel("Allele Frequency")
ax_lin.set_ylabel("Number of Variants")
ax_lin.legend()

# Bottom: separate log-scale subplot per gene
for i, gene in enumerate(GENES):
    ax = fig.add_subplot(2, 3, 4 + i)
    data = overall[(overall["Gene"] == gene) & (overall["AF"] > 0)]["AF"]
    log_bins = np.logspace(np.log10(data.min()), np.log10(data.max()), 50)
    ax.hist(data, bins=log_bins, color=GENE_COLORS[gene], alpha=0.85, edgecolor="white")
    ax.set_xscale("log")
    ax.set_title(f"{gene} — Log Scale", fontsize=11, fontweight="bold")
    ax.set_xlabel("Allele Frequency (log scale)")
    ax.set_ylabel("Number of Variants")

plt.tight_layout()
plt.savefig(f"{PLOT_DIR}/01_overall_AF_distribution.png", dpi=150)
plt.close()
print(f"\nSaved → {PLOT_DIR}/01_overall_AF_distribution.png")

# ── 3. Plot 2: AF by population — boxplots ─────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(18, 6))
fig.suptitle("AF Distribution by Population and Gene (Log Scale)", fontsize=14, fontweight="bold")

for ax, gene in zip(axes, GENES):
    gene_data = by_pop[(by_pop["Gene"] == gene) & (by_pop["AF"] > 0)]
    pop_order = (
        gene_data.groupby("Population")["AF"]
        .median().sort_values(ascending=False).index.tolist()
    )
    pop_data     = [np.log10(gene_data[gene_data["Population"] == p]["AF"].values) for p in pop_order]
    short_labels = [POP_SHORT.get(p, p) for p in pop_order]

    ax.boxplot(pop_data, tick_labels=short_labels, patch_artist=True,
               boxprops=dict(facecolor=GENE_COLORS[gene], alpha=0.6),
               medianprops=dict(color="black", linewidth=2))
    ax.set_title(gene, fontsize=13, fontweight="bold")
    ax.set_xlabel("Population")
    ax.set_ylabel("log₁₀(AF)")
    ax.tick_params(axis="x", rotation=45)

plt.tight_layout()
plt.savefig(f"{PLOT_DIR}/02_AF_by_population.png", dpi=150)
plt.close()
print(f"Saved → {PLOT_DIR}/02_AF_by_population.png")

# ── 4. Plot 3: frequency class breakdown ───────────────────────────────────────
bins   = [0, 0.001, 0.01, 0.05, 0.5, 1.0]
labels = ["<0.1%", "0.1–1%", "1–5%", "5–50%", ">50%"]

fig, axes = plt.subplots(1, 3, figsize=(15, 5))
fig.suptitle("Variant Frequency Class Breakdown per Gene", fontsize=14, fontweight="bold")

for ax, gene in zip(axes, GENES):
    data   = overall[overall["Gene"] == gene]["AF"].dropna()
    data   = data[data > 0]
    cats   = pd.cut(data, bins=bins, labels=labels, include_lowest=True)
    counts = cats.value_counts().reindex(labels).fillna(0).astype(int)

    bars = ax.bar(labels, counts.values, color=GENE_COLORS[gene], alpha=0.8, edgecolor="white")
    ax.set_yscale("log")
    ax.set_ylim(bottom=0.5)
    ax.set_title(gene, fontsize=13, fontweight="bold")
    ax.set_xlabel("Frequency Class")
    ax.set_ylabel("Number of Variants (log scale)")
    ax.tick_params(axis="x", rotation=30)

    for bar, count in zip(bars, counts.values):
        if count > 0:
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() * 1.2,
                    str(count), ha="center", va="bottom", fontsize=9, fontweight="bold")

plt.tight_layout()
plt.savefig(f"{PLOT_DIR}/03_frequency_class_breakdown.png", dpi=150)
plt.close()
print(f"Saved → {PLOT_DIR}/03_frequency_class_breakdown.png")

print("\nDone. All outputs saved.")