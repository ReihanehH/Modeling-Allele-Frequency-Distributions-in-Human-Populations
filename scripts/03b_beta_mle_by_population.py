import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import warnings
from scipy import stats
import os

# ── Setup ──────────────────────────────────────────────────────────────────────
DATA_PATH  = "data/processed/all_genes_tidy.csv"
PLOT_DIR   = "data/plots/modeling"
OUTPUT_DIR = "data/output"
os.makedirs(PLOT_DIR,   exist_ok=True)
os.makedirs(OUTPUT_DIR, exist_ok=True)

GENES = ["APOE", "CYP2C19", "HLA-B"]
GENE_COLORS = {"APOE": "#E63946", "CYP2C19": "#457B9D", "HLA-B": "#2A9D8F"}
POP_SHORT = {
    "African/African American": "AFR", "Admixed American": "AMR",
    "Ashkenazi Jewish": "ASJ",         "East Asian": "EAS",
    "European (Finnish)": "FIN",       "Middle Eastern": "MID",
    "European (non-Finnish)": "NFE",   "Amish": "AMI",
    "South Asian": "SAS",              "Remaining": "REM",
}
POPULATIONS = list(POP_SHORT.keys())

MIN_VARIANTS = 10   # minimum variants required for a reliable MLE fit

# ── Load data ──────────────────────────────────────────────────────────────────
df     = pd.read_csv(DATA_PATH)
by_pop = df[df["Population"] != "Overall"].dropna(subset=["AF"])
by_pop = by_pop[by_pop["AF"] > 0]

# ── 1. Fit Beta MLE for every gene × population ────────────────────────────────
print("=" * 60)
print("BETA MLE FITTING — per gene × population")
print("=" * 60)

records = []

for gene in GENES:
    for pop in POPULATIONS:
        af = by_pop[(by_pop["Gene"] == gene) &
                    (by_pop["Population"] == pop)]["AF"].values
        af = af[(af > 0) & (af < 1)]   # strict interior of [0,1] required by Beta

        n = len(af)

        if n < MIN_VARIANTS:
            alpha, beta_par, ok = np.nan, np.nan, False
        else:
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    a, b, _, _ = stats.beta.fit(af, floc=0, fscale=1)
                ok = bool((a > 0) and (b > 0) and np.isfinite(a) and np.isfinite(b))
                alpha, beta_par = (a, b) if ok else (np.nan, np.nan)
            except Exception:
                alpha, beta_par, ok = np.nan, np.nan, False

        records.append({
            "Gene":       gene,
            "Population": pop,
            "Pop_Short":  POP_SHORT[pop],
            "n_variants": n,
            "mean_AF":    round(float(af.mean()), 8) if n > 0 else np.nan,
            "alpha":      round(alpha,    6) if ok else np.nan,
            "beta":       round(beta_par, 6) if ok else np.nan,
            "converged":  ok,
        })

        status = f"α={alpha:.4f}, β={beta_par:.4f}" if ok else "FAILED (too few variants)"
        print(f"  {gene:8s} | {POP_SHORT[pop]:3s} | n={n:5d} | {status}")

results = pd.DataFrame(records)

# ── Save table ─────────────────────────────────────────────────────────────────
results.to_csv(f"{OUTPUT_DIR}/beta_mle_by_population.csv", index=False)
print(f"\nSaved → {OUTPUT_DIR}/beta_mle_by_population.csv")

# ── Print summary ──────────────────────────────────────────────────────────────
print("\n" + "=" * 60)
print("FITTED PARAMETERS SUMMARY (converged fits only)")
print("=" * 60)
fmt = results[results["converged"]][["Gene", "Pop_Short", "n_variants", "mean_AF", "alpha", "beta"]]
print(fmt.to_string(index=False))

# ── 2. Plot 4: Heatmap of fitted α per gene × population ──────────────────────
def param_heatmap(results, param, title, cmap, filename):
    pop_shorts = [POP_SHORT[p] for p in POPULATIONS]
    matrix = np.full((len(GENES), len(POPULATIONS)), np.nan)
    for i, gene in enumerate(GENES):
        for j, pop in enumerate(POPULATIONS):
            row = results[(results["Gene"] == gene) & (results["Population"] == pop)]
            if not row.empty and row[param].notna().any():
                matrix[i, j] = row[param].values[0]

    log_mat = np.where(matrix > 0, np.log10(matrix), np.nan)

    fig, ax = plt.subplots(figsize=(14, 4))
    vmin, vmax = np.nanmin(log_mat), np.nanmax(log_mat)
    im = ax.imshow(log_mat, cmap=cmap, aspect="auto", vmin=vmin, vmax=vmax)

    ax.set_xticks(range(len(POPULATIONS)))
    ax.set_xticklabels(pop_shorts, fontsize=11)
    ax.set_yticks(range(len(GENES)))
    ax.set_yticklabels(GENES, fontsize=12, fontweight="bold")
    ax.set_title(title, fontsize=13, fontweight="bold")

    for i in range(len(GENES)):
        for j in range(len(POPULATIONS)):
            val = matrix[i, j]
            if not np.isnan(val):
                txt = f"{val:.3f}" if val < 10 else f"{val:.1f}"
                ax.text(j, i, txt, ha="center", va="center",
                        fontsize=8.5, color="black", fontweight="bold")

    cbar = plt.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label(f"log₁₀({param})", fontsize=10)
    tick_locs = np.linspace(vmin, vmax, 5)
    cbar.set_ticks(tick_locs)
    cbar.set_ticklabels([f"{10**t:.3g}" for t in tick_locs])

    plt.tight_layout()
    plt.savefig(f"{PLOT_DIR}/{filename}", dpi=150)
    plt.close()
    print(f"Saved → {PLOT_DIR}/{filename}")

param_heatmap(
    results, "alpha",
    title="Fitted Beta α per Gene × Population\n"
          "(log₁₀ colour scale; α < 1 → J-shaped / rare-skewed; α > 1 → bell-shaped)",
    cmap="YlOrRd",
    filename="04_beta_alpha_heatmap.png"
)

# ── 3. Plot 5: Heatmap of fitted β per gene × population ──────────────────────
param_heatmap(
    results, "beta",
    title="Fitted Beta β per Gene × Population\n"
          "(log₁₀ colour scale; higher β → more density concentrated near zero)",
    cmap="Blues",
    filename="05_beta_beta_heatmap.png"
)

# ── 4. Plot 6: α vs β scatter — one panel per gene ────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(16, 5))
fig.suptitle(
    "Fitted Beta Parameters: α vs β per Gene and Population\n"
    "(log–log scale; each point = one population; dashed line = α = β symmetric Beta)",
    fontsize=13, fontweight="bold"
)

for ax, gene in zip(axes, GENES):
    sub = results[(results["Gene"] == gene) & (results["converged"])]

    ax.scatter(sub["alpha"], sub["beta"],
               color=GENE_COLORS[gene], s=80, zorder=3,
               edgecolors="white", linewidths=0.5)

    for _, row in sub.iterrows():
        ax.annotate(row["Pop_Short"],
                    (row["alpha"], row["beta"]),
                    textcoords="offset points", xytext=(5, 4),
                    fontsize=8.5, color="black")

    if len(sub) > 0:
        lim_min = min(sub["alpha"].min(), sub["beta"].min()) * 0.5
        lim_max = max(sub["alpha"].max(), sub["beta"].max()) * 2
        diag = np.linspace(lim_min, lim_max, 200)
        ax.plot(diag, diag, color="grey", linestyle="--", linewidth=0.9, label="α = β")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("α  (shape — skewed toward 0 when α < 1)", fontsize=10)
    ax.set_ylabel("β  (shape — tail weight)", fontsize=10)
    ax.set_title(gene, fontsize=13, fontweight="bold")
    ax.legend(fontsize=9)
    ax.grid(True, which="both", alpha=0.3)

plt.tight_layout()
plt.savefig(f"{PLOT_DIR}/06_alpha_vs_beta_scatter.png", dpi=150)
plt.close()
print(f"Saved → {PLOT_DIR}/06_alpha_vs_beta_scatter.png")

print("\nDone.")