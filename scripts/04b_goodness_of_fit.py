import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import warnings
from scipy import stats
import os

# ── Setup ──────────────────────────────────────────────────────────────────────
DATA_PATH   = "data/processed/all_genes_tidy.csv"
PARAMS_PATH = "data/output/beta_mle_by_population.csv"
PLOT_DIR    = "data/plots/evaluation"
OUTPUT_DIR  = "data/output"
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

# ── Load data ──────────────────────────────────────────────────────────────────
df     = pd.read_csv(DATA_PATH)
params = pd.read_csv(PARAMS_PATH)

by_pop = df[df["Population"] != "Overall"].dropna(subset=["AF"])
by_pop = by_pop[by_pop["AF"] > 0]

# ── 1. Kolmogorov–Smirnov test for every converged fit ────────────────────────
# H0: the observed AF values follow Beta(α, β) with the MLE-fitted parameters.
# We use the one-sample KS test comparing the empirical CDF to the fitted Beta CDF.
# Note: because parameters are estimated from the same data, the true p-value is
# conservative (the test is slightly anti-conservative), but the D-statistic
# remains a valid measure of the maximum discrepancy between the two CDFs.

print("=" * 65)
print("KOLMOGOROV–SMIRNOV GOODNESS-OF-FIT TEST")
print("H0: data follow Beta(α, β) with MLE-fitted parameters")
print("=" * 65)
print(f"\n{'Gene':8s} {'Pop':3s} | {'n':>5s} | {'D-stat':>8s} | {'p-value':>12s} | Reject H0?")
print("-" * 65)

ks_records = []

for _, row in params[params["converged"]].iterrows():
    gene, pop = row["Gene"], row["Population"]
    a, b      = row["alpha"], row["beta"]

    af = by_pop[(by_pop["Gene"] == gene) &
                (by_pop["Population"] == pop)]["AF"].values
    af = af[(af > 0) & (af < 1)]

    # One-sample KS test against fitted Beta CDF
    ks_stat, p_val = stats.kstest(af, lambda x: stats.beta.cdf(x, a, b))
    reject = p_val < 0.05

    ks_records.append({
        "Gene":       gene,
        "Population": pop,
        "Pop_Short":  row["Pop_Short"],
        "n":          len(af),
        "alpha":      a,
        "beta":       b,
        "KS_stat":    round(ks_stat, 6),
        "p_value":    p_val,
        "reject_H0":  reject,
    })

    flag = "YES ***" if reject else "NO"
    print(f"{gene:8s} {row['Pop_Short']:3s} | {len(af):5d} | "
          f"{ks_stat:8.4f} | {p_val:12.4e} | {flag}")

ks_df = pd.DataFrame(ks_records)
ks_df.to_csv(f"{OUTPUT_DIR}/ks_test_results.csv", index=False)
print(f"\nSaved → {OUTPUT_DIR}/ks_test_results.csv")

n_total   = len(ks_df)
n_reject  = ks_df["reject_H0"].sum()
print(f"\nSummary: {n_reject}/{n_total} fits reject H0 at α=0.05")
print("\nKS statistic range by gene:")
print(ks_df.groupby("Gene")["KS_stat"].agg(["min", "max", "mean"]).round(4))

# ── 2. Plot 10: KS statistic heatmap — gene × population ──────────────────────
print("\n── Generating plots ──")

pop_shorts = [POP_SHORT[p] for p in POPULATIONS]
matrix_ks  = np.full((len(GENES), len(POPULATIONS)), np.nan)

for i, gene in enumerate(GENES):
    for j, pop in enumerate(POPULATIONS):
        row = ks_df[(ks_df["Gene"] == gene) & (ks_df["Population"] == pop)]
        if not row.empty:
            matrix_ks[i, j] = row["KS_stat"].values[0]

fig, ax = plt.subplots(figsize=(14, 4))
im = ax.imshow(matrix_ks, cmap="RdYlGn_r", aspect="auto", vmin=0, vmax=0.55)

ax.set_xticks(range(len(POPULATIONS)))
ax.set_xticklabels(pop_shorts, fontsize=11)
ax.set_yticks(range(len(GENES)))
ax.set_yticklabels(GENES, fontsize=12, fontweight="bold")
ax.set_title(
    "Kolmogorov–Smirnov Statistic per Gene × Population\n"
    "(higher D → larger discrepancy between observed and fitted Beta; "
    "grey = insufficient data)",
    fontsize=13, fontweight="bold"
)

for i in range(len(GENES)):
    for j in range(len(POPULATIONS)):
        val = matrix_ks[i, j]
        if not np.isnan(val):
            # Mark the one non-rejected fit with a circle
            row = ks_df[(ks_df["Gene"] == GENES[i]) &
                        (ks_df["Population"] == POPULATIONS[j])]
            if not row.empty and not row["reject_H0"].values[0]:
                ax.text(j, i, f"{val:.3f}\n✓", ha="center", va="center",
                        fontsize=8.5, color="black", fontweight="bold")
            else:
                ax.text(j, i, f"{val:.3f}", ha="center", va="center",
                        fontsize=8.5, color="black", fontweight="bold")
        else:
            ax.text(j, i, "n/a", ha="center", va="center",
                    fontsize=8, color="grey")

cbar = plt.colorbar(im, ax=ax, shrink=0.8)
cbar.set_label("KS D-statistic", fontsize=10)
tick_locs = np.linspace(0, 0.55, 6)
cbar.set_ticks(tick_locs)
cbar.set_ticklabels([f"{t:.2f}" for t in tick_locs])

plt.tight_layout()
plt.savefig(f"{PLOT_DIR}/10_ks_statistic_heatmap.png", dpi=150)
plt.close()
print(f"Saved → {PLOT_DIR}/10_ks_statistic_heatmap.png")

# ── 3. Plot 11: QQ plots — one figure per gene, one panel per population ───────
# QQ plot: observed quantiles vs theoretical Beta quantiles.
# Points on the diagonal y = x indicate perfect fit.
# Systematic deviations reveal where the model fails.

for gene in GENES:
    fig, axes = plt.subplots(2, 5, figsize=(21, 8))
    fig.suptitle(
        f"QQ Plots: Observed vs Beta({gene}) Quantiles — {gene}\n"
        f"(points on the diagonal = perfect fit; "
        f"deviations show where the Beta model fails)",
        fontsize=13, fontweight="bold"
    )
    axes_flat = axes.flatten()

    for idx, pop in enumerate(POPULATIONS):
        ax = axes_flat[idx]
        pop_short = POP_SHORT[pop]

        row = ks_df[(ks_df["Gene"] == gene) & (ks_df["Population"] == pop)]

        af = by_pop[(by_pop["Gene"] == gene) &
                    (by_pop["Population"] == pop)]["AF"].values
        af = np.sort(af[(af > 0) & (af < 1)])

        if row.empty or len(af) < 3:
            ax.set_title(pop_short, fontsize=10, fontweight="bold")
            ax.text(0.5, 0.5, "insufficient\ndata", transform=ax.transAxes,
                    ha="center", va="center", fontsize=9, color="grey")
            ax.axis("off")
            continue

        a    = row["alpha"].values[0]
        b    = row["beta"].values[0]
        ks_d = row["KS_stat"].values[0]
        p_v  = row["p_value"].values[0]
        rej  = row["reject_H0"].values[0]

        # Theoretical quantiles: Beta ppf at the same probability points
        probs       = (np.arange(1, len(af) + 1)) / (len(af) + 1)
        theoretical = stats.beta.ppf(probs, a, b)

        ax.scatter(theoretical, af, s=6, alpha=0.5,
                   color=GENE_COLORS[gene], linewidths=0)

        # y = x diagonal — perfect fit line
        lim = max(theoretical.max(), af.max()) * 1.05
        ax.plot([0, lim], [0, lim], color="black", linewidth=1,
                linestyle="--", label="y = x (perfect fit)")

        # Annotate with KS result
        p_label = f"p={p_v:.2e}" if p_v >= 1e-4 else f"p<10⁻⁴"
        reject_label = "reject H₀" if rej else "fail to reject H₀"
        color_label  = "#C0392B" if rej else "#27AE60"
        ax.set_title(pop_short, fontsize=10, fontweight="bold")
        ax.text(0.03, 0.97,
                f"D={ks_d:.3f}\n{p_label}\n{reject_label}",
                transform=ax.transAxes, ha="left", va="top",
                fontsize=7.5, color=color_label,
                bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.7))

        ax.set_xlabel("Theoretical Beta quantiles", fontsize=8)
        ax.set_ylabel("Observed AF quantiles",      fontsize=8)
        ax.tick_params(labelsize=7)
        ax.legend(fontsize=7, loc="lower right")

    plt.tight_layout()
    fname = f"11_qq_plot_{gene}.png"
    plt.savefig(f"{PLOT_DIR}/{fname}", dpi=150)
    plt.close()
    print(f"Saved → {PLOT_DIR}/{fname}")

print("\nDone.")
