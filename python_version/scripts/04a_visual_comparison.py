import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import warnings
from scipy import stats
import os

# ── Setup ──────────────────────────────────────────────────────────────────────
DATA_PATH  = "data/processed/all_genes_tidy.csv"
PARAMS_PATH = "data/output/beta_mle_by_population.csv"
PLOT_DIR   = "data/plots/evaluation"
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
POPULATIONS = list(POP_SHORT.keys())

# ── Load data ──────────────────────────────────────────────────────────────────
df     = pd.read_csv(DATA_PATH)
params = pd.read_csv(PARAMS_PATH)

by_pop = df[df["Population"] != "Overall"].dropna(subset=["AF"])
by_pop = by_pop[by_pop["AF"] > 0]

# ── Helper: one panel — observed histogram + Beta PDF overlay ──────────────────
def plot_panel(ax, af, alpha, beta, pop_short, gene_color, converged):
    """
    Draw one subplot: log₁₀(AF) histogram with Beta PDF overlaid.
    af          : 1-D array of raw AF values (already filtered > 0)
    alpha, beta : fitted Beta parameters (NaN if not converged)
    pop_short   : population label string
    converged   : bool — whether MLE succeeded
    """
    log_af = np.log10(af)

    # Histogram of log₁₀(AF) — linear bins in log space, density=True
    ax.hist(log_af, bins=40, density=True, alpha=0.45,
            color=gene_color, label="Observed AF")

    if converged and np.isfinite(alpha) and np.isfinite(beta):
        # Fine grid in log₁₀ space for smooth PDF curve
        y_fine = np.linspace(log_af.min() - 0.2, 0, 1000)
        x_fine = 10 ** y_fine                       # back to AF scale

        # Beta PDF in log₁₀ space via Jacobian: f_log(y) = f_AF(x) * x * ln(10)
        pdf_log = stats.beta.pdf(x_fine, alpha, beta) * x_fine * np.log(10)

        ax.plot(y_fine, pdf_log, color="black", linewidth=1.8,
                label=f"Beta(α={alpha:.3f}, β={beta:.2f})")
    else:
        ax.text(0.5, 0.5, "insufficient\ndata", transform=ax.transAxes,
                ha="center", va="center", fontsize=9, color="grey")

    ax.set_title(pop_short, fontsize=10, fontweight="bold")
    ax.set_xlabel("AF (log₁₀ scale)", fontsize=8)
    ax.set_ylabel("Density", fontsize=8)
    ax.tick_params(labelsize=7)

    # x-axis: show actual AF values
    x_ticks = [t for t in range(-6, 1) if t >= log_af.min() - 0.5]
    ax.set_xticks(x_ticks)
    ax.set_xticklabels([f"$10^{{{t}}}$" for t in x_ticks], fontsize=6)
    ax.legend(fontsize=7, loc="upper right")


# ── Main plotting loop — one figure per gene ───────────────────────────────────
GRID = {
    "APOE":    (2, 5),   # 10 panels — AMI shown as "insufficient data"
    "CYP2C19": (2, 5),   # 10 panels — AMI shown as "insufficient data"
    "HLA-B":   (2, 5),   # 10 panels — all converge
}
FILE_NUM = {"APOE": "07", "CYP2C19": "08", "HLA-B": "09"}

for gene in GENES:
    nrows, ncols = GRID[gene]
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 4.2, nrows * 3.5))
    fig.suptitle(
        f"Observed AF vs Fitted Beta Distribution — {gene}\n"
        f"(one panel per population, log₁₀ space)",
        fontsize=13, fontweight="bold"
    )
    axes_flat = axes.flatten()

    for idx, pop in enumerate(POPULATIONS):
        ax = axes_flat[idx]

        # Raw AF values for this gene × population
        af = by_pop[(by_pop["Gene"] == gene) &
                    (by_pop["Population"] == pop)]["AF"].values
        af = af[(af > 0) & (af < 1)]

        # Fitted parameters
        row = params[(params["Gene"] == gene) & (params["Population"] == pop)]
        if not row.empty and row["converged"].values[0]:
            alpha = row["alpha"].values[0]
            beta  = row["beta"].values[0]
            conv  = True
        else:
            alpha, beta, conv = np.nan, np.nan, False

        if len(af) >= 3:
            plot_panel(ax, af, alpha, beta, POP_SHORT[pop],
                       GENE_COLORS[gene], conv)
        else:
            ax.set_title(POP_SHORT[pop], fontsize=10, fontweight="bold")
            ax.text(0.5, 0.5, "insufficient\ndata", transform=ax.transAxes,
                    ha="center", va="center", fontsize=9, color="grey")
            ax.axis("off")

    # Hide any unused axes (only possible for APOE/CYP2C19 3×3 with 10 pops)
    for idx in range(len(POPULATIONS), len(axes_flat)):
        axes_flat[idx].set_visible(False)

    plt.tight_layout()
    fname = f"{FILE_NUM[gene]}_fit_overlay_{gene}.png"
    plt.savefig(f"{PLOT_DIR}/{fname}", dpi=150)
    plt.close()
    print(f"Saved → {PLOT_DIR}/{fname}")

print("\nDone.")
