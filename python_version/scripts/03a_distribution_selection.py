import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import os

# ── Setup ──────────────────────────────────────────────────────────────────────
PLOT_DIR = "data/plots/modeling"
os.makedirs(PLOT_DIR, exist_ok=True)

GENE_COLORS = {"APOE": "#E63946", "CYP2C19": "#457B9D", "HLA-B": "#2A9D8F"}

# ── Plot 1: Beta distribution shapes for different parameter regimes ───────────
# Show how flexible Beta is: J-shaped, U-shaped, bell-shaped
x = np.linspace(0.001, 0.999, 1000)

params = [
    (0.1, 1.0,  "α=0.1, β=1.0  → strongly J-shaped (ultra-rare variants)"),
    (0.5, 2.0,  "α=0.5, β=2.0  → J-shaped (rare-skewed)"),
    (0.5, 0.5,  "α=0.5, β=0.5  → U-shaped (two peaks)"),
    (2.0, 5.0,  "α=2.0, β=5.0  → bell-shaped, skewed right"),
    (1.0, 1.0,  "α=1.0, β=1.0  → Uniform (flat)"),
]
colors = ["#E63946", "#457B9D", "#2A9D8F", "#F4A261", "#6A0572"]

fig, ax = plt.subplots(figsize=(11, 6))
for (a, b, label), color in zip(params, colors):
    pdf = stats.beta.pdf(x, a, b)
    ax.plot(x, pdf, linewidth=2.5, label=label, color=color)

ax.set_xlim(0, 1)
ax.set_ylim(0, 10)
ax.set_xlabel("Allele Frequency", fontsize=12)
ax.set_ylabel("Probability Density", fontsize=12)
ax.set_title("Beta Distribution PDF for Different Parameter Regimes\n"
             "Beta(α, β) is defined on [0,1] — naturally suited for allele frequencies",
             fontsize=13, fontweight="bold")
ax.legend(fontsize=10)
ax.axvline(0.01, color="grey", linestyle="--", linewidth=0.8, alpha=0.6)
ax.text(0.012, 9, "AF = 1%\n(rare threshold)", color="grey", fontsize=9)
plt.tight_layout()
plt.savefig(f"{PLOT_DIR}/01_beta_distribution_shapes.png", dpi=150)
plt.close()
print(f"Saved → {PLOT_DIR}/01_beta_distribution_shapes.png")



# ── Load data (needed for Plots 2 and 3) ──────────────────────────────────────
import pandas as pd
df      = pd.read_csv("data/processed/all_genes_tidy.csv")
overall = df[df["Population"] == "Overall"].dropna(subset=["AF"])
overall = overall[overall["AF"] > 0]

# ── Plot 2: Why Beta fits AF data ...


# ── Plot 2: Why Beta fits AF data — overlay on observed histogram ──────────────
fig, axes = plt.subplots(1, 3, figsize=(16, 5))
fig.suptitle("Observed AF Histogram vs. Candidate Distributions\n"
             "(working in log₁₀ space for proper density comparison)",
             fontsize=13, fontweight="bold")

y_fine = np.linspace(-6, 0, 2000)   # log₁₀(AF) values on linear axis
x_fine = 10 ** y_fine               # back to AF for PDF evaluation

for ax, gene in zip(axes, ["APOE", "CYP2C19", "HLA-B"]):
    data     = overall[overall["Gene"] == gene]["AF"].values
    log_data = np.log10(data)        # transform data to log space

    # Histogram of log₁₀(AF) — linear bins in log space, density in log space
    ax.hist(log_data, bins=50, density=True, alpha=0.4,
            color=GENE_COLORS[gene], label="Observed AF")

    # Fit Beta MLE and plot in log space (Jacobian: f(x)*x*ln(10))
    try:
        a_hat, b_hat, loc, scale = stats.beta.fit(data, floc=0, fscale=1)
        beta_pdf = stats.beta.pdf(x_fine, a_hat, b_hat) * x_fine * np.log(10)
        ax.plot(y_fine, beta_pdf, color="black", linewidth=2,
                label=f"Beta MLE\n(α={a_hat:.3f}, β={b_hat:.3f})")
    except Exception as e:
        print(f"Beta fit failed for {gene}: {e}")

    # Fit Log-normal MLE and plot in log space
    try:
        shape, loc_ln, scale_ln = stats.lognorm.fit(data, floc=0)
        lognorm_pdf = stats.lognorm.pdf(x_fine, shape, loc=loc_ln, scale=scale_ln) * x_fine * np.log(10)
        ax.plot(y_fine, lognorm_pdf, color="orange", linewidth=2,
                linestyle="--", label="Log-normal MLE")
    except Exception as e:
        print(f"Log-normal fit failed for {gene}: {e}")

    # Format x-axis ticks to show actual AF values
    ax.set_xticks(range(-6, 1))
    ax.set_xticklabels([f'$10^{{{i}}}$' for i in range(-6, 1)], fontsize=8)
    ax.set_title(gene, fontsize=13, fontweight="bold")
    ax.set_xlabel("Allele Frequency")
    ax.set_ylabel("Density (log space)")
    ax.legend(fontsize=9)

plt.tight_layout()
plt.savefig(f"{PLOT_DIR}/02_distribution_candidates.png", dpi=150)
plt.close()
print(f"Saved → {PLOT_DIR}/02_distribution_candidates.png")

# ── Plot 3: Why Beta is theoretically motivated ────────────────────────────────
# Simulate Wright-Fisher stationary distribution and compare to Beta
fig, axes = plt.subplots(1, 3, figsize=(16, 5))
fig.suptitle("Wright-Fisher Stationary Distribution vs. Beta Distribution\n"
             "(theoretical motivation for using Beta to model AF)",
             fontsize=13, fontweight="bold")

np.random.seed(42)

for ax, (gene, Ne, mu) in zip(axes, [
    ("APOE",    10000, 1e-5),
    ("CYP2C19", 10000, 1e-5),
    ("HLA-B",   10000, 5e-5),   # higher mutation rate → more diversity
]):
    # Stationary Beta distribution: α = 4*Ne*mu, β = 4*Ne*mu
    # For neutral allele: both params equal 4*Ne*mu
    theta = 4 * Ne * mu
    x_plot = np.linspace(0.001, 0.999, 500)
    pdf = stats.beta.pdf(x_plot, theta, theta)

    ax.plot(x_plot, pdf, color=GENE_COLORS[gene], linewidth=2.5,
            label=f"Beta(θ={theta:.2f}, θ)\nNe={Ne}, μ={mu:.0e}")
    ax.fill_between(x_plot, pdf, alpha=0.2, color=GENE_COLORS[gene])
    ax.set_title(gene, fontsize=13, fontweight="bold")
    ax.set_xlabel("Allele Frequency")
    ax.set_ylabel("Probability Density")
    ax.legend(fontsize=10)

plt.tight_layout()
plt.savefig(f"{PLOT_DIR}/03_wright_fisher_beta.png", dpi=150)
plt.close()
print(f"Saved → {PLOT_DIR}/03_wright_fisher_beta.png")

print("\nDone.")