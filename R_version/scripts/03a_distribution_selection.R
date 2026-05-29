library(tidyverse)
library(fitdistrplus) # contains flexible standard ML fitting alternatives
library(gridExtra)
library(grid)

# ── Setup ──────────────────────────────────────────────────────────────────────
PLOT_DIR <- "data/plots/modeling"
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)
GENE_COLORS <- c("APOE" = "#E63946", "CYP2C19" = "#457B9D", "HLA-B" = "#2A9D8F")

# ── Plot 1: Beta Distribution Regimes ──────────────────────────────────────────
x_vals <- seq(0.001, 0.999, length.out = 1000)
regimes <- list(
  list(a=0.1, b=1.0, l="α=0.1, β=1.0  → strongly J-shaped (ultra-rare variants)", c="#E63946"),
  list(a=0.5, b=2.0, l="α=0.5, β=2.0  → J-shaped (rare-skewed)", c="#457B9D"),
  list(a=0.5, b=0.5, l="α=0.5, β=0.5  → U-shaped (two peaks)", c="#2A9D8F"),
  list(a=2.0, b=5.0, l="α=2.0, β=5.0  → bell-shaped, skewed right", c="#F4A261"),
  list(a=1.0, b=1.0, l="α=1.0, β=1.0  → Uniform (flat)", c="#6A0572")
)

beta_df <- map_df(regimes, function(r) {
  tibble(x = x_vals, pdf = dbeta(x_vals, r$a, r$b), label = r$l, color = r$c)
})

p1 <- ggplot(beta_df, aes(x = x, y = pdf, color = label)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = setNames(map_chr(regimes, ~.$c), map_chr(regimes, ~.$l))) +
  ylim(0, 10) +
  geom_vline(xintercept = 0.01, linetype = "dashed", color = "grey") +
  annotate("text", x = 0.012, y = 9, label = "AF = 1%\n(rare threshold)", color = "grey", size = 3, hjust = 0) +
  labs(title = "Beta Distribution PDF for Different Parameter Regimes", 
       subtitle = "Beta(α, β) is defined on [0,1] — naturally suited for allele frequencies",
       x = "Allele Frequency", y = "Probability Density", color = "") +
  theme_minimal() + theme(legend.position = "bottom", legend.direction = "vertical")

ggsave(file.path(PLOT_DIR, "01_beta_distribution_shapes.png"), plot = p1, width = 11, height = 6, dpi = 150)

# ── Plot 2 & 3 Dependencies ───────────────────────────────────────────────────
df      <- read_csv("data/processed/all_genes_tidy.csv")
overall <- df %>% filter(Population == "Overall") %>% drop_na(AF) %>% filter(AF > 0)

y_fine  <- seq(-6, 0, length.out = 2000)
x_fine  <- 10^y_fine

# Custom panels generation
library(gridExtra)
plots_cands <- list()

for (gene_name in names(GENE_COLORS)) {
  sub_data <- overall %>% filter(Gene == gene_name) %>% pull(AF)
  log_data <- log10(sub_data)
  
  # Fit Beta via fitdistrplus
  fit_b <- tryCatch(fitdist(sub_data, "beta", method="mle"), error = function(e) NULL)
  # Fit Lognormal
  fit_ln <- tryCatch(fitdist(sub_data, "lnorm", method="mle"), error = function(e) NULL)
  
  p_frame <- tibble(y = y_fine, x = x_fine)
  
  if (!is.null(fit_b)) {
    p_frame$beta_pdf <- dbeta(x_fine, fit_b$estimate["shape1"], fit_b$estimate["shape2"]) * x_fine * log(10)
  }
  if (!is.null(fit_ln)) {
    p_frame$ln_pdf <- dlnorm(x_fine, fit_ln$estimate["meanlog"], fit_ln$estimate["sdlog"]) * x_fine * log(10)
  }
  
  p_curr <- ggplot() +
    geom_histogram(data = tibble(y = log_data), aes(x = y, y = ..density..), bins = 50, fill = GENE_COLORS[gene_name], alpha = 0.4)
  
  if ("beta_pdf" %names_in% p_frame) {
    lbl <- sprintf("Beta MLE\n(α=%.3f, β=%.3f)", fit_b$estimate["shape1"], fit_b$estimate["shape2"])
    p_curr <- p_curr + geom_line(data = p_frame, aes(x = y, y = beta_pdf, color = "Beta"), linewidth = 1)
  }
  if ("ln_pdf" %names_in% p_frame) {
    p_curr <- p_curr + geom_line(data = p_frame, aes(x = y, y = ln_pdf, color = "Lognormal"), linetype = "dashed", linewidth = 1)
  }
  
  p_curr <- p_curr +
    scale_x_continuous(breaks = -6:0, labels = parse(text = paste0("10^", -6:0))) +
    scale_color_manual(values = c("Beta" = "black", "Lognormal" = "orange")) +
    labs(title = gene_name, x = "Allele Frequency", y = "Density (log space)", color = "") +
    theme_minimal()
  
  plots_cands[[gene_name]] <- p_curr
}
`%names_in%` <- function(x, y) x %in% names(y)

png(file.path(PLOT_DIR, "02_distribution_candidates.png"), width = 16, height = 5, units = "in", res = 150)
grid.arrange(do.call(arrangeGrob, c(plots_cands, ncol = 3)), 
             top = textGrob("Observed AF Histogram vs. Candidate Distributions\n(working in log10 space for proper density comparison)", gp = gpar(fontsize = 13, fontface = "bold")))
dev.off()

# ── Plot 3: Wright-Fisher Stationary Distribution ──────────────────────────────
wf_params <- list(
  list(g="APOE", Ne=10000, mu=1e-5),
  list(g="CYP2C19", Ne=10000, mu=1e-5),
  list(g="HLA-B", Ne=10000, mu=5e-5)
)
x_wf <- seq(0.001, 0.999, length.out = 500)

wf_plots <- map(wf_params, function(p) {
  theta <- 4 * p$Ne * p$mu
  wf_df <- tibble(x = x_wf, pdf = dbeta(x_wf, theta, theta))
  
  ggplot(wf_df, aes(x = x, y = pdf)) +
    geom_line(color = GENE_COLORS[p$g], linewidth = 1) +
    geom_area(fill = GENE_COLORS[p$g], alpha = 0.2) +
    labs(title = p$g, subtitle = sprintf("Beta(θ=%.2f, θ)\nNe=%d, μ=%.0e", theta, p$Ne, p$mu),
         x = "Allele Frequency", y = "Probability Density") +
    theme_minimal()
})

png(file.path(PLOT_DIR, "03_wright_fisher_beta.png"), width = 16, height = 5, units = "in", res = 150)
grid.arrange(do.call(arrangeGrob, c(wf_plots, ncol = 3)), 
             top = textGrob("Wright-Fisher Stationary Distribution vs. Beta Distribution\n(theoretical motivation for using Beta to model AF)", gp = gpar(fontsize = 13, fontface = "bold")))
dev.off()