library(tidyverse)
library(fitdistrplus)

# ── Setup ──────────────────────────────────────────────────────────────────────
DATA_PATH  <- "data/processed/all_genes_tidy.csv"
PLOT_DIR   <- "data/plots/modeling"
OUTPUT_DIR <- "data/output"

GENES <- c("APOE", "CYP2C19", "HLA-B")
GENE_COLORS <- c("APOE" = "#E63946", "CYP2C19" = "#457B9D", "HLA-B" = "#2A9D8F")
POP_SHORT <- c(
  "African/African American" = "AFR", "Admixed American" = "AMR",
  "Ashkenazi Jewish" = "ASJ",         "East Asian" = "EAS",
  "European (Finnish)" = "FIN",       "Middle Eastern" = "MID",
  "European (non-Finnish)" = "NFE",   "Amish" = "AMI",
  "South Asian" = "SAS",              "Remaining" = "REM"
)
POPULATIONS <- names(POP_SHORT)
MIN_VARIANTS <- 10

df <- read_csv(DATA_PATH)
by_pop <- df %>% filter(Population != "Overall") %>% drop_na(AF) %>% filter(AF > 0)

message("============================================================")
message("BETA MLE FITTING — per gene × population")
message("============================================================")

records <- list()

for (gene in GENES) {
  for (pop in POPULATIONS) {
    af <- by_pop %>% filter(Gene == gene, Population == pop) %>% pull(AF)
    af <- af[af > 0 & af < 1]
    n <- length(af)
    
    if (n < MIN_VARIANTS) {
      alpha <- NA; beta_par <- NA; ok <- FALSE
    } else {
      fit <- tryCatch(fitdist(af, "beta", method="mle"), error = function(e) NULL)
      if (!is.null(fit)) {
        alpha <- fit$estimate["shape1"]
        beta_par <- fit$estimate["shape2"]
        ok <- TRUE
      } else {
        alpha <- NA; beta_par <- NA; ok <- FALSE
      }
    }
    
    status_msg <- if(ok) sprintf("α=%.4f, β=%.4f", alpha, beta_par) else "FAILED (too few variants)"
    message(sprintf("  %-8s | %s | n=%5d | %s", gene, POP_SHORT[pop], n, status_msg))
    
    records[[length(records) + 1]] <- tibble(
      Gene = gene, Population = pop, Pop_Short = POP_SHORT[pop], n_variants = n,
      mean_AF = if(n > 0) mean(af) else NA, alpha = alpha, beta = beta_par, converged = ok
    )
  }
}

results <- bind_rows(records)
write_csv(results, file.path(OUTPUT_DIR, "beta_mle_by_population.csv"))

# Heatmap function
param_heatmap <- function(res, param, title, low_c, high_c, filename) {
  mat_data <- res %>% filter(converged == TRUE) %>% rename(val = !!sym(param))
  
  p <- ggplot(mat_data, aes(x = Pop_Short, y = Gene, fill = log10(val))) +
    geom_tile() +
    geom_text(aes(label = if_else(val < 10, sprintf("%.3f", val), sprintf("%.1f", val))), fontface="bold", size=3) +
    scale_fill_gradient(low = low_c, high = high_c, labels = function(x) sprintf("%.3g", 10^x)) +
    labs(title = title, x = "Population", fill = paste0("log10(", param, ")")) +
    theme_minimal()
  
  ggsave(file.path(PLOT_DIR, filename), plot = p, width = 14, height = 4, dpi = 150)
}

param_heatmap(results, "alpha", "Fitted Beta α per Gene × Population\n(α < 1 → J-shaped; α > 1 → bell-shaped)", "#FFFFE5", "#FE9929", "04_beta_alpha_heatmap.png")
param_heatmap(results, "beta", "Fitted Beta β per Gene × Population\n(higher β → density concentrated near zero)", "#F7FBFF", "#08519C", "05_beta_beta_heatmap.png")

# Scatter Plot
p_scatter <- results %>% 
  filter(converged == TRUE) %>%
  ggplot(aes(x = alpha, y = beta, color = Gene)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_text(aes(label = Pop_Short), hjust = -0.2, vjust = -0.2, size = 3, color = "black") +
  scale_x_log10() + scale_y_log10() +
  scale_color_manual(values = GENE_COLORS) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey") +
  facet_wrap(~Gene, scales = "free") +
  labs(title = "Fitted Beta Parameters: α vs β per Gene and Population", x = "α (shape)", y = "β (tail weight)") +
  theme_minimal() + theme(legend.position = "none")

ggsave(file.path(PLOT_DIR, "06_alpha_vs_beta_scatter.png"), plot = p_scatter, width = 16, height = 5, dpi = 150)