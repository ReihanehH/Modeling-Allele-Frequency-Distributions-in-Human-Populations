library(tidyverse)
# CRITICAL: We load 'dgof' because the base ks.test step fails or defaults warnings 
# when assessing against user-defined step functions or closures in standard ties.
library(dgof) 
library(gridExtra)
library(grid)

# ── Setup ──────────────────────────────────────────────────────────────────────
DATA_PATH   <- "data/processed/all_genes_tidy.csv"
PARAMS_PATH <- "data/output/beta_mle_by_population.csv"
PLOT_DIR    <- "data/plots/evaluation"
OUTPUT_DIR  <- "data/output"

GENES <- c("APOE", "CYP2C19", "HLA-B")
POP_SHORT <- c(
  "African/African American" = "AFR", "Admixed American" = "AMR",
  "Ashkenazi Jewish" = "ASJ",         "East Asian" = "EAS",
  "European (Finnish)" = "FIN",       "Middle Eastern" = "MID",
  "European (non-Finnish)" = "NFE",   "Amish" = "AMI",
  "South Asian" = "SAS",              "Remaining" = "REM"
)

df     <- read_csv(DATA_PATH)
params <- read_csv(PARAMS_PATH)
by_pop <- df %>% filter(Population != "Overall") %>% drop_na(AF) %>% filter(AF > 0)

message("============================================================")
message("KOLMOGOROV–SMIRNOV GOODNESS-OF-FIT TEST")
message("============================================================")

ks_records <- list()

converged_fits <- params %>% filter(converged == TRUE)

for (i in seq_len(nrow(converged_fits))) {
  row  <- converged_fits[i, ]
  gene <- row$Gene; pop <- row$Population
  a    <- row$alpha; b <- row$beta
  
  af <- by_pop %>% filter(Gene == gene, Population == pop) %>% pull(AF)
  af <- af[af > 0 & af < 1]
  
  # Perform KS Test against theoretical target
  ks_res <- dgof::ks.test(af, "pbeta", a, b)
  
  reject <- ks_res$p.value < 0.05
  
  ks_records[[length(ks_records) + 1]] <- tibble(
    Gene = gene, Population = pop, Pop_Short = row$Pop_Short, n = length(af),
    alpha = a, beta = b, KS_stat = round(ks_res$statistic, 6), p_value = ks_res$p.value, reject_H0 = reject
  )
  
  flag <- if(reject) "YES ***" else "NO"
  message(sprintf("%-8s %s | %5d | %8.4f | %.4e | %s", gene, row$Pop_Short, length(af), ks_res$statistic, ks_res$p.value, flag))
}

ks_df <- bind_rows(ks_records)
write_csv(ks_df, file.path(OUTPUT_DIR, "ks_test_results.csv"))

# Heatmap
p_ks_heat <- ggplot(ks_df, aes(x = Pop_Short, y = Gene, fill = KS_stat)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.3f%s", KS_stat, if_else(!reject_H0, "\n✓", ""))), fontface = "bold", size = 3) +
  scale_fill_gradient(low = "#2A9D8F", high = "#E63946") +
  labs(title = "Kolmogorov–Smirnov Statistic per Gene × Population", x = "Population", fill = "KS D-stat") +
  theme_minimal()

ggsave(file.path(PLOT_DIR, "10_ks_statistic_heatmap.png"), plot = p_ks_heat, width = 14, height = 4, dpi = 150)

# QQ Plots 
for (gene_name in GENES) {
  qq_panels <- list()
  
  for (pop_name in names(POP_SHORT)) {
    af <- by_pop %>% filter(Gene == gene_name, Population == pop_name) %>% pull(AF)
    af <- sort(af[af > 0 & af < 1])
    
    ks_row <- ks_df %>% filter(Gene == gene_name, Population == pop_name)
    
    if (nrow(ks_row) == 0 || length(af) < 3) {
      qq_panels[[POP_SHORT[pop_name]]] <- ggplot() + annotate("text", x = 0.5, y = 0.5, label = "insufficient\ndata") + theme_void()
      next
    }
    
    probs <- seq_along(af) / (length(af) + 1)
    theoretical <- qbeta(probs, ks_row$alpha[1], ks_row$beta[1])
    
    lbl_color <- if(ks_row$reject_H0[1]) "#C0392B" else "#27AE60"
    anno_text <- sprintf("D=%.3f\np=%.2e\n%s", ks_row$KS_stat[1], ks_row$p_value[1], if(ks_row$reject_H0[1]) "reject H₀" else "fail to reject H₀")
    
    p_qq <- tibble(th = theoretical, ob = af) %>%
      ggplot(aes(x = th, y = ob)) +
      geom_point(alpha = 0.5, size = 0.8, color = "#457B9D") +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
      annotate("label", x = 0, y = max(af), label = list(anno_text), parse=FALSE, hjust=0, vjust=1, color=lbl_color, size=2.5, fill="white", alpha=0.7) +
      labs(title = POP_SHORT[pop_name], x = "Theoretical", y = "Observed") +
      theme_minimal() + theme(axis.text = element_text(size = 7))
    
    qq_panels[[POP_SHORT[pop_name]]] <- p_qq
  }
  
  png(file.path(PLOT_DIR, sprintf("11_qq_plot_%s.png", gene_name)), width = 21, height = 8, units = "in", res = 150)
  grid.arrange(do.call(arrangeGrob, c(qq_panels, ncol = 5)), 
               top = textGrob(sprintf("QQ Plots: Observed vs Beta Quantiles — %s", gene_name), gp = gpar(fontsize = 13, fontface = "bold")))
  dev.off()
}