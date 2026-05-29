library(tidyverse)
library(gridExtra)
library(grid)


# ── Setup ──────────────────────────────────────────────────────────────────────
DATA_PATH   <- "data/processed/all_genes_tidy.csv"
PARAMS_PATH <- "data/output/beta_mle_by_population.csv"
PLOT_DIR    <- "data/plots/evaluation"
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

GENE_COLORS <- c("APOE" = "#E63946", "CYP2C19" = "#457B9D", "HLA-B" = "#2A9D8F")
POP_SHORT <- c(
  "African/African American" = "AFR", "Admixed American" = "AMR",
  "Ashkenazi Jewish" = "ASJ",         "East Asian" = "EAS",
  "European (Finnish)" = "FIN",       "Middle Eastern" = "MID",
  "European (non-Finnish)" = "NFE",   "Amish" = "AMI",
  "South Asian" = "SAS",              "Remaining" = "REM"
)
FILE_NUM <- c("APOE" = "07", "CYP2C19" = "08", "HLA-B" = "09")

df     <- read_csv(DATA_PATH)
params <- read_csv(PARAMS_PATH)
by_pop <- df %>% filter(Population != "Overall") %>% drop_na(AF) %>% filter(AF > 0)

for (gene_name in names(GENE_COLORS)) {
  panels <- list()
  
  for (pop_name in names(POP_SHORT)) {
    af <- by_pop %>% filter(Gene == gene_name, Population == pop_name) %>% pull(AF)
    af <- af[af > 0 & af < 1]
    
    p_row <- params %>% filter(Gene == gene_name, Population == pop_name)
    converged <- nrow(p_row) > 0 && p_row$converged[1]
    
    p <- ggplot()
    
    if (length(af) >= 3) {
      log_af <- log10(af)
      p <- p + geom_histogram(data = tibble(x = log_af), aes(x = x, y = ..density..), bins = 40, fill = GENE_COLORS[gene_name], alpha = 0.45)
      
      if (converged) {
        a <- p_row$alpha[1]; b <- p_row$beta[1]
        y_fine <- seq(min(log_af) - 0.2, 0, length.out = 1000)
        x_fine <- 10^y_fine
        pdf_log <- dbeta(x_fine, a, b) * x_fine * log(10)
        
        p <- p + geom_line(data = tibble(x = y_fine, y = pdf_log), aes(x = x, y = y), color = "black", linewidth = 0.8)
      }
      p <- p + scale_x_continuous(labels = function(x) parse(text = paste0("10^", round(x))))
    } else {
      p <- p + annotate("text", x = 0.5, y = 0.5, label = "insufficient\ndata", color = "grey", size = 3) + theme_void()
    }
    
    p <- p + labs(title = POP_SHORT[pop_name], x = "AF (log10)", y = "Density") + theme_minimal() +
      theme(plot.title = element_text(size = 10, face = "bold"), axis.text = element_text(size = 6))
    
    panels[[POP_SHORT[pop_name]]] <- p
  }
  
  png(file.path(PLOT_DIR, sprintf("%s_fit_overlay_%s.png", FILE_NUM[gene_name], gene_name)), width = 21, height = 7, units = "in", res = 150)
  grid.arrange(do.call(arrangeGrob, c(panels, ncol = 5)), 
               top = textGrob(sprintf("Observed AF vs Fitted Beta Distribution — %s\n(one panel per population, log10 space)", gene_name), gp = gpar(fontsize = 13, fontface = "bold")))
  dev.off()
}