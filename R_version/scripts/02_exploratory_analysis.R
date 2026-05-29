library(tidyverse)
library(gridExtra)
library(grid)

# ── Setup ──────────────────────────────────────────────────────────────────────
DATA_PATH  <- "data/processed/all_genes_tidy.csv"
PLOT_DIR   <- "data/plots/exploratory"
OUTPUT_DIR <- "data/output"
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

GENE_COLORS <- c("APOE" = "#E63946", "CYP2C19" = "#457B9D", "HLA-B" = "#2A9D8F")
POP_SHORT <- c(
  "African/African American" = "AFR", "Admixed American" = "AMR",
  "Ashkenazi Jewish" = "ASJ",         "East Asian" = "EAS",
  "European (Finnish)" = "FIN",       "Middle Eastern" = "MID",
  "European (non-Finnish)" = "NFE",   "Amish" = "AMI",
  "South Asian" = "SAS",              "Remaining" = "REM"
)

df      <- read_csv(DATA_PATH)
overall <- df %>% filter(Population == "Overall") %>% drop_na(AF)
by_pop  <- df %>% filter(Population != "Overall") %>% drop_na(AF)

# ── 1. Summary Statistics ──────────────────────────────────────────────────────
summary_overall <- overall %>%
  group_by(Gene) %>%
  summarise(
    count = n(), mean = mean(AF), std = sd(AF), min = min(AF),
    `1%` = quantile(AF, 0.01), `5%` = quantile(AF, 0.05), `25%` = quantile(AF, 0.25),
    `50%` = quantile(AF, 0.50), `75%` = quantile(AF, 0.75), `95%` = quantile(AF, 0.95),
    `99%` = quantile(AF, 0.99), max = max(AF)
  )
write_csv(summary_overall, file.path(OUTPUT_DIR, "summary_statistics_overall.csv"))

summary_pop <- by_pop %>%
  group_by(Gene, Population) %>%
  summarise(
    count = n(), mean = mean(AF), std = sd(AF), min = min(AF),
    `25%` = quantile(AF, 0.25), `50%` = quantile(AF, 0.50), `75%` = quantile(AF, 0.75), max = max(AF)
  )
write_csv(summary_pop, file.path(OUTPUT_DIR, "summary_statistics_by_population.csv"))

message("\n============================================================")
message("RARE VARIANT PREVALENCE (AF < 0.01)")
message("============================================================")
overall %>%
  group_by(Gene) %>%
  summarise(rare = sum(AF < 0.01), total = n()) %>%
  purrr::pwalk(function(Gene, rare, total) {
    cat(sprintf("  %-8s: %4d / %d variants are rare  (%.1f%%)\n", Gene, rare, total, 100 * rare / total))
  })

# ── 2. Plot 1: Linear overview + Log Histograms ─────────────────────────────────
library(gridExtra)
p_lin <- ggplot(overall, aes(x = AF, fill = Gene)) +
  geom_histogram(bins = 80, alpha = 0.6, position = "identity") +
  scale_fill_manual(values = GENE_COLORS) +
  labs(title = "Linear Scale — all three genes", x = "Allele Frequency", y = "Number of Variants") +
  theme_minimal()

log_plots <- overall %>%
  filter(AF > 0) %>%
  group_split(Gene) %>%
  map(function(sub) {
    gene_name <- unique(sub$Gene)
    ggplot(sub, aes(x = AF)) +
      geom_histogram(bins = 50, fill = GENE_COLORS[gene_name], color = "white", alpha = 0.85) +
      scale_x_log10() +
      labs(title = paste(gene_name, "— Log Scale"), x = "Allele Frequency (log scale)", y = "Number of Variants") +
      theme_minimal() +
      theme(plot.title = element_text(face = "bold"))
  })

png(file.path(PLOT_DIR, "01_overall_AF_distribution.png"), width = 16, height = 10, units = "in", res = 150)
grid.arrange(p_lin, do.call(arrangeGrob, c(log_plots, ncol = 3)), nrow = 2, 
             top = textGrob("Overall AF Distribution per Gene", gp = gpar(fontsize = 14, fontface = "bold")))
dev.off()

# ── 3. Plot 2: Boxplots by Population ───────────────────────────────────────────
box_plots <- by_pop %>%
  filter(AF > 0) %>%
  group_split(Gene) %>%
  map(function(sub) {
    gene_name <- unique(sub$Gene)
    
    # Calculate median sorting logic matching your Python workflow
    pop_order <- sub %>%
      group_by(Population) %>%
      summarise(med = median(AF)) %>%
      arrange(desc(med)) %>%
      pull(Population)
    
    sub_ordered <- sub %>% 
      mutate(Population = factor(Population, levels = pop_order))
    
    ggplot(sub_ordered, aes(x = Population, y = log10(AF))) +
      geom_boxplot(fill = GENE_COLORS[gene_name], alpha = 0.6, outlier.size = 0.5) +
      scale_x_discrete(labels = POP_SHORT) +
      labs(title = gene_name, x = "Population", y = "log10(AF)") +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1), # <-- Standard clean R text rotation
        plot.title = element_text(face = "bold")
      )
  })

png(file.path(PLOT_DIR, "02_AF_by_population.png"), width = 18, height = 6, units = "in", res = 150)
grid.arrange(do.call(arrangeGrob, c(box_plots, ncol = 3)), 
             top = textGrob("AF Distribution by Population and Gene (Log Scale)", gp = gpar(fontsize = 14, fontface = "bold")))
dev.off()

# ── 4. Plot 3: Frequency Class Breakdown ────────────────────────────────────────
bins_classes <- c(0, 0.001, 0.01, 0.05, 0.5, 1.0)
labels_classes <- c("<0.1%", "0.1–1%", "1–5%", "5–50%", ">50%")

class_plots <- overall %>%
  filter(AF > 0) %>%
  mutate(Class = cut(AF, breaks = bins_classes, labels = labels_classes, include.lowest = TRUE)) %>%
  group_split(Gene) %>%
  map(function(sub) {
    gene_name <- unique(sub$Gene)
    counts <- sub %>% group_by(Class) %>% tally()
    
    ggplot(counts, aes(x = Class, y = n)) +
      geom_bar(stat = "identity", fill = GENE_COLORS[gene_name], alpha = 0.8, color = "white") +
      geom_text(aes(label = n), vjust = -0.5, fontface = "bold", size = 3) +
      scale_y_log10(limits = c(0.5, max(counts$n) * 5)) +
      labs(title = gene_name, x = "Frequency Class", y = "Number of Variants (log scale)") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 30, hjust = 1), plot.title = element_text(face = "bold"))
  })

png(file.path(PLOT_DIR, "03_frequency_class_breakdown.png"), width = 15, height = 5, units = "in", res = 150)
grid.arrange(do.call(arrangeGrob, c(class_plots, ncol = 3)), 
             top = textGrob("Variant Frequency Class Breakdown per Gene", gp = gpar(fontsize = 14, fontface = "bold")))
dev.off()