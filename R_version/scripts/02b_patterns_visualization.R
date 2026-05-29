library(tidyverse)
library(scales)

# ── Setup ──────────────────────────────────────────────────────────────────────
DATA_PATH <- "data/processed/all_genes_tidy.csv"
PLOT_DIR  <- "data/plots/exploratory"

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

# ── Plot 4: ECDF ───────────────────────────────────────────────────────────────
p_ecdf <- ggplot(overall, aes(x = AF, color = Gene)) +
  stat_ecdf(geom = "line", linewidth = 1) +
  scale_x_log10() +
  scale_y_continuous(labels = percent) +
  scale_color_manual(values = GENE_COLORS) +
  geom_vline(xintercept = c(0.001, 0.01, 0.05), linetype = "dashed", color = "grey", alpha = 0.7) +
  annotate("text", x = c(0.0011, 0.011, 0.055), y = 0.05, label = c("0.1%", "1%", "5%"), color = "grey", size = 3, hjust = 0) +
  labs(title = "ECDF of Allele Frequency per Gene", 
       subtitle = "(shows what fraction of variants fall below each AF threshold)",
       x = "Allele Frequency (log scale)", y = "Cumulative Proportion of Variants") +
  theme_minimal()

ggsave(file.path(PLOT_DIR, "04_ecdf_AF.png"), plot = p_ecdf, width = 10, height = 6, dpi = 150)

# ── Plot 5: Heatmap Matrix ─────────────────────────────────────────────────────
matrix_data <- by_pop %>%
  filter(AF > 0) %>%
  group_by(Gene, Population) %>%
  summarise(median_log = log10(median(AF)), .groups = "drop") %>%
  mutate(Pop_Short = POP_SHORT[Population])

p_heat <- ggplot(matrix_data, aes(x = Pop_Short, y = Gene, fill = median_log)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.1f", median_log)), fontface = "bold", size = 3) +
  scale_fill_gradient2(low = "#E63946", mid = "#FFFFBF", high = "#2A9D8F", midpoint = -3) +
  labs(title = "Median log10(AF) per Gene and Population",
       subtitle = "(greener = higher AF, redder = lower AF)", x = "Population", fill = "Median log10(AF)") +
  theme_minimal()

ggsave(file.path(PLOT_DIR, "05_heatmap_median_AF.png"), plot = p_heat, width = 13, height = 4, dpi = 150)

# ── Plot 6: Rare variant proportion ────────────────────────────────────────────
p_rare <- by_pop %>%
  group_by(Gene, Population) %>%
  summarise(prop_rare = sum(AF < 0.01) / n() * 100, .groups = "drop") %>%
  mutate(Pop_Short = POP_SHORT[Population]) %>%
  ggplot(aes(x = Pop_Short, y = prop_rare, fill = Gene)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7, color = "white", alpha = 0.85) +
  scale_fill_manual(values = GENE_COLORS) +
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits = c(0, 105)) +
  geom_hline(yintercept = 100, linetype = "dashed", color = "grey", alpha = 0.5) +
  labs(title = "Proportion of Rare Variants (AF < 1%) by Population and Gene",
       x = "Population", y = "% of Variants with AF < 1%") +
  theme_minimal()

ggsave(file.path(PLOT_DIR, "06_rare_variant_proportion_by_population.png"), plot = p_rare, width = 13, height = 6, dpi = 150)