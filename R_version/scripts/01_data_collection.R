library(tidyverse)

# ── Setup ──────────────────────────────────────────────────────────────────────
RAW_DIR  <- "data/raw"
OUT_DIR  <- "data/processed"
dir.create(RAW_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

FILES <- list(
  "APOE"    = "data/raw/gnomad_APOE.csv",
  "CYP2C19" = "data/raw/gnomad_CYP2C19.csv",
  "HLA-B"   = "data/raw/gnomad_HLA-B.csv"
)

POPULATIONS <- c(
  "African/African American", "Admixed American", "Ashkenazi Jewish",
  "East Asian", "European (Finnish)", "Middle Eastern",
  "European (non-Finnish)", "Amish", "South Asian", "Remaining"
)

tidy_gene <- function(gene, filepath) {
  df <- read_csv(filepath, col_types = cols(.default = "c")) # read as character to mimic low_memory=False flexibility
  
  # Ensure necessary numeric conversions
  df <- df %>% mutate(
    across(c(`Allele Count`, `Allele Number`, `Allele Frequency`, `Homozygote Count`), as.numeric)
  )
  
  base_cols <- c("gnomAD ID", "Chromosome", "Position", "Reference", "Alternate", "VEP Annotation", "Gene")
  df$Gene <- gene
  
  # 1. Overall rows
  overall_df <- df %>% 
    select(any_of(c(base_cols, "Allele Count", "Allele Number", "Allele Frequency", "Homozygote Count"))) %>%
    rename(AC = `Allele Count`, AN = `Allele Number`, AF = `Allele Frequency`) %>%
    mutate(Population = "Overall")
  
  # 2. Per-Population rows
  pop_list <- list()
  for (pop in POPULATIONS) {
    ac_col <- paste0("Allele Count ", pop)
    an_col <- paste0("Allele Number ", pop)
    hom_col <- paste0("Homozygote Count ", pop)
    
    pop_df <- df %>%
      select(any_of(base_cols), all_of(c(ac_col, an_col))) %>%
      mutate(
        AC = as.numeric(.[[ac_col]]),
        AN = as.numeric(.[[an_col]]),
        Homozygote_Count = if(hom_col %in% names(df)) as.numeric(df[[hom_col]]) else NA_real_,
        Population = pop,
        AF = if_else(!is.na(AN) & AN > 0 & !is.na(AC), AC / AN, NA_real_)
      ) %>%
      select(-all_of(c(ac_col, an_col)))
    
    pop_list[[pop]] <- pop_df
  }
  
  combined_gene <- bind_rows(overall_df, bind_rows(pop_list)) %>%
    filter(!is.na(AC) | !is.na(AN))
  
  return(combined_gene)
}

main <- function() {
  all_dfs <- list()
  
  for (gene in names(FILES)) {
    message("Processing ", gene, "...")
    tidy_df <- tidy_gene(gene, FILES[[gene]])
    
    out_path <- file.path(OUT_DIR, paste0(gene, "_tidy.csv"))
    write_csv(tidy_df, out_path)
    message("  ", nrow(tidy_df), " rows -> ", out_path)
    
    all_dfs[[gene]] <- tidy_df
  }
  
  combined <- bind_rows(all_dfs)
  write_csv(combined, file.path(OUT_DIR, "all_genes_tidy.csv"))
  message("\nCombined: ", nrow(combined), " rows -> ", file.path(OUT_DIR, "all_genes_tidy.csv"))
  
  message("\n── Variant counts by gene and population ──")
  summary_counts <- combined %>%
    filter(Population != "Overall") %>%
    group_by(Gene, Population) %>%
    tally() %>%
    pivot_wider(names_from = Population, values_from = n, values_fill = 0)
  
  print(as.data.frame(summary_counts))
}

main()