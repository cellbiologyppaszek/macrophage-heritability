# =============================================================================
# scRNA-seq CV fold-change analysis: clonal vs mixed populations
#
# For each gene:
#   1. Import logcount matrices for JM01–JM12
#   2. Calculate mean expression per gene within each sample
#   3. Keep genes with mean expression > 1 in at least one clonal sample
#   4. Calculate CV across clonal sample means and mixed sample means
#   5. Calculate FCcv = CV_clonal / CV_mixed
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(readxl)
  library(writexl)
})

# ---- Set directory -----------------------------------------------------------
project_dir <- getwd()
data_dir <- file.path(project_dir, "data")
results_dir <- file.path(project_dir, "results")
output_file <- file.path(results_dir, "scRNAseq_CV_fold_change_clonal_vs_mixed.xlsx")


files <- list.files(
  path = data_dir,
  pattern = "^combined_logcounts_JM[0-9]{2}\\.xlsx$",
  full.names = TRUE
)

sample_names <- sprintf("JM%02d", 1:12)

logcount_list <- files %>%
  set_names(
    basename(.) %>%
      str_extract("JM[0-9]{2}")
  ) %>%
  map(~ read_excel(.x, sheet = 1))

logcount_list <- logcount_list[sample_names]

# ---- Define populations ------------------------------------------------------

clonal_samples <- sprintf("JM%02d", 1:8)
mixed_samples  <- sprintf("JM%02d", 9:12)

# ---- Calculate mean expression per gene per sample ---------------------------

get_sample_gene_means <- function(df, sample_name) {
  
  expr_cols <- 3:ncol(df)   # assumes columns 1–2 are gene metadata
  
  df %>%
    mutate(
      mean_sample = rowMeans(across(all_of(expr_cols)), na.rm = TRUE)
    ) %>%
    transmute(
      gene = Symbol,
      sample = sample_name,
      mean_sample = mean_sample
    )
}

sample_gene_means <- imap_dfr(
  logcount_list,
  ~ get_sample_gene_means(.x, .y)
) %>%
  mutate(
    population = case_when(
      sample %in% clonal_samples ~ "clonal",
      sample %in% mixed_samples  ~ "mixed",
      TRUE ~ NA_character_
    )
  )

# ---- Filter expressed genes --------------------------------------------------
# Keep genes with mean expression > 1 in at least one clonal sample

expressed_genes <- sample_gene_means %>%
  filter(population == "clonal") %>%
  group_by(gene) %>%
  summarise(
    expressed_in_clonal = any(mean_sample > 1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(expressed_in_clonal) %>%
  pull(gene)

sample_gene_means_filtered <- sample_gene_means %>%
  filter(gene %in% expressed_genes)

# ---- CV helper ---------------------------------------------------------------

safe_cv <- function(x) {
  m <- mean(x, na.rm = TRUE)
  s <- sd(x, na.rm = TRUE)
  
  if (is.na(m) || m == 0) {
    return(NA_real_)
  }
  
  s / m
}

# ---- CV fold-change analysis -------------------------------------------------

cv_fc_results <- sample_gene_means_filtered %>%
  group_by(gene, population) %>%
  summarise(
    mean_of_sample_means = mean(mean_sample, na.rm = TRUE),
    sd_between_samples   = sd(mean_sample, na.rm = TRUE),
    cv                   = safe_cv(mean_sample),
    n_samples            = n(),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = population,
    values_from = c(
      mean_of_sample_means,
      sd_between_samples,
      cv,
      n_samples
    )
  ) %>%
  mutate(
    FCcv = cv_clonal / cv_mixed,
    heritable = FCcv >= 2.5
  ) %>%
  arrange(desc(FCcv))

# ---- Save outputs ------------------------------------------------------------

write_xlsx(
  cv_fc_results,
  path = output_file
)

# ---- Quick checks ------------------------------------------------------------

n_expressed_genes <- length(expressed_genes)
n_final_genes <- nrow(cv_fc_results)

print(paste("Number of expressed genes:", n_expressed_genes))
print(paste("Number of genes in CV analysis:", n_final_genes))

head(cv_fc_results, 20)
