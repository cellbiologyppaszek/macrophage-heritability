#!/usr/bin/env Rscript

# =============================================================================
# CV fold-change analysis for MemorySeq hereditary gene sets
# =============================================================================
#
# Purpose:
#   Calculate the coefficient of variation (CV) for each gene clonal and mixed 
#   populations and report the CV fold change between stimulated clonal and
#   stimulated mixed samples:
#
#       CV fold change = CV_stimulated_clonal / CV_stimulated_mixed
#     
#   * CV is computed on log2(normalised read count + 1).
#
#   Genes with CV fold change >= 2.5 are classified as heritable (= 1), otherwise
#   non-heritable (= 0).
#
# Input workbook:
#   data/MemorySeq_data.xlsx
#
# Expected sheets:
#   1. JM_20210805_minus_JM15_JM18
#      - contains gene-level expression values
#      - contains a gene_name column
#      - sample columns follow the naming pattern:
#          1SC.y, 2SC.y, ..., 43SC.y    = stimulated clonal
#          44SM.y, ..., 86SM.y          = stimulated mixed
#          87UM.y, ..., 96UM.y          = unstimulated mixed
#
#   2. regulated genes heretability
#      - contains gene lists in the columns:
#          upregulated genes
#          downregulated genes
#          TLR-independent genes
#
# Output workbook:
#   results/CV_fold_change_heritability.xlsx
#
# Output sheets:
#   1. up_regulated
#   2. downregulated
#   3. TLR_independent
#

# =============================================================================


# ---- 1. Package checks -------------------------------------------------------

required_packages <- c(
  "readxl",
  "dplyr",
  "stringr",
  "tibble",
  "writexl"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "The following R packages are required but not installed: ",
    paste(missing_packages, collapse = ", "),
    "\nInstall them before running the script. Do not install packages inside ",
    "the analysis script submitted with the manuscript."
  )
}

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(writexl)
})


# ---- 2. User-editable configuration -----------------------------------------

project_dir <- getwd()
data_dir <- file.path(project_dir, "data")
results_dir <- file.path(project_dir, "results")

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

input_file <- file.path(data_dir, "MemorySeq data.xlsx")

expression_sheet <- "JM_20210805_minus_JM15_JM18"
gene_list_sheet  <- "genes list"

output_file <- file.path(results_dir, "CV_fold_change_heritability_MemorySeq.xlsx")

heritability_threshold <- 2.5
pseudocount <- 1

# Gene-list columns in the `genes list` sheet.
upregulated_gene_column   <- "upregulated genes"
downregulated_gene_column <- "downregulated genes"
tlr_independent_column    <- "TLR-independent genes"


# ---- 3. Helper functions -----------------------------------------------------

check_file_exists <- function(path) {
  if (!file.exists(path)) {
    stop("Required file not found: ", path)
  }
  invisible(path)
}

normalise_gene_names <- function(x) {
  # Some gene names in the workbook are stored as strings such as "'Cxcl2'".
  # This removes leading/trailing spaces and leading/trailing apostrophes.
  x <- as.character(x)
  x <- stringr::str_trim(x)
  x <- stringr::str_replace_all(x, "^'+|'+$", "")
  x[x == "" | is.na(x)] <- NA_character_
  x
}

find_first_matching_column <- function(column_names, pattern, preferred = NULL) {
  if (!is.null(preferred) && preferred %in% column_names) {
    return(preferred)
  }

  matches <- column_names[stringr::str_detect(column_names, pattern)]

  if (length(matches) == 0) {
    stop(
      "Could not find a column matching pattern: ", pattern,
      "\nAvailable columns are:\n",
      paste(column_names, collapse = ", ")
    )
  }

  matches[1]
}

read_gene_list <- function(df, column_name) {
  if (!column_name %in% colnames(df)) {
    stop(
      "Column `", column_name, "` was not found in the gene-list sheet.\n",
      "Available columns are:\n",
      paste(colnames(df), collapse = ", ")
    )
  }

  genes <- normalise_gene_names(df[[column_name]])
  unique(stats::na.omit(genes))
}

detect_sample_columns <- function(expression_df) {
  column_info <- tibble::tibble(column_name = colnames(expression_df)) %>%
    mutate(
      sample_id = suppressWarnings(
        as.integer(stringr::str_match(column_name, "^(\\d+)(SC|SM|UM)\\.y$")[, 2])
      ),
      sample_group = stringr::str_match(column_name, "^(\\d+)(SC|SM|UM)\\.y$")[, 3],
      condition = dplyr::case_when(
        sample_group == "SC" & sample_id >= 1  & sample_id <= 43 ~ "stimulated_clonal",
        sample_group == "SM" & sample_id >= 44 & sample_id <= 86 ~ "stimulated_mixed",
        sample_group == "UM" & sample_id >= 87 & sample_id <= 96 ~ "unstimulated_mixed",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(condition))

  required_conditions <- c(
    "stimulated_clonal",
    "stimulated_mixed",
    "unstimulated_mixed"
  )

  missing_conditions <- setdiff(required_conditions, unique(column_info$condition))

  if (length(missing_conditions) > 0) {
    stop(
      "No sample columns were detected for the following condition(s): ",
      paste(missing_conditions, collapse = ", "),
      "\nExpected column names like 1SC.y, 44SM.y and 87UM.y."
    )
  }

  column_info
}

as_numeric_matrix <- function(df) {
  mat <- as.matrix(df)
  suppressWarnings(storage.mode(mat) <- "numeric")
  mat
}

row_n_nonmissing <- function(mat) {
  rowSums(!is.na(mat))
}


row_mean_log2_safe <- function(mat, pseudocount = 1) {
  apply(mat, 1, function(x) {
    x <- as.numeric(x)
    x <- x[!is.na(x)]
    
    if (length(x) == 0) {
      return(NA_real_)
    }
    
    mean(log2(x + pseudocount))
  })
}

row_sd_log2_safe <- function(mat) {
  apply(mat, 1, function(x) {
    x <- x[!is.na(x)]
    if (length(x) < 2) {
      return(NA_real_)
    }
    stats::sd(log2(x + pseudocount))
  })
}


row_cv_safe <- function(mat) {
  mu <- row_mean_log2_safe(mat)
  sigma <- row_sd_log2_safe(mat)

  cv <- sigma / mu
  cv[is.na(mu) | mu == 0] <- NA_real_
  cv
}

calculate_cv_metrics <- function(expression_df, sample_column_info) {
  stim_clonal_cols <- sample_column_info %>%
    filter(condition == "stimulated_clonal") %>%
    pull(column_name)

  stim_mixed_cols <- sample_column_info %>%
    filter(condition == "stimulated_mixed") %>%
    pull(column_name)

  unstim_mixed_cols <- sample_column_info %>%
    filter(condition == "unstimulated_mixed") %>%
    pull(column_name)

  stim_clonal_mat <- as_numeric_matrix(expression_df[, stim_clonal_cols, drop = FALSE])
  stim_mixed_mat <- as_numeric_matrix(expression_df[, stim_mixed_cols, drop = FALSE])
  unstim_mixed_mat <- as_numeric_matrix(expression_df[, unstim_mixed_cols, drop = FALSE])

  cv_stim_clonal <- row_cv_safe(stim_clonal_mat)
  cv_stim_mixed <- row_cv_safe(stim_mixed_mat)
  cv_unstim_mixed <- row_cv_safe(unstim_mixed_mat)

  cv_fold_change <- cv_stim_clonal / cv_stim_mixed
  cv_fold_change[is.infinite(cv_fold_change)] <- NA_real_

  tibble::tibble(
    mean_stimulated_clonal_log = row_mean_log2_safe (stim_clonal_mat, pseudocount),
    mean_stimulated_mixed_log = row_mean_log2_safe (stim_mixed_mat, pseudocount),
    mean_unstimulated_mixed_log = row_mean_log2_safe (unstim_mixed_mat, pseudocount),


    sd_stimulated_clonal = row_sd_log2_safe(stim_clonal_mat),
    sd_stimulated_mixed = row_sd_log2_safe(stim_mixed_mat),
    sd_unstimulated_mixed = row_sd_log2_safe(unstim_mixed_mat),

    cv_stimulated_clonal = cv_stim_clonal,
    cv_stimulated_mixed = cv_stim_mixed,
    cv_unstimulated_mixed = cv_unstim_mixed,

    cv_fold_change_clonal_over_mixed = cv_fold_change,
    heritability = dplyr::if_else(
      !is.na(cv_fold_change) & cv_fold_change >= heritability_threshold,
      1L,
      0L
    ),

    n_stimulated_clonal = row_n_nonmissing(stim_clonal_mat),
    n_stimulated_mixed = row_n_nonmissing(stim_mixed_mat),
    n_unstimulated_mixed = row_n_nonmissing(unstim_mixed_mat)
  )
}

make_output_sheet <- function(metrics_df, gene_list, sheet_label) {
  missing_genes <- setdiff(gene_list, metrics_df$gene_name)

  if (length(missing_genes) > 0) {
    writeLines(
      missing_genes,
      file.path(results_dir, paste0("missing_genes_", sheet_label, ".txt"))
    )
  }

  metrics_df %>%
    filter(gene_name %in% gene_list) %>%
    arrange(desc(cv_fold_change_clonal_over_mixed), gene_name)
}


# ---- 4. Read input workbook --------------------------------------------------

check_file_exists(input_file)

expression_df <- readxl::read_excel(
  input_file,
  sheet = expression_sheet,
  .name_repair = "unique"
)

gene_lists_df <- readxl::read_excel(
  input_file,
  sheet = gene_list_sheet,
  .name_repair = "unique"
)


# ---- 5. Identify gene and sample columns ------------------------------------

# There are two gene_name columns in the workbook after automatic name repair.
# The first gene_name column corresponds to the main gene identifier used here.
gene_name_column <- find_first_matching_column(
  colnames(expression_df),
  pattern = "^gene_name(\\.\\.\\.[0-9]+)?$",
  preferred = "gene_name...10"
)

# Optional description column, retained if present.
gene_description_column <- find_first_matching_column(
  colnames(expression_df),
  pattern = "^gene_description(\\.\\.\\.[0-9]+)?$",
  preferred = "gene_description"
)

expression_df <- expression_df %>%
  mutate(gene_name = normalise_gene_names(.data[[gene_name_column]]))

sample_column_info <- detect_sample_columns(expression_df)

message("Detected sample columns:")
message("  stimulated clonal:   ", sum(sample_column_info$condition == "stimulated_clonal"))
message("  stimulated mixed:    ", sum(sample_column_info$condition == "stimulated_mixed"))
message("  unstimulated mixed:  ", sum(sample_column_info$condition == "unstimulated_mixed"))


# ---- 6. Calculate CV metrics -------------------------------------------------

cv_metrics <- calculate_cv_metrics(expression_df, sample_column_info)

metrics_df <- bind_cols(
  expression_df %>%
    transmute(
      gene_name = gene_name,
      gene_description = .data[[gene_description_column]]
    ),
  cv_metrics
) %>%
  filter(!is.na(gene_name))

# If duplicate gene names are present, keep the first occurrence and report them.
duplicate_genes <- metrics_df$gene_name[duplicated(metrics_df$gene_name)]

if (length(duplicate_genes) > 0) {
  writeLines(
    unique(duplicate_genes),
    file.path(results_dir, "duplicate_gene_names_in_expression_sheet.txt")
  )

  metrics_df <- metrics_df %>%
    distinct(gene_name, .keep_all = TRUE)
}


# ---- 7. Read gene lists ------------------------------------------------------

upregulated_genes <- read_gene_list(gene_lists_df, upregulated_gene_column)
downregulated_genes <- read_gene_list(gene_lists_df, downregulated_gene_column)
tlr_independent_genes <- read_gene_list(gene_lists_df, tlr_independent_column)


# ---- 8. Create output sheets -------------------------------------------------

upregulated_output <- make_output_sheet(
  metrics_df,
  upregulated_genes,
  sheet_label = "up_regulated"
)

downregulated_output <- make_output_sheet(
  metrics_df,
  downregulated_genes,
  sheet_label = "down_regulated"
)

TLR_independent_output <- make_output_sheet(
  metrics_df,
  tlr_independent_genes,
  sheet_label = "TLR_independent"
)


# ---- 9. Write output workbook ------------------------------------------------

writexl::write_xlsx(
  list(
    up_regulated = upregulated_output,
    downregulated = downregulated_output,
    TLR_independent = TLR_independent_output
  ),
  path = output_file
)


# ---- 10. Save run information ------------------------------------------------

sink(file.path(results_dir, "sessionInfo.txt"))
print(sessionInfo())
sink()

message("Analysis complete.")
message("Output workbook saved to: ", output_file)
message("Genes classified as heritable if CV_stimulated_clonal / CV_stimulated_mixed >= ", heritability_threshold)

