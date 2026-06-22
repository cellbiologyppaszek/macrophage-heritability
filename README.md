# Publication code and processed data

This repository contains the analysis code and small processed data tables used for the publication.

Large raw microscopy images, large intermediate image-analysis objects, and raw sequencing input files are not stored in GitHub. They are available from the BioImage Archive and should be downloaded by readers who want to rerun the full analysis.

## Data availability

BioImage Archive accession: `S-BIAD2517`

BioImage Archive study page:

https://www.ebi.ac.uk/biostudies/bioimages/studies/S-BIAD2517

Direct download link for the external data archive:

https://ftp.ebi.ac.uk/biostudies/fire/S-BIAD/517/S-BIAD2517/Files/external_data.zip

After downloading and unzipping `external_data.zip`, place the external files locally under:

```text
data/external/
```

Do not commit `data/external/` to GitHub.

## 1. Recommended GitHub structure

Organize the repository like this:

```text
official_code/
|
|-- README.md
|-- .gitignore
|-- .gitattributes
|
|-- data/
|   |
|   |-- processed/
|   |   |-- Cv_7_signals_FCV.csv
|   |   |-- cv_experiment_vs_theory_p65_relb_cd36.csv
|   |   |-- cv_experiment_vs_theory_tnfa_il1b_combine.csv
|   |   |-- cv_experiment_vs_theory_tnfr2_f480.csv
|   |
|   |-- external/
|       |-- README_external_data.md
|       |
|       | Large raw or intermediate input files downloaded from
|       | BioImage Archive should be placed here locally.
|       | This folder should not be committed to GitHub.
|
|-- notebooks/
|   |-- 7gens_plot.ipynb
|   |-- cd36_listeria.ipynb
|   |-- Density-CD36 IL1-programmes.ipynb
|   |-- NMF_CD36_IL1.ipynb
|   |-- p65_relb_cd36_costainning.ipynb
|   |-- tnfa_il1b_costainning.ipynb
|   |-- tnfr2_f480_costainning.ipynb
|
|-- R_scripts/
|   |-- Constant_Switching_colony_tracking_Simulation_kmax_beta.R
|   |-- Constant_Switching_Parameter_Simulation.R
|   |-- Constant_Switching_Rates_Calculation.Rmd
|   |-- Density_Dependent_colony_tracking_Simulation_kmax_beta.R
|   |-- Density_Dependent_Parameter_Simulation.R
|   |-- Density_Dependent_Parameter_Simulation_kmax_beta.R
|   |-- Optimized_Parameter_Estimation.R
|   |-- cv_fold_change_MemorySeq.R
|   |-- cv_fold_change_scRNAseq.R
|
|-- results/
|   | Generated tables, pickles, and intermediate output files.
|   | This folder should not be committed to GitHub.
|
|-- figures/
    | Generated SVG, PNG, or PDF figures.
    | This folder should not be committed to GitHub.
```

## 2. What belongs in GitHub

| File type | Store in GitHub? | Reason |
|---|---:|---|
| `README.md` | Yes | Explains how to use the repository. |
| Jupyter notebooks, `.ipynb` | Yes | Main Python analysis workflows. |
| R scripts, `.R` and `.Rmd` | Yes | Simulation, parameter estimation, and sequencing CV analyses. |
| Small processed CSV tables | Yes | Allow readers to reproduce summary plots quickly. |
| Raw microscopy images | No | Too large for GitHub; deposit in BioImage Archive. |
| Large `.pkl`, `.pkl.gz`, `.h5`, `.tif`, `.czi`, `.nd2` files | No | Store in BioImage Archive or another data archive. |
| Temporary files, checkpoints, logs | No | Not needed for reproducibility. |

## 3. Quick start for readers

### 3.1 Clone the GitHub repository

```bash
git clone https://github.com/cellbiologyppaszek/macrophage-heritability.git
cd macrophage-heritability
```

### 3.2 Download external data from BioImage Archive

Download the external data archive:

```text
https://ftp.ebi.ac.uk/biostudies/fire/S-BIAD/517/S-BIAD2517/Files/external_data.zip
```

Unzip it locally and place the required files under:

```text
data/external/
```

The small processed CSV files are already included in:

```text
data/processed/
```

### 3.3 Install Python packages

A standard Python environment with Jupyter is sufficient for most notebooks.

```bash
python -m venv .venv

# Windows PowerShell:
.venv\Scripts\Activate.ps1

# macOS/Linux:
# source .venv/bin/activate

pip install pandas numpy matplotlib seaborn scipy scikit-learn statsmodels jupyter
pip install mygene statannotations libpysal esda splot openpyxl
```

The spatial-analysis packages `libpysal`, `esda`, and `splot` are mainly needed for `notebooks/cd36_listeria.ipynb`.

### 3.4 Install R packages

Use R version 4.0 or newer.

```r
install.packages(c(
  "readxl",
  "ggplot2",
  "dplyr",
  "purrr",
  "tidyr",
  "patchwork",
  "rmarkdown",
  "knitr",
  "BB",
  "readr",
  "stringr",
  "tibble",
  "writexl"
))
```

## 4. Important path note

The notebooks were originally written as interactive analysis notebooks. Several notebooks read files using simple relative filenames such as:

```python
pd.read_pickle("p65_relb_cd36.pkl.gz", compression="gzip")
pd.read_csv("cv_experiment_vs_theory_p65_relb_cd36.csv")
```

If you use the recommended GitHub structure, update the first data-loading cells so that files are read from `data/processed/` or `data/external/`.

For example, if the notebook is inside the `notebooks/` folder:

```python
from pathlib import Path

PROCESSED_DIR = Path("../data/processed")
EXTERNAL_DIR = Path("../data/external")
RESULTS_DIR = Path("../results")
FIGURES_DIR = Path("../figures")

CSV_P65_RELB_CD36 = PROCESSED_DIR / "cv_experiment_vs_theory_p65_relb_cd36.csv"
CSV_TNFR2_F480 = PROCESSED_DIR / "cv_experiment_vs_theory_tnfr2_f480.csv"
CSV_TNFA_IL1B = PROCESSED_DIR / "cv_experiment_vs_theory_tnfa_il1b_combine.csv"
```

For raw image-analysis notebooks, update pickle paths similarly:

```python
file_path = EXTERNAL_DIR / "pkl_intermediate_files" / "p65_relb_cd36.pkl.gz"
quant_df = pd.read_pickle(file_path, compression="gzip")
```

## 5. Included processed CSV data

These files are small processed summary tables. They are suitable for GitHub and allow readers to reproduce the final coefficient-of-variation comparison plots without downloading the full raw microscopy dataset.

| File | Rows x columns | Role | Main columns | Used by |
|---|---:|---|---|---|
| `data/processed/Cv_7_signals_FCV.csv` | 112 x 11 | Combined seven-signal summary table containing experiment/theory CV summaries and fold-change CV values from MemorySeq and scRNA-seq. | `Cell_Bin`, `Lipid A`, `Mean`, `SD`, `Source`, `Signal`, `Sampling`, `gene_name`, `FCV MemSeq`, `FCV scRNA-sq` | Final summary/reporting table |
| `data/processed/cv_experiment_vs_theory_p65_relb_cd36.csv` | 48 x 6 | Experiment-versus-theory CV summary for p65, RelB, and CD36. | `Cell_Bin`, `Lipid A`, `Mean`, `SD`, `Source`, `Signal` | `notebooks/7gens_plot.ipynb` |
| `data/processed/cv_experiment_vs_theory_tnfa_il1b_combine.csv` | 32 x 6 | Experiment-versus-theory CV summary for TNF-alpha and IL-1-beta. | `Cell_Bin`, `Lipid A`, `Mean`, `SD`, `Source`, `Signal` | `notebooks/7gens_plot.ipynb` |
| `data/processed/cv_experiment_vs_theory_tnfr2_f480.csv` | 32 x 6 | Experiment-versus-theory CV summary for TNFR2 and F4/80. | `Cell_Bin`, `Lipid A`, `Mean`, `SD`, `Source`, `Signal` | `notebooks/7gens_plot.ipynb` |

### Column definitions

| Column | Meaning |
|---|---|
| `Cell_Bin` | Colony/cell-number bin. Values include `<4`, `4-6`, `6-8`, and `>8`. |
| `Lipid A` | Stimulation condition. Values include `Mock` and `500 ng/mL`. |
| `Mean` | Mean coefficient-of-variation summary value for the group. |
| `SD` | Standard deviation of the coefficient-of-variation summary value for the group. |
| `Source` | Indicates whether the value comes from `Experiment` or `Theory`. |
| `Signal` | Measured marker/signal, for example `per_cd36_mean`, `per_p65_mean`, `per_relb_mean`, `per_il1b_mean`, `per_tnfa_mean`, `per_tnfr2_mean`, or `per_f480_mean`. |
| `Sampling` | Sampling scheme used in the combined seven-signal FCV table. |
| `gene_name` | Gene name corresponding to the signal. |
| `FCV MemSeq` | Fold-change CV value calculated from MemorySeq data. |
| `FCV scRNA-sq` | Fold-change CV value calculated from scRNA-seq data. |

## 6. External data expected from BioImage Archive

The following files are required only if readers want to rerun the full raw-data analysis. Download them from the BioImage Archive and place them locally in `data/external/`.

### 6.1 Image-analysis intermediate files

Recommended location:

```text
data/external/pkl_intermediate_files/
```

| Expected file | Used by | Purpose |
|---|---|---|
| `p65_relb_cd36.pkl.gz` | `notebooks/p65_relb_cd36_costainning.ipynb` | Raw or intermediate single-cell image-analysis table for p65, RelB, and CD36 co-staining. |
| `WellsInfo1.pkl` | `notebooks/p65_relb_cd36_costainning.ipynb` | Plate/well layout metadata for the p65/RelB/CD36 analysis. |
| `p65_relb_cd36_clone_threshold.csv` | `notebooks/p65_relb_cd36_costainning.ipynb` | Clone threshold table used to classify or filter p65/RelB/CD36 clone populations. |
| `tnfr2_f480.pkl.gz` | `notebooks/tnfr2_f480_costainning.ipynb` | Raw or intermediate single-cell image-analysis table for TNFR2 and F4/80 co-staining. |
| `WellsInfo2.pkl` | `notebooks/tnfr2_f480_costainning.ipynb` | Plate/well layout metadata for the TNFR2/F4/80 analysis. |
| `IF22_co.pkl.gz` | `notebooks/tnfa_il1b_costainning.ipynb` | Raw or intermediate single-cell image-analysis table for one TNF-alpha/IL-1-beta experiment. |
| `IF24_co.pkl.gz` | `notebooks/tnfa_il1b_costainning.ipynb` | Raw or intermediate single-cell image-analysis table for one TNF-alpha/IL-1-beta experiment. |
| `IF48_listeria.pkl.gz` | `notebooks/cd36_listeria.ipynb` | Raw or intermediate single-cell image-analysis table for one CD36/Listeria experiment. |
| `IF50_listeria.pkl.gz` | `notebooks/cd36_listeria.ipynb` | Raw or intermediate single-cell image-analysis table for one CD36/Listeria experiment. |
| `IF54_listeria.pkl.gz` | `notebooks/cd36_listeria.ipynb` | Raw or intermediate single-cell image-analysis table for one CD36/Listeria experiment. |

### 6.2 Sequencing input files

Recommended location:

```text
data/external/sequencing_inputs/
```

| Expected file or file pattern | Used by | Purpose |
|---|---|---|
| `*_quant.sf` | `notebooks/Density-CD36 IL1-programmes.ipynb` | Transcript-level quantification files used to build gene-level expression summaries. |
| `scRNA_filtered_genes.pkl` | `notebooks/NMF_CD36_IL1.ipynb` | Filtered scRNA-seq gene-expression matrix used for NMF/program analysis. |
| `MemorySeq data.xlsx` | `R_scripts/cv_fold_change_MemorySeq.R` | MemorySeq expression workbook and gene-list workbook sheets. |
| `combined_logcounts_JM01.xlsx` to `combined_logcounts_JM12.xlsx` | `R_scripts/cv_fold_change_scRNAseq.R` | scRNA-seq logcount matrices for clonal and mixed samples. |

### 6.3 Excel input files for simulation and model fitting

Recommended location:

```text
data/external/excel_inputs/
```

| Expected file | Used by | Purpose |
|---|---|---|
| `IL! and CD36 fractions.xlsx` | `R_scripts/Constant_Switching_Parameter_Simulation.R`, `R_scripts/Density_Dependent_Parameter_Simulation.R`, `R_scripts/Density_Dependent_Parameter_Simulation_kmax_beta.R`, `R_scripts/Optimized_Parameter_Estimation.R`, `R_scripts/Constant_Switching_Rates_Calculation.Rmd` | Experimental fraction data used for parameter estimation and simulation benchmarking. |
| `Corrected_odds_ratio_CD36_tnfa_il1b.xlsx` | `R_scripts/Constant_Switching_Rates_Calculation.Rmd` | Odds-ratio workbook used to estimate switching rates. |

If filenames are changed during BioImage Archive download or extraction, update the corresponding path variables in the notebook or R script.

## 7. How to reproduce the final seven-signal CV comparison figure

This is the simplest analysis to rerun because it uses only the processed CSV files included in GitHub.

1. Make sure the processed CSV files are in:

   ```text
   data/processed/
   ```

2. Start Jupyter:

   ```bash
   jupyter lab
   ```

3. Open:

   ```text
   notebooks/7gens_plot.ipynb
   ```

4. Update the CSV path variables if needed:

   ```python
   from pathlib import Path

   PROCESSED_DIR = Path("../data/processed")

   CSV_P65_RELB_CD36 = PROCESSED_DIR / "cv_experiment_vs_theory_p65_relb_cd36.csv"
   CSV_TNFR2_F480 = PROCESSED_DIR / "cv_experiment_vs_theory_tnfr2_f480.csv"
   CSV_TNFA_IL1B = PROCESSED_DIR / "cv_experiment_vs_theory_tnfa_il1b_combine.csv"

   SVG_PATH_GRID = "../figures/cv_grid.svg"
   SVG_PATH_SINGLE = "../figures/cv_single.svg"
   ```

5. Run all cells.

Expected output:

```text
figures/cv_grid.svg
```

This figure compares experimental and theoretical CV summaries across the seven measured signals.

## 8. Notebook guide

| Notebook | Main role | Main input data | Main outputs |
|---|---|---|---|
| `notebooks/7gens_plot.ipynb` | Generates the final seven-signal experiment-versus-theory CV plot. | Processed CSV files in `data/processed/`. | `cv_grid.svg`, optional single-signal SVG plots. |
| `notebooks/p65_relb_cd36_costainning.ipynb` | Processes p65, RelB, and CD36 co-staining image-analysis data; performs colony/clone summaries and CV comparison. | `p65_relb_cd36.pkl.gz`, `WellsInfo1.pkl`, `p65_relb_cd36_clone_threshold.csv`. | Cleaned pickle files, clone morphology tables, average signal tables, `cv_experiment_vs_theory_p65_relb_cd36.csv`, and SVG figures. |
| `notebooks/tnfa_il1b_costainning.ipynb` | Processes TNF-alpha and IL-1-beta co-staining image-analysis data. | `IF22_co.pkl.gz`, `IF24_co.pkl.gz`. | Cleaned/intermediate pickle files for TNF-alpha/IL-1-beta analysis. |
| `notebooks/tnfr2_f480_costainning.ipynb` | Processes TNFR2 and F4/80 co-staining image-analysis data. | `tnfr2_f480.pkl.gz`, `WellsInfo2.pkl`. | Cleaned/intermediate pickle files and plate overview figures. |
| `notebooks/cd36_listeria.ipynb` | Processes CD36/Listeria perinuclear signal data and spatial statistics. | `IF48_listeria.pkl.gz`, `IF50_listeria.pkl.gz`, `IF54_listeria.pkl.gz`. | `cd36_Listeria_clean.pkl` and spatial/statistical plots. |
| `notebooks/Density-CD36 IL1-programmes.ipynb` | Builds gene-level expression tables from transcript quantification files. | Files matching `*_quant.sf`. | `df_gene.pkl`, `df_gene.csv`. |
| `notebooks/NMF_CD36_IL1.ipynb` | Performs latent program/NMF analysis for CD36 and IL1B-related gene programs. | `scRNA_filtered_genes.pkl`. | NMF/program annotation figures including SVG, PNG, and PDF outputs. |

## 9. R script guide

| Script | Main role | Main input data | Main outputs |
|---|---|---|---|
| `R_scripts/Constant_Switching_Rates_Calculation.Rmd` | Estimates constant switching rates from odds-ratio and fraction workbooks. | `Corrected_odds_ratio_CD36_tnfa_il1b.xlsx`, `IL! and CD36 fractions.xlsx`. | Rendered HTML report and estimated rate summaries printed in the report. |
| `R_scripts/Optimized_Parameter_Estimation.R` | Optimizes density-dependent switching parameters by matching simulated and experimental mean/CV profiles. | `IL! and CD36 fractions.xlsx`. | Optimized parameter values printed to the console. |
| `R_scripts/Constant_Switching_Parameter_Simulation.R` | Runs Gillespie simulations with constant ON/OFF switching rates and compares with experimental fraction data. | `IL! and CD36 fractions.xlsx`. | Simulation plots; optional CSV export if the commented `write.csv()` line is enabled. |
| `R_scripts/Density_Dependent_Parameter_Simulation.R` | Runs density-dependent switching simulations with a Hill-type activation function. | `IL! and CD36 fractions.xlsx`. | `CD36_simulated_density_constant_rates.csv`, `simulated_loess_trend_data_kON_constant_wCI.csv`. |
| `R_scripts/Density_Dependent_Parameter_Simulation_kmax_beta.R` | Runs density-dependent simulations with Beta-distributed maximum activation rates. | `IL! and CD36 fractions.xlsx`. | `CD36_simulated_density_beta_rates_doublerate14_4.csv`, `simulated_loess_trend_data_kON_constant_wCI.csv`, `experimental_loess_trend_data_kON_constant_wCI.csv`. |
| `R_scripts/Constant_Switching_colony_tracking_Simulation_kmax_beta.R` | Simulates lineage-resolved colony growth with constant phenotypic switching rates. | No external file required; parameters are defined inside the script. | Lineage tracking table in memory and plots printed to the R graphics device. |
| `R_scripts/Density_Dependent_colony_tracking_Simulation_kmax_beta.R` | Simulates lineage-resolved colony growth with density-dependent switching rates. | No external file required; parameters are defined inside the script. | Lineage tracking table in memory and plots printed to the R graphics device; optional CSV export if the commented `write.csv()` line is enabled. |
| `R_scripts/cv_fold_change_MemorySeq.R` | Calculates MemorySeq CV fold-change and heritability classification for gene sets. | `data/external/sequencing_inputs/MemorySeq data.xlsx`. | `results/CV_fold_change_heritability_MemorySeq.xlsx`, missing-gene text files, and `sessionInfo.txt`. |
| `R_scripts/cv_fold_change_scRNAseq.R` | Calculates scRNA-seq CV fold-change between clonal and mixed populations. | `combined_logcounts_JM01.xlsx` to `combined_logcounts_JM12.xlsx`. | `results/scRNAseq_CV_fold_change_clonal_vs_mixed.xlsx`. |

## 10. Running the R scripts

Run R scripts from the repository root unless otherwise noted.

### Example: MemorySeq CV fold-change analysis

Expected input:

```text
data/external/sequencing_inputs/MemorySeq data.xlsx
```

Run:

```bash
Rscript R_scripts/cv_fold_change_MemorySeq.R
```

Expected output:

```text
results/CV_fold_change_heritability_MemorySeq.xlsx
results/sessionInfo.txt
```

### Example: scRNA-seq CV fold-change analysis

Expected inputs:

```text
data/external/sequencing_inputs/combined_logcounts_JM01.xlsx
data/external/sequencing_inputs/combined_logcounts_JM02.xlsx
...
data/external/sequencing_inputs/combined_logcounts_JM12.xlsx
```

Run:

```bash
Rscript R_scripts/cv_fold_change_scRNAseq.R
```

Expected output:

```text
results/scRNAseq_CV_fold_change_clonal_vs_mixed.xlsx
```

### Example: render the rate-calculation R Markdown file

```bash
Rscript -e "rmarkdown::render('R_scripts/Constant_Switching_Rates_Calculation.Rmd')"
```

### Important R path note

Some R scripts currently contain local Windows paths, for example:

```r
DATA_PATH <- "C:/Users/apurv/OneDrive/Desktop/Data Analysis/rate switch/IL! and CD36 fractions.xlsx"
```

For publication, replace these with relative paths, for example:

```r
DATA_PATH <- file.path("data", "external", "excel_inputs", "IL! and CD36 fractions.xlsx")
```

## 11. Suggested order for full reproducibility

The exact order depends on which part of the paper the reader wants to reproduce.

### A. Reproduce final summary plots only

Use this route if the reader only wants to recreate the final seven-signal CV plot.

```text
data/processed/*.csv
        |
        v
notebooks/7gens_plot.ipynb
        |
        v
figures/cv_grid.svg
```

### B. Rerun image-analysis notebooks from large BioImage Archive files

Use this route if the reader wants to regenerate cleaned single-cell/clone tables from image-analysis outputs.

```text
data/external/pkl_intermediate_files/*.pkl.gz
data/external/pkl_intermediate_files/WellsInfo*.pkl
        |
        v
notebooks/*costainning.ipynb and notebooks/cd36_listeria.ipynb
        |
        v
results/*.pkl, results/*.csv, figures/*.svg
```

### C. Rerun stochastic switching simulations

Use this route if the reader wants to reproduce model simulations and parameter-estimation workflows.

```text
data/external/excel_inputs/IL! and CD36 fractions.xlsx
data/external/excel_inputs/Corrected_odds_ratio_CD36_tnfa_il1b.xlsx
        |
        v
R_scripts/Constant_Switching_Rates_Calculation.Rmd
R_scripts/Optimized_Parameter_Estimation.R
R_scripts/*Simulation*.R
        |
        v
simulation summary CSVs and model-comparison figures
```

### D. Rerun MemorySeq/scRNA-seq CV fold-change analyses

Use this route if the reader wants to reproduce the FCV values reported in the combined seven-signal table.

```text
data/external/sequencing_inputs/MemorySeq data.xlsx
data/external/sequencing_inputs/combined_logcounts_JM*.xlsx
        |
        v
R_scripts/cv_fold_change_MemorySeq.R
R_scripts/cv_fold_change_scRNAseq.R
        |
        v
results/*.xlsx
```

## 12. Suggested `.gitignore`

Use a `.gitignore` file to prevent raw or temporary files from being committed accidentally:

```text
# External large files downloaded from BioImage Archive
data/external/

# Generated outputs
results/
figures/

# Jupyter temporary files
.ipynb_checkpoints/

# Python temporary files
__pycache__/
*.pyc

# R temporary files
.Rhistory
.RData
.Rproj.user/

# System files
.DS_Store
Thumbs.db

# Logs
*.log

# Large raw or intermediate data
*.pkl
*.pkl.gz
*.h5
*.hdf5
*.tif
*.tiff
*.czi
*.nd2
*.lif
*.zip
*.tar.gz
```

The four processed CSV files in `data/processed/` should remain tracked by Git.

## 13. Notes for users

1. The processed CSV files in `data/processed/` are sufficient to reproduce the final CV comparison plot.
2. The full raw-data workflow requires external files from BioImage Archive.
3. Some notebooks and R scripts may require minor path edits because they were originally developed interactively.
4. For exact reproducibility, keep the downloaded BioImage Archive filenames unchanged, or update the path variables consistently.
5. Generated output files should be written to `results/` or `figures/`, not mixed with source code.

## 14. Citation

If you use this code, processed data, or the associated external data, please cite the associated publication:

```text
DOI: 10.64898/2026.06.17.732820
```

Publication DOI:

https://doi.org/10.64898/2026.06.17.732820

Data availability:

```text
Raw microscopy and large analysis files are available from the BioImage Archive:
https://ftp.ebi.ac.uk/biostudies/fire/S-BIAD/517/S-BIAD2517/Files/external_data.zip
```
