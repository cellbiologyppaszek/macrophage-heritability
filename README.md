# Heritable single-cell gene expression analysis

This repository contains code and processed data used for the analysis of heritable single-cell gene expression programs in macrophage populations.

The repository is intended to help readers understand, inspect, and reproduce the computational analysis associated with the publication.

Large raw image files, large intermediate `.pkl` / `.pkl.gz` files, and sequencing input files are not stored in this GitHub repository. These files should be downloaded separately from the BioImage Archive and placed locally in `data/external/`.

BioImage Archive accession: `TO_BE_ADDED`

## Repository structure

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
|       | Large raw or intermediate input files downloaded from the
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
|   |-- Constant_Switching_Parameter_Simulation.R
|   |-- Constant_Switching_Rates_Calculation.Rmd
|   |-- Density_Dependent_Parameter_Simulation.R
|   |-- Density_Dependent_Parameter_Simulation_kmax_beta.R
|   |-- Optimized_Parameter_Estimation.R
|   |-- cv_fold_change_MemorySeq.R
|   |-- cv_fold_change_scRNAseq.R
```

## What is included in GitHub

This GitHub repository includes:

| Folder | Content | Purpose |
|---|---|---|
| `data/processed/` | Small processed CSV files | Processed summary data used for plotting and comparison |
| `notebooks/` | Jupyter notebooks | Python analysis and figure-generation workflows |
| `R_scripts/` | R and R Markdown scripts | Simulation, parameter estimation, and CV fold-change analysis |
| `README.md` | Main documentation | Instructions for readers |
| `.gitignore` | Git exclusion rules | Prevents raw and large files from being uploaded |
| `.gitattributes` | Line-ending rules | Keeps text files consistent across operating systems |

## What is not included in GitHub

The following files should not be committed to GitHub:

| File type | Example | Where to store |
|---|---|---|
| Raw microscopy images | `.tif`, `.tiff`, `.czi`, `.nd2`, `.lif` | BioImage Archive |
| Large intermediate image-analysis files | `.pkl`, `.pkl.gz` | BioImage Archive |
| Large sequencing input files | `*_quant.sf`, large `.xlsx`, large `.pkl` | BioImage Archive |
| Temporary files | `.ipynb_checkpoints/`, `.Rhistory`, `__pycache__/` | Do not archive |
| Generated local output | `results/`, `figures/` | Regenerate locally |

Locally, these external files should be placed under:

```text
data/external/
```

## Processed data files included in GitHub

The following processed CSV files are included in `data/processed/`.

| File | Role |
|---|---|
| `Cv_7_signals_FCV.csv` | Processed coefficient-of-variation and fold-change data for multiple measured signals |
| `cv_experiment_vs_theory_p65_relb_cd36.csv` | Processed comparison between experimental and theoretical CV values for p65, RelB, and CD36 |
| `cv_experiment_vs_theory_tnfa_il1b_combine.csv` | Processed comparison between experimental and theoretical CV values for TNFA and IL1B |
| `cv_experiment_vs_theory_tnfr2_f480.csv` | Processed comparison between experimental and theoretical CV values for TNFR2 and F4/80 |

These files are small enough to be stored directly in GitHub and can be used by readers without downloading the full raw BioImage Archive dataset.

## External data expected from BioImage Archive

To rerun the complete analysis from large raw or intermediate files, download the required files from the BioImage Archive and place them under `data/external/`.

Recommended local structure:

```text
data/external/
|
|-- raw_image_files/
|   | Raw microscopy image files
|
|-- pkl_intermediate_files/
|   | Image-analysis intermediate tables
|
|-- sequencing_inputs/
|   | Sequencing or transcriptomic input files
|
|-- excel_inputs/
    | Excel files used by R scripts
```

### Image-analysis intermediate files

Place these files in:

```text
data/external/pkl_intermediate_files/
```

| File | Used by |
|---|---|
| `p65_relb_cd36.pkl.gz` | `notebooks/p65_relb_cd36_costainning.ipynb` |
| `WellsInfo1.pkl` | `notebooks/p65_relb_cd36_costainning.ipynb` |
| `tnfr2_f480.pkl.gz` | `notebooks/tnfr2_f480_costainning.ipynb` |
| `WellsInfo2.pkl` | `notebooks/tnfr2_f480_costainning.ipynb` |
| `IF22_co.pkl.gz` | `notebooks/tnfa_il1b_costainning.ipynb` |
| `IF24_co.pkl.gz` | `notebooks/tnfa_il1b_costainning.ipynb` |
| `IF48_listeria.pkl.gz` | `notebooks/cd36_listeria.ipynb` |
| `IF50_listeria.pkl.gz` | `notebooks/cd36_listeria.ipynb` |
| `IF54_listeria.pkl.gz` | `notebooks/cd36_listeria.ipynb` |

### Sequencing input files

Place these files in:

```text
data/external/sequencing_inputs/
```

| File or file pattern | Used by |
|---|---|
| `*_quant.sf` | `notebooks/Density-CD36 IL1-programmes.ipynb` |
| `scRNA_filtered_genes.pkl` | `notebooks/NMF_CD36_IL1.ipynb` |
| `MemorySeq data.xlsx` | `R_scripts/cv_fold_change_MemorySeq.R` |
| `combined_logcounts_JM01.xlsx` to `combined_logcounts_JM12.xlsx` | `R_scripts/cv_fold_change_scRNAseq.R` |

### Excel input files for model fitting and simulation

Place these files in:

```text
data/external/excel_inputs/
```

| File | Used by |
|---|---|
| `IL1 and CD36 fractions.xlsx` | `R_scripts/Constant_Switching_Parameter_Simulation.R`, `R_scripts/Density_Dependent_Parameter_Simulation.R`, `R_scripts/Density_Dependent_Parameter_Simulation_kmax_beta.R`, `R_scripts/Optimized_Parameter_Estimation.R`, `R_scripts/Constant_Switching_Rates_Calculation.Rmd` |
| `odds_ratio_CD36_tnfa_il1b.xlsx` | `R_scripts/Constant_Switching_Rates_Calculation.Rmd` |

Note: Please check whether the filename `IL1 and CD36 fractions.xlsx` is intentional. If the intended name is `IL1 and CD36 fractions.xlsx` or `IL-1 and CD36 fractions.xlsx`, update the filename and the scripts consistently.

## Software requirements

### Python

Recommended Python version:

```text
Python 3.9 or newer
```

Common Python packages used by the notebooks include:

```text
numpy
pandas
matplotlib
seaborn
scipy
scikit-learn
jupyter
openpyxl
```

Install them using:

```bash
pip install numpy pandas matplotlib seaborn scipy scikit-learn jupyter openpyxl
```

Or, using conda:

```bash
conda create -n macrophage_expression python=3.10
conda activate macrophage_expression
conda install numpy pandas matplotlib seaborn scipy scikit-learn jupyter openpyxl
```

### R

Recommended R version:

```text
R 4.0 or newer
```

Common R packages used by the scripts include:

```text
readxl
dplyr
ggplot2
tidyr
data.table
```

Install them in R using:

```r
install.packages(c("readxl", "dplyr", "ggplot2", "tidyr", "data.table"))
```

## How to run the analysis

### 1. Clone the repository

```bash
git clone https://github.com/cellbiologyppaszek/macrophage-heritability.git
cd macrophage-heritability
```

### 2. Download external files

Download the external raw and intermediate files from the BioImage Archive.

Place them locally under:

```text
data/external/
```

For example:

```text
data/external/pkl_intermediate_files/
data/external/sequencing_inputs/
data/external/excel_inputs/
```

### 3. Run the Jupyter notebooks

Start Jupyter from the repository root:

```bash
jupyter notebook
```

Open notebooks from the `notebooks/` folder.

Recommended order:

| Step | Notebook | Purpose |
|---|---|---|
| 1 | `notebooks/7gens_plot.ipynb` | Plot and summarize processed CV and fold-change data |
| 2 | `notebooks/p65_relb_cd36_costainning.ipynb` | Analyze p65, RelB, and CD36 co-staining data |
| 3 | `notebooks/tnfa_il1b_costainning.ipynb` | Analyze TNFA and IL1B co-staining data |
| 4 | `notebooks/tnfr2_f480_costainning.ipynb` | Analyze TNFR2 and F4/80 co-staining data |
| 5 | `notebooks/cd36_listeria.ipynb` | Analyze CD36-related Listeria infection data |
| 6 | `notebooks/Density-CD36 IL1-programmes.ipynb` | Analyze CD36 and IL1-related transcriptional programs |
| 7 | `notebooks/NMF_CD36_IL1.ipynb` | Perform NMF analysis of CD36 and IL1-associated programs |

Some notebooks may need path updates depending on where the BioImage Archive files are placed locally.

For example, if a notebook currently reads:

```python
pd.read_pickle("p65_relb_cd36.pkl.gz")
```

change it to:

```python
pd.read_pickle("../data/external/pkl_intermediate_files/p65_relb_cd36.pkl.gz")
```

If a notebook reads a processed CSV file from the current folder, change it to:

```python
pd.read_csv("../data/processed/Cv_7_signals_FCV.csv")
```

when running from the `notebooks/` folder.

### 4. Run the R scripts

Run R scripts from the repository root, or update paths inside the scripts.

Example from command line:

```bash
Rscript R_scripts/Constant_Switching_Parameter_Simulation.R
Rscript R_scripts/Density_Dependent_Parameter_Simulation.R
Rscript R_scripts/Optimized_Parameter_Estimation.R
Rscript R_scripts/cv_fold_change_MemorySeq.R
Rscript R_scripts/cv_fold_change_scRNAseq.R
```

For the R Markdown file:

```bash
Rscript -e "rmarkdown::render('R_scripts/Constant_Switching_Rates_Calculation.Rmd')"
```

Some R scripts may contain local absolute paths from the original analysis computer. Replace them with relative paths such as:

```r
data_path <- "data/external/excel_inputs/IL1 and CD36 fractions.xlsx"
```

or:

```r
data_path <- "data/processed/Cv_7_signals_FCV.csv"
```


## Reproducibility notes

1. The processed CSV files in `data/processed/` are included so that readers can inspect key processed data without downloading the full raw dataset.
2. Full rerunning of the image-analysis workflows requires external files from the BioImage Archive.
3. Some notebooks and R scripts may require path adjustments after reorganizing the repository.
4. Large raw data and large intermediate files should remain in BioImage Archive, not GitHub.
5. GitHub should contain only code, documentation, and small processed data files.

## Citation

If you use this repository, please cite the associated publication:

```text
TO_BE_ADDED
```

## Contact

For questions about the code or data organization, please contact the corresponding author of the associated publication.
