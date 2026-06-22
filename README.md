# Publication code and processed data

This repository contains the analysis code and small processed data tables used for the publication.  
Large raw microscopy images, large intermediate image-analysis objects, and raw sequencing input files should **not** be stored in GitHub. They should be deposited in the BioImage Archive and downloaded by readers when they want to rerun the full analysis.

```text
BioImage Archive accession: S-BIAD2517
BioImage Archive download link: https://www.ebi.ac.uk/biostudies/bioimages/studies/S-BIAD2517

Raw microscopy and large analysis files are available from the BioImage Archive:
https://ftp.ebi.ac.uk/biostudies/fire/S-BIAD/517/S-BIAD2517/Files/external\_data.zip

\---

## 1\. Recommended GitHub structure

Organize the repository like this:

```text
official\_code/
├── README.md
├── .gitignore
│
├── data/
│   ├── processed/
│   │   ├── Cv\_7\_signals\_FCV.csv
│   │   ├── cv\_experiment\_vs\_theory\_p65\_relb\_cd36.csv
│   │   ├── cv\_experiment\_vs\_theory\_tnfa\_il1b\_combine.csv
│   │   └── cv\_experiment\_vs\_theory\_tnfr2\_f480.csv
│   │
│   └── external/
│       └── README\_external\_data.md
│       # Raw or large input files downloaded from BioImage Archive go here.
│       # This folder should not be committed to GitHub.
│
├── notebooks/
│   ├── 7gens\_plot.ipynb
│   ├── cd36\_listeria.ipynb
│   ├── Density-CD36 IL1-programmes.ipynb
│   ├── NMF\_CD36\_IL1.ipynb
│   ├── p65\_relb\_cd36\_costainning.ipynb
│   ├── tnfa\_il1b\_costainning.ipynb
│   └── tnfr2\_f480\_costainning.ipynb
│
├── R\_scripts/
│   ├── Constant\_Switching\_colony\_tracking\_Simulation\_kmax\_beta.R
│   ├── Constant\_Switching\_Parameter\_Simulation.R
│   ├── Constant\_Switching\_Rates\_Calculation.Rmd
│   ├── Density\_Dependent\_colony\_tracking\_Simulation\_kmax\_beta.R
│   ├── Density\_Dependent\_Parameter\_Simulation.R
│   ├── Density\_Dependent\_Parameter\_Simulation\_kmax\_beta.R
│   ├── Optimized\_Parameter\_Estimation.R
│   ├── cv\_fold\_change\_MemorySeq.R
│   └── cv\_fold\_change\_scRNAseq.R
│
├── results/
│   # Generated tables, pickles, and intermediate output files.
│
└── figures/
    # Generated SVG, PNG, or PDF figures.
```

### What belongs in GitHub?

|File type|Store in GitHub?|Reason|
|-|-:|-|
|`README.md`|Yes|Explains how to use the repository.|
|Jupyter notebooks, `.ipynb`|Yes|Main Python analysis workflows.|
|R scripts, `.R` / `.Rmd`|Yes|Simulation, parameter estimation, and sequencing CV analyses.|
|Small processed CSV tables|Yes|Allow readers to reproduce summary plots quickly.|
|Raw microscopy images|No|Too large for GitHub; deposit in BioImage Archive.|
|Large `.pkl`, `.pkl.gz`, `.h5`, `.tif`, `.czi`, `.nd2` files|No|Deposit in BioImage Archive or another data archive.|
|Temporary files, checkpoints, logs|No|Not needed for reproducibility.|

\---

## 2\. Quick start for readers

### 2.1 Clone the GitHub repository

```bash
git clone https://github.com/\[USER]/\[REPOSITORY].git
cd \[REPOSITORY]
```

### 2.2 Download external data from BioImage Archive

After the BioImage Archive accession is available, download the raw and large processed input files and place them in:

```text
data/external/
```

The small processed CSV files are already included in:

```text
data/processed/
```

### 2.3 Install Python packages

A standard Python environment with Jupyter is sufficient for most notebooks.

```bash
python -m venv .venv
# Windows PowerShell:
.\\.venv\\Scripts\\Activate.ps1

# macOS/Linux:
# source .venv/bin/activate

pip install pandas numpy matplotlib seaborn scipy scikit-learn statsmodels jupyter
pip install mygene statannotations libpysal esda splot
```

The spatial-analysis packages `libpysal`, `esda`, and `splot` are mainly needed for the `cd36\_listeria.ipynb` notebook.

### 2.4 Install R packages

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

\---

## 3\. Important path note

The notebooks were originally written as interactive analysis notebooks. Several notebooks read files using simple relative filenames such as:

```python
pd.read\_pickle("p65\_relb\_cd36.pkl.gz", compression="gzip")
pd.read\_csv("cv\_experiment\_vs\_theory\_p65\_relb\_cd36.csv")
```

If you use the recommended GitHub structure, update the first data-loading cells so that files are read from `data/processed/` or `data/external/`.

For example, if the notebook is inside the `notebooks/` folder:

```python
from pathlib import Path

PROCESSED\_DIR = Path("../data/processed")
EXTERNAL\_DIR = Path("../data/external")
RESULTS\_DIR = Path("../results")
FIGURES\_DIR = Path("../figures")

CSV\_P65\_RELB\_CD36 = PROCESSED\_DIR / "cv\_experiment\_vs\_theory\_p65\_relb\_cd36.csv"
CSV\_TNFR2\_F480 = PROCESSED\_DIR / "cv\_experiment\_vs\_theory\_tnfr2\_f480.csv"
CSV\_TNFA\_IL1B = PROCESSED\_DIR / "cv\_experiment\_vs\_theory\_tnfa\_il1b\_combine.csv"
```

For raw image-analysis notebooks, update pickle paths similarly:

```python
file\_path = EXTERNAL\_DIR / "p65\_relb\_cd36.pkl.gz"
quant\_df = pd.read\_pickle(file\_path, compression="gzip")
```

\---

## 4\. Included processed CSV data

These files are small processed summary tables. They are suitable for GitHub and allow readers to reproduce the final coefficient-of-variation comparison plots without downloading the full raw microscopy dataset.

|File|Rows × columns|Role|Main columns|Used by|
|-|-:|-|-|-|
|`data/processed/Cv\_7\_signals\_FCV.csv`|112 × 11|Combined seven-signal summary table containing experiment/theory CV summaries and fold-change CV values from MemorySeq and scRNA-seq.|`Cell\_Bin`, `Lipid A`, `Mean`, `SD`, `Source`, `Signal`, `Sampling`, `gene\_name`, `FCV MemSeq`, `FCV scRNA-sq`|Final summary/reporting table|
|`data/processed/cv\_experiment\_vs\_theory\_p65\_relb\_cd36.csv`|48 × 6|Experiment-versus-theory CV summary for p65, RelB, and CD36.|`Cell\_Bin`, `Lipid A`, `Mean`, `SD`, `Source`, `Signal`|`notebooks/7gens\_plot.ipynb`|
|`data/processed/cv\_experiment\_vs\_theory\_tnfa\_il1b\_combine.csv`|32 × 6|Experiment-versus-theory CV summary for TNFα and IL-1β.|`Cell\_Bin`, `Lipid A`, `Mean`, `SD`, `Source`, `Signal`|`notebooks/7gens\_plot.ipynb`|
|`data/processed/cv\_experiment\_vs\_theory\_tnfr2\_f480.csv`|32 × 6|Experiment-versus-theory CV summary for TNFR2 and F4/80.|`Cell\_Bin`, `Lipid A`, `Mean`, `SD`, `Source`, `Signal`|`notebooks/7gens\_plot.ipynb`|

### Column definitions

|Column|Meaning|
|-|-|
|`Cell\_Bin`|Colony/cell-number bin. Values include `<4`, `4-6`, `6-8`, and `>8`.|
|`Lipid A`|Stimulation condition. Values include `Mock` and `500 ng/mL`.|
|`Mean`|Mean coefficient-of-variation summary value for the group.|
|`SD`|Standard deviation of the coefficient-of-variation summary value for the group.|
|`Source`|Indicates whether the value comes from `Experiment` or `Theory`.|
|`Signal`|Measured marker/signal, for example `per\_cd36\_mean`, `per\_p65\_mean`, `per\_relb\_mean`, `per\_il1b\_mean`, `per\_tnfa\_mean`, `per\_tnfr2\_mean`, or `per\_f480\_mean`.|
|`Sampling`|Sampling scheme used in the combined seven-signal FCV table.|
|`gene\_name`|Gene name corresponding to the signal.|
|`FCV MemSeq`|Fold-change CV value calculated from MemorySeq data.|
|`FCV scRNA-sq`|Fold-change CV value calculated from scRNA-seq data.|

\---

## 5\. External data expected from BioImage Archive

The following files are required only if readers want to rerun the full raw-data analysis. These files should be deposited in the BioImage Archive and placed locally in `data/external/` after download.

|Expected file or file pattern|Used by|Purpose|
|-|-|-|
|`p65\_relb\_cd36.pkl.gz`|`notebooks/p65\_relb\_cd36\_costainning.ipynb`|Raw or intermediate single-cell image-analysis table for p65, RelB, and CD36 co-staining.|
|`WellsInfo1.pkl`|`notebooks/p65\_relb\_cd36\_costainning.ipynb`|Plate/well layout metadata for the p65/RelB/CD36 analysis.|
|`p65\_relb\_cd36\_clone\_threshold.csv`|`notebooks/p65\_relb\_cd36\_costainning.ipynb`|Clone threshold table used to classify or filter p65/RelB/CD36 clone populations.|
|`tnfr2\_f480.pkl.gz`|`notebooks/tnfr2\_f480\_costainning.ipynb`|Raw or intermediate single-cell image-analysis table for TNFR2 and F4/80 co-staining.|
|`WellsInfo2.pkl`|`notebooks/tnfr2\_f480\_costainning.ipynb`|Plate/well layout metadata for the TNFR2/F4/80 analysis.|
|`IF22\_co.pkl.gz`|`notebooks/tnfa\_il1b\_costainning.ipynb`|Raw or intermediate single-cell image-analysis table for one TNFα/IL-1β experiment.|
|`IF24\_co.pkl.gz`|`notebooks/tnfa\_il1b\_costainning.ipynb`|Raw or intermediate single-cell image-analysis table for one TNFα/IL-1β experiment.|
|`IF48\_listeria.pkl.gz`|`notebooks/cd36\_listeria.ipynb`|Raw or intermediate single-cell image-analysis table for one CD36/Listeria experiment.|
|`IF50\_listeria.pkl.gz`|`notebooks/cd36\_listeria.ipynb`|Raw or intermediate single-cell image-analysis table for one CD36/Listeria experiment.|
|`IF54\_listeria.pkl.gz`|`notebooks/cd36\_listeria.ipynb`|Raw or intermediate single-cell image-analysis table for one CD36/Listeria experiment.|
|`\*\_quant.sf`|`notebooks/Density-CD36 IL1-programmes.ipynb`|Transcript-level quantification files used to build gene-level expression summaries.|
|`scRNA\_filtered\_genes.pkl`|`notebooks/NMF\_CD36\_IL1.ipynb`|Filtered scRNA-seq gene-expression matrix used for NMF/program analysis.|
|`MemorySeq data.xlsx`|`R\_scripts/cv\_fold\_change\_MemorySeq.R`|MemorySeq expression workbook and gene-list workbook sheets.|
|`combined\_logcounts\_JM01.xlsx` to `combined\_logcounts\_JM12.xlsx`|`R\_scripts/cv\_fold\_change\_scRNAseq.R`|scRNA-seq logcount matrices for clonal and mixed samples.|
|`IL! and CD36 fractions.xlsx`|`R\_scripts/Constant\_Switching\_Parameter\_Simulation.R`, `R\_scripts/Density\_Dependent\_Parameter\_Simulation.R`, `R\_scripts/Density\_Dependent\_Parameter\_Simulation\_kmax\_beta.R`, `R\_scripts/Optimized\_Parameter\_Estimation.R`, `R\_scripts/Constant\_Switching\_Rates\_Calculation.Rmd`|Experimental fraction data used for parameter estimation and simulation benchmarking.|
|`Corrected\_odds\_ratio\_CD36\_tnfa\_il1b.xlsx`|`R\_scripts/Constant\_Switching\_Rates\_Calculation.Rmd`|Odds-ratio workbook used to estimate switching rates.|

If filenames are changed during BioImage Archive upload or download, update the corresponding path variables in the notebook or R script.

\---

## 6\. How to reproduce the final seven-signal CV comparison figure

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
   notebooks/7gens\_plot.ipynb
   ```

4. Update the CSV path variables if needed:

```python
   from pathlib import Path

   PROCESSED\_DIR = Path("../data/processed")

   CSV\_P65\_RELB\_CD36 = PROCESSED\_DIR / "cv\_experiment\_vs\_theory\_p65\_relb\_cd36.csv"
   CSV\_TNFR2\_F480 = PROCESSED\_DIR / "cv\_experiment\_vs\_theory\_tnfr2\_f480.csv"
   CSV\_TNFA\_IL1B = PROCESSED\_DIR / "cv\_experiment\_vs\_theory\_tnfa\_il1b\_combine.csv"

   SVG\_PATH\_GRID = "../figures/cv\_grid.svg"
   SVG\_PATH\_SINGLE = "../figures/cv\_single.svg"
   ```

5. Run all cells.

Expected output:

```text
figures/cv\_grid.svg
```

This figure compares experimental and theoretical CV summaries across the seven measured signals.

\---

## 7\. Notebook guide

|Notebook|Main role|Main input data|Main outputs|
|-|-|-|-|
|`notebooks/7gens\_plot.ipynb`|Generates the final seven-signal experiment-versus-theory CV plot.|Processed CSV files in `data/processed/`.|`cv\_grid.svg`, optional single-signal SVG plots.|
|`notebooks/p65\_relb\_cd36\_costainning.ipynb`|Processes p65, RelB, and CD36 co-staining image-analysis data; performs colony/clone summaries and CV comparison.|`p65\_relb\_cd36.pkl.gz`, `WellsInfo1.pkl`, `p65\_relb\_cd36\_clone\_threshold.csv`.|Cleaned pickle files, clone morphology tables, average signal tables, `cv\_experiment\_vs\_theory\_p65\_relb\_cd36.csv`, and SVG figures.|
|`notebooks/tnfa\_il1b\_costainning.ipynb`|Processes TNFα and IL-1β co-staining image-analysis data.|`IF22\_co.pkl.gz`, `IF24\_co.pkl.gz`.|Cleaned/intermediate pickle files for TNFα/IL-1β analysis.|
|`notebooks/tnfr2\_f480\_costainning.ipynb`|Processes TNFR2 and F4/80 co-staining image-analysis data.|`tnfr2\_f480.pkl.gz`, `WellsInfo2.pkl`.|Cleaned/intermediate pickle files and plate overview figures.|
|`notebooks/cd36\_listeria.ipynb`|Processes CD36/Listeria perinuclear signal data and spatial statistics.|`IF48\_listeria.pkl.gz`, `IF50\_listeria.pkl.gz`, `IF54\_listeria.pkl.gz`.|`cd36\_Listeria\_clean.pkl` and spatial/statistical plots.|
|`notebooks/Density-CD36 IL1-programmes.ipynb`|Builds gene-level expression tables from transcript quantification files.|Files matching `\*\_quant.sf`.|`df\_gene.pkl`, `df\_gene.csv`.|
|`notebooks/NMF\_CD36\_IL1.ipynb`|Performs latent program/NMF analysis for CD36 and IL1B-related gene programs.|`scRNA\_filtered\_genes.pkl`.|NMF/program annotation figures including SVG, PNG, and PDF outputs.|

\---

## 8\. R script guide

|Script|Main role|Main input data|Main outputs|
|-|-|-|-|
|`R\_scripts/Constant\_Switching\_Rates\_Calculation.Rmd`|Estimates constant switching rates from odds-ratio and fraction workbooks.|`Corrected\_odds\_ratio\_CD36\_tnfa\_il1b.xlsx`, `IL! and CD36 fractions.xlsx`.|Rendered HTML report and estimated rate summaries printed in the report.|
|`R\_scripts/Optimized\_Parameter\_Estimation.R`|Optimizes density-dependent switching parameters by matching simulated and experimental mean/CV profiles.|`IL! and CD36 fractions.xlsx`.|Optimized parameter values printed to the console.|
|`R\_scripts/Constant\_Switching\_Parameter\_Simulation.R`|Runs Gillespie simulations with constant ON/OFF switching rates and compares with experimental fraction data.|`IL! and CD36 fractions.xlsx`.|Simulation plots; optional CSV export if the commented `write.csv()` line is enabled.|
|`R\_scripts/Density\_Dependent\_Parameter\_Simulation.R`|Runs density-dependent switching simulations with a Hill-type activation function.|`IL! and CD36 fractions.xlsx`.|`CD36\_simulated\_density\_constant\_rates.csv`, `simulated\_loess\_trend\_data\_kON\_constant\_wCI.csv`.|
|`R\_scripts/Density\_Dependent\_Parameter\_Simulation\_kmax\_beta.R`|Runs density-dependent simulations with Beta-distributed maximum activation rates.|`IL! and CD36 fractions.xlsx`.|`CD36\_simulated\_density\_beta\_rates\_doublerate14\_4.csv`, `simulated\_loess\_trend\_data\_kON\_constant\_wCI.csv`, `experimental\_loess\_trend\_data\_kON\_constant\_wCI.csv`.|
|`R\_scripts/Constant\_Switching\_colony\_tracking\_Simulation\_kmax\_beta.R`|Simulates lineage-resolved colony growth with constant phenotypic switching rates.|No external file required; parameters are defined inside the script.|Lineage tracking table in memory and plots printed to the R graphics device.|
|`R\_scripts/Density\_Dependent\_colony\_tracking\_Simulation\_kmax\_beta.R`|Simulates lineage-resolved colony growth with density-dependent switching rates.|No external file required; parameters are defined inside the script.|Lineage tracking table in memory and plots printed to the R graphics device; optional CSV export if the commented `write.csv()` line is enabled.|
|`R\_scripts/cv\_fold\_change\_MemorySeq.R`|Calculates MemorySeq CV fold-change and heritability classification for gene sets.|`data/MemorySeq data.xlsx` or `data/external/MemorySeq data.xlsx`, depending on the chosen path convention.|`results/CV\_fold\_change\_heritability\_MemorySeq.xlsx`, missing-gene text files, and `sessionInfo.txt`.|
|`R\_scripts/cv\_fold\_change\_scRNAseq.R`|Calculates scRNA-seq CV fold-change between clonal and mixed populations.|`combined\_logcounts\_JM01.xlsx` to `combined\_logcounts\_JM12.xlsx`.|`results/scRNAseq\_CV\_fold\_change\_clonal\_vs\_mixed.xlsx`.|

\---

## 9\. Running the R scripts

Run R scripts from the repository root unless otherwise noted.

### Example: MemorySeq CV fold-change analysis

Expected input:

```text
data/MemorySeq data.xlsx
```

Run:

```bash
Rscript R\_scripts/cv\_fold\_change\_MemorySeq.R
```

Expected output:

```text
results/CV\_fold\_change\_heritability\_MemorySeq.xlsx
results/sessionInfo.txt
```

### Example: scRNA-seq CV fold-change analysis

Expected inputs:

```text
data/combined\_logcounts\_JM01.xlsx
data/combined\_logcounts\_JM02.xlsx
...
data/combined\_logcounts\_JM12.xlsx
```

Run:

```bash
Rscript R\_scripts/cv\_fold\_change\_scRNAseq.R
```

Expected output:

```text
results/scRNAseq\_CV\_fold\_change\_clonal\_vs\_mixed.xlsx
```

### Example: render the rate-calculation R Markdown file

```bash
Rscript -e "rmarkdown::render('R\_scripts/Constant\_Switching\_Rates\_Calculation.Rmd')"
```

### Important R path note

Some R scripts currently contain local Windows paths, for example:

```r
DATA\_PATH <- "C:/Users/apurv/OneDrive/Desktop/Data Analysis/rate switch/IL! and CD36 fractions.xlsx"
```

For publication, replace these with relative paths, for example:

```r
DATA\_PATH <- file.path("data", "external", "IL! and CD36 fractions.xlsx")
```

or, if you decide to keep Excel inputs directly in `data/`:

```r
DATA\_PATH <- file.path("data", "IL! and CD36 fractions.xlsx")
```

\---

## 10\. Suggested order for full reproducibility

The exact order depends on which part of the paper the reader wants to reproduce.

### A. Reproduce final summary plots only

Use this route if the reader only wants to recreate the final seven-signal CV plot.

```text
data/processed/\*.csv
        ↓
notebooks/7gens\_plot.ipynb
        ↓
figures/cv\_grid.svg
```

### B. Rerun image-analysis notebooks from large BioImage Archive files

Use this route if the reader wants to regenerate cleaned single-cell/clone tables from image-analysis outputs.

```text
data/external/\*.pkl.gz and WellsInfo\*.pkl
        ↓
notebooks/\*costainning.ipynb and notebooks/cd36\_listeria.ipynb
        ↓
results/\*.pkl, results/\*.csv, figures/\*.svg
```

### C. Rerun stochastic switching simulations

Use this route if the reader wants to reproduce model simulations and parameter-estimation workflows.

```text
data/external/IL! and CD36 fractions.xlsx
data/external/Corrected\_odds\_ratio\_CD36\_tnfa\_il1b.xlsx
        ↓
R\_scripts/Constant\_Switching\_Rates\_Calculation.Rmd
R\_scripts/Optimized\_Parameter\_Estimation.R
R\_scripts/\*Simulation\*.R
        ↓
simulation summary CSVs and model-comparison figures
```

### D. Rerun MemorySeq/scRNA-seq CV fold-change analyses

Use this route if the reader wants to reproduce the FCV values reported in the combined seven-signal table.

```text
data/external/MemorySeq data.xlsx
data/external/combined\_logcounts\_JM\*.xlsx
        ↓
R\_scripts/cv\_fold\_change\_MemorySeq.R
R\_scripts/cv\_fold\_change\_scRNAseq.R
        ↓
results/\*.xlsx
```

\---

## 11\. Suggested `.gitignore`

Use a `.gitignore` file to prevent raw or temporary files from being committed accidentally:

```text
# External large files downloaded from BioImage Archive
data/external/

# Generated outputs
results/
figures/

# Jupyter temporary files
.ipynb\_checkpoints/

# Python temporary files
\_\_pycache\_\_/
\*.pyc

# R temporary files
.Rhistory
.RData
.Rproj.user/

# System files
.DS\_Store
Thumbs.db

# Logs
\*.log

# Large raw or intermediate data
\*.pkl
\*.pkl.gz
\*.h5
\*.hdf5
\*.tif
\*.tiff
\*.czi
\*.nd2
\*.lif
\*.zip
\*.tar.gz
```

The four processed CSV files in `data/processed/` should remain tracked by Git.

\---

## 12\. Notes for users

1. The processed CSV files in `data/processed/` are sufficient to reproduce the final CV comparison plot.
2. The full raw-data workflow requires external files from BioImage Archive.
3. Some notebooks and R scripts may require minor path edits because they were originally developed interactively.
4. For exact reproducibility, keep the downloaded BioImage Archive filenames unchanged, or update the path variables consistently.
5. Generated output files should be written to `results/` or `figures/`, not mixed with source code.

\---

## 13\. Citation

If you use this code, processed data, or the associated external data, please cite the associated publication:

```text
DOI: 10.64898/2026.06.17.732820
```

Publication DOI: [https://doi.org/10.64898/2026.06.17.732820](https://doi.org/10.64898/2026.06.17.732820)

Data availability:

```text
Raw microscopy and large analysis files are available from the BioImage Archive:
https://ftp.ebi.ac.uk/biostudies/fire/S-BIAD/517/S-BIAD2517/Files/external\_data.zip
```

