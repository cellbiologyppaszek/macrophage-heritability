# Heritable single-cell gene expression analysis

This repository contains code and processed data used for the analysis of heritable single-cell gene expression programs in macrophage populations.

The repository is intended to help readers understand, inspect, and reproduce the computational analysis associated with the publication.

Large raw image files, large intermediate `.pkl` / `.pkl.gz` files, and sequencing input files are not stored in this GitHub repository. These files should be downloaded separately from the BioImage Archive and placed locally in `data/external/`.

BioImage Archive accession: `TO\_BE\_ADDED`

## Repository structure

```text
official\_code/
|
|-- README.md
|-- .gitignore
|-- .gitattributes
|
|-- data/
|   |
|   |-- processed/
|   |   |-- Cv\_7\_signals\_FCV.csv
|   |   |-- cv\_experiment\_vs\_theory\_p65\_relb\_cd36.csv
|   |   |-- cv\_experiment\_vs\_theory\_tnfa\_il1b\_combine.csv
|   |   |-- cv\_experiment\_vs\_theory\_tnfr2\_f480.csv
|   |
|   |-- external/
|       |-- README\_external\_data.md
|       |
|       | Large raw or intermediate input files downloaded from the
|       | BioImage Archive should be placed here locally.
|       | This folder should not be committed to GitHub.
|
|-- notebooks/
|   |-- 7gens\_plot.ipynb
|   |-- cd36\_listeria.ipynb
|   |-- Density-CD36 IL1-programmes.ipynb
|   |-- NMF\_CD36\_IL1.ipynb
|   |-- p65\_relb\_cd36\_costainning.ipynb
|   |-- tnfa\_il1b\_costainning.ipynb
|   |-- tnfr2\_f480\_costainning.ipynb
|
|-- R\_scripts/
|   |-- Constant\_Switching\_colony\_tracking\_Simulation\_kmax\_beta.R
|   |-- Constant\_Switching\_Parameter\_Simulation.R
|   |-- Constant\_Switching\_Rates\_Calculation.Rmd
|   |-- Density\_Dependent\_colony\_tracking\_Simulation\_kmax\_beta.R
|   |-- Density\_Dependent\_Parameter\_Simulation.R
|   |-- Density\_Dependent\_Parameter\_Simulation\_kmax\_beta.R
|   |-- Optimized\_Parameter\_Estimation.R
|   |-- cv\_fold\_change\_MemorySeq.R
|   |-- cv\_fold\_change\_scRNAseq.R
|
|-- results/
|   | Generated tables and intermediate output files.
|   | This folder is not committed to GitHub.
|
|-- figures/
    | Generated SVG, PNG, or PDF figures.
    | This folder is not committed to GitHub.
```

## What is included in GitHub

This GitHub repository includes:

|Folder|Content|Purpose|
|-|-|-|
|`data/processed/`|Small processed CSV files|Processed summary data used for plotting and comparison|
|`notebooks/`|Jupyter notebooks|Python analysis and figure-generation workflows|
|`R\_scripts/`|R and R Markdown scripts|Simulation, parameter estimation, and CV fold-change analysis|
|`README.md`|Main documentation|Instructions for readers|
|`.gitignore`|Git exclusion rules|Prevents raw and large files from being uploaded|
|`.gitattributes`|Line-ending rules|Keeps text files consistent across operating systems|

## What is not included in GitHub

The following files should not be committed to GitHub:

|File type|Example|Where to store|
|-|-|-|
|Raw microscopy images|`.tif`, `.tiff`, `.czi`, `.nd2`, `.lif`|BioImage Archive|
|Large intermediate image-analysis files|`.pkl`, `.pkl.gz`|BioImage Archive|
|Large sequencing input files|`\*\_quant.sf`, large `.xlsx`, large `.pkl`|BioImage Archive|
|Temporary files|`.ipynb\_checkpoints/`, `.Rhistory`, `\_\_pycache\_\_/`|Do not archive|
|Generated local output|`results/`, `figures/`|Regenerate locally|

Locally, these external files should be placed under:

```text
data/external/
```

## Processed data files included in GitHub

The following processed CSV files are included in `data/processed/`.

|File|Role|
|-|-|
|`Cv\_7\_signals\_FCV.csv`|Processed coefficient-of-variation and fold-change data for multiple measured signals|
|`cv\_experiment\_vs\_theory\_p65\_relb\_cd36.csv`|Processed comparison between experimental and theoretical CV values for p65, RelB, and CD36|
|`cv\_experiment\_vs\_theory\_tnfa\_il1b\_combine.csv`|Processed comparison between experimental and theoretical CV values for TNFA and IL1B|
|`cv\_experiment\_vs\_theory\_tnfr2\_f480.csv`|Processed comparison between experimental and theoretical CV values for TNFR2 and F4/80|

These files are small enough to be stored directly in GitHub and can be used by readers without downloading the full raw BioImage Archive dataset.

## External data expected from BioImage Archive

To rerun the complete analysis from large raw or intermediate files, download the required files from the BioImage Archive and place them under `data/external/`.

Recommended local structure:

```text
data/external/
|
|-- raw\_image\_files/
|   | Raw microscopy image files
|
|-- pkl\_intermediate\_files/
|   | Image-analysis intermediate tables
|
|-- sequencing\_inputs/
|   | Sequencing or transcriptomic input files
|
|-- excel\_inputs/
    | Excel files used by R scripts
```

### Image-analysis intermediate files

Place these files in:

```text
data/external/pkl\_intermediate\_files/
```

|File|Used by|
|-|-|
|`p65\_relb\_cd36.pkl.gz`|`notebooks/p65\_relb\_cd36\_costainning.ipynb`|
|`WellsInfo1.pkl`|`notebooks/p65\_relb\_cd36\_costainning.ipynb`|
|`tnfr2\_f480.pkl.gz`|`notebooks/tnfr2\_f480\_costainning.ipynb`|
|`WellsInfo2.pkl`|`notebooks/tnfr2\_f480\_costainning.ipynb`|
|`IF22\_co.pkl.gz`|`notebooks/tnfa\_il1b\_costainning.ipynb`|
|`IF24\_co.pkl.gz`|`notebooks/tnfa\_il1b\_costainning.ipynb`|
|`IF48\_listeria.pkl.gz`|`notebooks/cd36\_listeria.ipynb`|
|`IF50\_listeria.pkl.gz`|`notebooks/cd36\_listeria.ipynb`|
|`IF54\_listeria.pkl.gz`|`notebooks/cd36\_listeria.ipynb`|

### Sequencing input files

Place these files in:

```text
data/external/sequencing\_inputs/
```

|File or file pattern|Used by|
|-|-|
|`\*\_quant.sf`|`notebooks/Density-CD36 IL1-programmes.ipynb`|
|`scRNA\_filtered\_genes.pkl`|`notebooks/NMF\_CD36\_IL1.ipynb`|
|`MemorySeq data.xlsx`|`R\_scripts/cv\_fold\_change\_MemorySeq.R`|
|`combined\_logcounts\_JM01.xlsx` to `combined\_logcounts\_JM12.xlsx`|`R\_scripts/cv\_fold\_change\_scRNAseq.R`|

### Excel input files for model fitting and simulation

Place these files in:

```text
data/external/excel\_inputs/
```

|File|Used by|
|-|-|
|`IL1 and CD36 fractions.xlsx`|`R\_scripts/Constant\_Switching\_Parameter\_Simulation.R`, `R\_scripts/Density\_Dependent\_Parameter\_Simulation.R`, `R\_scripts/Density\_Dependent\_Parameter\_Simulation\_kmax\_beta.R`, `R\_scripts/Optimized\_Parameter\_Estimation.R`, `R\_scripts/Constant\_Switching\_Rates\_Calculation.Rmd`|
|`Corrected\_odds\_ratio\_CD36\_tnfa\_il1b.xlsx`|`R\_scripts/Constant\_Switching\_Rates\_Calculation.Rmd`|

Note: Please check whether the filename `IL! and CD36 fractions.xlsx` is intentional. If the intended name is `IL1 and CD36 fractions.xlsx` or `IL-1 and CD36 fractions.xlsx`, update the filename and the scripts consistently.

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
conda create -n macrophage\_expression python=3.10
conda activate macrophage\_expression
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

### 1\. Clone the repository

```bash
git clone https://github.com/cellbiologyppaszek/Heritable-single-cell-gene-expression-ppaszek.git
cd Heritable-single-cell-gene-expression-ppaszek
```

### 2\. Download external files

Download the external raw and intermediate files from the BioImage Archive.

Place them locally under:

```text
data/external/
```

For example:

```text
data/external/pkl\_intermediate\_files/
data/external/sequencing\_inputs/
data/external/excel\_inputs/
```

Do not commit these external files to GitHub.

### 3\. Run the Jupyter notebooks

Start Jupyter from the repository root:

```bash
jupyter notebook
```

Open notebooks from the `notebooks/` folder.

Recommended order:

|Step|Notebook|Purpose|
|-|-|-|
|1|`notebooks/7gens\_plot.ipynb`|Plot and summarize processed CV and fold-change data|
|2|`notebooks/p65\_relb\_cd36\_costainning.ipynb`|Analyze p65, RelB, and CD36 co-staining data|
|3|`notebooks/tnfa\_il1b\_costainning.ipynb`|Analyze TNFA and IL1B co-staining data|
|4|`notebooks/tnfr2\_f480\_costainning.ipynb`|Analyze TNFR2 and F4/80 co-staining data|
|5|`notebooks/cd36\_listeria.ipynb`|Analyze CD36-related Listeria infection data|
|6|`notebooks/Density-CD36 IL1-programmes.ipynb`|Analyze CD36 and IL1-related transcriptional programs|
|7|`notebooks/NMF\_CD36\_IL1.ipynb`|Perform NMF analysis of CD36 and IL1-associated programs|

Some notebooks may need path updates depending on where the BioImage Archive files are placed locally.

For example, if a notebook currently reads:

```python
pd.read\_pickle("p65\_relb\_cd36.pkl.gz")
```

change it to:

```python
pd.read\_pickle("../data/external/pkl\_intermediate\_files/p65\_relb\_cd36.pkl.gz")
```

If a notebook reads a processed CSV file from the current folder, change it to:

```python
pd.read\_csv("../data/processed/Cv\_7\_signals\_FCV.csv")
```

when running from the `notebooks/` folder.

### 4\. Run the R scripts

Run R scripts from the repository root, or update paths inside the scripts.

Example from command line:

```bash
Rscript R\_scripts/Constant\_Switching\_Parameter\_Simulation.R
Rscript R\_scripts/Density\_Dependent\_Parameter\_Simulation.R
Rscript R\_scripts/Optimized\_Parameter\_Estimation.R
Rscript R\_scripts/cv\_fold\_change\_MemorySeq.R
Rscript R\_scripts/cv\_fold\_change\_scRNAseq.R
```

For the R Markdown file:

```bash
Rscript -e "rmarkdown::render('R\_scripts/Constant\_Switching\_Rates\_Calculation.Rmd')"
```

Some R scripts may contain local absolute paths from the original analysis computer. Replace them with relative paths such as:

```r
data\_path <- "data/external/excel\_inputs/IL! and CD36 fractions.xlsx"
```

or:

```r
data\_path <- "data/processed/Cv\_7\_signals\_FCV.csv"
```

## Expected outputs

Generated outputs should be written locally to:

```text
results/
figures/
```

These folders are for local results only and should not be committed unless specific final figures or tables are intentionally included.

|Folder|Content|
|-|-|
|`results/`|Generated tables, model outputs, intermediate processed files|
|`figures/`|Generated plots in SVG, PNG, or PDF format|

## Recommended `.gitignore`

The repository should include a `.gitignore` file to prevent external or large files from being committed.

Recommended content:

```text
# External BioImage Archive data
data/external/

# Raw microscopy files
\*.tif
\*.tiff
\*.czi
\*.nd2
\*.lif

# Large intermediate files
\*.pkl
\*.pkl.gz
\*.h5
\*.hdf5

# Sequencing and compressed files
\*.zip
\*.tar.gz
\*.rar
\*.7z

# Local results
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
TO\_BE\_ADDED
```

## Contact

For questions about the code or data organization, please contact the corresponding author of the associated publication.

