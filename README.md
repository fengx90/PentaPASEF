# PentaPASEF Identification Filtering and Performance Scoring

This repository contains two R scripts for analyzing **DIA-NN** results from a **PentaPASEF Classic DIA** experiment (example shown for the **45 min gradient**) and for computing **1–5 performance scores** for profiling depth and quantitative accuracy.

## System requirements

### Operating system
- Tested/assumed: **Windows 11** (based on the original absolute paths)
- Should also work on **macOS/Linux** if you replace `setwd()` / file paths with platform-appropriate paths.

### Packages used (Script 1)
- `arrow` (read parquet)
- `diann` (DIA-NN helper functions: `diann_matrix()`, `diann_maxlfq()`)
- `dplyr`, `tidyr`, `stringr` (data manipulation)
- `readxl`, `readr` (I/O)
- `ggplot2`, `ggpubr`, `ggbreak` (plots)
- `dplyr` (scoring calculations)

### Install packages
Install CRAN packages:
```r
install.packages(c(
  "arrow",
  "dplyr",
  "tidyr",
  "stringr",
  "readxl",
  "readr",
  "ggplot2",
  "ggpubr",
  "ggbreak"
))
### R environment
- **R >= 4.1** recommended
- RStudio optional but recommended

## Software / external tools

- **DIA-NN** to generate `report.parquet` (export your DIA-NN report as parquet).
  - DIA-NN repository: https://github.com/vdemichev/diann
---

# DIA-PASEF (PentaPASEF Classic DIA) — 45 min Gradient (DIA-NN → R analysis)
This code contains an R workflow to post-process **DIA-NN** results for a **PentaPASEF Classic DIA** experiment (focused here on the **45 min gradient**). The script:

- Reads DIA-NN output (`report.parquet`)
- Renames DIA-NN run names to real sample names using an Excel mapping
- Removes contaminants and decoys
- Annotates peptides/proteins into **L. murinus**, **S. ruber**, or **background**
- Excludes false assignments based on negative control (`NSCtrl`) samples
- Generates quantitative matrices at **1% FDR**:
  - precursor-level matrix
  - peptide-level (Modified.Sequence) **MaxLFQ** matrix
  - protein group **MaxLFQ** matrix
- Summarizes IDs (precursor/peptide counts)
- Computes CVs (injection replicate and dilution replicate)
- Produces dynamic range and quantitative accuracy benchmarking plots
- Exports multiple CSV/TSV outputs, including a **MetaLab-ready** peptide intensity table

## How to run

Export DIA-NN report as parquet: report.parquet
Create/update mapping file 45min.xlsx (metadata, FileName, NewSampleName)
Open R/RStudio
Set working directory (or update the script to use project-relative paths)
Run the script top-to-bottom

## Outputs
1: False-assignment exclusion
precursors_to_exclude_basedOn_NGCtrl.csv

2: Precursor-level
Filtered_Precursors_all_FDR0.01.csv
Precursor_IDs.csv

3: Peptide-level (Modified.Sequence MaxLFQ)
Filtered_peptides_Quant_all_FDR0.01_45min.csv
ModifiedPeptides_IDs_45min.csv
Filtered_ModifiedPeptides_IDs_LM_45min-FAST.csv
Filtered_ModifiedPeptides_IDs_SR_45min.csv
Filtered_combined_peptide_count_LM_SR_45min.tsv

4: CV benchmarking
45min-FAST_peptide_averaged_intensity_by_injection_CV.csv
45min-FAST_peptide_averaged_intensity_by_Dilution_CV.csv

5: Quant accuracy
DIA_45min_QuantitativeAccuracy_Dil1-Dil2.csv

6: Taxa annotation / Unipept prep
45min_Injection_Averaged_Intensity.csv
45min_Injection_Averaged_Intensity_WO_HumanPeptide.csv
45min_Injection_Averaged_Intensity_For_MetaLab.csv

7: Protein group
Filtered_ProteinGroup_Quant_FDR0.01.csv
ProteinGroup_45min_Injection_Averaged_Intensity_CV.csv
ProteinGroup_45min_Injection_Averaged_Intensity.csv

---

# Scoring: Profiling Depth & Quantitative Accuracy (Peptides / Proteins)

This code computes **performance scores (1–5 scale)** for different DIA methods/gradients based on:

1. **Profiling depth** (mean IDs) penalized by variability (CV)
   - peptide-level scoring from `total_summary`
   - protein-level scoring from `total_summary_protein`
   - spike-in peptide depth scoring from `spike_df`

2. **Quantitative accuracy scoring** for spike-in species (non-background) using:
   - absolute log2-ratio error vs expected ratio
   - stability penalty (CV of error)
   - **coverage penalty** (fewer peptides → lower score)

All scoring is normalized **within each Gradient** (and in some cases within `Gradient × Group × Dilution_Replicate`) to produce a comparable 1–5 score.

---

## Inputs (required objects)

This code assumes the following data frames already exist in the R environment:

### 1) `total_summary` (peptide depth summary)
Must contain at least:
- `Gradient`
- `mean_peptides`  (mean peptide IDs or mean peptide counts across replicates)
- `sd_peptides`    (SD across replicates)

### 2) `total_summary_protein` (protein depth summary)
Must contain at least:
- `Gradient`
- `Origin`         (e.g., background vs spike-in, or source category)
- `mean_protein`
- `sd_protein`

### 3) `spike_df` (spike-in peptide ID counts)
Must contain at least:
- `Method`
- `Gradient`
- `Group`              (e.g., species group like *L.murinus* / *S.ruber*)
- `Dilution_Replicate`
- `Peptide_Count`

### 4) `data1` (quant accuracy input table; spike-in only)
Used for quantitative accuracy scoring; expected columns:
- `Method`
- `Gradient`
- `species`            (must include `background` and spike-in species labels)
- `ModifiedPeptide`    (for peptide universe / distinct counting)
- `log2ratio`          (observed log2 ratio)
- `ExpectedRatio`      (expected log2 ratio)

---
##
scoring.csv
scoring_protein.csv
spike_scoring.csv
spike_mean_scores.csv
accuracy_scoring.csv
---
## Interpretation of scores
- Scores range from 1 (worst) to 5 (best) after normalization.
- For depth scoring, higher mean IDs and lower variability (CV) improve score.
- For accuracy scoring, smaller deviation from expected ratios improves score, but the score is penalized if only a small fraction of spike-in peptides are quantified (coverage penalty).
---
