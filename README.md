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

---



# Scoring: Profiling Depth & Quantitative Accuracy (Peptides / Proteins)

This module computes **performance scores (1–5 scale)** for different DIA methods/gradients based on:

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

> The code filters out background with:
> `data1[!grepl("background", data1$species), ]`

---
