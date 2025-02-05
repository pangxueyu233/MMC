# MMC

This page recorded the codes and data used and mentioned in [*xxx*](XXX). And you could downloaded this paper by clicking [here](pdf/XXX)



abstract .

To effectively demonstrate our step-by-step analysis of single-cell RNA sequencing (scRNA-seq), single-cell T-cell receptor sequencing (scTCR-seq), and single-cell spatial transcriptomics (scSpatial), we have meticulously compiled and stored detailed procedural information in Markdown files. This extensive documentation encompasses records of quality control measures, batch effect reduction, dimensionality reduction, cell clustering, pseudotime construction, identification of dynamically expressed genes, and pathway enrichment.

# **1. Codes of analyzing and visualization**

**Introduction to Our Script Compilation for Analysis**

In our comprehensive analysis of antifungal treatment in the MMC cohort, we have organized our scripts into six distinct chapters. Each chapter focuses on a specific aspect of the analysis, enabling a detailed examination of different facets of antifungal treatment effects in IBD. Below is a guide to the content of each chapter:

- **[Chapter 0](Chapter0.md): Pre-processing of Microbiome and Metabolome Data**
  
  - Content: Initial processing steps for **multi-omics data (ITS, metagenomics, metabolomics)**, including data cleaning, normalization, and quality control.
  
- **[Chapter 1](Chapter1.md): Microbiome and Metabolome Data Integration and Processing**

  - **Content:** This chapter integrates **multi-omics data (ITS, metagenomics, metabolomics)** with **longitudinal clinical metadata** to examine **microbiome dynamics in response to treatment**.

  **[Chapter 2](Chapter2.md): Fungal Abundance Data Analysis**

  - **Content:** This section integrates **fungal abundance data (ITS sequencing) with clinical metadata** to assess **treatment-induced shifts in fungal communities**.

  **[Chapter 3](Chapter3.md): Bacterial Analysis – Data Processing, Visualization, and Statistical Testing**

  - **Content:** This analysis incorporates **ITS sequencing and metagenomic profiling** from the **MMC cohort**, focusing on **treatment responses and gut microbial shifts**. It provides additional evidence on **how different antifungal treatments influence the gut bacterial community**.

  **[Chapter 4](Chapter4.md): Metabolomic Data Analysis**

  - **Content:** This section investigates **metabolomic shifts under antifungal treatments** (Nystatin, Clotrimazole, and Fluconazole).

  **[Chapter 5](Chapter5.md): Clinical Data Analysis**

  - **Content:** This chapter presents **clinical outcome analyses**, including **disease scores (DAI, UC, CD), patient metadata, and their association with antifungal treatments**.


Each chapter contains detailed scripts, methodologies, and analyses relevant to the specific aspect of antifungal treatment it addresses. This structured approach allows researchers to navigate our comprehensive analysis with ease, enhancing their understanding of antifungal treatment effects in IBD.

# **2. Raw data download**

- **Description**: This section includes all the raw FASTQ files from our study. These files are crucial for in-depth data analysis and understanding the sequencing results from single-cell RNA, T-cell receptor, and spatial transcriptomics.

- **Download**: You can access and download these files from the [GEO database](https://chat.openai.com/c/link-to-download).

Below is a detailed annotation of the file structure and contents:

```shell
[4.0K]  .
├── [4.0K]  scRNA
│   ├── [ 44G]  ACC10_RNA_S1_L001_R1_001.fastq.gz
3 directories, 42 files
```

# **3. Processed Data Download**

## 3.1. CellRanger and SpaceRanger Output

- **Description**: This section includes the output files from Cell Ranger and Space Ranger, essential for the initial data processing and analysis of single-cell RNA, T-cell receptor, and spatial transcriptomics data.
- **Download**: These files are available for access and download from the [GEO database](https://chat.openai.com/c/link-to-download).

Below is a detailed annotation of the file structure and contents:

```shell
plaintextCopy code[4.0K]  .
├── [4.0K]  cellranger_output
│   ├── [424M]  ACC10_RNA.tar.gz
```

- Contents
  - Each `_RNA.tar.gz` file includes the filtered_feature_bc_matrix output and loupe file from the Cell Ranger count model.
  - Each `_TCR.tar.gz` file contains the filtered_contig_annotations.csv, clonotypes.csv, and loupe.vloupe files generated from the Cell Ranger TCR model.
  - Each `_spatial.tar.gz` file includes the filtered_feature_bc_matrix output and loupe.cloupe file from the Space Ranger count model.

## 3.2. R Data Files Generated in This Study

- **Description**: All R data files (.rds) related to single-cell RNA sequencing (scRNA-seq), single-cell T-cell receptor sequencing (scTCR-seq), single-cell spatial transcriptomics (scSpatial), and The Cancer Genome Atlas Adrenocortical Carcinoma (TCGA-ACC) data are available. These files encompass a comprehensive range of analyses and findings from our study.
- **Download**: You can download these files from [Zenodo Zenodo10416598](https://zenodo.org/records/10416598).

Here's an annotation for each file to give you :

~~~shell
tree -lh
[4.0K]  .
├── [ 16M]  CCI.cor_genes.GSEA.rds
CCI.cor_genes.GSEA.rds: Gene Set Enrichment Analysis (GSEA) results for genes correlated with Confused Cell Identity (CCI) in ACC.

├── [791K]  CCI.cor_genes.rds
CCI.cor_genes.rds: List or data frame of genes correlated with Confused Cell Identity in ACC.

0 directories, 33 files
~~~

Each file seems to contain specific data subsets or analysis results, crucial for a comprehensive understanding of ACC and normal adrenal tissues at the single-cell level.

# **Citation**

Our paper has been published in [*XXX Journal*](https://chat.openai.com/c/xxxx). For further reference and details, you can access the publication at the provided link.

The raw data supporting the findings of this study can be downloaded from the following repositories:

- **GEO Database**: Access our dataset by visiting [GSEXXX](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSEXXX). This link will take you directly to the dataset's page.
- **Zenodo**: Additional data files are available on Zenodo. Download them at [Zenodo10416598](https://zenodo.org/records/10416598).