# MMC

This page recorded the codes and data used and mentioned in [*xxx*](XXX). And you could downloaded this paper by clicking [here](pdf/XXX)



abstract .

Gut fungal dysbiosis and the presence of pathogenic fungal strains is increasingly recognized as a contributing factor to inflammatory bowel disease (IBD) pathogenesis. The subpopulation of IBD patients as well as clinical strategies if targeting gut fungi in IBD remain underexplored. Here, we conducted a longitudinal, multi-omics analysis in a prospective, head-to-head comparative study involving 40 IBD patients presenting with oral thrush and mild to moderate IBD disease activity. Notably, *C. albicans* strains detected in these patients were shared between the oral cavity and the gut, enabling a therapeutic strategy that distinguishes between oral and orogastrointestinal fungal targeting. Based on this, patients were assigned to standard antifungal treatments with site-specific activity: oral nystatin targeting *Candida* in the oral cavity (ORNT), or systemic fluconazole targeting both oral and gastrointestinal *Candida* (GIFT). By integrating fungal culturing, ITS and metagenomic sequencing with untargeted metabolomics we found that GIFT, but not ORNT, effectively targets gut *C. albicans* and other pathogenic fungi in IBD patients. Cross-kingdom analyses revealed that in the GIFT group reductions in fungal diversity correlated with increased bacterial diversity and the emergence of beneficial microbial metabolites including butyrate and lithocholic acid. Metabolomic modules negatively associated with fungal burden were enriched in bile acid and fatty acid biosynthesis pathways, while modules linked to persistence of *C. albicans* were enriched in pro-inflammatory mediators.  These findings highlight an oral-gut axis of fungal dissemination and establish a mechanistic rationale for testing antifungal therapy as personalized co-treatment strategies for IBD patients suffering pathogenic gut fungal expansion.

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

- **Description**: This section includes all the raw FASTQ files from our study. These files are crucial for in-depth data analysis and understanding the sequencing results from multiple-Omics.

- **Download**: You can access and download these files from the [GEO database](https://chat.openai.com/c/link-to-download).

Below is a detailed annotation of the file structure and contents:

```shell
[4.0K]  .
```

# **3. Processed Data Download**

## 3.1. XXX and XXX  Output

- **Description**: This section includes the output files fromXXXX, essential for the initial data processing and analysis.
- **Download**: These files are available for access and download from the [GEO database](https://chat.openai.com/c/link-to-download).

Below is a detailed annotation of the file structure and contents:

```shell
code[4.0K]  .
```

- Contents
  - Each `XXX.tar.gz` file includes the XX file from XXX count model.
  
    

## 3.2. R Data Files Generated in This Study

- **Description**: All R data files (.rds) related to XXX data are available. These files encompass a comprehensive range of analyses and findings from our study.
- **Download**: You can download these files from [Zenodo Zenodo XXXX](https://zenodo.org/records/XXX).

Here's an annotation for each file to give you :

~~~shell
tree -lh
[4.0K]  .
~~~

Each file seems to contain specific data subsets or analysis results, crucial for a comprehensive understanding of ACC and normal adrenal tissues at the single-cell level.

# **Citation**

Our paper has been published in [*XXX Journal*](https://chat.openai.com/c/xxxx). For further reference and details, you can access the publication at the provided link.

The raw data supporting the findings of this study can be downloaded from the following repositories:

- **GEO Database**: Access our dataset by visiting [GSEXXX](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSEXXX). This link will take you directly to the dataset's page.
- **Zenodo**: Additional data files are available on Zenodo. Download them at [Zenodo XXX](https://zenodo.org/records/XXX).