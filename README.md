# MMC

This page recorded the codes and data used and mentioned in [*xxx*](XXX). And you could downloaded this paper by clicking [here](pdf/XXX)

Gut fungal dysbiosis has been implicated in inflammatory bowel disease (IBD), yet strategies for targeting the gut mycobiota in IBD remain unexplored. Here we leveraged the observation that *Candida albicans* strains are shared between the oral cavity and gut in IBD patients with mild oral thrush to design a prospective observational study where we compared oral antifungal therapy (swish‑and‑spit nystatin, ORNT) versus oro‑gastrointestinal antifungal therapy (fluconazole, GIFT). In 53 patients with mild‑to‑moderate ulcerative colitis (UC) or Crohn’s disease (CD), fluconazole, but not nystatin, effectively reduced intestinal *Candida* burden and reshaped gut fungal composition. Intestinal fungal targeting in the GIFT group was associated with increased bacterial diversity, expansion of short‑chain fatty acid-producing taxa, restoration of anti‑inflammatory microbial metabolites, and durable shifts in cross‑kingdom microbial networks. This microbiome and metabolomic changes in the GIFT group were associated with improved disease activity indices and a decreased risk of disease progression over the 8‑week follow‑up period. This study demonstrates the feasibility of mycobiome‑based patient stratification and establishes a framework for investigating antifungal co-therapy for IBD patients with fungal manifestations.

# **1. Codes of analyzing and visualization**

**Introduction to Our Script Compilation for Analysis**

In our comprehensive analysis of antifungal treatment in the MMC cohort, we have organized our scripts into five analysis chapters (Chapters 2–6). Each chapter focuses on a specific aspect of the analysis, enabling a detailed examination of different facets of antifungal treatment effects in IBD. Below is a guide to the content of each chapter and the script file that contains it:

**[Chapter 1 — Fungal Abundance Data Analysis](Fungi.md)** (`Fungi.md`)

- **Content:** Integrates **fungal abundance data (ITS sequencing) with clinical metadata** to assess **treatment-induced shifts in fungal communities**, including alpha/beta diversity, *Candida* burden, taxonomic composition, and visualization of fungal responses across treatment groups and time points.

**[Chapter 2 — Bacterial Analysis: Data Processing, Visualization, and Statistical Testing](Bacterial.md)** (`Bacterial.md`)

- **Content:** Incorporates **ITS sequencing and metagenomic profiling** from the **MMC cohort**, focusing on **treatment responses and gut microbial shifts**. It provides additional evidence on **how different antifungal treatments influence the gut bacterial community**, covering diversity metrics, differential abundance, and SCFA-producing taxa.

**[Chapter 3 — Metabolomic Data Analysis](Metabolomics.md)** (`Metabolomics.md`)

- **Content:** Investigates **metabolomic shifts under antifungal treatments** (Nystatin and Fluconazole), including data normalization, differential metabolite testing, Mfuzz temporal clustering, and pathway enrichment.

**[Chapter 4 — Clinical Data Analysis](Clincal.md)** (`Clincal.md`)

- **Content:** Presents **clinical outcome analyses**, including **disease scores (DAI, UC, CD), patient metadata, medications, disease progression, and their association with antifungal treatments**.

**[Chapter 5 — Cross-kingdom and Clinical Correlation Analysis](Cor_Fungi_Metabolites.md)**

- **Content:** Integrates the multi-omics layers above to characterize **cross-kingdom microbial networks** and their links to clinical outcomes. This chapter is split across three scripts:
  - **[5.1 Fungal features vs metabolite modules](Cor_Fungi_Metabolites.md)** (`Cor_Fungi_Metabolites.md`) — correlates fungal diversity / *Candida* burden with metabolite modules.
  - **[5.2 Bacterial taxa vs metabolite modules](Cor_Bac_Metabolites.md)** (`Cor_Bac_Metabolites.md`) — correlates Fluconazole-responsive bacterial taxa with metabolite modules.
  - **[5.3 Clinical indices vs microbiome/metabolome features](Cor_Clincal_micro.md)** (`Cor_Clincal_micro.md`) — relates SCFA-producing/probiotic taxa and metabolite features to clinical disease indices.

Each chapter contains detailed scripts, methodologies, and analyses relevant to the specific aspect of antifungal treatment it addresses. This structured approach allows researchers to navigate our comprehensive analysis with ease, enhancing their understanding of antifungal treatment effects in IBD.

# **2. Raw data download**

- **Description**: This section includes all the raw FASTQ files from our study. These files are crucial for in-depth data analysis and understanding the sequencing results from multiple-Omics.

- **Download**: You can access and download these files from the [PRJNA1449485](https://dataview.ncbi.nlm.nih.gov/object/PRJNA1449485).

# **Citation**

Our paper has been published in [*Nature Medicine*](https://c). For further reference and details, you can access the publication at the provided link.

The raw data supporting the findings of this study can be downloaded from the following repositories:

- **Database**: Access our dataset by visiting [PRJNA1449485](https://dataview.ncbi.nlm.nih.gov/object/PRJNA1449485). 