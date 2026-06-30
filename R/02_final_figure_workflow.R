# Minimal final-figure workflow skeleton.
# This script assumes controlled, de-identified derived tables are available.

source("R/00_cohort_definition.R")
source("R/01_stats_utils.R")

metadata_path <- "data/derived/sample_metadata_public_columns.csv"
feature_path <- "data/derived/final_feature_table.csv"
output_dir <- "outputs"

if (!file.exists(metadata_path) || !file.exists(feature_path)) {
  stop(
    "Required derived input tables are not present in the public package. ",
    "Prepare de-identified derived tables inside the controlled workspace first."
  )
}

metadata <- read.csv(metadata_path, check.names = FALSE)
features <- read.csv(feature_path, check.names = FALSE)

cohort <- define_main_cohort(metadata)
counts <- cohort_counts(cohort)
write_cohort_counts(counts, file.path(output_dir, "cohort_counts"))

analysis_data <- merge(
  cohort,
  features,
  by = "sample_id",
  all.x = TRUE,
  all.y = FALSE
)

# Example primary association family. Replace these with the manuscript's final
# prespecified feature names in the controlled workspace.
microbiome_features <- intersect(
  c("Candida_albicans", "Fungal_shannon", "Bacterial_shannon"),
  colnames(analysis_data)
)
clinical_features <- intersect(
  c("DAI", "UC_score", "CD_score"),
  colnames(analysis_data)
)

if (length(microbiome_features) > 0 && length(clinical_features) > 0) {
  correlation_results <- spearman_grid(
    analysis_data,
    x_vars = microbiome_features,
    y_vars = clinical_features
  )
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(correlation_results, file.path(output_dir, "primary_spearman_results.csv"), row.names = FALSE)
}

# Plotting code should read final figure-specific summary tables and write
# figures directly from script. No manual figure edits are part of this workflow.
