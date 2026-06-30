# Cohort definition utilities for the minimal public workflow.
# These functions use metadata fields only; no patient/sample IDs are hard-coded.

required_columns <- function(data, columns) {
  missing <- setdiff(columns, colnames(data))
  if (length(missing) > 0) {
    stop("Missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

define_main_cohort <- function(metadata,
                               participant_col = "participant_id",
                               sample_col = "sample_id",
                               treatment_col = "treatment",
                               diagnosis_col = "diagnosis",
                               time_col = "time",
                               include_col = "include_main_analysis",
                               allowed_treatments = c("Nystatin", "Fluconazole"),
                               allowed_times = c("W0", "W1", "W2", "W4", "W8")) {
  required_columns(
    metadata,
    c(participant_col, sample_col, treatment_col, diagnosis_col, time_col, include_col)
  )

  cohort <- metadata
  cohort <- cohort[cohort[[include_col]] %in% c(TRUE, 1, "TRUE", "true", "Yes", "yes"), , drop = FALSE]
  cohort <- cohort[cohort[[treatment_col]] %in% allowed_treatments, , drop = FALSE]
  cohort <- cohort[cohort[[time_col]] %in% allowed_times, , drop = FALSE]

  cohort[[treatment_col]] <- factor(cohort[[treatment_col]], levels = allowed_treatments)
  cohort[[time_col]] <- factor(cohort[[time_col]], levels = allowed_times)
  cohort
}

define_extended_cohort <- function(metadata,
                                   cohort_phase_col = "cohort_phase",
                                   extended_level = "extended",
                                   ...) {
  required_columns(metadata, cohort_phase_col)
  main <- define_main_cohort(metadata, ...)
  main[main[[cohort_phase_col]] %in% c("main", extended_level), , drop = FALSE]
}

cohort_counts <- function(metadata,
                          participant_col = "participant_id",
                          sample_col = "sample_id",
                          treatment_col = "treatment",
                          diagnosis_col = "diagnosis",
                          time_col = "time") {
  required_columns(metadata, c(participant_col, sample_col, treatment_col, diagnosis_col, time_col))

  list(
    n_participants = length(unique(metadata[[participant_col]])),
    n_samples = length(unique(metadata[[sample_col]])),
    by_treatment = as.data.frame(table(metadata[[treatment_col]], useNA = "ifany")),
    by_diagnosis = as.data.frame(table(metadata[[diagnosis_col]], useNA = "ifany")),
    by_time = as.data.frame(table(metadata[[time_col]], useNA = "ifany"))
  )
}

write_cohort_counts <- function(counts, outdir) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  write.csv(data.frame(metric = "n_participants", value = counts$n_participants),
            file.path(outdir, "n_participants.csv"), row.names = FALSE)
  write.csv(data.frame(metric = "n_samples", value = counts$n_samples),
            file.path(outdir, "n_samples.csv"), row.names = FALSE)
  write.csv(counts$by_treatment, file.path(outdir, "counts_by_treatment.csv"), row.names = FALSE)
  write.csv(counts$by_diagnosis, file.path(outdir, "counts_by_diagnosis.csv"), row.names = FALSE)
  write.csv(counts$by_time, file.path(outdir, "counts_by_time.csv"), row.names = FALSE)
  invisible(outdir)
}
