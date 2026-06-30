# Data Notes

This directory intentionally does not contain controlled clinical metadata or derived individual-level multi-omics tables.

Expected inputs for the public scripts are de-identified derived tables prepared inside the controlled analysis workspace:

- sample metadata with one row per sample,
- participant-level metadata with one row per participant,
- feature matrices for fungal, bacterial, metabolomic, and clinical variables,
- final figure-specific summary tables where controlled metadata cannot be redistributed.

Minimum recommended metadata columns:

- `participant_id`: de-identified participant code used only inside the controlled workspace,
- `sample_id`: de-identified sample code,
- `treatment`: treatment group,
- `diagnosis`: diagnosis group,
- `time`: study time point,
- `include_main_analysis`: Boolean or 0/1 flag for the primary analysis cohort,
- `cohort_phase`: optional field for prespecified main vs extended cohort membership.

Do not commit files containing dates of birth, medical record numbers, names, contact information, full visit dates, or unrestricted individual-level clinical records.
