#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(broom)
  library(dplyr)
  library(survival)
  library(tidyr)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) > 0) {
  normalizePath(sub("^--file=", "", file_arg[1]))
} else {
  normalizePath(getwd())
}

analysis_dir <- dirname(script_path)
combined_df_path <- Sys.getenv("COMBINED_COHORT_DF", unset = file.path(analysis_dir, "combined_df.csv"))
out_dir <- Sys.getenv("CNA_SURVIVAL_OUT_DIR", unset = analysis_dir)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(combined_df_path)) {
  stop("Missing combined cohort input at ", combined_df_path, call. = FALSE)
}

combined_df <- read.csv(combined_df_path, stringsAsFactors = FALSE, check.names = FALSE)
required_cols <- c(
  "patient_id", "cohort", "present_1q", "present_8", "present_12",
  "present_14q", "present_16q", "time", "event"
)
missing_cols <- setdiff(required_cols, names(combined_df))
if (length(missing_cols) > 0) {
  stop("combined cohort input is missing columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
}

cox_pool <- coxph(
  Surv(time, event) ~ present_1q + present_8 + present_12 + present_14q + present_16q +
    strata(cohort),
  data = combined_df
)

counts <- combined_df %>%
  summarise(
    n_1q = sum(present_1q),
    n_8 = sum(present_8),
    n_12 = sum(present_12),
    n_14q = sum(present_14q),
    n_16q = sum(present_16q)
  ) %>%
  pivot_longer(everything(), names_to = "key", values_to = "n") %>%
  mutate(
    term = recode(
      key,
      n_8 = "present_8",
      n_1q = "present_1q",
      n_16q = "present_16q",
      n_14q = "present_14q",
      n_12 = "present_12"
    )
  )

model_summary <- broom::tidy(cox_pool, exponentiate = TRUE, conf.int = TRUE) %>%
  filter(grepl("^present_", term)) %>%
  mutate(term = sub("TRUE$", "", term)) %>%
  left_join(counts, by = "term") %>%
  mutate(
    CNA = recode(
      term,
      present_8 = "chr8 gain",
      present_1q = "chr1q gain",
      present_16q = "chr16q loss",
      present_14q = "chr14q gain",
      present_12 = "chr12 gain"
    )
  )

out_summary <- file.path(out_dir, "CNA_survival_cox_model_summary.csv")
write.csv(model_summary, out_summary, row.names = FALSE)
message("Wrote computational output: ", out_summary)
