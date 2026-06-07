#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

imputed_dir <- Sys.getenv(
  "IMPUTED_DIR",
  unset = "dep_shctrl_adjusted_no_rdes_strong_kd_plus_mhhes1_cellline_weighted_78_decontam_results"
)
nonimputed_dir <- Sys.getenv(
  "NONIMPUTED_DIR",
  unset = "dep_shctrl_adjusted_no_rdes_strong_kd_plus_mhhes1_cellline_weighted_78_nonimputed_decontam_results"
)
output_dir <- Sys.getenv(
  "OUTPUT_DIR",
  unset = file.path(imputed_dir, "tables", "publication")
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

read_hits <- function(base_dir, prefix) {
  readr::read_csv(
    file.path(base_dir, "tables", "pooled", "stringent_shctrl_adjusted_responsive_concordant_proteins.csv"),
    show_col_types = FALSE
  ) %>%
    select(
      feature_id,
      gene_symbol,
      response_direction,
      Protein.Group,
      Protein.Names,
      Genes,
      First.Protein.Description,
      logFC,
      P.Value,
      adj.P.Val,
      imputed,
      num_NAs,
      n_clone_adj_sig_down,
      n_clone_adj_sig_up,
      min_clone_logFC,
      max_clone_logFC,
      min_loo_logFC,
      max_loo_logFC,
      max_loo_adj_p
    ) %>%
    rename_with(
      ~ paste(prefix, .x, sep = "_"),
      c(
        logFC,
        P.Value,
        adj.P.Val,
        imputed,
        num_NAs,
        n_clone_adj_sig_down,
        n_clone_adj_sig_up,
        min_clone_logFC,
        max_clone_logFC,
        min_loo_logFC,
        max_loo_logFC,
        max_loo_adj_p
      )
    )
}

imputed_hits <- read_hits(imputed_dir, "imputed")
nonimputed_hits <- read_hits(nonimputed_dir, "nonimputed")

publication_hits <- full_join(
  imputed_hits,
  nonimputed_hits,
  by = c(
    "feature_id",
    "gene_symbol",
    "response_direction",
    "Protein.Group",
    "Protein.Names",
    "Genes",
    "First.Protein.Description"
  )
) %>%
  mutate(
    detected_in_imputed = !is.na(imputed_logFC),
    detected_in_nonimputed = !is.na(nonimputed_logFC),
    confidence_class = case_when(
      detected_in_imputed & detected_in_nonimputed ~ "high_confidence_overlap",
      detected_in_imputed ~ "imputed_only",
      detected_in_nonimputed ~ "nonimputed_only",
      TRUE ~ "unclassified"
    ),
    primary_logFC = coalesce(imputed_logFC, nonimputed_logFC),
    primary_adj.P.Val = coalesce(imputed_adj.P.Val, nonimputed_adj.P.Val)
  ) %>%
  arrange(
    factor(confidence_class, levels = c("high_confidence_overlap", "imputed_only", "nonimputed_only")),
    factor(response_direction, levels = c("down", "up")),
    primary_adj.P.Val,
    desc(abs(primary_logFC))
  )

readr::write_csv(
  publication_hits,
  file.path(output_dir, "sf3b4_shctrl_adjusted_publication_responsive_hits.csv")
)

readr::write_csv(
  publication_hits %>% filter(confidence_class == "high_confidence_overlap"),
  file.path(output_dir, "sf3b4_shctrl_adjusted_high_confidence_overlap_responsive_hits.csv")
)

summary_tbl <- publication_hits %>%
  count(confidence_class, response_direction, name = "n_hits") %>%
  arrange(confidence_class, response_direction)

readr::write_csv(
  summary_tbl,
  file.path(output_dir, "sf3b4_shctrl_adjusted_publication_responsive_hit_summary.csv")
)

message("Wrote publication responsive hit tables to: ", normalizePath(output_dir))
