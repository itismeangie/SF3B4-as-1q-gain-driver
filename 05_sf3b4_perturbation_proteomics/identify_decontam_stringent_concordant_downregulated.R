#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(limma)
  library(readr)
  library(tibble)
  library(tidyr)
})

base_dir <- Sys.getenv("BASE_DIR", unset = "dep_no_shctrl_decontam_results")
tables_dir <- file.path(base_dir, "tables")
pooled_tables_dir <- file.path(tables_dir, "pooled")

pooled_fdr_cutoff <- as.numeric(Sys.getenv("POOLED_FDR_CUTOFF", unset = "0.01"))
pooled_logfc_cutoff <- as.numeric(Sys.getenv("POOLED_LOGFC_CUTOFF", unset = "-0.30"))
clone_logfc_cutoff <- as.numeric(Sys.getenv("CLONE_LOGFC_CUTOFF", unset = "-0.20"))
clone_adj_p_cutoff <- as.numeric(Sys.getenv("CLONE_ADJ_P_CUTOFF", unset = "0.05"))

stop_if_missing <- function(path) {
  if (!file.exists(path)) {
    stop("Required input not found: ", path, call. = FALSE)
  }

  TRUE
}

input_paths <- c(
  imputed_matrix = file.path(tables_dir, "imputed_matrix.csv"),
  clone_limma = file.path(tables_dir, "limma_all_contrasts_long.csv"),
  retained_clones = file.path(pooled_tables_dir, "retained_clones.csv"),
  pooled_metadata = file.path(pooled_tables_dir, "pooled_sample_metadata.csv"),
  previous_pooled_limma = file.path(pooled_tables_dir, "limma_pooled_DOX_vs_CTRL.csv")
)

invisible(vapply(input_paths, stop_if_missing, logical(1)))

message("Reading decontaminated DEP outputs from: ", base_dir)
imputed_tbl <- readr::read_csv(input_paths[["imputed_matrix"]], show_col_types = FALSE)
clone_tbl <- readr::read_csv(input_paths[["clone_limma"]], show_col_types = FALSE)
retained_clones <- readr::read_csv(input_paths[["retained_clones"]], show_col_types = FALSE)
pooled_metadata <- readr::read_csv(input_paths[["pooled_metadata"]], show_col_types = FALSE)
previous_pooled_tbl <- readr::read_csv(input_paths[["previous_pooled_limma"]], show_col_types = FALSE)

if (anyDuplicated(imputed_tbl$feature_id) > 0) {
  stop("feature_id must be unique in imputed_matrix.csv", call. = FALSE)
}

if (anyDuplicated(previous_pooled_tbl$feature_id) > 0) {
  stop("feature_id must be unique in limma_pooled_DOX_vs_CTRL.csv", call. = FALSE)
}

retained_clones <- retained_clones %>%
  mutate(clone_id = paste(cell_line, clone, sep = "__"))

pooled_metadata <- pooled_metadata %>%
  semi_join(retained_clones, by = c("cell_line", "clone")) %>%
  mutate(
    treatment = factor(treatment, levels = c("CTRL", "DOX")),
    clone_id = paste(cell_line, clone, sep = "__"),
    clone_id = factor(clone_id)
  )

expected_clone_contrasts <- nrow(retained_clones)
min_clone_adj_support <- floor(expected_clone_contrasts / 2) + 1L

if (expected_clone_contrasts < 2) {
  stop("At least two retained clones are required for concordance analysis.", call. = FALSE)
}

if (nrow(pooled_metadata) == 0) {
  stop("No pooled samples remained after joining retained clones.", call. = FALSE)
}

if (any(is.na(pooled_metadata$ID))) {
  stop("pooled_sample_metadata.csv contains missing DEP sample IDs.", call. = FALSE)
}

missing_sample_ids <- setdiff(pooled_metadata$ID, colnames(imputed_tbl))
if (length(missing_sample_ids) > 0) {
  stop(
    "Imputed matrix is missing pooled sample IDs: ",
    paste(missing_sample_ids, collapse = ", "),
    call. = FALSE
  )
}

feature_annotation <- imputed_tbl %>%
  select(any_of(c(
    "feature_id",
    "Protein.Group",
    "Protein.Names",
    "Genes",
    "First.Protein.Description",
    "name",
    "ID",
    "imputed",
    "num_NAs"
  )))

expression_matrix <- imputed_tbl %>%
  select(all_of(pooled_metadata$ID)) %>%
  as.matrix()
storage.mode(expression_matrix) <- "numeric"
rownames(expression_matrix) <- imputed_tbl$feature_id

fit_clone_adjusted_pooled_limma <- function(sample_metadata, expression_matrix) {
  sample_ids <- sample_metadata$ID
  missing_sample_ids <- setdiff(sample_ids, colnames(expression_matrix))
  if (length(missing_sample_ids) > 0) {
    stop(
      "Expression matrix is missing sample IDs: ",
      paste(missing_sample_ids, collapse = ", "),
      call. = FALSE
    )
  }

  model_metadata <- sample_metadata %>%
    mutate(
      treatment = factor(as.character(treatment), levels = c("CTRL", "DOX")),
      clone_id = factor(as.character(clone_id))
    )

  if (nlevels(model_metadata$treatment) != 2 || any(table(model_metadata$treatment) == 0)) {
    stop("Both CTRL and DOX samples are required for pooled limma fitting.", call. = FALSE)
  }

  if (nlevels(model_metadata$clone_id) < 2) {
    stop("At least two clone_id levels are required for pooled limma fitting.", call. = FALSE)
  }

  design <- model.matrix(~ treatment + clone_id, data = model_metadata)
  if (!"treatmentDOX" %in% colnames(design)) {
    stop("Could not find treatmentDOX coefficient in pooled design.", call. = FALSE)
  }

  fit <- limma::lmFit(expression_matrix[, sample_ids, drop = FALSE], design)
  fit <- limma::eBayes(fit)

  limma::topTable(
    fit,
    coef = "treatmentDOX",
    number = Inf,
    sort.by = "P"
  ) %>%
    rownames_to_column("feature_id")
}

message("Fitting clone-adjusted pooled limma model: ~ treatment + clone_id")
clone_adjusted_pooled_tbl <- fit_clone_adjusted_pooled_limma(
  pooled_metadata,
  expression_matrix
) %>%
  left_join(feature_annotation, by = "feature_id") %>%
  mutate(
    gene_symbol = dplyr::coalesce(name, Genes),
    contrast = "pooled_DOX_vs_CTRL",
    model = "treatment + clone_id"
  )

message("Fitting leave-one-clone-out pooled models")
loo_tbl <- lapply(seq_len(nrow(retained_clones)), function(idx) {
  held_out <- retained_clones[idx, ]
  loo_metadata <- pooled_metadata %>%
    filter(clone_id != held_out$clone_id)

  fit_clone_adjusted_pooled_limma(loo_metadata, expression_matrix) %>%
    transmute(
      feature_id,
      held_out_cell_line = held_out$cell_line,
      held_out_clone = held_out$clone,
      held_out_clone_id = held_out$clone_id,
      loo_logFC = logFC,
      loo_P.Value = P.Value,
      loo_adj.P.Val = adj.P.Val
    )
}) %>%
  bind_rows()

validated_clone_tbl <- clone_tbl %>%
  semi_join(retained_clones, by = c("cell_line", "clone")) %>%
  transmute(
    feature_id,
    gene_symbol,
    contrast,
    cell_line,
    clone,
    clone_id = paste(cell_line, clone, sep = "__"),
    clone_label = paste(cell_line, clone, sep = " | "),
    clone_logFC = logFC,
    clone_p_value = P.Value,
    clone_adj_p = adj.P.Val,
    clone_imputed = imputed,
    clone_num_NAs = num_NAs,
    clone_down = clone_logFC < 0,
    clone_down_effect_floor = clone_logFC <= clone_logfc_cutoff,
    clone_down_adj_sig = clone_logFC < 0 & clone_adj_p < clone_adj_p_cutoff
  )

clone_validation_summary <- validated_clone_tbl %>%
  mutate(
    down_evidence_z = dplyr::if_else(
      clone_logFC < 0,
      stats::qnorm(1 - pmax(clone_p_value, .Machine$double.xmin) / 2),
      -stats::qnorm(1 - pmax(clone_p_value, .Machine$double.xmin) / 2)
    )
  ) %>%
  group_by(feature_id, gene_symbol) %>%
  summarise(
    n_clone_contrasts = n(),
    expected_clone_contrasts = expected_clone_contrasts,
    has_all_retained_clones = n_clone_contrasts == expected_clone_contrasts,
    n_clone_down = sum(clone_down, na.rm = TRUE),
    n_clone_down_effect_floor = sum(clone_down_effect_floor, na.rm = TRUE),
    n_clone_adj_sig_down = sum(clone_down_adj_sig, na.rm = TRUE),
    all_clones_down = has_all_retained_clones & n_clone_down == expected_clone_contrasts,
    all_clones_down_effect_floor =
      has_all_retained_clones & n_clone_down_effect_floor == expected_clone_contrasts,
    clone_adj_sig_majority = n_clone_adj_sig_down >= min_clone_adj_support,
    all_clone_adj_sig_down = n_clone_adj_sig_down == expected_clone_contrasts,
    min_clone_logFC = min(clone_logFC, na.rm = TRUE),
    max_clone_logFC = max(clone_logFC, na.rm = TRUE),
    directional_stouffer_down_z = sum(down_evidence_z, na.rm = TRUE) / sqrt(n()),
    directional_stouffer_down_p =
      stats::pnorm(directional_stouffer_down_z, lower.tail = FALSE),
    .groups = "drop"
  )

loo_summary <- loo_tbl %>%
  group_by(feature_id) %>%
  summarise(
    n_loo_models = n(),
    expected_loo_models = expected_clone_contrasts,
    all_loo_logFC_negative = n_loo_models == expected_loo_models &
      all(loo_logFC < 0, na.rm = TRUE),
    n_loo_adj_sig = sum(loo_logFC < 0 & loo_adj.P.Val < clone_adj_p_cutoff, na.rm = TRUE),
    min_loo_logFC = min(loo_logFC, na.rm = TRUE),
    max_loo_logFC = max(loo_logFC, na.rm = TRUE),
    max_loo_adj_p = max(loo_adj.P.Val, na.rm = TRUE),
    .groups = "drop"
  )

previous_pooled_summary <- previous_pooled_tbl %>%
  transmute(
    feature_id,
    previous_cell_line_model_logFC = logFC,
    previous_cell_line_model_adj_p = adj.P.Val
  )

stringent_tbl <- clone_adjusted_pooled_tbl %>%
  left_join(previous_pooled_summary, by = "feature_id") %>%
  left_join(
    clone_validation_summary %>% select(-gene_symbol),
    by = "feature_id"
  ) %>%
  left_join(loo_summary, by = "feature_id") %>%
  mutate(
    clone_adjusted_pooled_fdr_lt_cutoff = adj.P.Val < pooled_fdr_cutoff,
    clone_adjusted_pooled_down_effect_floor = logFC <= pooled_logfc_cutoff,
    previous_pooled_fdr_lt_cutoff =
      previous_cell_line_model_adj_p < pooled_fdr_cutoff,
    previous_pooled_down_effect_floor =
      previous_cell_line_model_logFC <= pooled_logfc_cutoff,
    both_pooled_models_supported =
      clone_adjusted_pooled_fdr_lt_cutoff &
      clone_adjusted_pooled_down_effect_floor &
      previous_pooled_fdr_lt_cutoff &
      previous_pooled_down_effect_floor,
    stringent_down_concordant =
      both_pooled_models_supported &
      all_clones_down &
      all_clones_down_effect_floor &
      clone_adj_sig_majority &
      all_loo_logFC_negative
  ) %>%
  arrange(adj.P.Val, desc(abs(logFC)))

stringent_hits <- stringent_tbl %>%
  filter(stringent_down_concordant) %>%
  arrange(adj.P.Val, desc(abs(logFC)))

stringent_clone_validation_long <- validated_clone_tbl %>%
  semi_join(stringent_hits, by = c("feature_id", "gene_symbol")) %>%
  arrange(match(feature_id, stringent_hits$feature_id), cell_line, clone)

stringent_loo_long <- loo_tbl %>%
  semi_join(stringent_hits, by = "feature_id") %>%
  left_join(
    stringent_hits %>% select(feature_id, gene_symbol),
    by = "feature_id"
  ) %>%
  arrange(match(feature_id, stringent_hits$feature_id), held_out_cell_line, held_out_clone)

stringent_heatmap_matrix <- stringent_clone_validation_long %>%
  select(gene_symbol, clone_label, clone_logFC, clone_adj_p) %>%
  pivot_wider(
    names_from = clone_label,
    values_from = c(clone_logFC, clone_adj_p)
  )

filter_summary <- stringent_tbl %>%
  summarise(
    retained_clones = max(.data$expected_clone_contrasts, na.rm = TRUE),
    step_1_both_pooled_models_adj_p_lt_cutoff_and_logFC_lt_0 = sum(
      clone_adjusted_pooled_fdr_lt_cutoff &
        previous_pooled_fdr_lt_cutoff &
        logFC < 0 &
        previous_cell_line_model_logFC < 0,
      na.rm = TRUE
    ),
    step_2_plus_both_pooled_models_logFC_effect_floor = sum(
      clone_adjusted_pooled_fdr_lt_cutoff &
        previous_pooled_fdr_lt_cutoff &
        clone_adjusted_pooled_down_effect_floor &
        previous_pooled_down_effect_floor,
      na.rm = TRUE
    ),
    step_3_plus_negative_logFC_in_all_retained_clones = sum(
      both_pooled_models_supported &
        all_clones_down,
      na.rm = TRUE
    ),
    step_4_plus_clone_logFC_effect_floor_in_all_retained_clones = sum(
      both_pooled_models_supported &
        all_clones_down &
        all_clones_down_effect_floor,
      na.rm = TRUE
    ),
    step_5_plus_clone_adj_p_support_in_majority = sum(
      both_pooled_models_supported &
        all_clones_down &
        all_clones_down_effect_floor &
        clone_adj_sig_majority,
      na.rm = TRUE
    ),
    step_6_plus_negative_in_all_leave_one_clone_out_models = sum(
      stringent_down_concordant,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = everything(),
    names_to = "filter_step",
    values_to = "remaining_proteins"
  )

rules_tbl <- tibble(
  pooled_model_primary = "treatment + cell_line",
  pooled_model_sensitivity = "treatment + clone_id",
  pooled_fdr_cutoff = pooled_fdr_cutoff,
  pooled_logfc_cutoff = pooled_logfc_cutoff,
  clone_logfc_cutoff = clone_logfc_cutoff,
  clone_adj_p_cutoff = clone_adj_p_cutoff,
  expected_clone_contrasts = expected_clone_contrasts,
  min_clone_adj_support = min_clone_adj_support,
  leave_one_clone_out_requirement = "pooled logFC < 0 after holding out each retained clone",
  definition = paste(
    "Stringent down-concordant proteins require pooled adj.P.Val <",
    pooled_fdr_cutoff,
    "and pooled logFC <=",
    pooled_logfc_cutoff,
    "in both pooled models (treatment + cell_line and treatment + clone_id),",
    "negative clone-level logFC in all retained clones, clone-level logFC <=",
    clone_logfc_cutoff,
    "in all retained clones, clone-level adj.P.Val <",
    clone_adj_p_cutoff,
    "in a strict majority of retained clones, and negative pooled logFC in every leave-one-clone-out refit."
  )
)

readr::write_csv(
  clone_adjusted_pooled_tbl,
  file.path(pooled_tables_dir, "limma_pooled_DOX_vs_CTRL_clone_adjusted.csv")
)

readr::write_csv(
  loo_tbl,
  file.path(pooled_tables_dir, "leave_one_clone_out_pooled_limma_long.csv")
)

readr::write_csv(
  stringent_tbl,
  file.path(pooled_tables_dir, "stringent_concordance_all_proteins.csv")
)

readr::write_csv(
  stringent_hits,
  file.path(pooled_tables_dir, "stringent_downregulated_clone_concordant_proteins.csv")
)

readr::write_csv(
  stringent_clone_validation_long,
  file.path(pooled_tables_dir, "stringent_downregulated_clone_validation_long.csv")
)

readr::write_csv(
  stringent_loo_long,
  file.path(pooled_tables_dir, "stringent_downregulated_leave_one_clone_out_long.csv")
)

readr::write_csv(
  stringent_heatmap_matrix,
  file.path(pooled_tables_dir, "stringent_downregulated_clone_concordant_heatmap_matrix.csv")
)

readr::write_csv(
  filter_summary,
  file.path(pooled_tables_dir, "stringent_concordance_filter_summary.csv")
)

readr::write_csv(
  rules_tbl,
  file.path(pooled_tables_dir, "stringent_concordance_rules.csv")
)

violations <- stringent_hits %>%
  filter(
    !both_pooled_models_supported |
      !all_clones_down |
      !all_clones_down_effect_floor |
      !clone_adj_sig_majority |
      !all_loo_logFC_negative
  )

if (nrow(violations) > 0) {
  stop("Internal validation failed: stringent hits include rule violations.", call. = FALSE)
}

message("Wrote stringent concordance tables to: ", pooled_tables_dir)
message("Stringent down-concordant proteins: ", nrow(stringent_hits))
if (nrow(stringent_hits) > 0) {
  message(paste(stringent_hits$gene_symbol, collapse = ", "))
}
