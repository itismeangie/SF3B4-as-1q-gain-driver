#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)

local_r_lib <- Sys.getenv("LOCAL_R_LIB", unset = "")
if (nzchar(local_r_lib) && dir.exists(local_r_lib)) {
  .libPaths(c(local_r_lib, .libPaths()))
}

suppressPackageStartupMessages({
  library(NetBID2)
  library(edgeR)
  library(tidyverse)
  library(Biobase)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(key, default = NULL) {
  hit <- grep(paste0("^--", key, "="), args, value = TRUE)
  if (length(hit) == 0) return(default)
  sub(paste0("^--", key, "="), "", hit[1])
}

tcga_base_dir <- normalizePath(get_arg("tcga_base_dir", "."), mustWork = TRUE)
tpm_dir <- file.path(tcga_base_dir, "tpm_matrices")
copy_number_path <- file.path(tcga_base_dir, "PANCAN_ArmCallsAndAneuploidyScore_092817.txt")
project_main_dir <- file.path(tcga_base_dir, get_arg("project_main_dir", "netbid2_runs_sensitivity_forced_1q_genes"))
dir.create(project_main_dir, recursive = TRUE, showWarnings = FALSE)

selected_projects_arg <- get_arg("projects", "TCGA-UCEC,TCGA-READ,TCGA-LGG,TCGA-THCA,TCGA-KIRC")
selected_projects <- trimws(strsplit(selected_projects_arg, ",", fixed = TRUE)[[1]])
selected_projects <- selected_projects[nzchar(selected_projects)]
if (length(selected_projects) == 0L) stop("No projects provided.")

forced_gene_map <- data.table(
  gene_symbol = c("SF3B4", "TPR", "TFB2M"),
  ensembl_gene_id = c("ENSG00000143368", "ENSG00000047410", "ENSG00000162851")
)
current_date <- format(Sys.time(), "%Y%m%d")

classify_copy_number <- function(arm_vector) {
  if (is.null(arm_vector["1q"]) || is.na(arm_vector["1q"])) return("other")
  if (any(is.na(arm_vector))) return("other")
  if (all(arm_vector == 0)) return("diploid")
  other_arms <- arm_vector[names(arm_vector) != "1q"]
  if (arm_vector["1q"] > 0 && all(other_arms == 0)) return("1q_gain_only")
  "other"
}

load_copy_number_annotations <- function(path) {
  raw <- readr::read_tsv(path, show_col_types = FALSE, progress = FALSE)
  arm_cols <- setdiff(colnames(raw), c("Sample", "Type", "Aneuploidy Score"))
  raw %>%
    mutate(across(all_of(arm_cols), as.numeric)) %>%
    rowwise() %>%
    mutate(copy_number_category = classify_copy_number(setNames(c_across(all_of(arm_cols)), arm_cols))) %>%
    ungroup() %>%
    mutate(
      cancer_project = paste0("TCGA-", Type),
      copy_number_category = factor(copy_number_category, levels = c("1q_gain_only", "diploid", "other"))
    ) %>%
    select(sample_barcode = Sample, cancer_project, copy_number_category)
}

construct_expression_set <- function(tpm_path, project_id, copy_number_annotations, forced_genes) {
  expr_tbl <- readr::read_tsv(tpm_path, show_col_types = FALSE, progress = FALSE)
  stopifnot("gene_id" %in% colnames(expr_tbl))

  gene_ids_with_version <- expr_tbl$gene_id
  gene_ids <- stringr::str_remove(gene_ids_with_version, "\\.\\d+$")
  if (anyDuplicated(gene_ids) > 0) {
    warning("Duplicated Ensembl IDs detected after removing version numbers. Keeping first occurrence.")
  }

  expr_mat <- as.matrix(expr_tbl[, -1, drop = FALSE])
  rownames(expr_mat) <- gene_ids
  colnames(expr_mat) <- colnames(expr_tbl)[-1]
  expr_mat[expr_mat == 0] <- 0.1

  log2_mat <- log2(expr_mat)
  dge <- DGEList(counts = expr_mat)
  keep <- filterByExpr(dge)
  keep <- keep | (rownames(expr_mat) %in% forced_genes)
  filtered_mat <- log2_mat[keep, , drop = FALSE]

  feature_df <- tibble(
    ensembl_gene_id = gene_ids,
    ensembl_gene_id_version = gene_ids_with_version
  ) %>%
    distinct(ensembl_gene_id, .keep_all = TRUE) %>%
    filter(ensembl_gene_id %in% rownames(filtered_mat)) %>%
    arrange(match(ensembl_gene_id, rownames(filtered_mat)))
  rownames(feature_df) <- feature_df$ensembl_gene_id

  sample_annotations <- tibble(aliquot_barcode = colnames(filtered_mat)) %>%
    mutate(
      sample_barcode = substr(aliquot_barcode, 1, 15),
      patient_barcode = substr(aliquot_barcode, 1, 12),
      cancer_project = project_id
    ) %>%
    left_join(copy_number_annotations, by = c("sample_barcode", "cancer_project")) %>%
    mutate(
      copy_number_category = replace_na(as.character(copy_number_category), "other"),
      copy_number_category = factor(copy_number_category, levels = c("1q_gain_only", "diploid", "other"))
    )
  rownames(sample_annotations) <- sample_annotations$aliquot_barcode

  ExpressionSet(
    assayData = filtered_mat,
    phenoData = AnnotatedDataFrame(sample_annotations),
    featureData = AnnotatedDataFrame(feature_df)
  )
}

run_netbid_for_project <- function(project_id, tpm_path, copy_number_annotations, forced_gene_map) {
  message("Preparing sensitivity NetBID2 inputs for ", project_id)
  eset <- construct_expression_set(
    tpm_path = tpm_path,
    project_id = project_id,
    copy_number_annotations = copy_number_annotations,
    forced_genes = forced_gene_map$ensembl_gene_id
  )

  project_name <- sprintf("%s_forced1q_%s", project_id, current_date)
  project_dir <- file.path(project_main_dir, project_name)
  network.par <- NetBID.network.dir.create(
    project_main_dir = project_main_dir,
    project_name = project_name
  )

  network.par$net.eset <- eset
  NetBID.saveRData(network.par = network.par, step = "exp-load")

  use_gene_type <- "ensembl_gene_id"
  use_genes <- fData(network.par$net.eset)$ensembl_gene_id
  use_list <- get.TF_SIG.list(use_genes, use_gene_type = use_gene_type)
  forced_present <- forced_gene_map$ensembl_gene_id[forced_gene_map$ensembl_gene_id %in% use_genes]
  use_list$sig <- sort(unique(c(use_list$sig, forced_present)))

  phe <- pData(network.par$net.eset)
  use.samples <- rownames(phe)
  prj.name <- network.par$project.name

  SJAracne.prepare(
    eset = network.par$net.eset,
    use.samples = use.samples,
    TF_list = use_list$tf,
    SIG_list = use_list$sig,
    IQR.thre = 0.5,
    IQR.loose_thre = 0.1,
    SJAR.project_name = prj.name,
    SJAR.main_dir = network.par$out.dir.SJAR
  )

  fwrite(
    data.table(
      project = project_id,
      project_name = project_name,
      forced_gene_symbol = forced_gene_map$gene_symbol,
      forced_ensembl_gene_id = forced_gene_map$ensembl_gene_id,
      present_in_eset = forced_gene_map$ensembl_gene_id %in% use_genes,
      forced_into_sig = forced_gene_map$ensembl_gene_id %in% use_list$sig
    ),
    file.path(project_dir, "forced_gene_manifest.tsv"),
    sep = "\t"
  )
}

stopifnot(dir.exists(tpm_dir), file.exists(copy_number_path))
invisible(db.preload(use_level = "gene", use_spe = "human", update = FALSE))
copy_number_annotations <- load_copy_number_annotations(copy_number_path)

available_projects <- tibble(
  tpm_path = list.files(tpm_dir, pattern = "_tpm_matrix\\.tsv$", full.names = TRUE)
) %>%
  mutate(project_id = sub("_tpm_matrix\\.tsv$", "", basename(tpm_path))) %>%
  arrange(project_id)

unknown_projects <- setdiff(selected_projects, available_projects$project_id)
if (length(unknown_projects) > 0L) {
  stop("No TPM matrices found for: ", paste(unknown_projects, collapse = ", "))
}

project_plan <- available_projects %>% filter(project_id %in% selected_projects)
purrr::walk2(
  project_plan$project_id,
  project_plan$tpm_path,
  ~ run_netbid_for_project(.x, .y, copy_number_annotations, forced_gene_map)
)
