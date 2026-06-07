#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)

local_r_lib <- Sys.getenv("LOCAL_R_LIB", unset = "")
if (nzchar(local_r_lib) && dir.exists(local_r_lib)) {
  .libPaths(c(local_r_lib, .libPaths()))
}

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(Biobase)
  library(limma)
  library(NetBID2)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(key, default = NULL) {
  hit <- grep(paste0("^--", key, "="), args, value = TRUE)
  if (length(hit) == 0) return(default)
  sub(paste0("^--", key, "="), "", hit[1])
}

project_arg <- get_arg("projects", "TCGA-UCEC,TCGA-READ,TCGA-LGG,TCGA-THCA,TCGA-KIRC")
project_prefixes <- trimws(strsplit(project_arg, ",", fixed = TRUE)[[1]])
project_prefixes <- project_prefixes[nzchar(project_prefixes)]
if (length(project_prefixes) == 0L) stop("No projects provided via --projects.")

tcga_base_dir <- normalizePath(get_arg("tcga_base_dir", "."), mustWork = TRUE)
netbid_dir <- normalizePath(get_arg("netbid_dir", file.path(tcga_base_dir, "netbid2_runs")), mustWork = TRUE)
summary_file <- normalizePath(
  get_arg(
    "summary_file",
    Sys.getenv("NETBID_SUMMARY_FILE", unset = file.path(tcga_base_dir, "per_cancer_v5_project_summary.tsv"))
  ),
  mustWork = TRUE
)
cohort_file <- normalizePath(
  get_arg(
    "cohort_file",
    Sys.getenv("NETBID_COHORT_FILE", unset = file.path(tcga_base_dir, "panTCGA_1q_v2_patient_covariates_maxpanel.tsv.gz"))
  ),
  mustWork = TRUE
)
clinical_xlsx <- normalizePath(
  get_arg(
    "clinical_xlsx",
    Sys.getenv("TCGA_CLINICAL_XLSX", unset = file.path(tcga_base_dir, "TCGA-CDR-SupplementalTableS1.xlsx"))
  ),
  mustWork = TRUE
)
arm_calls_file <- normalizePath(
  get_arg(
    "arm_calls_file",
    Sys.getenv("TCGA_ARM_CALLS_FILE", unset = file.path(tcga_base_dir, "PANCAN_ArmCallsAndAneuploidyScore_092817.txt"))
  ),
  mustWork = TRUE
)
output_label <- get_arg("output_label", "1q_gain_vs_neutral_best_v5_model")
min_group_size <- as.integer(get_arg("min_group_size", "10"))
min_regulon_size <- as.integer(get_arg("min_regulon_size", "5"))
max_regulon_size <- as.integer(get_arg("max_regulon_size", "1000"))

group1_label <- "1q_gain"
group0_label <- "no_1q_gain"
comp_name <- paste0(group1_label, "_vs_", group0_label)

normalize_tcga_sample <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x[x == "" | x == "NA"] <- NA_character_
  x <- sub("^(TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-[0-9]{2}).*$", "\\1", x)
  x[!grepl("^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-[0-9]{2}$", x)] <- NA_character_
  x
}

normalize_tcga_patient <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x[x == "" | x == "NA"] <- NA_character_
  x <- sub("^(TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}).*$", "\\1", x)
  x[!grepl("^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}$", x)] <- NA_character_
  x
}

safe_numeric <- function(x) {
  x <- as.character(x)
  x[x %in% c("NaN", "nan", "NA", "", "NULL", "[Not Available]", "[Not Applicable]")] <- NA_character_
  suppressWarnings(as.numeric(x))
}

zscore <- function(x) {
  x <- as.numeric(x)
  s <- stats::sd(x, na.rm = TRUE)
  m <- mean(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) {
    out <- rep(0, length(x))
    out[is.na(x)] <- NA_real_
    return(out)
  }
  (x - m) / s
}

clean_histology <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "[Not Available]")] <- NA_character_
  x
}

clean_grade <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x[x %in% c("", "NA", "[NOT AVAILABLE]", "[DISCREPANCY]")] <- NA_character_
  x[x == "HIGH GRADE"] <- "G3"
  x
}

clean_race <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x[x %in% c("", "NA", "[NOT AVAILABLE]", "[NOT APPLICABLE]")] <- NA_character_
  ifelse(
    grepl("^WHITE", x),
    "White",
    ifelse(
      grepl("BLACK|AFRICAN", x),
      "Black",
      ifelse(
        grepl("ASIAN", x),
        "Asian",
        ifelse(is.na(x), NA_character_, "Other")
      )
    )
  )
}

collapse_factor <- function(
    x,
    max_levels = 3L,
    min_count = 15L,
    other_label = "Other",
    missing_label = "Missing"
) {
  x <- trimws(as.character(x))
  x[x == "" | x == "NA"] <- NA_character_

  counts <- sort(table(x), decreasing = TRUE)
  keep <- names(counts)[counts >= min_count]
  keep <- head(keep, max_levels)

  y <- rep(NA_character_, length(x))
  nonmiss <- !is.na(x)

  if (length(keep) == 0L) {
    y[nonmiss] <- other_label
  } else {
    y[nonmiss] <- ifelse(x[nonmiss] %in% keep, x[nonmiss], other_label)
  }

  y[is.na(y)] <- missing_label
  factor(y)
}

sanitize_dt_names <- function(dt) {
  setnames(dt, names(dt), trimws(names(dt)))
  dt
}

prepare_extra_clinical <- function(path) {
  raw <- suppressWarnings(as.data.table(read_excel(path, .name_repair = "minimal")))
  sanitize_dt_names(raw)
  out <- raw[, .(
    patient_id = normalize_tcga_patient(bcr_patient_barcode),
    cancer_type = as.character(type),
    race_raw = as.character(race),
    diagnosis_year = safe_numeric(initial_pathologic_dx_year)
  )]
  out <- out[!is.na(patient_id) & !is.na(cancer_type)]
  unique(out, by = c("patient_id", "cancer_type"))
}

build_arm_metrics <- function(path, target_arm = "1q", target_state = "gain") {
  arm <- fread(path, sep = "\t", header = TRUE, check.names = FALSE)
  sanitize_dt_names(arm)
  setnames(arm, c("Sample", "Type"), c("sample_id", "cancer_type"), skip_absent = TRUE)
  if (!target_arm %in% names(arm)) stop("Target arm column not found: ", target_arm)
  arm[, sample_id := normalize_tcga_sample(sample_id)]
  arm <- arm[!is.na(sample_id)]
  target_multiplier <- if (target_state == "gain") 1 else -1
  arm_cols <- setdiff(names(arm), c("sample_id", "cancer_type", "Aneuploidy Score"))
  arm_cols <- arm_cols[arm_cols != target_arm]
  for (col in arm_cols) arm[, (col) := safe_numeric(get(col))]
  arm[, target_raw := safe_numeric(get(target_arm))]
  arm[, oneq_num := as.numeric(target_multiplier * target_raw)]
  arm[, oneq_state := fifelse(oneq_num == 1, "gain", fifelse(oneq_num == -1, "loss", "neutral"))]
  arm[, gain_1q := as.integer(oneq_num == 1)]
  arm[, loss_1q := as.integer(oneq_num == -1)]
  arm[, arm_event_count_excl_1q := rowSums(abs(.SD) > 0, na.rm = TRUE), .SDcols = arm_cols]
  arm[, onep_loss := as.integer(`1p` < 0)]
  arm[, nineteenq_loss := as.integer(`19q` < 0)]
  arm[, onep19q_codel := as.integer(onep_loss == 1 & nineteenq_loss == 1)]
  arm[, .(
    sample_id,
    oneq_state,
    oneq_num,
    gain_1q,
    loss_1q,
    arm_event_count_excl_1q,
    onep19q_codel
  )]
}

parse_var_list <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x)] <- ""
  x <- gsub('"', "", x, fixed = TRUE)
  if (!nzchar(x)) return(character(0))
  vals <- trimws(strsplit(x, ",", fixed = TRUE)[[1]])
  vals[nzchar(vals)]
}

is_numeric_var <- function(var_name) {
  var_name %in% c(
    "age_at_dx", "diagnosis_year", "purity", "ImmuneScore",
    "StromaScore", "tmb_log1p", "arm_event_count_excl_1q"
  )
}

is_binary_var <- function(var_name) {
  startsWith(var_name, "drv_") ||
    startsWith(var_name, "feat_") ||
    var_name %in% c("gender_male", "onep19q_codel")
}

find_latest_run <- function(prefix) {
  runs <- list.files(netbid_dir, pattern = paste0("^", prefix, "_"), full.names = TRUE)
  if (length(runs) == 0L) return(NA_character_)
  sort(runs)[length(runs)]
}

find_consensus_network <- function(run_dir, tag) {
  primary <- file.path(run_dir, "SJAR", sprintf("%s%s", basename(run_dir), tag), "consensus_network_ncol_.txt")
  if (file.exists(primary)) return(primary)
  alt <- list.files(run_dir, pattern = "consensus_network_ncol_.txt", recursive = TRUE, full.names = TRUE)
  alt <- alt[grepl(tag, alt, fixed = TRUE)]
  if (length(alt) == 0L) stop("No consensus network found for ", tag, " in ", run_dir)
  alt[[1]]
}

prepare_analysis_objects <- function(run_dir) {
  network_par_path <- file.path(run_dir, "DATA", "network.par.Step.exp-load.RData")
  if (!file.exists(network_par_path)) stop("Missing network.par at ", network_par_path)
  tf_network_file <- find_consensus_network(run_dir, "_TF")
  sig_network_file <- find_consensus_network(run_dir, "_SIG")

  load(network_par_path)
  tf_network <- get.SJAracne.network(network_file = tf_network_file)
  sig_network <- get.SJAracne.network(network_file = sig_network_file)
  merge_network <- merge_TF_SIG.network(TF_network = tf_network, SIG_network = sig_network)

  phe_df <- as.data.frame(Biobase::pData(network.par$net.eset))
  rownames(phe_df) <- sampleNames(network.par$net.eset)
  ac_mat <- cal.Activity(
    target_list = merge_network$target_list,
    cal_mat = Biobase::exprs(network.par$net.eset),
    es.method = "weightedmean"
  )
  ac_eset <- generate.eset(
    exp_mat = ac_mat,
    phenotype_info = phe_df[colnames(ac_mat), , drop = FALSE],
    feature_info = NULL,
    annotation_info = "activity in net-dataset"
  )

  list(
    net_eset = network.par$net.eset,
    tf_network = tf_network,
    sig_network = sig_network,
    merge_network = merge_network,
    merge_ac_eset = ac_eset
  )
}

make_symbol_map <- function(network_dat) {
  map <- setNames(network_dat$target.symbol, network_dat$target)
  if ("source" %in% colnames(network_dat) && "source.symbol" %in% colnames(network_dat)) {
    base_id <- sub("_(TF|SIG)$", "", network_dat$source)
    map <- c(map, setNames(network_dat$source.symbol, base_id))
  }
  map
}

apply_symbol_names <- function(df, map) {
  sym <- map[rownames(df)]
  sym[is.na(sym)] <- rownames(df)[is.na(sym)]
  rownames(df) <- make.unique(sym)
  df
}

signed_z_from_p <- function(logfc, pval) {
  pval <- pmax(as.numeric(pval), .Machine$double.xmin)
  sign(as.numeric(logfc)) * stats::qnorm(pval / 2, lower.tail = FALSE)
}

run_adjusted_limma <- function(mat, pheno, coef_name = "gain_1q", g1_label = group1_label, g0_label = group0_label) {
  use_samples <- rownames(pheno$design)
  mat <- mat[, use_samples, drop = FALSE]
  fit <- limma::lmFit(mat, pheno$design)
  fit <- limma::eBayes(fit, trend = TRUE)
  tt <- limma::topTable(fit, coef = coef_name, number = Inf, sort.by = "none")
  tt <- tt[rownames(mat), , drop = FALSE]
  tt <- cbind(ID = rownames(tt), tt, stringsAsFactors = FALSE)
  tt$`Z-statistics` <- signed_z_from_p(tt$logFC, tt$P.Value)

  gain_samples <- rownames(pheno$sample_info)[pheno$sample_info$gain_1q == 1]
  neutral_samples <- rownames(pheno$sample_info)[pheno$sample_info$gain_1q == 0]
  tt[[paste0("Ave.", g0_label)]] <- rowMeans(mat[, neutral_samples, drop = FALSE])
  tt[[paste0("Ave.", g1_label)]] <- rowMeans(mat[, gain_samples, drop = FALSE])
  tt <- tt[order(tt$P.Value, decreasing = FALSE), , drop = FALSE]
  tt
}

build_design_for_project <- function(project_dt, covariates_used, strata_used) {
  sample_info <- as.data.frame(project_dt)
  rownames(sample_info) <- sample_info$aliquot_barcode

  design_df <- data.frame(row.names = sample_info$aliquot_barcode)
  design_df$gain_1q <- as.numeric(sample_info$gain_1q)

  all_adjusters <- unique(c(covariates_used, strata_used))
  used_vars <- character(0)
  dropped_vars <- character(0)
  added_terms <- character(0)

  for (var_name in all_adjusters) {
    if (!var_name %in% colnames(sample_info)) {
      dropped_vars <- c(dropped_vars, var_name)
      next
    }
    x <- sample_info[[var_name]]

    if (is_numeric_var(var_name)) {
      if (all(is.na(x)) || uniqueN(x[!is.na(x)]) <= 1L) {
        dropped_vars <- c(dropped_vars, var_name)
        next
      }
      miss <- is.na(x)
      x_imp <- x
      x_imp[miss] <- stats::median(x[!miss])
      z_col <- paste0(var_name, "_z")
      design_df[[z_col]] <- zscore(x_imp)
      added_terms <- c(added_terms, z_col)
      if (any(miss)) {
        miss_col <- paste0(var_name, "_missing")
        if (uniqueN(as.integer(miss)) > 1L) {
          design_df[[miss_col]] <- as.integer(miss)
          added_terms <- c(added_terms, miss_col)
        }
      }
      used_vars <- c(used_vars, var_name)
      next
    }

    if (is_binary_var(var_name)) {
      x_num <- suppressWarnings(as.numeric(x))
      if (all(is.na(x_num)) || uniqueN(x_num[!is.na(x_num)]) <= 1L) {
        dropped_vars <- c(dropped_vars, var_name)
        next
      }
      miss <- is.na(x_num)
      x_imp <- x_num
      x_imp[miss] <- as.integer(stats::median(x_num[!miss]) >= 0.5)
      design_df[[var_name]] <- as.integer(x_imp)
      added_terms <- c(added_terms, var_name)
      if (any(miss)) {
        miss_col <- paste0(var_name, "_missing")
        if (uniqueN(as.integer(miss)) > 1L) {
          design_df[[miss_col]] <- as.integer(miss)
          added_terms <- c(added_terms, miss_col)
        }
      }
      used_vars <- c(used_vars, var_name)
      next
    }

    x_chr <- trimws(as.character(x))
    x_chr[is.na(x_chr) | x_chr == ""] <- "Missing"
    if (uniqueN(x_chr) <= 1L) {
      dropped_vars <- c(dropped_vars, var_name)
      next
    }
    design_df[[var_name]] <- factor(x_chr)
    added_terms <- c(added_terms, var_name)
    used_vars <- c(used_vars, var_name)
  }

  rhs <- c("gain_1q", added_terms)
  if (length(rhs) == 0L) stop("No model terms available after preprocessing.")
  formula_txt <- paste("~", paste(rhs, collapse = " + "))
  design <- model.matrix(stats::as.formula(formula_txt), data = design_df)
  non_est <- limma::nonEstimable(design)
  if (!is.null(non_est) && length(non_est) > 0L) {
    if ("gain_1q" %in% non_est) stop("gain_1q became non-estimable in the design matrix.")
    design <- design[, !colnames(design) %in% non_est, drop = FALSE]
  }

  list(
    sample_info = sample_info,
    design_data = design_df,
    design = design,
    formula = formula_txt,
    used_vars = unique(used_vars),
    dropped_vars = unique(dropped_vars),
    dropped_design_columns = if (is.null(non_est)) character(0) else non_est
  )
}

get_project_subset <- function(net_eset, cohort_dt, project_prefix) {
  phe <- as.data.table(as.data.frame(Biobase::pData(net_eset)), keep.rownames = "aliquot_rowname")
  if (!"aliquot_barcode" %in% names(phe)) {
    setnames(phe, "aliquot_rowname", "aliquot_barcode")
  } else {
    phe[, aliquot_rowname := NULL]
  }
  phe[, sample_barcode := normalize_tcga_sample(sample_barcode)]
  phe[, patient_barcode := normalize_tcga_patient(patient_barcode)]

  project_cohort <- copy(cohort_dt[project == project_prefix])
  merged <- merge(
    phe,
    project_cohort,
    by.x = c("sample_barcode", "patient_barcode"),
    by.y = c("sample_id", "patient_id"),
    all.x = FALSE,
    all.y = FALSE,
    suffixes = c(".netbid", ".cov")
  )

  if (nrow(merged) == 0L) stop("No matched samples after merging NetBID and covariate cohort for ", project_prefix)
  merged <- merged[oneq_state %in% c("gain", "neutral")]
  merged[, histology_simple := collapse_factor(clean_histology(histological_type), max_levels = 3L, min_count = 15L)]
  merged[, grade_simple := collapse_factor(clean_grade(histological_grade), max_levels = 3L, min_count = 15L)]
  merged[, race_simple := collapse_factor(clean_race(race_raw), max_levels = 4L, min_count = 15L)]
  merged[, stage_factor := ifelse(is.na(stage_num), "Missing", as.character(as.integer(stage_num)))]
  merged[, stage_factor := factor(stage_factor)]
  merged
}

write_lines_file <- function(lines, path) {
  writeLines(as.character(lines), con = path, useBytes = TRUE)
}

invisible(db.preload(use_level = "gene", use_spe = "human", update = FALSE))
summary_dt <- fread(summary_file, sep = "\t")
summary_dt <- summary_dt[project %in% project_prefixes]
if (nrow(summary_dt) != length(project_prefixes)) {
  missing_projects <- setdiff(project_prefixes, summary_dt$project)
  stop("Missing project summary rows for: ", paste(missing_projects, collapse = ", "))
}

cohort_dt <- fread(cmd = paste("zcat -f", shQuote(cohort_file)))
invisible(sanitize_dt_names(cohort_dt))
cohort_dt[, sample_id := normalize_tcga_sample(sample_id)]
cohort_dt[, patient_id := normalize_tcga_patient(patient_id)]
cohort_dt[, project := toupper(trimws(as.character(project)))]
invisible(cohort_dt[!grepl("^TCGA-", project), project := paste0("TCGA-", project)])
cohort_dt <- cohort_dt[!is.na(sample_id) & !is.na(patient_id) & !is.na(project)]

arm_dt <- build_arm_metrics(arm_calls_file)
extra_clinical_dt <- prepare_extra_clinical(clinical_xlsx)

drop_cols <- intersect(c("oneq_state", "oneq_num", "gain_1q", "loss_1q", "arm_event_count_excl_1q", "onep19q_codel"), names(cohort_dt))
if (length(drop_cols) > 0L) cohort_dt[, (drop_cols) := NULL]
cohort_dt <- merge(cohort_dt, arm_dt, by = "sample_id", all.x = TRUE)
cohort_dt <- merge(cohort_dt, extra_clinical_dt, by = c("patient_id", "cancer_type"), all.x = TRUE)
cohort_dt <- unique(cohort_dt, by = c("project", "sample_id"))

project_qc_rows <- list()
qc_idx <- 1L

for (project_prefix in project_prefixes) {
  message("Processing adjusted NetBID2 analysis for ", project_prefix)
  run_dir <- find_latest_run(project_prefix)
  if (is.na(run_dir)) {
    warning("Skipping ", project_prefix, ": no NetBID run directory found.")
    next
  }

  summary_row <- summary_dt[project == project_prefix][1]
  covariates_used <- parse_var_list(summary_row$covariates_used)
  strata_used <- parse_var_list(summary_row$strata_used)
  analysis_objects <- prepare_analysis_objects(run_dir)
  project_dt <- get_project_subset(analysis_objects$net_eset, cohort_dt, project_prefix)

  n_gain <- sum(project_dt$gain_1q == 1, na.rm = TRUE)
  n_neutral <- sum(project_dt$gain_1q == 0, na.rm = TRUE)
  n_loss_excluded <- sum(cohort_dt[project == project_prefix]$oneq_state == "loss", na.rm = TRUE)
  if (n_gain < min_group_size || n_neutral < min_group_size) {
    warning(
      "Skipping ", project_prefix,
      ": insufficient matched samples after excluding losses (gain=", n_gain,
      ", neutral=", n_neutral, ")."
    )
    next
  }

  project_dt <- project_dt[order(aliquot_barcode)]
  pheno <- build_design_for_project(project_dt, covariates_used, strata_used)

  out_dir <- file.path(run_dir, "analysis", output_label)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  expr_mat <- Biobase::exprs(analysis_objects$net_eset)[, rownames(pheno$design), drop = FALSE]
  ac_mat <- Biobase::exprs(analysis_objects$merge_ac_eset)[, rownames(pheno$design), drop = FALSE]
  de_gene <- run_adjusted_limma(expr_mat, pheno)
  da_driver <- run_adjusted_limma(ac_mat, pheno)

  de_list <- setNames(list(de_gene), comp_name)
  da_list <- setNames(list(da_driver), comp_name)

  ms_tab <- generate.masterTable(
    use_comp = names(de_list),
    DE = de_list,
    DA = da_list,
    target_list = analysis_objects$merge_network$target_list,
    tf_sigs = tf_sigs,
    z_col = "Z-statistics",
    display_col = c("logFC", "P.Value"),
    main_id_type = "ensembl_gene_id"
  )
  ms_tab <- as.data.table(ms_tab)
  ms_tab <- ms_tab[Size >= min_regulon_size & Size <= max_regulon_size]
  fwrite(ms_tab, file.path(out_dir, "orig_driver_table.csv"))

  da_p_col <- paste0("P.Value.", comp_name, "_DA")
  da_ave_col <- paste0("Ave.", group1_label, ".", comp_name, "_DA")
  da_logfc_col <- paste0("logFC.", comp_name, "_DA")
  sig_pos <- ms_tab[get(da_p_col) < 0.05 & get(da_ave_col) > 0][order(get(da_p_col), -abs(get(da_logfc_col)))]
  fwrite(sig_pos, file.path(out_dir, "sig_pos_driver_table.csv"))

  fwrite(as.data.table(pheno$sample_info), file.path(out_dir, "sample_info.tsv"), sep = "\t")
  fwrite(as.data.table(pheno$design_data, keep.rownames = "aliquot_barcode"), file.path(out_dir, "design_data.tsv"), sep = "\t")
  fwrite(as.data.table(pheno$design, keep.rownames = "aliquot_barcode"), file.path(out_dir, "design_matrix.tsv"), sep = "\t")
  fwrite(as.data.table(de_gene), file.path(out_dir, "de_gene_limma.tsv"), sep = "\t")
  fwrite(as.data.table(da_driver), file.path(out_dir, "da_driver_limma.tsv"), sep = "\t")

  metadata_dt <- data.table(
    project = project_prefix,
    run_dir = run_dir,
    output_dir = out_dir,
    best_model = summary_row$best_model,
    summary_rule = summary_row$summary_rule,
    covariates_used_requested = paste(covariates_used, collapse = ","),
    strata_used_requested = paste(strata_used, collapse = ","),
    covariates_used_final = paste(pheno$used_vars, collapse = ","),
    covariates_dropped_final = paste(unique(c(pheno$dropped_vars, pheno$dropped_design_columns)), collapse = ","),
    model_formula = pheno$formula,
    min_regulon_size = min_regulon_size,
    max_regulon_size = max_regulon_size,
    n_matched_samples = nrow(project_dt),
    n_gain = n_gain,
    n_neutral = n_neutral,
    n_loss_excluded = n_loss_excluded
  )
  fwrite(metadata_dt, file.path(out_dir, "analysis_metadata.tsv"), sep = "\t")
  write_lines_file(pheno$formula, file.path(out_dir, "model_formula.txt"))

  save(
    analysis_objects,
    project_dt,
    pheno,
    de_gene,
    da_driver,
    ms_tab,
    sig_pos,
    metadata_dt,
    file = file.path(out_dir, "all_objects.RData")
  )

  project_qc_rows[[qc_idx]] <- metadata_dt
  qc_idx <- qc_idx + 1L
}

if (length(project_qc_rows) > 0L) {
  qc_dt <- rbindlist(project_qc_rows, use.names = TRUE, fill = TRUE)
  fwrite(qc_dt, file.path(netbid_dir, paste0(output_label, "_project_summary.tsv")), sep = "\t")
}
