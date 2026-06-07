#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)

local_r_lib <- Sys.getenv("LOCAL_R_LIB", unset = "")
if (nzchar(local_r_lib) && dir.exists(local_r_lib)) {
  .libPaths(c(local_r_lib, .libPaths()))
}

suppressPackageStartupMessages({
  library(data.table)
})

# ------------------------------ Args ------------------------------
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(key, default = NULL) {
  hit <- grep(paste0("^--", key, "="), args, value = TRUE)
  if (length(hit) == 0) return(default)
  sub(paste0("^--", key, "="), "", hit[1])
}

flag_arg <- function(key, default = FALSE) {
  x <- tolower(get_arg(key, if (default) "true" else "false"))
  x %in% c("true", "1", "yes", "y")
}

root <- get_arg("root", Sys.getenv("INFORM_ALL_ROOT", unset = file.path("external", "inform_all")))
mode <- tolower(get_arg("mode", "sample_genes_results"))
valid_modes <- c("sample_genes_results", "specimen_matrix")
if (!mode %in% valid_modes) {
  stop("Unsupported mode: ", mode, ". Supported modes: ", paste(valid_modes, collapse = ", "))
}

weights_file <- get_arg("weights", file.path(root, "reference_data", "mRNAsi_SC_PCBC_stemSig_weights.tsv"))
meta_file <- get_arg("meta", file.path(root, "wgcna_analysis2_reactome", "sample_info.tsv"))
covar_file <- get_arg("covar", file.path(root, "confounders_full", "covariates_full.tsv"))

expr_matrix_file <- get_arg(
  "expr_matrix",
  Sys.getenv("INFORM_TPM_SYMBOL_CSV", unset = file.path(root, "INFORM_EWS_BCOR_CIC_merged_TPM_SYMBOL.csv"))
)
specimen_meta_file <- get_arg("specimen_meta", file.path(root, "tp53_cdkn2a_specimen_status_resolved.tsv"))
specimen_cnv_file <- get_arg("specimen_cnv", file.path(root, "cnv_mapping_loss_breakdown.csv"))
manual_arm_file <- get_arg("manual_arm_file", file.path(root, "manual_final_inform.csv"))
processed_inform_dir <- get_arg(
  "processed_inform_dir",
  Sys.getenv("PROCESSED_INFORM_DIR", unset = file.path(root, "processed_inform"))
)
default_immune_candidates <- c(
  Sys.getenv("INFORM_IMMUNE_FILE", unset = ""),
  file.path(root, "results_immune_deconv_full_matrix", "deconvolution_long.csv"),
  file.path(root, "results_immune_deconv", "deconvolution_long.csv")
)
default_immune_candidates <- default_immune_candidates[nzchar(default_immune_candidates)]
default_immune_file <- default_immune_candidates[file.exists(default_immune_candidates)][1]
if (is.na(default_immune_file)) default_immune_file <- default_immune_candidates[1]
immune_file <- get_arg(
  "immune_file",
  default_immune_file
)
immune_method <- get_arg("immune_method", "xcell")
cell_panel_file <- get_arg("cell_panel_file", file.path(root, "cell_marker_custom_panels.tsv"))
allow_legacy_covariate_fallback <- flag_arg("allow_legacy_covariate_fallback", TRUE)
include_cellcycle <- flag_arg("include_cellcycle", TRUE)

gtf_paths <- strsplit(get_arg(
  "gtf_paths",
  paste(
    c(
      Sys.getenv("GTF_PATH", unset = ""),
      file.path(root, "1q_nfcore_rnaseq_results_grch38/genome/Homo_sapiens.GRCh38.dna.primary_assembly.filtered.gtf")
    ),
    collapse = ","
  )
), ",")[[1]]
gtf_paths <- gtf_paths[nzchar(gtf_paths)]

out_dir <- get_arg(
  "out_dir",
  Sys.getenv(
    "MRNASI_OUT_DIR",
    unset = if (mode == "specimen_matrix") "stemness_outputs_specimen_matrix" else "stemness_outputs"
  )
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_scores <- file.path(out_dir, "mRNAsi_scores.tsv")
out_score_meta <- file.path(out_dir, "mRNAsi_score_meta.tsv")
out_summary <- file.path(out_dir, "mRNAsi_cohort_summary.tsv")
out_models <- file.path(out_dir, "mRNAsi_models.tsv")
out_coef_table <- file.path(out_dir, "mRNAsi_adjusted_coefficients.tsv")
out_oneq_sens_tsv <- file.path(out_dir, "mRNAsi_1q_sensitivity.tsv")

covar_na_policy <- tolower(get_arg("covar_na_policy", "complete_case"))
if (!covar_na_policy %in% c("complete_case", "median_impute")) {
  message("Unknown covar_na_policy=", covar_na_policy, "; using complete_case")
  covar_na_policy <- "complete_case"
}

# ------------------------------ Helpers ------------------------------
stop_if_missing <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing)) stop("Missing required input(s):\n", paste(missing, collapse = "\n"))
}

coerce_binary <- function(x) {
  if (is.logical(x)) return(as.integer(x))
  if (is.numeric(x)) return(as.integer(x))
  x <- trimws(tolower(as.character(x)))
  out <- rep(NA_integer_, length(x))
  out[x %in% c("1", "true", "t", "yes", "y")] <- 1L
  out[x %in% c("0", "false", "f", "no", "n")] <- 0L
  out
}

get_gtf <- function(paths) {
  p <- paths[file.exists(paths)][1]
  if (is.na(p)) stop("No GTF found in paths.")
  p
}

extract_gene_map <- function(gtf_path) {
  con <- file(gtf_path, "r")
  on.exit(close(con), add = TRUE)

  gene_id <- character()
  gene_symbol <- character()
  gene_type <- character()

  while (TRUE) {
    line <- readLines(con, n = 1)
    if (length(line) == 0) break
    if (startsWith(line, "#")) next
    if (!grepl("\tgene\t", line, fixed = FALSE)) next

    gid <- sub('.*gene_id "([^"]+)".*', "\\1", line)
    gname <- sub('.*gene_name "([^"]+)".*', "\\1", line)

    gt <- NA_character_
    if (grepl('gene_type "', line, fixed = TRUE)) {
      gt <- sub('.*gene_type "([^"]+)".*', "\\1", line)
    } else if (grepl('gene_biotype "', line, fixed = TRUE)) {
      gt <- sub('.*gene_biotype "([^"]+)".*', "\\1", line)
    }

    if (gid != line && gname != line) {
      gene_id <- c(gene_id, gid)
      gene_symbol <- c(gene_symbol, gname)
      gene_type <- c(gene_type, gt)
    }
  }

  data.table(gene_id = gene_id, gene_symbol = gene_symbol, gene_type = gene_type)
}

load_or_build_gene_map <- function(gtf_paths, cache_path = "", protein_coding_only = TRUE) {
  if (nzchar(cache_path) && file.exists(cache_path)) {
    m <- readRDS(cache_path)
    if (is.data.frame(m) && all(c("gene_id", "gene_symbol") %in% names(m))) {
      m <- as.data.table(m)
      m[, gene_id := sub("\\..*$", "", gene_id)]
      if (protein_coding_only && "gene_type" %in% names(m)) {
        m <- m[is.na(gene_type) | gene_type == "protein_coding"]
      }
      return(m[, .(gene_id, gene_symbol)])
    }
  }

  gtf <- get_gtf(gtf_paths)
  m <- extract_gene_map(gtf)
  m[, gene_id := sub("\\..*$", "", gene_id)]
  m <- m[!is.na(gene_id) & gene_id != "" & !is.na(gene_symbol) & gene_symbol != ""]
  m <- unique(m, by = c("gene_id", "gene_symbol"))
  if (protein_coding_only) {
    m <- m[is.na(gene_type) | gene_type == "protein_coding"]
  }
  if (nzchar(cache_path)) saveRDS(m, cache_path)
  m[, .(gene_id, gene_symbol)]
}

read_gene_tpm <- function(path) {
  header <- strsplit(readLines(path, n = 1), "\t")[[1]]
  col_classes <- rep("NULL", length(header))
  col_classes[header == "gene_id"] <- "character"
  col_classes[header == "TPM"] <- "numeric"
  read.delim(path, colClasses = col_classes, check.names = FALSE)
}

load_bulk_tpm <- function(root_dir, sample_ids, gene_map) {
  gene_files <- list.files(root_dir, pattern = "\\.genes\\.results$", recursive = TRUE, full.names = TRUE)
  gene_samples <- sub("\\.genes\\.results$", "", basename(gene_files))
  keep_idx <- gene_samples %in% sample_ids
  gene_files <- gene_files[keep_idx]
  gene_samples <- gene_samples[keep_idx]
  if (anyDuplicated(gene_samples)) {
    dup_ids <- unique(gene_samples[duplicated(gene_samples)])
    warning(
      "Duplicate genes.results sample IDs detected; keeping first occurrence for: ",
      paste(dup_ids, collapse = ", ")
    )
    keep_first <- !duplicated(gene_samples)
    gene_files <- gene_files[keep_first]
    gene_samples <- gene_samples[keep_first]
  }
  if (length(gene_files) == 0) stop("No genes.results files for selected samples.")

  first <- read_gene_tpm(gene_files[1])
  genes <- first$gene_id
  tpm <- matrix(NA_real_, nrow = length(genes), ncol = length(gene_files))
  colnames(tpm) <- gene_samples
  tpm[, 1] <- first$TPM

  if (length(gene_files) > 1) {
    for (i in 2:length(gene_files)) {
      dt <- read_gene_tpm(gene_files[i])
      if (!identical(dt$gene_id, genes)) {
        dt <- dt[match(genes, dt$gene_id), ]
      }
      tpm[, i] <- dt$TPM
    }
  }
  rownames(tpm) <- genes

  ens_ids <- sub("\\..*$", "", rownames(tpm))
  map_idx <- match(ens_ids, gene_map$gene_id)
  map_syms <- gene_map$gene_symbol[map_idx]
  keep <- !is.na(map_syms) & map_syms != ""
  expr_mat <- tpm[keep, , drop = FALSE]
  expr_mat <- rowsum(expr_mat, group = map_syms[keep])
  expr_mat
}

load_symbol_tpm_matrix <- function(path) {
  if (!file.exists(path)) stop("Missing expression matrix: ", path)
  dt <- fread(path)
  if (ncol(dt) < 2) stop("Expression matrix needs at least one sample column: ", path)
  gene_col <- names(dt)[1]
  genes <- as.character(dt[[gene_col]])
  keep <- !is.na(genes) & genes != ""
  dt <- dt[keep]
  genes <- genes[keep]
  mat <- as.matrix(dt[, -1, with = FALSE])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- genes
  if (anyDuplicated(rownames(mat))) {
    mat <- rowsum(mat, group = rownames(mat), reorder = FALSE)
  }
  mat
}

load_weights <- function(path) {
  if (!file.exists(path)) stop("Missing weights file: ", path)
  dt <- fread(path)
  if (ncol(dt) < 2) stop("Weights file needs >=2 columns")
  gene_candidates <- c("gene", "Gene", "symbol", "Symbol", "gene_symbol", "genes")
  weight_candidates <- c("weight", "Weight", "coef", "Coef", "beta", "Beta", "coeff", "Coeff")
  gene_col <- gene_candidates[gene_candidates %in% names(dt)][1]
  weight_col <- weight_candidates[weight_candidates %in% names(dt)][1]
  if (is.na(gene_col) || is.na(weight_col)) {
    gene_col <- names(dt)[1]
    weight_col <- names(dt)[2]
  }
  w <- dt[, .(gene = as.character(get(gene_col)),
              weight = suppressWarnings(as.numeric(get(weight_col))))]
  w <- w[!is.na(gene) & gene != "" & !is.na(weight)]
  w <- w[!duplicated(gene)]
  w
}

compute_oclr_rankcorr <- function(expr_log2, weights) {
  idx <- match(weights$gene, rownames(expr_log2))
  keep <- !is.na(idx)
  if (sum(keep) < 10) stop("Insufficient overlap between weights and expression.")
  w_vec <- weights$weight[keep]
  expr_sub <- expr_log2[idx[keep], , drop = FALSE]
  w_rank <- rank(w_vec, ties.method = "average")
  raw <- apply(expr_sub, 2, function(x) {
    cor(w_rank, rank(x, ties.method = "average"), method = "pearson", use = "complete.obs")
  })
  min_v <- min(raw, na.rm = TRUE)
  max_v <- max(raw, na.rm = TRUE)
  scaled <- if (is.finite(min_v) && is.finite(max_v) && max_v > min_v) {
    (raw - min_v) / (max_v - min_v)
  } else {
    rep(NA_real_, length(raw))
  }
  list(raw = raw, scaled = scaled)
}

get_coef_table <- function(fit, se_type = "HC3") {
  se_type <- toupper(se_type)
  if (se_type %in% c("HC0", "HC1", "HC2", "HC3", "HC4", "HC5") &&
      requireNamespace("sandwich", quietly = TRUE) &&
      requireNamespace("lmtest", quietly = TRUE)) {
    robust <- suppressWarnings(tryCatch({
      vc <- sandwich::vcovHC(fit, type = se_type)
      lmtest::coeftest(fit, vcov. = vc)
    }, error = function(e) NULL))
    if (!is.null(robust)) {
      return(list(ct = robust, se_type = se_type))
    }
  }
  list(ct = summary(fit)$coefficients, se_type = "CLASSICAL")
}

extract_coef <- function(ct, term, se_type_used, df_resid) {
  if (!term %in% rownames(ct)) return(NULL)
  est <- ct[term, 1]
  se <- ct[term, 2]
  p <- ct[term, 4]
  crit <- qt(0.975, df = df_resid)
  ci_low <- est - crit * se
  ci_high <- est + crit * se
  data.frame(
    beta = est,
    se = se,
    p = p,
    ci_low = ci_low,
    ci_high = ci_high,
    se_type = se_type_used,
    stringsAsFactors = FALSE
  )
}

extract_all_coefs <- function(ct, se_type_used, df_resid, score_term, n_model,
                              score_label = NA_character_, model_variant = NA_character_,
                              variant_label = NA_character_) {
  if (is.null(ct) || nrow(ct) == 0) return(data.table())
  est <- as.numeric(ct[, 1])
  se <- as.numeric(ct[, 2])
  p <- as.numeric(ct[, 4])
  crit <- qt(0.975, df = df_resid)
  data.table(
    score = score_term,
    score_label = if (is.na(score_label)) score_term else score_label,
    model_variant = if (is.na(model_variant)) "default" else model_variant,
    variant_label = if (is.na(variant_label)) "Default adjustment" else variant_label,
    term = rownames(ct),
    beta = est,
    se = se,
    p = p,
    ci_low = est - crit * se,
    ci_high = est + crit * se,
    se_type = se_type_used,
    n = n_model
  )
}

scale_rows <- function(mat) {
  z <- t(scale(t(mat)))
  z[is.na(z)] <- 0
  z
}

standardize_to_reference <- function(x, ref_idx) {
  x <- as.numeric(x)
  ref_idx <- as.logical(ref_idx)
  ref <- x[ref_idx & !is.na(x)]
  if (length(ref) == 0) return(rep(NA_real_, length(x)))
  mu <- mean(ref)
  sdv <- stats::sd(ref)
  if (!is.finite(sdv) || sdv == 0) return(rep(0, length(x)))
  (x - mu) / sdv
}

term_label <- function(term) {
  labels <- c(
    "(Intercept)" = "Intercept",
    "gain_1q" = "1q gain",
    "gain_8p" = "8p gain",
    "gain_12p" = "12p gain",
    "gain_12q" = "12q gain",
    "loss_16q" = "16q loss",
    "ImmuneScore" = "Immune score",
    "StromaScore" = "Stroma score",
    "cellcycle_score" = "Cell-cycle score",
    "metastasis_specimen" = "Metastatic specimen",
    "tp53_mutated" = "TP53 mutated",
    "cdkn2a_mutated" = "CDKN2A mutated"
  )
  out <- labels[term]
  out[is.na(out)] <- term[is.na(out)]
  unname(out)
}

model_variant_label <- function(x) {
  labels <- c(
    "base_adjusted" = "Base adjusted",
    "base_plus_cellcycle" = "Base adjusted + cell-cycle",
    "base_plus_tp53" = "Base adjusted + TP53",
    "base_plus_tp53_cellcycle" = "Base adjusted + TP53 + cell-cycle",
    "single_model" = "Single adjusted model"
  )
  out <- labels[x]
  out[is.na(out)] <- x[is.na(out)]
  unname(out)
}

format_covariate_caption <- function(covars, prefix) {
  covars <- unique(covars)
  if (length(covars) == 0) {
    return(paste0(prefix, " no additional covariates"))
  }
  txt <- paste(term_label(covars), collapse = ", ")
  paste(strwrap(paste0(prefix, " ", txt), width = 88), collapse = "\n")
}

arm_status_to_binary <- function(x, event = c("gain", "loss")) {
  event <- match.arg(event)
  x <- trimws(tolower(as.character(x)))
  out <- rep(NA_integer_, length(x))
  out[x == paste("complete", event)] <- 1L
  out[x %in% c("clean", "neutral")] <- 0L
  out
}

specimen_metastasis_binary <- function(x) {
  x <- trimws(tolower(as.character(x)))
  out <- rep(NA_integer_, length(x))
  out[x == ""] <- NA_integer_
  out[grepl("^metastasis", x)] <- 1L
  out[!is.na(x) & x != "" & is.na(out)] <- 0L
  out
}

load_panel_genes <- function(path, panel_name) {
  if (!file.exists(path)) return(character())
  dt <- fread(path)
  if (!all(c("panel", "genes") %in% names(dt))) return(character())
  hit <- dt[panel == panel_name]
  if (nrow(hit) == 0) return(character())
  unique(trimws(unlist(strsplit(hit$genes[1], ",", fixed = TRUE))))
}

compute_panel_score <- function(expr_log2, genes) {
  genes <- unique(genes)
  genes <- intersect(genes, rownames(expr_log2))
  out <- rep(NA_real_, ncol(expr_log2))
  names(out) <- colnames(expr_log2)
  if (length(genes) < 2) return(out)
  z <- scale_rows(expr_log2[genes, , drop = FALSE])
  score <- colMeans(z, na.rm = TRUE)
  score <- as.numeric(scale(score))
  names(score) <- colnames(expr_log2)
  score
}

load_specimen_immune_scores <- function(path, method = "xcell") {
  if (!file.exists(path)) stop("Missing immune file: ", path)
  dt <- fread(path)
  req <- c("cell_type", "sample_id", "score", "method")
  if (!all(req %in% names(dt))) {
    stop("Immune file missing required columns: ", paste(setdiff(req, names(dt)), collapse = ", "))
  }
  keep_types <- c(
    "immune score",
    "stroma score",
    "ImmuneScore",
    "StromaScore",
    "Endothelial cell",
    "Endothelial cells",
    "ly Endothelial cells",
    "mv Endothelial cells"
  )
  method_name <- method
  sub <- dt[get("method") == method_name & cell_type %in% keep_types, .(sample_id, cell_type, score)]
  if (nrow(sub) == 0) {
    return(data.table(
      sample_id = character(),
      ImmuneScore = numeric(),
      StromaScore = numeric(),
      EndothelialScore = numeric()
    ))
  }
  wide <- dcast(
    sub,
    sample_id ~ cell_type,
    value.var = "score",
    fun.aggregate = function(x) x[1]
  )
  out <- data.table(sample_id = wide$sample_id)
  if ("ImmuneScore" %in% names(wide)) out[, ImmuneScore := wide[["ImmuneScore"]]]
  if (!"ImmuneScore" %in% names(out) && "immune score" %in% names(wide)) out[, ImmuneScore := wide[["immune score"]]]
  if ("StromaScore" %in% names(wide)) out[, StromaScore := wide[["StromaScore"]]]
  if (!"StromaScore" %in% names(out) && "stroma score" %in% names(wide)) out[, StromaScore := wide[["stroma score"]]]
  if ("EndothelialScore" %in% names(wide)) out[, EndothelialScore := wide[["EndothelialScore"]]]
  if (!"EndothelialScore" %in% names(out) && "Endothelial cell" %in% names(wide)) out[, EndothelialScore := wide[["Endothelial cell"]]]
  if (!"EndothelialScore" %in% names(out) && "Endothelial cells" %in% names(wide)) out[, EndothelialScore := wide[["Endothelial cells"]]]
  if ("ly Endothelial cells" %in% names(wide)) out[, LymphaticEndothelialScore := wide[["ly Endothelial cells"]]]
  if ("mv Endothelial cells" %in% names(wide)) out[, MicrovascularEndothelialScore := wide[["mv Endothelial cells"]]]
  out
}

load_manual_cnv_calls <- function(path) {
  if (!file.exists(path)) stop("Missing manual CNV file: ", path)
  dt <- fread(path)
  if (!"patient_ID" %in% names(dt)) stop("Manual CNV file missing patient_ID: ", path)

  pick_col <- function(options) {
    hit <- options[options %in% names(dt)][1]
    if (is.na(hit)) stop("Manual CNV file missing expected column(s): ", paste(options, collapse = " / "))
    hit
  }

  col_1q <- pick_col(c("1q", "X1q"))
  col_8p <- pick_col(c("8p", "X8p"))
  col_12p <- pick_col(c("12p", "X12p"))
  col_12q <- pick_col(c("12q", "X12q"))
  col_16q <- pick_col(c("16q", "X16q"))

  unique(
    dt[, .(
      sentrix = as.character(patient_ID),
      gain_1q = arm_status_to_binary(get(col_1q), event = "gain"),
      gain_8p = arm_status_to_binary(get(col_8p), event = "gain"),
      gain_12p = arm_status_to_binary(get(col_12p), event = "gain"),
      gain_12q = arm_status_to_binary(get(col_12q), event = "gain"),
      loss_16q = arm_status_to_binary(get(col_16q), event = "loss")
    )],
    by = "sentrix"
  )
}

classify_segment_arm_binary <- function(seg_dt, chromosome, arm_side, event = c("gain", "loss")) {
  event <- match.arg(event)
  sub <- seg_dt[Chromosome == chromosome & Arm == arm_side]
  if (nrow(sub) == 0) return(NA_integer_)

  status <- trimws(tolower(as.character(sub$CNA_status)))
  status[is.na(status) | status == ""] <- "unknown"
  pct <- suppressWarnings(as.numeric(sub$percent_affected))
  pct[!is.finite(pct)] <- 0

  pct_by_status <- tapply(pct, status, sum)
  target_pct <- unname(if (event %in% names(pct_by_status)) pct_by_status[[event]] else 0)
  neutral_pct <- unname(if ("neutral" %in% names(pct_by_status)) pct_by_status[["neutral"]] else 0)
  other_nonneutral <- sum(pct_by_status[setdiff(names(pct_by_status), c("neutral", event))], na.rm = TRUE)

  if (all(status == "neutral")) return(0L)
  if (target_pct >= 80 && other_nonneutral == 0) return(1L)
  if (target_pct == 0 && other_nonneutral == 0 && neutral_pct >= 80) return(0L)
  NA_integer_
}

load_segment_fallback_cnv <- function(sentrix_ids, processed_inform_dir) {
  sentrix_ids <- unique(as.character(sentrix_ids))
  sentrix_ids <- sentrix_ids[!is.na(sentrix_ids) & sentrix_ids != ""]
  if (length(sentrix_ids) == 0) {
    return(data.table(
      sentrix = character(),
      gain_1q = integer(),
      gain_8p = integer(),
      gain_12p = integer(),
      gain_12q = integer(),
      loss_16q = integer()
    ))
  }

  rows <- vector("list", length(sentrix_ids))
  keep <- 0L
  for (sid in sentrix_ids) {
    seg_path <- file.path(processed_inform_dir, sid, "CNVsegments_modified.seg")
    if (!file.exists(seg_path)) next
    seg_dt <- tryCatch(
      fread(seg_path, select = c("Chromosome", "Arm", "percent_affected", "CNA_status")),
      error = function(e) NULL
    )
    if (is.null(seg_dt) || nrow(seg_dt) == 0) next
    keep <- keep + 1L
    rows[[keep]] <- data.table(
      sentrix = sid,
      gain_1q = classify_segment_arm_binary(seg_dt, "chr1", "q", event = "gain"),
      gain_8p = classify_segment_arm_binary(seg_dt, "chr8", "p", event = "gain"),
      gain_12p = classify_segment_arm_binary(seg_dt, "chr12", "p", event = "gain"),
      gain_12q = classify_segment_arm_binary(seg_dt, "chr12", "q", event = "gain"),
      loss_16q = classify_segment_arm_binary(seg_dt, "chr16", "q", event = "loss")
    )
  }

  if (keep == 0L) {
    return(data.table(
      sentrix = character(),
      gain_1q = integer(),
      gain_8p = integer(),
      gain_12p = integer(),
      gain_12q = integer(),
      loss_16q = integer()
    ))
  }
  rbindlist(rows[seq_len(keep)], fill = TRUE)
}

count_stage <- function(dt, label) {
  data.table(
    stage = label,
    n = nrow(dt),
    n_gain = if ("gain_1q" %in% names(dt)) sum(dt$gain_1q == 1, na.rm = TRUE) else NA_integer_,
    n_nogain = if ("gain_1q" %in% names(dt)) sum(dt$gain_1q == 0, na.rm = TRUE) else NA_integer_,
    n_tp53_mut = if ("tp53_mutated" %in% names(dt)) sum(dt$tp53_mutated == 1, na.rm = TRUE) else NA_integer_,
    n_cdkn2a_mut = if ("cdkn2a_mutated" %in% names(dt)) sum(dt$cdkn2a_mutated == 1, na.rm = TRUE) else NA_integer_
  )
}

build_sample_meta <- function(meta_file, covar_file) {
  stop_if_missing(c(meta_file, covar_file))
  meta <- fread(meta_file)
  if (!"sample_id" %in% names(meta)) stop("sample_info.tsv missing sample_id")
  if (!"gain_1q" %in% names(meta)) stop("sample_info.tsv missing gain_1q")

  logical_cols <- setdiff(names(meta), c("inform_id", "sample_id", "group"))
  for (col in logical_cols) {
    meta[, (col) := coerce_binary(get(col))]
  }

  covar_dt <- fread(covar_file)
  if (!"sample_id" %in% names(covar_dt)) stop("covariates_full.tsv missing sample_id")
  meta <- merge(meta, covar_dt, by = "sample_id", all.x = TRUE)
  meta
}

build_specimen_meta <- function(specimen_meta_file,
                                specimen_cnv_file,
                                manual_arm_file,
                                processed_inform_dir,
                                immune_file,
                                immune_method,
                                legacy_meta_file,
                                legacy_covar_file,
                                expr_log2,
                                cell_panel_file,
                                include_cellcycle,
                                allow_legacy_covariate_fallback) {
  stop_if_missing(c(specimen_meta_file, immune_file, legacy_meta_file, manual_arm_file))
  meta <- fread(specimen_meta_file)
  req <- c(
    "inform_id",
    "unique_tumor_sample_id",
    "sentrix",
    "tp53_status",
    "cdkn2a_status",
    "inform_exome_specimen_count"
  )
  missing_req <- setdiff(req, names(meta))
  if (length(missing_req)) {
    stop("Specimen meta missing required columns: ", paste(missing_req, collapse = ", "))
  }

  meta <- unique(meta, by = "unique_tumor_sample_id")
  meta[, sample_id := unique_tumor_sample_id]
  meta[, tp53_mutated := as.integer(tp53_status == "somatic_functional_mutated")]
  meta[, cdkn2a_mutated := as.integer(cdkn2a_status == "somatic_functional_mutated")]
  metastasis_label <- if ("meth_sample_type_tumor" %in% names(meta)) {
    meta[["meth_sample_type_tumor"]]
  } else if ("exome_specimen_label" %in% names(meta)) {
    meta[["exome_specimen_label"]]
  } else {
    rep(NA_character_, nrow(meta))
  }
  meta[, metastasis_specimen := specimen_metastasis_binary(metastasis_label)]

  cnv_map <- data.table(
    sentrix = character(),
    gain_1q = integer(),
    gain_8p = integer(),
    gain_12p = integer(),
    gain_12q = integer(),
    loss_16q = integer()
  )
  if (file.exists(specimen_cnv_file)) {
    cnv <- fread(specimen_cnv_file)
    cnv_req <- c("sentrix", "oneq_gain", "gain_8p", "gain_12p", "gain_12q", "loss_16q")
    missing_cnv <- setdiff(cnv_req, names(cnv))
    if (length(missing_cnv)) {
      stop("Specimen CNV file missing required columns: ", paste(missing_cnv, collapse = ", "))
    }
    cnv_map <- unique(
      cnv[, .(
        sentrix,
        gain_1q = coerce_binary(oneq_gain),
        gain_8p = coerce_binary(gain_8p),
        gain_12p = coerce_binary(gain_12p),
        gain_12q = coerce_binary(gain_12q),
        loss_16q = coerce_binary(loss_16q)
      )],
      by = "sentrix"
    )
  }
  manual_cnv <- load_manual_cnv_calls(manual_arm_file)

  setnames(cnv_map, old = c("gain_1q", "gain_8p", "gain_12p", "gain_12q", "loss_16q"),
           new = c("map_gain_1q", "map_gain_8p", "map_gain_12p", "map_gain_12q", "map_loss_16q"))
  setnames(manual_cnv, old = c("gain_1q", "gain_8p", "gain_12p", "gain_12q", "loss_16q"),
           new = c("manual_gain_1q", "manual_gain_8p", "manual_gain_12p", "manual_gain_12q", "manual_loss_16q"))
  meta <- merge(meta, cnv_map, by = "sentrix", all.x = TRUE, sort = FALSE)
  meta <- merge(meta, manual_cnv, by = "sentrix", all.x = TRUE, sort = FALSE)

  seg_needed <- meta[
    is.na(map_gain_1q) | is.na(map_gain_8p) | is.na(map_gain_12p) | is.na(map_gain_12q) | is.na(map_loss_16q),
    unique(sentrix)
  ]
  seg_cnv <- load_segment_fallback_cnv(seg_needed, processed_inform_dir)
  setnames(seg_cnv, old = c("gain_1q", "gain_8p", "gain_12p", "gain_12q", "loss_16q"),
           new = c("segment_gain_1q", "segment_gain_8p", "segment_gain_12p", "segment_gain_12q", "segment_loss_16q"))
  meta <- merge(meta, seg_cnv, by = "sentrix", all.x = TRUE, sort = FALSE)
  meta <- copy(meta)

  cnv_pairs <- list(
    gain_1q = c("map_gain_1q", "manual_gain_1q", "segment_gain_1q"),
    gain_8p = c("map_gain_8p", "manual_gain_8p", "segment_gain_8p"),
    gain_12p = c("map_gain_12p", "manual_gain_12p", "segment_gain_12p"),
    gain_12q = c("map_gain_12q", "manual_gain_12q", "segment_gain_12q"),
    loss_16q = c("map_loss_16q", "manual_loss_16q", "segment_loss_16q")
  )
  cnv_source_labels <- c("specimen_cnv_mapping", "manual_sentrix", "segment_fallback")
  cnv_source_cols <- character()
  for (nm in names(cnv_pairs)) {
    src_col <- paste0(nm, "_source")
    set(meta, j = nm, value = rep(NA_integer_, nrow(meta)))
    set(meta, j = src_col, value = rep(NA_character_, nrow(meta)))
    for (i in seq_along(cnv_pairs[[nm]])) {
      source_col <- cnv_pairs[[nm]][i]
      idx <- is.na(meta[[nm]]) & !is.na(meta[[source_col]])
      if (any(idx)) {
        meta[[nm]][idx] <- meta[[source_col]][idx]
        meta[[src_col]][idx] <- cnv_source_labels[i]
      }
    }
    cnv_source_cols <- c(cnv_source_cols, src_col)
  }

  immune_dt <- load_specimen_immune_scores(immune_file, immune_method)
  meta <- merge(meta, immune_dt, by = "sample_id", all.x = TRUE, sort = FALSE)

  meta[, cnv_source := apply(.SD, 1, function(x) {
    x <- unique(stats::na.omit(x))
    if (length(x) == 0) return("missing")
    if (length(x) == 1) return(x)
    "mixed_sources"
  }), .SDcols = cnv_source_cols]
  meta[, immune_source := ifelse(
    complete.cases(.SD),
    paste0("specimen_", immune_method),
    "missing"
  ), .SDcols = c("ImmuneScore", "StromaScore")]

  if (allow_legacy_covariate_fallback) {
    stop_if_missing(c(legacy_covar_file))

    legacy_meta <- fread(legacy_meta_file)
    legacy_meta_req <- c("inform_id", "gain_1q", "gain_8p", "gain_12p", "gain_12q", "loss_16q")
    missing_legacy_meta <- setdiff(legacy_meta_req, names(legacy_meta))
    if (length(missing_legacy_meta)) {
      stop("Legacy sample_info missing required columns: ", paste(missing_legacy_meta, collapse = ", "))
    }
    legacy_meta <- unique(
      legacy_meta[, .(
        inform_id,
        legacy_gain_1q = coerce_binary(gain_1q),
        legacy_gain_8p = coerce_binary(gain_8p),
        legacy_gain_12p = coerce_binary(gain_12p),
        legacy_gain_12q = coerce_binary(gain_12q),
        legacy_loss_16q = coerce_binary(loss_16q)
      )],
      by = "inform_id"
    )
    meta <- merge(meta, legacy_meta, by = "inform_id", all.x = TRUE, sort = FALSE)

    legacy_cov <- fread(legacy_covar_file)
    legacy_cov_req <- c("sample_id", "ImmuneScore", "StromaScore")
    missing_legacy_cov <- setdiff(legacy_cov_req, names(legacy_cov))
    if (length(missing_legacy_cov)) {
      stop("Legacy covariate file missing required columns: ", paste(missing_legacy_cov, collapse = ", "))
    }
    keep_cov <- c("sample_id", "ImmuneScore", "StromaScore")
    if ("cellcycle_z" %in% names(legacy_cov)) keep_cov <- c(keep_cov, "cellcycle_z")
    legacy_cov <- unique(legacy_cov[, ..keep_cov], by = "sample_id")
    setnames(
      legacy_cov,
      old = keep_cov,
      new = c("rna_sample_id", "legacy_ImmuneScore", "legacy_StromaScore", if ("cellcycle_z" %in% keep_cov) "legacy_cellcycle_z"),
      skip_absent = TRUE
    )
    if ("rna_sample_id" %in% names(meta)) {
      meta <- merge(meta, legacy_cov, by = "rna_sample_id", all.x = TRUE, sort = FALSE)
    }

    single_specimen <- meta$inform_exome_specimen_count == 1
    legacy_pairs <- list(
      gain_1q = "legacy_gain_1q",
      gain_8p = "legacy_gain_8p",
      gain_12p = "legacy_gain_12p",
      gain_12q = "legacy_gain_12q",
      loss_16q = "legacy_loss_16q"
    )
    used_legacy_cnv <- rep(FALSE, nrow(meta))
    for (nm in names(legacy_pairs)) {
      legacy_nm <- legacy_pairs[[nm]]
      idx <- is.na(meta[[nm]]) & single_specimen & !is.na(meta[[legacy_nm]])
      if (any(idx)) {
        meta[[nm]][idx] <- meta[[legacy_nm]][idx]
        used_legacy_cnv[idx] <- TRUE
      }
    }
    meta[used_legacy_cnv & cnv_source == "missing", cnv_source := "legacy_single_specimen"]

    used_legacy_immune <- rep(FALSE, nrow(meta))
    for (nm in c("ImmuneScore", "StromaScore")) {
      legacy_nm <- paste0("legacy_", nm)
      idx <- is.na(meta[[nm]]) & !is.na(meta[[legacy_nm]])
      if (any(idx)) {
        meta[[nm]][idx] <- meta[[legacy_nm]][idx]
        used_legacy_immune[idx] <- TRUE
      }
    }
    meta[used_legacy_immune & immune_source == "missing", immune_source := "legacy_rna_covariates"]
  }

  if (include_cellcycle) {
    cellcycle_genes <- load_panel_genes(cell_panel_file, "CM_Progenitor_cycling")
    cellcycle_score <- compute_panel_score(expr_log2, cellcycle_genes)
    cellcycle_dt <- data.table(sample_id = names(cellcycle_score), cellcycle_score = as.numeric(cellcycle_score))
    meta <- merge(meta, cellcycle_dt, by = "sample_id", all.x = TRUE, sort = FALSE)
    meta[, cellcycle_source := ifelse(!is.na(cellcycle_score), "expr_panel", "missing")]

    if (allow_legacy_covariate_fallback && "legacy_cellcycle_z" %in% names(meta)) {
      idx <- is.na(meta$cellcycle_score) & !is.na(meta$legacy_cellcycle_z)
      if (any(idx)) {
        meta$cellcycle_score[idx] <- meta$legacy_cellcycle_z[idx]
        meta$cellcycle_source[idx] <- "legacy_rna_covariates"
      }
    }
  }

  meta
}

# ------------------------------ Build expression + metadata ------------------------------
stop_if_missing(weights_file)
weights <- load_weights(weights_file)

if (mode == "sample_genes_results") {
  meta <- build_sample_meta(meta_file, covar_file)
  cache_path <- file.path(root, "tcga_stemness_scores", "gene_map_GRCh38.rds")
  map_dt <- load_or_build_gene_map(gtf_paths, cache_path = cache_path, protein_coding_only = TRUE)
  expr_mat <- load_bulk_tpm(root, meta$sample_id, map_dt)
  expr_log2 <- log2(expr_mat + 1)

  oclr <- compute_oclr_rankcorr(expr_log2, weights)
  score_dt <- data.table(
    sample_id = colnames(expr_log2),
    mRNAsi_raw = oclr$raw[colnames(expr_log2)],
    mRNAsi = oclr$scaled[colnames(expr_log2)]
  )
  score_meta <- merge(score_dt, meta, by = "sample_id", all.x = TRUE)

  exclude_cols <- c("sample_id", "inform_id", "group", "gain_1q")
  candidate_covars <- setdiff(names(score_meta), exclude_cols)
  candidate_covars <- setdiff(candidate_covars, c("mRNAsi_raw", "mRNAsi", "mRNAsi_z"))

  summary_dt <- rbindlist(list(
    count_stage(meta, "meta_rows"),
    count_stage(score_meta, "expression_overlap")
  ), fill = TRUE)
} else {
  stop_if_missing(c(expr_matrix_file, specimen_meta_file, immune_file, manual_arm_file))
  expr_mat <- load_symbol_tpm_matrix(expr_matrix_file)
  expr_log2 <- log2(expr_mat + 1)

  meta <- build_specimen_meta(
    specimen_meta_file = specimen_meta_file,
    specimen_cnv_file = specimen_cnv_file,
    manual_arm_file = manual_arm_file,
    processed_inform_dir = processed_inform_dir,
    immune_file = immune_file,
    immune_method = immune_method,
    legacy_meta_file = meta_file,
    legacy_covar_file = covar_file,
    expr_log2 = expr_log2,
    cell_panel_file = cell_panel_file,
    include_cellcycle = include_cellcycle,
    allow_legacy_covariate_fallback = allow_legacy_covariate_fallback
  )
  meta <- meta[sample_id %in% colnames(expr_log2)]

  oclr <- compute_oclr_rankcorr(expr_log2, weights)
  score_dt <- data.table(
    sample_id = colnames(expr_log2),
    mRNAsi_raw = oclr$raw[colnames(expr_log2)],
    mRNAsi = oclr$scaled[colnames(expr_log2)]
  )
  score_meta <- merge(score_dt, meta, by = "sample_id", all.x = FALSE, sort = FALSE)

  candidate_covars <- c(
    "gain_8p",
    "gain_12p",
    "gain_12q",
    "loss_16q",
    "ImmuneScore",
    "StromaScore",
    if (include_cellcycle) "cellcycle_score" else NULL,
    "metastasis_specimen",
    "tp53_mutated",
    "cdkn2a_mutated"
  )

  summary_dt <- rbindlist(list(
    count_stage(meta, "expr_overlap"),
    count_stage(meta[complete.cases(meta[, .(ImmuneScore, StromaScore)])], "with_immune_stroma"),
    count_stage(meta[complete.cases(meta[, .(gain_1q, gain_8p, gain_12p, gain_12q, loss_16q)])], "with_full_cnv"),
    count_stage(meta[complete.cases(meta[, .(gain_1q, gain_8p, gain_12p, gain_12q, loss_16q, ImmuneScore, StromaScore)])], "with_cnv_and_immune")
  ), fill = TRUE)
}

# ------------------------------ Modeling ------------------------------
candidate_covars <- candidate_covars[candidate_covars %in% names(score_meta)]

for (cv in candidate_covars) {
  if (!is.numeric(score_meta[[cv]])) {
    score_meta[[cv]] <- suppressWarnings(as.numeric(score_meta[[cv]]))
  }
}
if (!is.numeric(score_meta$gain_1q)) {
  score_meta$gain_1q <- suppressWarnings(as.numeric(score_meta$gain_1q))
}

covars_all <- candidate_covars[vapply(candidate_covars, function(v) {
  x <- score_meta[[v]]
  if (all(is.na(x))) return(FALSE)
  ux <- unique(x[!is.na(x)])
  length(ux) > 1
}, logical(1))]

has_cellcycle_variant <- "cellcycle_score" %in% covars_all
has_tp53_variant <- "tp53_mutated" %in% covars_all
base_covars <- setdiff(covars_all, c("cellcycle_score", "tp53_mutated"))
variant_defs <- list(base_adjusted = base_covars)
if (has_cellcycle_variant) {
  variant_defs$base_plus_cellcycle <- c(base_covars, "cellcycle_score")
}
if (has_tp53_variant) {
  variant_defs$base_plus_tp53 <- c(base_covars, "tp53_mutated")
}
if (has_tp53_variant && has_cellcycle_variant) {
  variant_defs$base_plus_tp53_cellcycle <- c(base_covars, "tp53_mutated", "cellcycle_score")
}
if (length(variant_defs) == 0) {
  variant_defs <- list(single_model = covars_all)
}
variant_names <- names(variant_defs)
default_variant <- if ("base_plus_tp53_cellcycle" %in% variant_names) {
  "base_plus_tp53_cellcycle"
} else if ("base_plus_tp53" %in% variant_names) {
  "base_plus_tp53"
} else if ("base_plus_cellcycle" %in% variant_names) {
  "base_plus_cellcycle"
} else {
  variant_names[1]
}
base_adjusted_caption <- format_covariate_caption(base_covars, "Base adjusted =")

base_keep_cols <- unique(c("sample_id", "mRNAsi", "gain_1q", unlist(variant_defs, use.names = FALSE)))
base_model_df <- as.data.frame(score_meta[, base_keep_cols, with = FALSE])

if (covar_na_policy == "median_impute") {
  for (cv in setdiff(base_keep_cols, c("sample_id", "mRNAsi", "gain_1q"))) {
    x <- base_model_df[[cv]]
    if (!is.numeric(x)) next
    med <- suppressWarnings(stats::median(x, na.rm = TRUE))
    if (is.finite(med)) base_model_df[[cv]][is.na(x)] <- med
  }
}

model_rows <- list()
coef_rows <- list()
score_terms <- c("mRNAsi", "mRNAsi_z")
summary_variant_rows <- list()
for (variant_name in variant_names) {
  variant_covars <- variant_defs[[variant_name]]
  flag_col <- paste0("model_complete_case_", variant_name)
  z_col <- paste0("mRNAsi_z_", variant_name)
  required_cols <- c("mRNAsi", "gain_1q", variant_covars)
  keep_variant <- complete.cases(base_model_df[, required_cols, drop = FALSE])
  ref_ids <- base_model_df$sample_id[keep_variant]
  score_meta[, (flag_col) := sample_id %in% ref_ids]
  score_meta[, (z_col) := standardize_to_reference(mRNAsi, get(flag_col))]
  summary_variant_rows[[length(summary_variant_rows) + 1L]] <- count_stage(
    score_meta[get(flag_col) == TRUE],
    paste0("model_complete_case_", variant_name)
  )

  variant_df <- base_model_df[keep_variant, c("sample_id", "mRNAsi", "gain_1q", variant_covars), drop = FALSE]
  if (nrow(variant_df) < 5) next
  variant_df$mRNAsi_z <- score_meta[[z_col]][match(variant_df$sample_id, score_meta$sample_id)]
  model_subset_meta <- score_meta[get(flag_col) == TRUE]

  for (score_term in score_terms) {
    form <- as.formula(paste(score_term, "~", paste(c("gain_1q", variant_covars), collapse = " + ")))
    fit <- lm(form, data = variant_df)
    sm <- summary(fit)
    coef_info <- get_coef_table(fit, se_type = "HC3")
    gain_row <- extract_coef(coef_info$ct, "gain_1q", coef_info$se_type, fit$df.residual)
    if (!is.null(gain_row)) {
      if (!is.finite(gain_row$se) || !is.finite(gain_row$p)) {
        coef_info <- get_coef_table(fit, se_type = "CLASSICAL")
        gain_row <- extract_coef(coef_info$ct, "gain_1q", coef_info$se_type, fit$df.residual)
      }
    }
    coef_rows[[length(coef_rows) + 1L]] <- extract_all_coefs(
      coef_info$ct,
      se_type_used = coef_info$se_type,
      df_resid = fit$df.residual,
      score_term = score_term,
      n_model = nrow(variant_df),
      score_label = if (score_term == "mRNAsi_z") "mRNAsi (z-score)" else "mRNAsi",
      model_variant = variant_name,
      variant_label = model_variant_label(variant_name)
    )
    if (!is.null(gain_row)) {
      model_rows[[length(model_rows) + 1L]] <- cbind(
        data.frame(
          mode = mode,
          model_variant = variant_name,
          variant_label = model_variant_label(variant_name),
          score = score_term,
          n = nrow(variant_df),
          n_gain = sum(variant_df$gain_1q == 1, na.rm = TRUE),
          n_nogain = sum(variant_df$gain_1q == 0, na.rm = TRUE),
          n_tp53_mut = if ("tp53_mutated" %in% names(model_subset_meta)) sum(model_subset_meta$tp53_mutated == 1, na.rm = TRUE) else NA_integer_,
          n_cdkn2a_mut = if ("cdkn2a_mutated" %in% names(model_subset_meta)) sum(model_subset_meta$cdkn2a_mutated == 1, na.rm = TRUE) else NA_integer_,
          covariates = paste(variant_covars, collapse = ","),
          adj_r2 = sm$adj.r.squared,
          na_policy = covar_na_policy,
          stringsAsFactors = FALSE
        ),
        gain_row
      )
    }
  }
}

score_meta[, model_complete_case := get(paste0("model_complete_case_", default_variant))]
score_meta[, mRNAsi_z := get(paste0("mRNAsi_z_", default_variant))]
score_dt_cols <- unique(c("sample_id", "mRNAsi_raw", "mRNAsi", "mRNAsi_z", paste0("mRNAsi_z_", variant_names)))
score_dt <- score_meta[, score_dt_cols, with = FALSE]
fwrite(score_dt, out_scores, sep = "\t")
fwrite(score_meta, out_score_meta, sep = "\t")

summary_dt <- rbindlist(c(list(summary_dt), summary_variant_rows), fill = TRUE)
fwrite(summary_dt, out_summary, sep = "\t")

model_dt <- if (length(model_rows) > 0) rbindlist(model_rows, fill = TRUE) else data.table()
fwrite(model_dt, out_models, sep = "\t")
coef_dt <- if (length(coef_rows) > 0) rbindlist(coef_rows, fill = TRUE) else data.table()
if (nrow(coef_dt) > 0) {
  coef_dt[, term_label := term_label(term)]
  coef_dt[, fdr_bh := p.adjust(p, method = "BH"), by = .(score, model_variant)]
}
fwrite(coef_dt, out_coef_table, sep = "\t")
oneq_sens_dt <- copy(model_dt)
if (nrow(oneq_sens_dt) > 0) {
  oneq_sens_dt[, variant_order := match(model_variant, variant_names)]
  setorder(oneq_sens_dt, score, variant_order)
}
fwrite(oneq_sens_dt, out_oneq_sens_tsv, sep = "\t")

message(
  "Mode: ", mode, "\n",
  "Wrote computational outputs:\n",
  out_scores, "\n",
  out_score_meta, "\n",
  out_summary, "\n",
  out_models, "\n",
  out_coef_table, "\n",
  out_oneq_sens_tsv
)
