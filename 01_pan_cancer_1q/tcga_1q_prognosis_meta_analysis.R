#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)

local_r_lib <- Sys.getenv("LOCAL_R_LIB", unset = "")
if (nzchar(local_r_lib) && dir.exists(local_r_lib)) {
  .libPaths(c(local_r_lib, .libPaths()))
}

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(survival)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(key, default = NULL) {
  hit <- grep(paste0("^--", key, "="), args, value = TRUE)
  if (length(hit) == 0) return(default)
  sub(paste0("^--", key, "="), "", hit[1])
}

projects_arg <- get_arg("projects", "")
projects <- if (nzchar(projects_arg)) trimws(strsplit(projects_arg, ",")[[1]]) else character(0)
oneq_mode <- tolower(get_arg("oneq_mode", "binary_exclude_loss"))
if (!oneq_mode %in% c("binary", "binary_exclude_loss")) {
  stop("Unsupported --oneq_mode. Use 'binary' or 'binary_exclude_loss'.")
}
min_group_size <- as.integer(get_arg("min_group_size", "10"))
min_events <- as.integer(get_arg("min_events", "10"))
tumor_sample_types_arg <- get_arg("tumor_sample_types", "01,02,03,04,05,06,07,08")
tumor_sample_types <- trimws(strsplit(tumor_sample_types_arg, ",")[[1]])
include_oncokb_models <- tolower(get_arg("include_oncokb_models", "true")) %in% c("1", "true", "t", "yes", "y")
purity_metric_order <- trimws(strsplit(get_arg("purity_metric_order", "CPE,ABSOLUTE,ESTIMATE,IHC,LUMP"), ",")[[1]])

tcga_root <- get_arg("tcga_root", Sys.getenv("TCGA_ROOT", unset = file.path("external", "TCGA")))
reference_root <- get_arg("reference_root", Sys.getenv("TCGA_REFERENCE_ROOT", unset = file.path("external", "tcga_clinical")))
out_root <- get_arg("out_root", Sys.getenv("TCGA_1Q_PROGNOSIS_OUT", unset = file.path("results", "tcga_1q_prognosis")))
arm_calls_file <- get_arg("arm_calls_file", file.path(tcga_root, "PANCAN_ArmCallsAndAneuploidyScore_092817.txt"))
infiltration_file <- get_arg("infiltration_file", file.path(tcga_root, "infiltration_estimation_for_tcga.csv"))
survival_xlsx <- get_arg("survival_xlsx", file.path(reference_root, "TCGA-CDR-SupplementalTableS1.xlsx"))
purity_xlsx <- get_arg("purity_xlsx", file.path(reference_root, "Aran_TCGA_purity_Supplementary_Data1.xlsx"))
mutation_root <- get_arg("mutation_root", Sys.getenv("TCGA_MUTATION_ROOT", unset = file.path("external", "tcga_mutation_lists")))
mutation_manifest_file <- get_arg("mutation_manifest", file.path(mutation_root, "tcga_mutation_list_manifest.tsv"))
mc3_cache_file <- get_arg("mc3_cache", Sys.getenv("TCGA_MC3_CACHE", unset = file.path("external", "tcga_stemness_scores", "mc3_mutation_covariates_fixed.tsv")))

dir.create(reference_root, recursive = TRUE, showWarnings = FALSE)
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

out_survival_clean <- file.path(reference_root, "tcga_cdr_os.tsv")
out_purity_clean <- file.path(reference_root, "tcga_aran_purity.tsv")
out_cohort <- file.path(out_root, "panTCGA_1q_patient_covariates.tsv.gz")
out_per_project <- file.path(out_root, "panTCGA_1q_prognosis_per_project.tsv")
out_meta <- file.path(out_root, "panTCGA_1q_prognosis_meta_summary.tsv")
out_pooled <- file.path(out_root, "panTCGA_1q_prognosis_pooled_models.tsv")

stop_if_missing <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0) {
    stop("Missing required input(s):\n", paste(missing, collapse = "\n"))
  }
}

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

parse_sample_type <- function(x) {
  sub("^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-([0-9]{2}).*$", "\\1", toupper(as.character(x)))
}

find_id_column <- function(dt) {
  candidates <- c("sample_id", "sample", "Sample", "SampleID", "sampleID",
                  "cell_type", "aliquot_barcode", "barcode")
  hit <- candidates[candidates %in% names(dt)]
  if (length(hit) > 0) return(hit[1])
  NA_character_
}

find_aneuploidy_col <- function(dt) {
  candidates <- c("Aneuploidy Score", "AneuploidyScore", "aneuploidy_score",
                  "aneuploidy score", "Aneuploidy", "aneuploidy")
  hit <- candidates[candidates %in% names(dt)]
  if (length(hit) > 0) return(hit[1])
  idx <- grep("aneuploidy", names(dt), ignore.case = TRUE, value = TRUE)
  if (length(idx) > 0) return(idx[1])
  NA_character_
}

safe_numeric <- function(x) {
  x <- as.character(x)
  x[x %in% c("NaN", "nan", "NA", "", "NULL")] <- NA_character_
  suppressWarnings(as.numeric(x))
}

mean_or_na <- function(x) {
  x <- as.numeric(x)
  if (all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

zscore <- function(x) {
  x <- as.numeric(x)
  if (all(is.na(x))) return(rep(NA_real_, length(x)))
  s <- stats::sd(x, na.rm = TRUE)
  m <- mean(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) {
    out <- rep(0, length(x))
    out[is.na(x)] <- NA_real_
    return(out)
  }
  (x - m) / s
}

meta_dl <- function(beta, se) {
  w <- 1 / (se ^ 2)
  mu_fixed <- sum(w * beta) / sum(w)
  q <- sum(w * (beta - mu_fixed) ^ 2)
  k <- length(beta)
  c_val <- sum(w) - sum(w ^ 2) / sum(w)
  tau2 <- if (k > 1) max(0, (q - (k - 1)) / c_val) else 0
  w_re <- 1 / (se ^ 2 + tau2)
  mu_re <- sum(w_re * beta) / sum(w_re)
  se_re <- sqrt(1 / sum(w_re))
  z <- mu_re / se_re
  p <- 2 * pnorm(-abs(z))
  i2 <- if (k > 1 && q > (k - 1)) max(0, (q - (k - 1)) / q) else 0
  data.table(
    beta = mu_re,
    se = se_re,
    hr = exp(mu_re),
    hr_low = exp(mu_re - 1.96 * se_re),
    hr_high = exp(mu_re + 1.96 * se_re),
    p_value = p,
    tau2 = tau2,
    Q = q,
    I2 = i2,
    k = k
  )
}

extract_gain_term <- function(fit, term = "gain_1q") {
  sm <- summary(fit)
  if (!term %in% rownames(sm$coefficients)) return(NULL)
  data.table(
    beta = unname(sm$coefficients[term, "coef"]),
    se = unname(sm$coefficients[term, "se(coef)"]),
    hr = unname(sm$coefficients[term, "exp(coef)"]),
    hr_low = unname(sm$conf.int[term, "lower .95"]),
    hr_high = unname(sm$conf.int[term, "upper .95"]),
    p_value = unname(sm$coefficients[term, "Pr(>|z|)"])
  )
}

prepare_survival_table <- function(path, out_tsv) {
  stop_if_missing(path)
  raw <- as.data.table(read_excel(path))
  out <- raw[, .(
    patient_id = normalize_tcga_patient(bcr_patient_barcode),
    cancer_type = as.character(type),
    os_event = safe_numeric(OS),
    os_time_days = safe_numeric(`OS.time`)
  )]
  out <- out[!is.na(patient_id) & !is.na(cancer_type)]
  out[, os_event := as.integer(os_event)]
  out <- out[!is.na(os_event) & !is.na(os_time_days) & os_time_days > 0]
  setorder(out, patient_id, -os_event, -os_time_days)
  out <- unique(out, by = "patient_id")
  fwrite(out, out_tsv, sep = "\t")
  out
}

prepare_purity_table <- function(path, out_tsv, metric_order) {
  stop_if_missing(path)
  raw <- as.data.table(read_excel(path, sheet = "Supp Data 1", skip = 3))
  drop_cols <- names(raw)[grepl("^\\.\\.\\.", names(raw))]
  if (length(drop_cols) > 0) raw[, (drop_cols) := NULL]
  for (col in intersect(metric_order, names(raw))) {
    raw[, (col) := safe_numeric(get(col))]
  }
  raw[, sample_id := normalize_tcga_sample(`Sample ID`)]
  raw[, cancer_type := as.character(`Cancer type`)]
  raw[, purity := NA_real_]
  raw[, purity_source := NA_character_]
  for (metric in metric_order) {
    if (!metric %in% names(raw)) next
    idx <- is.na(raw$purity) & !is.na(raw[[metric]])
    raw[idx, purity := raw[[metric]][idx]]
    raw[idx, purity_source := metric]
  }
  out <- raw[!is.na(sample_id), .(
    purity = mean_or_na(purity),
    purity_source = {
      src <- purity_source[!is.na(purity_source) & purity_source != ""]
      if (length(src) == 0) NA_character_ else src[1]
    }
  ), by = .(sample_id, cancer_type)]
  fwrite(out, out_tsv, sep = "\t")
  out
}

load_immune_stromal <- function(path) {
  stop_if_missing(path)
  infl <- fread(path)
  id_col <- find_id_column(infl)
  if (is.na(id_col)) stop("Could not find sample identifier column in infiltration file.")
  if (!all(c("immune score_XCELL", "stroma score_XCELL") %in% names(infl))) {
    stop("Infiltration file is missing xCell immune/stromal scores.")
  }
  out <- infl[, .(
    sample_id = normalize_tcga_sample(get(id_col)),
    ImmuneScore = safe_numeric(`immune score_XCELL`),
    StromaScore = safe_numeric(`stroma score_XCELL`)
  )]
  out <- out[!is.na(sample_id)]
  out <- out[, .(
    ImmuneScore = mean_or_na(ImmuneScore),
    StromaScore = mean_or_na(StromaScore)
  ), by = sample_id]
  out
}

load_arm_calls <- function(path, tumor_sample_types) {
  stop_if_missing(path)
  arm <- fread(path, sep = "\t", header = TRUE, check.names = FALSE)
  setnames(arm, c("Sample", "Type"), c("sample_id", "cancer_type"), skip_absent = TRUE)
  if (!all(c("sample_id", "cancer_type", "1q") %in% names(arm))) {
    stop("Arm-calls file is missing required columns.")
  }
  aneu_col <- find_aneuploidy_col(arm)
  arm[, sample_id := normalize_tcga_sample(sample_id)]
  arm[, patient_id := normalize_tcga_patient(sample_id)]
  arm[, sample_type := parse_sample_type(sample_id)]
  arm <- arm[!is.na(sample_id) & sample_type %in% tumor_sample_types]
  arm[, oneq_raw := safe_numeric(`1q`)]
  arm <- arm[!is.na(oneq_raw)]
  arm[, oneq_state := fifelse(oneq_raw == 1, "gain", fifelse(oneq_raw == -1, "loss", "neutral"))]
  arm[, gain_1q := as.integer(oneq_state == "gain")]
  if (!is.na(aneu_col)) {
    arm[, aneuploidy_score := safe_numeric(get(aneu_col))]
  } else {
    arm[, aneuploidy_score := NA_real_]
  }
  arm[, project := paste0("TCGA-", cancer_type)]
  arm[, .(project, cancer_type, sample_id, patient_id, sample_type, oneq_state, gain_1q, aneuploidy_score)]
}

load_mc3_cache <- function(path) {
  if (!nzchar(path) || !file.exists(path)) return(NULL)
  dt <- fread(path)
  if (!"sample_id" %in% names(dt)) return(NULL)
  dt[, sample_id := normalize_tcga_sample(sample_id)]
  dt <- dt[!is.na(sample_id)]
  keep <- intersect(c("sample_id", "tmb_nonsyn", "tmb_log1p"), names(dt))
  dt <- dt[, ..keep]
  for (col in setdiff(names(dt), "sample_id")) dt[, (col) := safe_numeric(get(col))]
  unique(dt, by = "sample_id")
}

load_mutation_summaries <- function(root, manifest_file = "") {
  if (!nzchar(root) || !dir.exists(root)) return(NULL)
  manifest <- NULL
  if (nzchar(manifest_file) && file.exists(manifest_file)) {
    manifest <- fread(manifest_file)
    if ("oncokb_annotated" %in% names(manifest)) {
      manifest[, oncokb_annotated := tolower(as.character(oncokb_annotated)) %in% c("true", "1", "yes")]
    }
  }
  sum_files <- list.files(root, pattern = "sample_mutation_summary\\.tsv$", recursive = TRUE, full.names = TRUE)
  if (length(sum_files) == 0) return(NULL)
  out_list <- vector("list", length(sum_files))
  for (i in seq_along(sum_files)) {
    f <- sum_files[i]
    project <- basename(dirname(f))
    dt <- fread(f)
    if (!"sample_id" %in% names(dt)) next
    keep <- intersect(
      c(
        "sample_id",
        "patient_id",
        "n_all_mutations",
        "n_nonsyn_mutations",
        "n_oncokb_annotated_mutations",
        "n_likely_driver_mutations",
        "n_likely_clinically_relevant_mutations",
        "n_therapeutically_actionable_mutations"
      ),
      names(dt)
    )
    dt <- dt[, ..keep]
    dt[, sample_id := normalize_tcga_sample(sample_id)]
    if ("patient_id" %in% names(dt)) {
      dt[, patient_id := normalize_tcga_patient(patient_id)]
    } else {
      dt[, patient_id := normalize_tcga_patient(sample_id)]
    }
    dt[, project := project]
    num_cols <- setdiff(names(dt), c("sample_id", "patient_id", "project"))
    for (col in num_cols) dt[, (col) := safe_numeric(get(col))]
    if (!is.null(manifest) && project %in% manifest$project && "oncokb_annotated" %in% names(manifest)) {
      ready <- manifest[["oncokb_annotated"]][match(project, manifest[["project"]])]
      dt[, oncokb_project_ready := as.logical(ready)]
    } else {
      dt[, oncokb_project_ready := any(n_oncokb_annotated_mutations > 0, na.rm = TRUE)]
    }
    dt[, driver_log1p := ifelse(oncokb_project_ready, log1p(pmax(0, n_likely_driver_mutations)), NA_real_)]
    out_list[[i]] <- dt
  }
  out <- rbindlist(out_list, use.names = TRUE, fill = TRUE)
  if (nrow(out) == 0) return(NULL)
  unique(out, by = c("project", "sample_id"))
}

select_one_sample_per_patient <- function(dt, tumor_sample_types) {
  type_rank <- setNames(seq_along(tumor_sample_types), tumor_sample_types)
  covar_cols <- intersect(
    c("purity", "ImmuneScore", "StromaScore", "aneuploidy_score", "tmb_log1p", "driver_log1p"),
    names(dt)
  )
  out <- copy(dt)
  out[, sample_type_rank := type_rank[sample_type]]
  out[is.na(sample_type_rank), sample_type_rank := length(tumor_sample_types) + 1L]
  if (length(covar_cols) > 0) {
    out[, covariate_nonmissing := rowSums(!is.na(.SD)), .SDcols = covar_cols]
  } else {
    out[, covariate_nonmissing := 0L]
  }
  out[, purity_present := as.integer(!is.na(purity))]
  setorder(out, project, patient_id, sample_type_rank, -covariate_nonmissing, -purity_present, -purity, sample_id)
  out <- out[!duplicated(out[, .(project, patient_id)])]
  out[, c("sample_type_rank", "covariate_nonmissing", "purity_present") := NULL]
  out
}

fit_project_model <- function(dt, model_name, covars, oncokb_only = FALSE) {
  res <- vector("list", length(unique(dt$project)))
  proj_list <- sort(unique(dt$project))
  for (i in seq_along(proj_list)) {
    proj_id <- proj_list[i]
    sub <- copy(dt[get("project") == proj_id])
    if (oneq_mode == "binary_exclude_loss") sub <- sub[oneq_state != "loss"]
    if (oncokb_only) sub <- sub[oncokb_project_ready == TRUE]
    status <- "ok"
    reason <- NA_character_
    beta <- se <- hr <- hr_low <- hr_high <- p_value <- NA_real_
    covars_used <- character(0)
    covars_dropped <- character(0)
    n_patients <- nrow(sub)
    events <- sum(sub$os_event == 1, na.rm = TRUE)
    n_gain <- sum(sub$gain_1q == 1, na.rm = TRUE)
    n_nogain <- sum(sub$gain_1q == 0, na.rm = TRUE)
    if (nrow(sub) == 0) {
      status <- "skipped"
      reason <- "no_samples_after_filters"
    } else {
      needed <- c("os_time_days", "os_event", "gain_1q", covars)
      sub <- sub[complete.cases(sub[, ..needed])]
      n_patients <- nrow(sub)
      events <- sum(sub$os_event == 1, na.rm = TRUE)
      n_gain <- sum(sub$gain_1q == 1, na.rm = TRUE)
      n_nogain <- sum(sub$gain_1q == 0, na.rm = TRUE)
      if (nrow(sub) == 0) {
        status <- "skipped"
        reason <- "no_complete_cases"
      } else if (n_gain < min_group_size || n_nogain < min_group_size) {
        status <- "skipped"
        reason <- "group_too_small"
      } else if (events < min_events) {
        status <- "skipped"
        reason <- "too_few_events"
      } else {
        for (col in covars) {
          if (uniqueN(sub[[col]]) <= 1) {
            covars_dropped <- c(covars_dropped, col)
          } else {
            z_col <- paste0(col, "_z")
            sub[, (z_col) := zscore(get(col))]
            covars_used <- c(covars_used, z_col)
          }
        }
        rhs <- c("gain_1q", covars_used)
        fit <- tryCatch(
          coxph(as.formula(paste("Surv(os_time_days, os_event) ~", paste(rhs, collapse = " + "))),
                data = sub, ties = "efron"),
          error = function(e) e
        )
        if (inherits(fit, "error")) {
          status <- "failed"
          reason <- conditionMessage(fit)
        } else {
          est <- extract_gain_term(fit, "gain_1q")
          if (is.null(est)) {
            status <- "failed"
            reason <- "gain_1q_term_missing"
          } else {
            beta <- est$beta
            se <- est$se
            hr <- est$hr
            hr_low <- est$hr_low
            hr_high <- est$hr_high
            p_value <- est$p_value
          }
        }
      }
    }
    res[[i]] <- data.table(
      project = proj_id,
      cancer_type = sub("^TCGA-", "", proj_id),
      model = model_name,
      oncokb_only = oncokb_only,
      oneq_mode = oneq_mode,
      status = status,
      reason = reason,
      n_patients = n_patients,
      n_events = events,
      n_gain = n_gain,
      n_nogain = n_nogain,
      beta_gain_1q = beta,
      se_gain_1q = se,
      hr_gain_1q = hr,
      hr_low_95 = hr_low,
      hr_high_95 = hr_high,
      p_value = p_value,
      covariates_requested = paste(covars, collapse = ","),
      covariates_used = paste(sub("_z$", "", covars_used), collapse = ","),
      covariates_dropped = paste(covars_dropped, collapse = ",")
    )
  }
  rbindlist(res, use.names = TRUE, fill = TRUE)
}

fit_pooled_model <- function(dt, model_name, covars, oncokb_only = FALSE) {
  sub <- copy(dt)
  if (oneq_mode == "binary_exclude_loss") sub <- sub[oneq_state != "loss"]
  if (oncokb_only) sub <- sub[oncokb_project_ready == TRUE]
  status <- "ok"
  reason <- NA_character_
  needed <- c("project", "cancer_type", "os_time_days", "os_event", "gain_1q", covars)
  sub <- sub[complete.cases(sub[, ..needed])]
  n_projects <- uniqueN(sub$project)
  n_patients <- nrow(sub)
  n_events <- sum(sub$os_event == 1, na.rm = TRUE)
  n_gain <- sum(sub$gain_1q == 1, na.rm = TRUE)
  n_nogain <- sum(sub$gain_1q == 0, na.rm = TRUE)
  beta <- se <- hr <- hr_low <- hr_high <- p_value <- NA_real_
  covars_used <- character(0)
  covars_dropped <- character(0)
  if (nrow(sub) == 0) {
    status <- "skipped"
    reason <- "no_complete_cases"
  } else if (n_projects < 2) {
    status <- "skipped"
    reason <- "too_few_projects"
  } else if (n_gain < min_group_size || n_nogain < min_group_size) {
    status <- "skipped"
    reason <- "group_too_small"
  } else if (n_events < min_events) {
    status <- "skipped"
    reason <- "too_few_events"
  } else {
    for (col in covars) {
      if (uniqueN(sub[[col]]) <= 1) {
        covars_dropped <- c(covars_dropped, col)
      } else {
        z_col <- paste0(col, "_z")
        sub[, (z_col) := zscore(get(col)), by = cancer_type]
        covars_used <- c(covars_used, z_col)
      }
    }
    rhs <- c("gain_1q", covars_used, "strata(cancer_type)")
    fit <- tryCatch(
      coxph(as.formula(paste("Surv(os_time_days, os_event) ~", paste(rhs, collapse = " + "))),
            data = sub, ties = "efron"),
      error = function(e) e
    )
    if (inherits(fit, "error")) {
      status <- "failed"
      reason <- conditionMessage(fit)
    } else {
      est <- extract_gain_term(fit, "gain_1q")
      if (is.null(est)) {
        status <- "failed"
        reason <- "gain_1q_term_missing"
      } else {
        beta <- est$beta
        se <- est$se
        hr <- est$hr
        hr_low <- est$hr_low
        hr_high <- est$hr_high
        p_value <- est$p_value
      }
    }
  }
  data.table(
    model = model_name,
    oncokb_only = oncokb_only,
    oneq_mode = oneq_mode,
    status = status,
    reason = reason,
    n_projects = n_projects,
    n_patients = n_patients,
    n_events = n_events,
    n_gain = n_gain,
    n_nogain = n_nogain,
    beta_gain_1q = beta,
    se_gain_1q = se,
    hr_gain_1q = hr,
    hr_low_95 = hr_low,
    hr_high_95 = hr_high,
    p_value = p_value,
    covariates_requested = paste(covars, collapse = ","),
    covariates_used = paste(sub("_z$", "", covars_used), collapse = ","),
    covariates_dropped = paste(covars_dropped, collapse = ",")
  )
}

stop_if_missing(c(arm_calls_file, infiltration_file, survival_xlsx, purity_xlsx, mc3_cache_file))

surv_dt <- prepare_survival_table(survival_xlsx, out_survival_clean)
purity_dt <- prepare_purity_table(purity_xlsx, out_purity_clean, purity_metric_order)
immune_dt <- load_immune_stromal(infiltration_file)
arm_dt <- load_arm_calls(arm_calls_file, tumor_sample_types)
mc3_dt <- load_mc3_cache(mc3_cache_file)
mut_dt <- load_mutation_summaries(mutation_root, mutation_manifest_file)

meta_dt <- copy(arm_dt)
meta_dt <- merge(meta_dt, immune_dt, by = "sample_id", all.x = TRUE)
meta_dt <- merge(meta_dt, purity_dt[, .(sample_id, purity, purity_source)], by = "sample_id", all.x = TRUE)
if (!is.null(mc3_dt)) {
  meta_dt <- merge(meta_dt, mc3_dt, by = "sample_id", all.x = TRUE)
}
if (!is.null(mut_dt)) {
  keep_mut <- intersect(
    c(
      "project", "sample_id", "n_all_mutations", "n_nonsyn_mutations",
      "n_oncokb_annotated_mutations", "n_likely_driver_mutations",
      "n_likely_clinically_relevant_mutations", "n_therapeutically_actionable_mutations",
      "driver_log1p", "oncokb_project_ready"
    ),
    names(mut_dt)
  )
  meta_dt <- merge(meta_dt, mut_dt[, ..keep_mut], by = c("project", "sample_id"), all.x = TRUE)
}
if (!"tmb_log1p" %in% names(meta_dt)) meta_dt[, tmb_log1p := NA_real_]
if (!"tmb_nonsyn" %in% names(meta_dt)) meta_dt[, tmb_nonsyn := NA_real_]
if (!"n_nonsyn_mutations" %in% names(meta_dt)) meta_dt[, n_nonsyn_mutations := NA_real_]
if (!"driver_log1p" %in% names(meta_dt)) meta_dt[, driver_log1p := NA_real_]
if (!"oncokb_project_ready" %in% names(meta_dt)) meta_dt[, oncokb_project_ready := FALSE]
meta_dt[is.na(oncokb_project_ready), oncokb_project_ready := FALSE]
meta_dt[is.na(tmb_log1p) & !is.na(n_nonsyn_mutations), tmb_log1p := log1p(pmax(0, n_nonsyn_mutations))]
meta_dt[is.na(tmb_nonsyn) & !is.na(n_nonsyn_mutations), tmb_nonsyn := pmax(0, n_nonsyn_mutations)]

meta_dt <- select_one_sample_per_patient(meta_dt, tumor_sample_types)
meta_dt <- merge(meta_dt, surv_dt, by = c("patient_id", "cancer_type"), all.x = TRUE)
meta_dt <- meta_dt[!is.na(os_time_days) & !is.na(os_event)]

if (length(projects) > 0) {
  meta_dt <- meta_dt[project %in% projects]
}
if (nrow(meta_dt) == 0) {
  stop("No samples remained after filtering and survival merge.")
}

setorder(meta_dt, project, patient_id)
fwrite(meta_dt, out_cohort, sep = "\t")

model_specs <- list(
  list(name = "core_no_purity", covars = c("tmb_log1p", "ImmuneScore", "StromaScore", "aneuploidy_score"), oncokb_only = FALSE),
  list(name = "core_with_purity", covars = c("tmb_log1p", "purity", "ImmuneScore", "StromaScore", "aneuploidy_score"), oncokb_only = FALSE)
)
if (include_oncokb_models) {
  model_specs <- c(
    model_specs,
    list(
      list(name = "oncokb_driver_no_purity", covars = c("tmb_log1p", "driver_log1p", "ImmuneScore", "StromaScore", "aneuploidy_score"), oncokb_only = TRUE),
      list(name = "oncokb_driver_with_purity", covars = c("tmb_log1p", "driver_log1p", "purity", "ImmuneScore", "StromaScore", "aneuploidy_score"), oncokb_only = TRUE)
    )
  )
}

per_project_res <- rbindlist(lapply(model_specs, function(spec) {
  fit_project_model(meta_dt, spec$name, spec$covars, spec$oncokb_only)
}), use.names = TRUE, fill = TRUE)
per_project_res[, fdr_gain_1q := NA_real_]
per_project_res[status == "ok" & is.finite(p_value),
                fdr_gain_1q := p.adjust(p_value, method = "BH"),
                by = model]

meta_res <- rbindlist(lapply(model_specs, function(spec) {
  sub <- per_project_res[model == spec$name & status == "ok" & is.finite(beta_gain_1q) & is.finite(se_gain_1q)]
  if (nrow(sub) < 2) {
    data.table(
      model = spec$name,
      oncokb_only = spec$oncokb_only,
      oneq_mode = oneq_mode,
      status = "skipped",
      reason = "too_few_projects_for_meta",
      k = nrow(sub),
      beta_gain_1q = NA_real_,
      se_gain_1q = NA_real_,
      hr_gain_1q = NA_real_,
      hr_low_95 = NA_real_,
      hr_high_95 = NA_real_,
      p_value = NA_real_,
      tau2 = NA_real_,
      Q = NA_real_,
      I2 = NA_real_
    )
  } else {
    ma <- meta_dl(sub$beta_gain_1q, sub$se_gain_1q)
    data.table(
      model = spec$name,
      oncokb_only = spec$oncokb_only,
      oneq_mode = oneq_mode,
      status = "ok",
      reason = NA_character_,
      k = ma$k,
      beta_gain_1q = ma$beta,
      se_gain_1q = ma$se,
      hr_gain_1q = ma$hr,
      hr_low_95 = ma$hr_low,
      hr_high_95 = ma$hr_high,
      p_value = ma$p_value,
      tau2 = ma$tau2,
      Q = ma$Q,
      I2 = ma$I2
    )
  }
}), use.names = TRUE, fill = TRUE)

pooled_res <- rbindlist(lapply(model_specs, function(spec) {
  fit_pooled_model(meta_dt, spec$name, spec$covars, spec$oncokb_only)
}), use.names = TRUE, fill = TRUE)

fwrite(per_project_res, out_per_project, sep = "\t")
fwrite(meta_res, out_meta, sep = "\t")
fwrite(pooled_res, out_pooled, sep = "\t")

message("Wrote cohort table: ", out_cohort)
message("Wrote per-project Cox results: ", out_per_project)
message("Wrote meta-analysis summary: ", out_meta)
message("Wrote pooled stratified Cox summary: ", out_pooled)
