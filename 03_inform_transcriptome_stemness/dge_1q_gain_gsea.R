#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(key, default = NULL) {
  hit <- grep(paste0("^--", key, "="), args, value = TRUE)
  if (length(hit) == 0) return(default)
  sub(paste0("^--", key, "="), "", hit[1])
}

# --------------------------- Install/load packages ---------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

install_if_missing <- function(pkgs, bioc = FALSE) {
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      if (bioc) {
        BiocManager::install(pkg, ask = FALSE, update = FALSE)
      } else {
        install.packages(pkg, repos = "https://cloud.r-project.org")
      }
    }
  }
}

install_if_missing(c("edgeR", "limma", "fgsea", "DESeq2", "BiocParallel"), bioc = TRUE)
install_if_missing(c("msigdbr"), bioc = FALSE)

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(DESeq2)
  library(fgsea)
  library(msigdbr)
  library(BiocParallel)
})

# Force serial to avoid BiocParallel port errors
BiocParallel::register(BiocParallel::SerialParam())
options(BiocParallel.progressbar = FALSE)

# ------------------------------ User settings ------------------------------
root <- get_arg("root", Sys.getenv("INFORM_ALL_ROOT", unset = file.path("external", "inform_all")))
out_root <- get_arg("out_root", file.path(root, "dge_1q_gain_gsea"))
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

min_gsea_size <- 10
max_gsea_size <- 500
gsea_padj_cutoff <- 0.05

# ------------------------------ Helpers ------------------------------
stop_if_missing <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing)) stop("Missing required input(s):\n", paste(missing, collapse = "\n"))
}

norm_sentrix <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("\\s+", "", x)
  x
}

is_gain_call <- function(x) {
  x <- trimws(tolower(as.character(x)))
  out <- rep(NA, length(x))
  out[grepl("gain|amp|amplification", x)] <- TRUE
  out[grepl("clean|normal|neutral|no\\s*gain|none|wt", x)] <- FALSE
  out
}

is_loss_call <- function(x) {
  x <- trimws(tolower(as.character(x)))
  out <- rep(NA, length(x))
  out[grepl("loss|del|deletion", x)] <- TRUE
  out[grepl("clean|normal|neutral|no\\s*loss|none|wt", x)] <- FALSE
  out
}

any_cna_in_row <- function(x) {
  x <- trimws(tolower(as.character(x)))
  x <- x[!is.na(x)]
  if (length(x) == 0) return(FALSE)
  any(grepl("gain|amp|amplification", x)) || any(grepl("loss|del|deletion", x))
}

all_clean_in_row <- function(x) {
  x <- trimws(tolower(as.character(x)))
  ok <- grepl("clean|normal|neutral|no\\s*gain|no\\s*loss|none|wt", x)
  ok[is.na(ok)] <- FALSE
  all(ok)
}

any_or_na <- function(x) {
  if (all(is.na(x))) NA else any(x, na.rm = TRUE)
}

all_or_na <- function(x) {
  if (all(is.na(x))) NA else all(x, na.rm = TRUE)
}

write_table <- function(dt, path) {
  fwrite(dt, path, sep = "\t")
}

make_pathways <- function(msig_dt) {
  msig_dt <- msig_dt[!is.na(ensembl_gene) & ensembl_gene != ""]
  msig_dt[, ensembl_gene := sub("\\..*$", "", ensembl_gene)]
  split(msig_dt$ensembl_gene, msig_dt$gs_name)
}

prep_stats <- function(stat, gene_id) {
  dt <- data.table(gene = sub("\\..*$", "", gene_id), stat = stat)
  dt <- dt[!is.na(gene) & gene != "" & is.finite(stat)]
  dt <- dt[, .SD[which.max(abs(stat))], by = gene]
  stats <- dt$stat
  names(stats) <- dt$gene
  sort(stats, decreasing = TRUE)
}

run_fgsea <- function(stats, pathways, out_path) {
  if (length(stats) < 10) {
    message("Skipping GSEA (too few stats): ", out_path)
    return(invisible(NULL))
  }
  fg <- fgseaMultilevel(
    pathways = pathways,
    stats = stats,
    minSize = min_gsea_size,
    maxSize = max_gsea_size,
    eps = 1e-10,
    BPPARAM = BiocParallel::SerialParam()
  )
  if (nrow(fg) == 0) {
    message("GSEA returned 0 rows: ", out_path)
    return(invisible(NULL))
  }
  fg_dt <- as.data.table(fg)
  fg_dt[, leadingEdge := vapply(leadingEdge, function(x) paste(x, collapse = ";"), character(1))]
  setorder(fg_dt, padj, pval)
  full_path <- sub("\\.tsv$", "_full.tsv", out_path)
  write_table(fg_dt, full_path)

  fg_dt <- fg_dt[!is.na(padj) & padj < gsea_padj_cutoff]
  write_table(fg_dt, out_path)
  invisible(fg_dt)
}

run_voom_de <- function(counts, sample_info, design_formula, coef_name, out_prefix) {
  sample_info <- as.data.table(sample_info)
  counts <- counts[, sample_info$sample_id, drop = FALSE]

  design <- model.matrix(design_formula, data = as.data.frame(sample_info))
  y <- DGEList(counts = counts)
  keep <- filterByExpr(y, design = design)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)

  v <- voom(y, design, plot = FALSE)
  fit <- lmFit(v, design)
  fit <- eBayes(fit)

  if (!coef_name %in% colnames(design)) {
    stop("Coefficient not in design: ", coef_name)
  }
  tt <- topTable(fit, coef = coef_name, number = Inf, sort.by = "P")
  res_dt <- data.table(gene_id = rownames(tt))
  res_dt <- cbind(res_dt, as.data.table(tt))

  write_table(res_dt, paste0(out_prefix, "_DE.tsv"))
  return(list(result = res_dt, fit = fit, design = design))
}

run_deseq2_de <- function(counts, sample_info, design_formula, contrast_vec, out_prefix) {
  sample_info <- as.data.table(sample_info)
  counts <- counts[, sample_info$sample_id, drop = FALSE]

  # DESeq2 requires integer counts; round RSEM expected_count
  counts_int <- round(counts)
  counts_int[counts_int < 0] <- 0

  dds <- DESeqDataSetFromMatrix(
    countData = counts_int,
    colData = as.data.frame(sample_info),
    design = design_formula
  )

  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep, ]
  dds <- DESeq(dds)

  res <- results(dds, contrast = contrast_vec)
  res_dt <- data.table(gene_id = rownames(res), as.data.table(res))
  setorder(res_dt, padj)
  write_table(res_dt, paste0(out_prefix, "_DESeq2.tsv"))

  return(list(result = res_dt, dds = dds))
}

# ------------------------------ Inputs check ------------------------------
required <- c(
  file.path(root, "inform_rnaseq_sample_map.csv"),
  file.path(root, "manual_final_inform.csv"),
  file.path(root, "meth_metadata.csv"),
  file.path(root, "my_metadata.csv")
)
stop_if_missing(required)

# ------------------------------ Sample map ------------------------------
map_dt <- fread(file.path(root, "inform_rnaseq_sample_map.csv"))
map_dt <- map_dt[toupper(in_nfcore_results) == "TRUE"]
map_dt[, inform_id := gsub("-", "_", inform_id)]
map_dt <- unique(map_dt[, .(inform_id, rnaseq_sample_id)])

dup_samples <- map_dt[, .N, by = rnaseq_sample_id][N > 1]
if (nrow(dup_samples) > 0) {
  message("Dropping rnaseq_sample_id mapping to multiple inform_id: ", nrow(dup_samples))
  map_dt <- map_dt[!rnaseq_sample_id %in% dup_samples$rnaseq_sample_id]
}

# ------------------------------ Read gene counts ------------------------------
gene_files <- list.files(root, pattern = "\\.genes\\.results$", recursive = TRUE, full.names = TRUE)
if (length(gene_files) == 0) stop("No *.genes.results files found under root.")

gene_samples <- sub("\\.genes\\.results$", "", basename(gene_files))
if (anyDuplicated(gene_samples)) {
  dupes <- unique(gene_samples[duplicated(gene_samples)])
  message("Duplicate genes.results for sample_id(s): ", paste(dupes, collapse = ", "), " -> keeping first.")
  keep_idx <- !duplicated(gene_samples)
  gene_files <- gene_files[keep_idx]
  gene_samples <- gene_samples[keep_idx]
}

keep_samples <- unique(map_dt$rnaseq_sample_id)
keep_idx <- gene_samples %in% keep_samples
gene_files <- gene_files[keep_idx]
gene_samples <- gene_samples[keep_idx]
if (length(gene_files) < 5) stop("Too few genes.results files after sample-map filtering.")

first <- fread(gene_files[1], select = c("gene_id", "expected_count"))
genes <- first$gene_id
counts <- matrix(NA_real_, nrow = length(genes), ncol = length(gene_files))
colnames(counts) <- gene_samples
counts[, 1] <- first$expected_count

for (i in 2:length(gene_files)) {
  dt <- fread(gene_files[i], select = c("gene_id", "expected_count"))
  if (!identical(dt$gene_id, genes)) {
    dt <- dt[match(genes, dt$gene_id)]
  }
  counts[, i] <- dt$expected_count
}
rownames(counts) <- genes

# ------------------------------ Choose one sample per patient ------------------------------
lib_size <- colSums(counts, na.rm = TRUE)
sample_lib <- data.table(sample_id = colnames(counts), lib_size = lib_size)
sample_lib <- merge(sample_lib, map_dt, by.x = "sample_id", by.y = "rnaseq_sample_id", all.x = TRUE)
sample_keep <- sample_lib[order(inform_id, -lib_size, sample_id), .SD[1], by = inform_id]
keep_sample_ids <- unique(sample_keep$sample_id)
counts <- counts[, keep_sample_ids, drop = FALSE]

# ------------------------------ CNA mapping ------------------------------
meth <- fread(file.path(root, "meth_metadata.csv"))
setnames(meth, c("INFORM_ID", "Sentrix ID (methylation)"), c("inform_id", "sentrix"), skip_absent = TRUE)
meth[, inform_id := gsub("-", "_", inform_id)]
meth[, sentrix := norm_sentrix(sentrix)]

my <- fread(file.path(root, "my_metadata.csv"))
setnames(my, c("INFORM_ID", "Sentrix_ID"), c("inform_id", "sentrix"), skip_absent = TRUE)
my[, inform_id := gsub("-", "_", inform_id)]
my[, sentrix := norm_sentrix(sentrix)]

map_sentrix <- unique(rbindlist(
  list(meth[, .(inform_id, sentrix)], my[, .(inform_id, sentrix)]),
  use.names = TRUE, fill = TRUE
))
map_sentrix <- map_sentrix[!is.na(sentrix) & sentrix != ""]

cna <- fread(file.path(root, "manual_final_inform.csv"))
setnames(cna, "patient_ID", "sentrix", skip_absent = TRUE)
cna[, sentrix := norm_sentrix(sentrix)]

arm_cols <- setdiff(names(cna), "sentrix")
arm_cols <- arm_cols[grepl("^[0-9]+[pq]$", arm_cols)]
if (!"1q" %in% arm_cols) stop("manual_final_inform.csv missing '1q' column.")

# Drop CNA rows with any "inadequate plot"
cna[, inadequate_plot := apply(.SD, 1, function(x) any(grepl("inadequate plot", x, ignore.case = TRUE))),
    .SDcols = arm_cols]
n_inadequate <- sum(cna$inadequate_plot, na.rm = TRUE)
if (n_inadequate > 0) message("Dropping CNA rows with inadequate plot: ", n_inadequate)
cna <- cna[inadequate_plot == FALSE][, inadequate_plot := NULL]

# Calls and clean/other CNA flags
other_arms <- setdiff(arm_cols, "1q")

cna[, gain_1q := is_gain_call(`1q`)]
cna[, loss_1q := is_loss_call(`1q`)]
cna[, gain_8p := is_gain_call(`8p`)]
cna[, gain_12p := is_gain_call(`12p`)]
cna[, gain_12q := is_gain_call(`12q`)]
cna[, loss_16q := is_loss_call(`16q`)]

cna[, other_cna := apply(.SD, 1, any_cna_in_row), .SDcols = other_arms]
cna[, clean_all_arms := apply(.SD, 1, all_clean_in_row), .SDcols = arm_cols]

map_cna <- merge(
  map_sentrix,
  cna[, c("sentrix", "gain_1q", "loss_1q", "other_cna", "clean_all_arms",
          "gain_8p", "gain_12p", "gain_12q", "loss_16q"), with = FALSE],
  by = "sentrix", all.x = TRUE
)

cnv_by_inform <- map_cna[, .(
  gain_1q = any_or_na(gain_1q),
  loss_1q = any_or_na(loss_1q),
  other_cna = any_or_na(other_cna),
  clean_all_arms = all_or_na(clean_all_arms),
  gain_8p = any_or_na(gain_8p),
  gain_12p = any_or_na(gain_12p),
  gain_12q = any_or_na(gain_12q),
  loss_16q = any_or_na(loss_16q)
), by = inform_id]

loss_excl <- cnv_by_inform[loss_1q == TRUE, .N]
if (loss_excl > 0) message("Excluding patients with 1q loss: ", loss_excl)
cnv_by_inform <- cnv_by_inform[is.na(loss_1q) | loss_1q == FALSE]

# ------------------------------ MSigDB gene sets ------------------------------
message("Loading MSigDB Hallmark, Reactome, and GO gene sets...")
m_hallmark <- as.data.table(msigdbr(species = "Homo sapiens", collection = "H"))
m_reactome <- as.data.table(msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME"))
m_go_bp <- as.data.table(msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP"))
m_go_cc <- as.data.table(msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:CC"))
m_go_mf <- as.data.table(msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:MF"))

pathways_hallmark <- make_pathways(m_hallmark)
pathways_reactome <- make_pathways(m_reactome)
pathways_go_bp <- make_pathways(m_go_bp)
pathways_go_cc <- make_pathways(m_go_cc)
pathways_go_mf <- make_pathways(m_go_mf)

# ------------------------------ Analysis 1 ------------------------------
# 1q gain only (no other CNA) vs clean everywhere
out1 <- file.path(out_root, "analysis1_1q_gain_only_vs_clean")
dir.create(out1, recursive = TRUE, showWarnings = FALSE)

meta1 <- merge(sample_keep[, .(inform_id, sample_id)], cnv_by_inform, by = "inform_id", all.x = TRUE)
meta1 <- meta1[complete.cases(meta1[, .(gain_1q, other_cna, clean_all_arms)])]

meta1[, group := NA_character_]
meta1[gain_1q == TRUE & (is.na(loss_1q) | loss_1q == FALSE) & other_cna == FALSE, group := "1q_gain_only"]
meta1[clean_all_arms == TRUE, group := "clean"]
meta1 <- meta1[!is.na(group)]

meta1[, group := factor(group, levels = c("clean", "1q_gain_only"))]

write_table(meta1, file.path(out1, "sample_info.tsv"))

n_clean <- sum(meta1$group == "clean")
n_gain <- sum(meta1$group == "1q_gain_only")
if (n_clean < 2 || n_gain < 2) stop("Analysis 1: need >=2 samples per group. clean=", n_clean, " gain_only=", n_gain)

design1 <- ~ group
res1 <- run_voom_de(
  counts = counts,
  sample_info = meta1,
  design_formula = design1,
  coef_name = "group1q_gain_only",
  out_prefix = file.path(out1, "analysis1")
)

stats1 <- prep_stats(res1$result$t, res1$result$gene_id)
run_fgsea(stats1, pathways_hallmark, file.path(out1, "gsea_hallmark.tsv"))
run_fgsea(stats1, pathways_reactome, file.path(out1, "gsea_reactome.tsv"))

res1_deseq <- run_deseq2_de(
  counts = counts,
  sample_info = meta1,
  design_formula = design1,
  contrast_vec = c("group", "1q_gain_only", "clean"),
  out_prefix = file.path(out1, "analysis1")
)

stats1_ds <- prep_stats(res1_deseq$result$stat, res1_deseq$result$gene_id)
run_fgsea(stats1_ds, pathways_hallmark, file.path(out1, "gsea_hallmark_deseq2.tsv"))
run_fgsea(stats1_ds, pathways_reactome, file.path(out1, "gsea_reactome_deseq2.tsv"))

# ------------------------------ Analysis 2 ------------------------------
# 1q gain vs no 1q gain, covariates: 8p gain, 12p gain, 12q gain, 16q loss
out2 <- file.path(out_root, "analysis2_1q_gain_vs_no_gain_covariates")
dir.create(out2, recursive = TRUE, showWarnings = FALSE)

meta2 <- merge(sample_keep[, .(inform_id, sample_id)], cnv_by_inform, by = "inform_id", all.x = TRUE)
covars <- c("gain_8p", "gain_12p", "gain_12q", "loss_16q")
meta2 <- meta2[complete.cases(meta2[, c("gain_1q", covars), with = FALSE])]
meta2[, group := ifelse(gain_1q, "gain_1q", "no_gain_1q")]
meta2[, group := factor(group, levels = c("no_gain_1q", "gain_1q"))]

for (cv in covars) meta2[[cv]] <- as.integer(meta2[[cv]])

write_table(meta2, file.path(out2, "sample_info.tsv"))

n_no <- sum(meta2$group == "no_gain_1q")
n_yes <- sum(meta2$group == "gain_1q")
if (n_no < 2 || n_yes < 2) stop("Analysis 2: need >=2 samples per group. no_gain=", n_no, " gain=", n_yes)

design2 <- as.formula(paste("~ group +", paste(covars, collapse = " + ")))
res2 <- run_voom_de(
  counts = counts,
  sample_info = meta2,
  design_formula = design2,
  coef_name = "groupgain_1q",
  out_prefix = file.path(out2, "analysis2")
)

stats2 <- prep_stats(res2$result$t, res2$result$gene_id)
run_fgsea(stats2, pathways_hallmark, file.path(out2, "gsea_hallmark.tsv"))
run_fgsea(stats2, pathways_reactome, file.path(out2, "gsea_reactome.tsv"))
run_fgsea(stats2, pathways_go_bp, file.path(out2, "gsea_go_bp.tsv"))
run_fgsea(stats2, pathways_go_cc, file.path(out2, "gsea_go_cc.tsv"))
run_fgsea(stats2, pathways_go_mf, file.path(out2, "gsea_go_mf.tsv"))

res2_deseq <- run_deseq2_de(
  counts = counts,
  sample_info = meta2,
  design_formula = design2,
  contrast_vec = c("group", "gain_1q", "no_gain_1q"),
  out_prefix = file.path(out2, "analysis2")
)

stats2_ds <- prep_stats(res2_deseq$result$stat, res2_deseq$result$gene_id)
run_fgsea(stats2_ds, pathways_hallmark, file.path(out2, "gsea_hallmark_deseq2.tsv"))
run_fgsea(stats2_ds, pathways_reactome, file.path(out2, "gsea_reactome_deseq2.tsv"))
run_fgsea(stats2_ds, pathways_go_bp, file.path(out2, "gsea_go_bp_deseq2.tsv"))
run_fgsea(stats2_ds, pathways_go_cc, file.path(out2, "gsea_go_cc_deseq2.tsv"))
run_fgsea(stats2_ds, pathways_go_mf, file.path(out2, "gsea_go_mf_deseq2.tsv"))

message("Done. Outputs in: ", out_root)
