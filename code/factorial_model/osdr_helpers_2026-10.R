# Shared inputs and model helpers for the 2026-10 OSDR refits.
# Sourced by refit_gene_level_2026-10.R and followups_2026-10.R; expects RAW (the
# raw_osdr directory) to be defined before the loaders are called.

rnaseq_files <- c(
  "OSD-37"  = "GLDS-37_rna_seq_Normalized_Counts_GLbulkRNAseq.csv",
  "OSD-120" = "GLDS-120_rna_seq_Normalized_Counts_GLbulkRNAseq.csv",
  "OSD-217" = "GLDS-217_rna_seq_Normalized_Counts.csv",
  "OSD-321" = "GLDS-321_rna_seq_Normalized_Counts_GLbulkRNAseq.csv",
  "OSD-314" = "GLDS-314_rna_seq_Normalized_Counts_GLbulkRNAseq.csv",
  "OSD-678" = "GLDS-612_rna_seq_Normalized_Counts_GLbulkRNAseq.csv",
  "OSD-38"  = "GLDS-38_rna_seq_Normalized_Counts_GLbulkRNAseq.csv")
rnaseq_raw_files <- c(
  "OSD-37"  = "GLDS-37_rna_seq_RSEM_Unnormalized_Counts_GLbulkRNAseq.csv",
  "OSD-120" = "GLDS-120_rna_seq_RSEM_Unnormalized_Counts_GLbulkRNAseq.csv",
  "OSD-217" = "GLDS-217_rna_seq_Unnormalized_Counts.csv",
  "OSD-321" = "GLDS-321_rna_seq_RSEM_Unnormalized_Counts_GLbulkRNAseq.csv",
  "OSD-314" = "GLDS-314_rna_seq_RSEM_Unnormalized_Counts_GLbulkRNAseq.csv",
  "OSD-678" = "GLDS-612_rna_seq_RSEM_Unnormalized_Counts_GLbulkRNAseq.csv",
  "OSD-38"  = "GLDS-38_rna_seq_RSEM_Unnormalized_Counts_GLbulkRNAseq.csv")
array_files <- c(
  "OSD-7"   = "GLDS-7_array_normalized_expression_probeset_GLmicroarray.csv",
  "OSD-17"  = "GLDS-17_array_normalized_expression_probeset_GLmicroarray.csv",
  "OSD-147" = "GLDS-147_array_normalized_expression_probeset.csv",
  "OSD-205" = "GLDS-205_array_normalized_expression_probeset.csv")
ANNOT <- c("SYMBOL", "GENENAME", "REFSEQ", "ENTREZID", "STRING_id", "GOSLIM_IDS", "ProbesetID", "TAIR")

load_rnaseq <- function(osd, fn) {
  d <- fread(file.path(RAW, osd, fn), data.table = FALSE)
  rownames(d) <- trimws(as.character(d[[1]])); d[[1]] <- NULL
  m <- as.matrix(d); colnames(m) <- paste0(osd, "__", colnames(m)); m
}
load_array <- function(osd, fn) {
  d <- fread(file.path(RAW, osd, fn), data.table = FALSE)
  d <- d[!is.na(d$TAIR), ]
  d$TAIR <- trimws(as.character(d$TAIR))
  d <- d[grepl("^AT[1-5]G[0-9]{5}", d$TAIR), ]   # pandas str.match anchors at the start only
  # same rule as the Python builder: every non-annotation column is kept
  # (count_ENSEMBL_mappings included; it is dropped later by the sample filter)
  cols <- setdiff(names(d), ANNOT)
  m <- as.matrix(d[, cols]); storage.mode(m) <- "double"
  m <- rowsum(m, d$TAIR) / as.vector(table(d$TAIR))   # both sorted by TAIR: probeset mean
  colnames(m) <- paste0(osd, "__", colnames(m)); m
}
outer_merge <- function(mats) {
  genes <- sort(unique(unlist(lapply(mats, rownames))))
  out <- matrix(NA_real_, length(genes), sum(sapply(mats, ncol)),
                dimnames = list(genes, unlist(lapply(mats, colnames))))
  for (m in mats) out[rownames(m), colnames(m)] <- m
  out
}

# ---- Shared model code (as in fit_factorial_model.R) -------------------------
# CPM > 1 in >= 10 % of samples, as in fit_factorial_model.R
cpm_keep <- function(M) { d <- DGEList(counts = M); rowSums(cpm(d) > 1) >= 0.10 * ncol(d) }

prep_meta <- function(meta) {
  meta$flight <- factor(meta$flight, levels = c("GC", "FLT"))
  meta$light <- factor(meta$light, levels = c("light", "dark", "red_light", "unspecified"))
  meta$osd_id <- factor(meta$osd_id)
  meta$g_bl <- as.numeric(meta$g_bl_mol_m2_s)
  meta$carbon_12h <- as.numeric(meta$carbon_12h_pct_earth)
  meta$cfd_hardware <- factor(meta$cfd_hardware, levels = c("BRIC", "CARA", "VEGGIE"))
  os <- as.character(meta$organ); os[os %in% "unknown"] <- "whole_seedling"
  meta$organ_simple <- factor(os)
  ec <- as.character(meta$ecotype); ec[ec %in% c("Ler-0", "Cvi-0")] <- "other"
  meta$ecotype_simple <- factor(ec, levels = c("Col-0", "WS", "other"))
  rownames(meta) <- meta$full_sample_id
  meta
}
drop_aliased <- function(X) {
  q <- qr(X); if (q$rank < ncol(X)) X <- X[, -q$pivot[(q$rank + 1):ncol(X)], drop = FALSE]
  colnames(X) <- make.names(colnames(X)); X
}
fit_voom <- function(dge, X, block) {
  v <- voom(dge, X)
  cf <- duplicateCorrelation(v, X, block = block)
  # single correlation pass, exactly as fit_factorial_model.R, so that only the
  # covariates change between the stored fit and this refit
  v <- voom(dge, X, correlation = cf$consensus.correlation, block = block)
  f <- eBayes(lmFit(v, X, block = block, correlation = cf$consensus.correlation))
  list(fit = f, rho = cf$consensus.correlation)
}
fit_limma <- function(E, X, block) {
  cf <- duplicateCorrelation(E, X, block = block)
  f <- eBayes(lmFit(E, X, block = block, correlation = cf$consensus.correlation), trend = TRUE)
  list(fit = f, rho = cf$consensus.correlation)
}
tt <- function(fit, coef = NULL, contrast = NULL) {
  if (!is.null(contrast)) { fit <- eBayes(contrasts.fit(fit, contrast)); coef <- 1 }
  r <- topTable(fit, coef = coef, number = Inf, sort.by = "none")
  r$gene_id <- rownames(r); r
}
n_sig <- function(r) c(fdr05 = sum(r$adj.P.Val < 0.05),
                       fdr05_lfc1 = sum(r$adj.P.Val < 0.05 & abs(r$logFC) > 1))
F1 <- ~ flight * cfd_hardware + organ_simple + light + ecotype_simple
F2 <- ~ flight + g_bl + carbon_12h + organ_simple + light + ecotype_simple

# Flight averaged over BRIC/CARA/VEGGIE (treatment coding, reference BRIC)
avg_flight <- function(X) {
  c <- setNames(numeric(ncol(X)), colnames(X))
  c["flightFLT"] <- 1
  ints <- grep("^flightFLT\\.cfd_hardware", colnames(X), value = TRUE)
  c[ints] <- 1 / 3
  matrix(c, ncol = 1, dimnames = list(colnames(X), "Flight_avg_over_HW"))
}
