#!/usr/bin/env Rscript
# =============================================================================
# Gene-level limma refit with the corrected CFD covariates (2026-10)
# =============================================================================
# 1. Rebuilds the harmonized matrix from the 11 OSDR processed files
#    (same inputs and rules as code/osdr_miner/build_harmonized_matrix.py).
# 2. Reconstructs the missing "harmonized_expression_matrix_filtered.tsv" step
#    (no script in the package produces it): samples restricted to the 307 in
#    the metadata, in metadata order; two gene rules are tried and the one that
#    reproduces the stored fit (25,453 genes, inter-study rho 0.7324) is kept.
# 3. Reproduction check: the original Model 1/2 with the original covariates,
#    compared with the stored tables/de_results_*.tsv.
# 4. Refit with tables/reanalysis_2026-10/expression_metadata_with_cfd_v2.tsv:
#      - Model 1 adds the flight effect averaged over hardware (the stored
#        "Flight_FLT_vs_GC" is the flight effect within BRIC only).
#      - Model 2 with the v2 g_bl and carbon_12h.
#    Plus a platform-separated sensitivity: RNA-seq studies with voom on counts,
#    microarray studies with limma on log-intensities. The original pipeline
#    passes microarray log-intensities through DGEList/voom as if they were
#    counts, which voom's mean-variance model does not support.
#
# Usage (paths are env-overridable):
#   OSDR_RAW=~/Documents/Airflow_omics_osdr_data/raw_osdr \
#   OSDR_DERIVED=~/Documents/Airflow_omics_osdr_data/derived \
#   Rscript code/factorial_model/refit_gene_level_2026-10.R
# Run from the repo root. Writes tables/reanalysis_2026-10/gene_level/.
# =============================================================================

suppressPackageStartupMessages({
  library(limma)
  library(edgeR)
  library(data.table)
})

RAW <- path.expand(Sys.getenv("OSDR_RAW", "~/Documents/Airflow_omics_osdr_data/raw_osdr"))
DERIVED <- path.expand(Sys.getenv("OSDR_DERIVED", "~/Documents/Airflow_omics_osdr_data/derived"))
OUT <- "tables/reanalysis_2026-10/gene_level"
dir.create(DERIVED, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
log_lines <- character()
say <- function(...) { s <- paste0(...); cat(s, "\n"); log_lines <<- c(log_lines, s) }

# ---- 1. Harmonized matrix ---------------------------------------------------
source("code/factorial_model/osdr_helpers_2026-10.R")

say("== 1. Harmonized matrix ==")
mats <- c(mapply(load_rnaseq, names(rnaseq_files), rnaseq_files, SIMPLIFY = FALSE),
          mapply(load_array, names(array_files), array_files, SIMPLIFY = FALSE))
for (n in names(mats)) say(sprintf("  %-8s %6d genes x %3d cols", n, nrow(mats[[n]]), ncol(mats[[n]])))
H <- outer_merge(mats)
say(sprintf("  merged (outer): %d genes x %d cols", nrow(H), ncol(H)))

meta1 <- fread("tables/expression_metadata_with_cfd.tsv", data.table = FALSE)
meta2 <- fread("tables/reanalysis_2026-10/expression_metadata_with_cfd_v2.tsv", data.table = FALSE)
stopifnot(identical(meta1$full_sample_id, meta2$full_sample_id))
stopifnot(all(meta1$full_sample_id %in% colnames(H)))
H <- H[, meta1$full_sample_id]
say(sprintf("  restricted to metadata samples: %d", ncol(H)))

# ---- 2. Reconstruct the filtered matrix --------------------------------------
variants <- list(
  complete_cases = H[rowSums(is.na(H)) == 0, ],
  na_to_zero = { Z <- H; Z[is.na(Z)] <- 0; Z })
for (n in names(variants)) say(sprintf("  variant %-15s %6d genes, %6d after CPM filter",
                                       n, nrow(variants[[n]]), sum(cpm_keep(variants[[n]]))))

# ---- 3. Reproduction check ----------------------------------------------------
say("\n== 3. Reproduction check (original covariates) ==")
m1 <- prep_meta(meta1)
stored <- function(nm) fread(file.path("tables", paste0("de_results_", nm, ".tsv")), data.table = FALSE)
pick <- NULL
repro_m1 <- list()
for (n in names(variants)) {
  M <- variants[[n]]; keep <- cpm_keep(M)
  dge <- calcNormFactors(DGEList(counts = M)[keep, , keep.lib.sizes = FALSE], method = "TMM")
  X1 <- drop_aliased(model.matrix(F1, m1))
  r1 <- fit_voom(dge, X1, m1$osd_id)
  s <- stored("Flight_FLT_vs_GC"); r <- tt(r1$fit, "flightFLT")
  g <- intersect(s$gene_id, r$gene_id)
  rt <- cor(s$t[match(g, s$gene_id)], r$t[match(g, r$gene_id)])
  say(sprintf("  %-15s genes=%d (stored 25453)  rho=%.4f (stored 0.7324)  t-corr vs stored=%.4f (n=%d)",
              n, nrow(dge), r1$rho, rt, length(g)))
  if (is.null(pick) && nrow(dge) == 25453) pick <- n
  repro_m1[[n]] <- c(genes = nrow(dge), rho = r1$rho, t_cor = rt)
}
if (is.null(pick)) {
  pick <- "complete_cases"
  say("  WARNING: neither variant reproduces 25,453 genes; continuing with complete_cases")
}
say("  using variant: ", pick)
M <- variants[[pick]]
fwrite(data.table(gene_id = rownames(M), M), file.path(DERIVED, "harmonized_expression_matrix_filtered.tsv"), sep = "\t")
keep <- cpm_keep(M)
dge <- calcNormFactors(DGEList(counts = M)[keep, , keep.lib.sizes = FALSE], method = "TMM")

X2o <- drop_aliased(model.matrix(F2, m1))
r2o <- fit_voom(dge, X2o, m1$osd_id)
repro <- rbindlist(lapply(c("CFD_g_bl" = "g_bl", "CFD_carbon_12h" = "carbon_12h",
                            "CFDmodel_Flight_FLT_vs_GC" = "flightFLT"), function(co) {
  nm <- names(which(c("CFD_g_bl" = "g_bl", "CFD_carbon_12h" = "carbon_12h",
                      "CFDmodel_Flight_FLT_vs_GC" = "flightFLT") == co))
  s <- stored(nm); r <- tt(r2o$fit, co); g <- intersect(s$gene_id, r$gene_id)
  data.table(contrast = nm, n_common = length(g),
             cor_t = cor(s$t[match(g, s$gene_id)], r$t[match(g, r$gene_id)]),
             stored_sig = sum(s$adj.P.Val < 0.05 & abs(s$logFC) > 1),
             refit_sig = n_sig(r)[["fdr05_lfc1"]])
}))
say(sprintf("  Model 2 rho=%.4f (stored 0.7782)", r2o$rho))
print(repro); fwrite(repro, file.path(OUT, "reproduction_check.tsv"), sep = "\t")

# ---- 4. Refit with v2 covariates --------------------------------------------
say("\n== 4. Refit, v2 covariates ==")
m2 <- prep_meta(meta2)
X1 <- drop_aliased(model.matrix(F1, m2))
r1 <- fit_voom(dge, X1, m2$osd_id)
X2 <- drop_aliased(model.matrix(F2, m2))
say("  Model 2 design columns: ", paste(colnames(X2), collapse = ", "))
r2 <- fit_voom(dge, X2, m2$osd_id)
say(sprintf("  rho: Model 1 %.4f, Model 2 %.4f", r1$rho, r2$rho))

# collinearity of the CFD covariates with study and hardware
vif <- function(X, col) { o <- setdiff(colnames(X), c(col, "X.Intercept.")); 1 / (1 - summary(lm(X[, col] ~ X[, o]))$r.squared) }
coll <- list()
for (co in intersect(c("g_bl", "carbon_12h"), colnames(X2))) {
  coll[[co]] <- c(vif = vif(X2, co), r2_study = summary(lm(X2[, co] ~ m2$osd_id))$r.squared)
  say(sprintf("  VIF %s = %.1f; R^2 on study dummies = %.3f", co, coll[[co]][["vif"]], coll[[co]][["r2_study"]]))
}

res <- list(
  Flight_within_BRIC     = tt(r1$fit, "flightFLT"),
  Flight_avg_over_HW     = tt(r1$fit, contrast = avg_flight(X1)),
  Flight_x_HW_CARA       = tt(r1$fit, "flightFLT.cfd_hardwareCARA"),
  Flight_x_HW_VEGGIE     = tt(r1$fit, "flightFLT.cfd_hardwareVEGGIE"),
  CFDv2_Flight           = tt(r2$fit, "flightFLT"),
  CFDv2_g_bl             = if ("g_bl" %in% colnames(X2)) tt(r2$fit, "g_bl"),
  CFDv2_carbon_12h       = if ("carbon_12h" %in% colnames(X2)) tt(r2$fit, "carbon_12h"))
res <- Filter(Negate(is.null), res)

# ---- 4b. Platform-separated sensitivity --------------------------------------
say("\n== 4b. Platform-separated sensitivity (v2) ==")
is_rna <- m2$omics_type == "RNAseq"
mr <- droplevels(m2[is_rna, ]); ma <- droplevels(m2[!is_rna, ])
say(sprintf("  RNA-seq: %d samples, %d studies; microarray: %d samples, %d studies",
            nrow(mr), nlevels(mr$osd_id), nrow(ma), nlevels(ma$osd_id)))
say("  RNA-seq hardware x flight:"); print(table(mr$cfd_hardware, mr$flight))
say("  microarray hardware x flight:"); print(table(ma$cfd_hardware, ma$flight))

# RNA-seq: voom on expected counts (RSEM/unnormalized), rounded
raw <- mapply(load_rnaseq, names(rnaseq_raw_files), rnaseq_raw_files, SIMPLIFY = FALSE)
C <- outer_merge(raw)
C <- round(C[rowSums(is.na(C)) == 0, rownames(mr)])
C <- C[grepl("^AT[1-5CM]G[0-9]{5}", rownames(C)), ]
dr <- DGEList(counts = C)
dr <- calcNormFactors(dr[filterByExpr(dr, group = interaction(mr$cfd_hardware, mr$flight)), , keep.lib.sizes = FALSE])
say(sprintf("  RNA-seq genes after filterByExpr: %d", nrow(dr)))
fit_platform <- function(f, meta, fitter, data) {
  # drop factors with one level in this subset (all microarray samples are BRIC)
  one <- Filter(function(v) is.factor(meta[[v]]) && nlevels(droplevels(meta[[v]])) < 2, all.vars(f))
  if (length(one)) {
    lab <- attr(terms(f), "term.labels")
    lab <- lab[!sapply(strsplit(lab, ":"), function(p) any(p %in% one))]
    f <- reformulate(lab)
  }
  X <- drop_aliased(model.matrix(f, meta)); r <- fitter(data, X, meta$osd_id); r$X <- X; r
}
rr1 <- fit_platform(F1, mr, fit_voom, dr)
rr2 <- fit_platform(F2, mr, fit_voom, dr)
# microarray: limma on log2 intensities (complete genes across the 4 arrays)
A <- H[, rownames(ma)]; A <- A[rowSums(is.na(A)) == 0, ]
say(sprintf("  microarray genes: %d; value range %.1f-%.1f (log2)", nrow(A), min(A), max(A)))
ra1 <- fit_platform(F1, ma, fit_limma, A)
ra2 <- fit_platform(F2, ma, fit_limma, A)
say(sprintf("  rho: RNA-seq M1 %.3f M2 %.3f; array M1 %.3f M2 %.3f", rr1$rho, rr2$rho, ra1$rho, ra2$rho))
say("  RNA-seq Model 2 columns: ", paste(colnames(rr2$X), collapse = ", "))
say("  array   Model 2 columns: ", paste(colnames(ra2$X), collapse = ", "))

plat <- function(r, nm, co = NULL, contrast = NULL) {
  if (!is.null(co) && !(co %in% colnames(r$X))) return(NULL)
  x <- tt(r$fit, co, contrast); x$contrast <- nm; x
}
res_p <- Filter(Negate(is.null), list(
  RNAseq_Flight_within_BRIC = plat(rr1, "RNAseq_Flight_within_BRIC", "flightFLT"),
  RNAseq_Flight_avg_over_HW = plat(rr1, "RNAseq_Flight_avg_over_HW", contrast = avg_flight(rr1$X)),
  RNAseq_CFDv2_g_bl         = plat(rr2, "RNAseq_CFDv2_g_bl", "g_bl"),
  RNAseq_CFDv2_carbon_12h   = plat(rr2, "RNAseq_CFDv2_carbon_12h", "carbon_12h"),
  Array_Flight_within_BRIC  = plat(ra1, "Array_Flight_within_BRIC", "flightFLT"),
  Array_CFDv2_g_bl          = plat(ra2, "Array_CFDv2_g_bl", "g_bl"),
  Array_CFDv2_carbon_12h    = plat(ra2, "Array_CFDv2_carbon_12h", "carbon_12h")))

# ---- Write ----------------------------------------------------------------------
all <- c(res, res_p)
summ <- rbindlist(lapply(names(all), function(n) {
  r <- all[[n]]; s <- n_sig(r)
  data.table(contrast = n, genes_tested = nrow(r), fdr05 = s[["fdr05"]], fdr05_absLFC1 = s[["fdr05_lfc1"]],
             up = sum(r$adj.P.Val < 0.05 & r$logFC > 0), down = sum(r$adj.P.Val < 0.05 & r$logFC < 0))
}))
print(summ)
fwrite(summ, file.path(OUT, "de_summary_v2.tsv"), sep = "\t")
for (n in names(all)) fwrite(all[[n]], file.path(OUT, paste0("de_v2_", n, ".tsv.gz")), sep = "\t")

# concordance of flight effects: pooled vs per-platform, within-BRIC vs averaged
conc <- function(a, b) { g <- intersect(a$gene_id, b$gene_id)
  cor(a$logFC[match(g, a$gene_id)], b$logFC[match(g, b$gene_id)]) }
cc <- data.table(
  pair = c("pooled within-BRIC vs pooled avg-over-HW", "pooled avg-over-HW vs RNA-seq avg-over-HW",
           "pooled within-BRIC vs array within-BRIC", "RNA-seq within-BRIC vs array within-BRIC"),
  logFC_cor = c(conc(res$Flight_within_BRIC, res$Flight_avg_over_HW),
                conc(res$Flight_avg_over_HW, res_p$RNAseq_Flight_avg_over_HW),
                conc(res$Flight_within_BRIC, res_p$Array_Flight_within_BRIC),
                conc(res_p$RNAseq_Flight_within_BRIC, res_p$Array_Flight_within_BRIC)))
print(cc); fwrite(cc, file.path(OUT, "flight_effect_concordance.tsv"), sep = "\t")
writeLines(c(log_lines, "", capture.output(sessionInfo())), file.path(OUT, "refit_log.txt"))

# ---- REANALYSIS.md section 9 (inserted before "## Not done here") -------------
v1n <- function(nm) { s <- stored(nm); sum(s$adj.P.Val < 0.05) }
S <- function(n, col = "fdr05") format(summ[contrast == n][[col]], big.mark = ",")
f3 <- function(x) sprintf("%.3f", x); f4 <- function(x) sprintf("%.4f", x)
rm1 <- repro_m1[[pick]]; cc_v <- setNames(cc$logFC_cor, cc$pair)
sec <- c(
  "## 9. Gene-level limma refit with the v2 covariates",
  "",
  "`code/factorial_model/refit_gene_level_2026-10.R`, from the 11 OSDR processed files (7 RNA-seq, 4 microarray; download manifest with SHA-256 in `Airflow_omics_osdr_data/raw_osdr/MANIFEST.json`). Gene tables: `gene_level/`.",
  "",
  "**Reconstructing the missing filtering step.** No script in the package writes `harmonized_expression_matrix_filtered.tsv`. Of two candidate rules, only one reproduces the stored fit:",
  "",
  "| rule | genes after CPM filter | inter-study ρ (Model 1) | t vs stored Flight_FLT_vs_GC |",
  "|---|---|---|---|",
  sprintf("| drop genes missing in any study | %s | %s | %s |", format(repro_m1$complete_cases[["genes"]], big.mark = ","),
          f4(repro_m1$complete_cases[["rho"]]), f4(repro_m1$complete_cases[["t_cor"]])),
  sprintf("| set missing values to 0 | %s | %s | %s |", format(repro_m1$na_to_zero[["genes"]], big.mark = ","),
          f4(repro_m1$na_to_zero[["rho"]]), f4(repro_m1$na_to_zero[["t_cor"]])),
  "| stored fit | 25,453 | 0.7324 | — |",
  "",
  sprintf("Model 2 reproduces too (ρ %s vs stored 0.7782; t-statistic correlation with the stored g_bl, carbon_12h and flight tables %s–%s).",
          f4(r2o$rho), f3(min(repro$cor_t)), f3(max(repro$cor_t))),
  sprintf("So the original pipeline **set unmeasured genes to 0**: of %s genes in the outer join, those measured only by RNA-seq enter every microarray sample as zero expression (and genes absent from one RNA-seq study are zero there). This builds study × platform differences into the matrix. It also passes microarray log2 intensities through `DGEList`/`voom` as if they were counts.",
          format(nrow(H), big.mark = ",")),
  "",
  "**Pooled refit (same pipeline, v2 covariates).** Genes at FDR < 0.05 (|logFC| > 1 in brackets):",
  "",
  "| contrast | v1 (stored) | v2 refit |",
  "|---|---|---|",
  sprintf("| flight within BRIC (stored name `Flight_FLT_vs_GC`) | %s | %s (%s) |", format(v1n("Flight_FLT_vs_GC"), big.mark = ","), S("Flight_within_BRIC"), S("Flight_within_BRIC", "fdr05_absLFC1")),
  sprintf("| **flight averaged over BRIC/CARA/VEGGIE** | not fitted | %s (%s) |", S("Flight_avg_over_HW"), S("Flight_avg_over_HW", "fdr05_absLFC1")),
  sprintf("| flight × CARA | %s | %s |", format(v1n("Flight_x_HW_CARA_vs_BRIC"), big.mark = ","), S("Flight_x_HW_CARA")),
  sprintf("| flight × VEGGIE | %s | %s |", format(v1n("Flight_x_HW_VEGGIE_vs_BRIC"), big.mark = ","), S("Flight_x_HW_VEGGIE")),
  sprintf("| Model 2: g_bl | %s | %s |", format(v1n("CFD_g_bl"), big.mark = ","), S("CFDv2_g_bl")),
  sprintf("| Model 2: carbon_12h | %s | %s |", format(v1n("CFD_carbon_12h"), big.mark = ","), S("CFDv2_carbon_12h")),
  "",
  sprintf("- The within-BRIC and averaged flight effects agree only weakly (logFC correlation %s). The stored headline is the BRIC-specific effect.", f3(cc_v[[1]])),
  sprintf("- g_bl: VIF %.1f in Model 2. It takes 4 values and recodes flight × scale (§4), so its gene count is not a boundary-layer effect.", coll$g_bl[["vif"]]),
  sprintf("- carbon_12h: study membership explains %.1f %% of its variance (v2 sets every BRIC sample to 1 %%). Its genes are a BRIC-studies-versus-the-rest contrast, not a carbon dose response. The |logFC| > 1 count is 0 because the unit is 1 %% of Earth carbon.", 100 * coll$carbon_12h[["r2_study"]]),
  "",
  "**Platform-separated sensitivity.** RNA-seq: voom on RSEM expected counts, `filterByExpr`, no zero-filling. Microarray: limma-trend on log2 intensities. All microarray samples are BRIC, so hardware terms and carbon_12h drop out there.",
  "",
  "| contrast | RNA-seq | microarray |",
  "|---|---|---|",
  sprintf("| samples / studies | %d / %d | %d / %d |", nrow(mr), nlevels(mr$osd_id), nrow(ma), nlevels(ma$osd_id)),
  sprintf("| flight within BRIC | %s | %s |", S("RNAseq_Flight_within_BRIC"), S("Array_Flight_within_BRIC")),
  sprintf("| flight averaged over hardware | %s | — |", S("RNAseq_Flight_avg_over_HW")),
  sprintf("| g_bl | %s | %s |", S("RNAseq_CFDv2_g_bl"), S("Array_CFDv2_g_bl")),
  sprintf("| carbon_12h | %s | not estimable |", S("RNAseq_CFDv2_carbon_12h")),
  "",
  sprintf("- Within BRIC, the RNA-seq and microarray flight effects are uncorrelated (logFC r = %s), and the pooled within-BRIC effect does not track the microarray one (r = %s). The pooled averaged effect does track the RNA-seq-only one (r = %s).",
          f3(cc_v[[4]]), f3(cc_v[[3]]), f3(cc_v[[2]])),
  "- **Reading:** the gene level agrees with §5 and §8. Neither CFD covariate is identifiable from study and hardware, so neither gene count supports a gas-exchange interpretation. The flight response differs by hardware and by platform/study era, and pooling with zero-filled genes hides that. Report the averaged flight effect, by platform, as the headline.",
  "- Caveats: the inter-study correlation uses one `duplicateCorrelation` pass, matching the original script. The variance-partition section of `fit_factorial_model.R` was not rerun (variancePartition not installed). Pathway (GSVA) scores were not recomputed from the rebuilt matrix.",
  "")
rpath <- "tables/reanalysis_2026-10/REANALYSIS.md"
txt <- paste(readLines(rpath), collapse = "\n")
txt <- sub("(?s)\n## 9\\. Gene-level limma refit.*?(?=\n## 10\\. |\n## Not done here|$)", "", txt, perl = TRUE)
# keep §9 ahead of §10 (followups_2026-10.R) when rerun
marker <- if (grepl("\n## 10. ", txt, fixed = TRUE)) "\n## 10. " else "\n## Not done here"
block <- paste(sec, collapse = "\n")
txt <- if (grepl(marker, txt, fixed = TRUE)) sub(marker, paste0("\n", block, marker), txt, fixed = TRUE) else paste0(txt, "\n", block)
writeLines(txt, rpath)
say("updated REANALYSIS.md section 9")
say("done")
