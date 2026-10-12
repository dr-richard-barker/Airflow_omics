#!/usr/bin/env Rscript
# =============================================================================
# Follow-ups to the gene-level refit (REANALYSIS §10)
# =============================================================================
# 4a. Pooled limma-voom fit without zero-filled genes: only genes measured in all
#     11 studies, same pipeline and v2 covariates as §9, compared with the
#     zero-filled fit.
# 4b. Pathway scores recomputed from the rebuilt matrix:
#     - reproduce the stored "GSVA-like" scores (hormone_carbon_model.py:
#       log2(x + 1) of the zero-filled matrix, genes z-scored across all 307
#       samples, mean per pathway);
#     - port the §5 within-study model (cfd_reanalysis.py) to R and check it
#       against the stored pathway_within_study_model.tsv;
#     - rescore per study: RNA-seq as log2-CPM from expected counts, microarray
#       as the published log2 values (no second log), genes z-scored within each
#       study, no zero-filling; rerun the within-study model.
#
# Usage (from the repo root; run after refit_gene_level_2026-10.R):
#   OSDR_RAW=~/Documents/Airflow_omics_osdr_data/raw_osdr \
#   Rscript code/factorial_model/followups_2026-10.R
# Writes tables/reanalysis_2026-10/followups/ and inserts §10 into REANALYSIS.md.
# =============================================================================

suppressPackageStartupMessages({
  library(limma)
  library(edgeR)
  library(data.table)
})

RAW <- path.expand(Sys.getenv("OSDR_RAW", "~/Documents/Airflow_omics_osdr_data/raw_osdr"))
OUT <- "tables/reanalysis_2026-10/followups"
GL <- "tables/reanalysis_2026-10/gene_level"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
log_lines <- character()
say <- function(...) { s <- paste0(...); cat(s, "\n"); log_lines <<- c(log_lines, s) }
source("code/factorial_model/osdr_helpers_2026-10.R")

meta2 <- fread("tables/reanalysis_2026-10/expression_metadata_with_cfd_v2.tsv", data.table = FALSE)
m2 <- prep_meta(meta2)
mats <- c(mapply(load_rnaseq, names(rnaseq_files), rnaseq_files, SIMPLIFY = FALSE),
          mapply(load_array, names(array_files), array_files, SIMPLIFY = FALSE))
H <- outer_merge(mats)[, m2$full_sample_id]

# =============================================================================
# 4a. Pooled fit on genes measured in every study
# =============================================================================
say("== 4a. Pooled fit without zero-filled genes ==")
CC <- H[rowSums(is.na(H)) == 0, ]
keep <- cpm_keep(CC)
dge <- calcNormFactors(DGEList(counts = CC)[keep, , keep.lib.sizes = FALSE], method = "TMM")
say(sprintf("  genes measured in all 11 studies: %d; after CPM filter: %d (zero-filled fit: 25,453)", nrow(CC), nrow(dge)))
X1 <- drop_aliased(model.matrix(F1, m2)); r1 <- fit_voom(dge, X1, m2$osd_id)
X2 <- drop_aliased(model.matrix(F2, m2)); r2 <- fit_voom(dge, X2, m2$osd_id)
say(sprintf("  rho: Model 1 %.4f, Model 2 %.4f (zero-filled: 0.7327, 0.7588)", r1$rho, r2$rho))
cc_res <- list(
  Flight_within_BRIC = tt(r1$fit, "flightFLT"),
  Flight_avg_over_HW = tt(r1$fit, contrast = avg_flight(X1)),
  Flight_x_HW_CARA   = tt(r1$fit, "flightFLT.cfd_hardwareCARA"),
  Flight_x_HW_VEGGIE = tt(r1$fit, "flightFLT.cfd_hardwareVEGGIE"),
  CFDv2_g_bl         = tt(r2$fit, "g_bl"),
  CFDv2_carbon_12h   = tt(r2$fit, "carbon_12h"))
zf <- function(n) read.delim(gzfile(file.path(GL, paste0("de_v2_", n, ".tsv.gz"))))   # fread needs R.utils for .gz
rna_avg <- zf("RNAseq_Flight_avg_over_HW")
cmp <- rbindlist(lapply(names(cc_res), function(n) {
  a <- cc_res[[n]]; z <- zf(n); g <- intersect(a$gene_id, z$gene_id)
  ia <- match(g, a$gene_id); iz <- match(g, z$gene_id)
  data.table(contrast = n,
             zero_filled_genes = nrow(z), zero_filled_fdr05 = sum(z$adj.P.Val < 0.05),
             complete_genes = nrow(a), complete_fdr05 = sum(a$adj.P.Val < 0.05),
             common_genes = length(g),
             logFC_cor_common = cor(a$logFC[ia], z$logFC[iz]),
             t_cor_common = cor(a$t[ia], z$t[iz]),
             sig_both = sum(a$adj.P.Val[ia] < 0.05 & z$adj.P.Val[iz] < 0.05))
}))
g <- intersect(cc_res$Flight_avg_over_HW$gene_id, rna_avg$gene_id)
cc_vs_rna <- cor(cc_res$Flight_avg_over_HW$logFC[match(g, cc_res$Flight_avg_over_HW$gene_id)],
                 rna_avg$logFC[match(g, rna_avg$gene_id)])
print(cmp); say(sprintf("  complete-case averaged flight vs RNA-seq-only averaged flight: logFC r = %.3f (n = %d)", cc_vs_rna, length(g)))
fwrite(cmp, file.path(OUT, "complete_case_vs_zero_filled.tsv"), sep = "\t")
for (n in names(cc_res)) fwrite(cc_res[[n]], file.path(OUT, paste0("de_v2_complete_", n, ".tsv.gz")), sep = "\t")

# =============================================================================
# 4b. Pathway scores
# =============================================================================
say("\n== 4b. Pathway scores ==")
gs_tab <- fread(file.path(OUT, "pathway_genesets.tsv"))   # written from hormone_carbon_model.py (see §10)
GS <- split(gs_tab$gene_id, factor(gs_tab$pathway, levels = unique(gs_tab$pathway)))
# curation fixes (pathway_geneset_corrections.tsv) apply to the rescoring only; the
# reproduction of the stored scores uses the original lists
corr <- fread(file.path(OUT, "pathway_geneset_corrections.tsv"))
GS_fixed <- GS
for (k in seq_len(nrow(corr))) if (corr$action[k] == "remove")
  GS_fixed[[corr$pathway[k]]] <- setdiff(GS_fixed[[corr$pathway[k]]], corr$gene_id[k])
n_dup <- sum(lengths(GS) - lengths(lapply(GS, unique))); n_dup_pw <- sum(lengths(GS) > lengths(lapply(GS, unique)))
rep_tab <- rbindlist(lapply(names(GS), function(pw) { t <- table(GS[[pw]]); data.table(pathway = pw, gene_id = names(t), times = as.integer(t), entries = length(GS[[pw]])) }))[times > 1][order(-times)]
fwrite(rep_tab, file.path(OUT, "pathway_geneset_repeats.tsv"), sep = "\t")

# genes x samples (already on log scale) -> samples x pathways. dedup = FALSE mirrors
# hormone_carbon_model.py, where a gene listed twice in a set is counted twice.
score_matrix <- function(E, dedup = TRUE, sets = if (dedup) GS_fixed else GS) {
  sapply(sets, function(genes) {
    v <- genes[genes %in% rownames(E)]; if (dedup) v <- unique(v)
    if (length(v) < 3) return(rep(NA_real_, ncol(E)))
    P <- E[v, , drop = FALSE]
    s <- apply(P, 1, sd); s[s == 0 | is.na(s)] <- 1
    colMeans((P - rowMeans(P)) / s, na.rm = TRUE)
  })
}

# reproduction of the stored scores
Z <- H; Z[is.na(Z)] <- 0
S_orig <- score_matrix(log2(Z + 1), dedup = FALSE); rownames(S_orig) <- colnames(H)
sd_ <- fread("tables/pathway_scores_gsva.tsv", data.table = FALSE)
stored <- as.matrix(sd_[, -1]); rownames(stored) <- sd_[[1]]
stored <- stored[rownames(S_orig), colnames(S_orig)]; storage.mode(stored) <- "double"
repro_max <- max(abs(S_orig - stored)); repro_cor <- min(diag(cor(S_orig, stored)))
say(sprintf("  rebuilt vs stored scores: max |diff| = %.2e, min per-pathway r = %.6f", repro_max, repro_cor))

# within-study model (port of cfd_reanalysis.py §5)
ws_design <- function(meta) {
  f <- ifelse(meta$flight == "FLT", 0.5, -0.5)
  s <- ifelse(meta$cfd_scale == "rosette", 0.5, -0.5)
  st <- model.matrix(~ factor(meta$osd_id))[, -1, drop = FALSE]
  X <- cbind(intercept = 1, flight = f, scale = s, flight_x_scale = f * s, st)
  X
}
ols_cluster <- function(y, X, groups) {
  n <- nrow(X); k <- ncol(X)
  XtXi <- MASS::ginv(crossprod(X))
  beta <- XtXi %*% crossprod(X, y); res <- as.vector(y - X %*% beta)
  meat <- matrix(0, k, k)
  for (g in unique(groups)) { i <- groups == g; s <- crossprod(X[i, , drop = FALSE], res[i]); meat <- meat + s %*% t(s) }
  G <- length(unique(groups)); cc <- (G / (G - 1)) * ((n - 1) / (n - k))
  V <- cc * XtXi %*% meat %*% XtXi
  se <- sqrt(diag(V)); tt_ <- as.vector(beta) / se
  list(beta = as.vector(beta), se = se, t = tt_, p = 2 * pt(abs(tt_), df = G - 1, lower.tail = FALSE))
}
within_model <- function(S, meta) {
  ok <- complete.cases(S)
  rbindlist(lapply(colnames(S), function(pw) {
    keep <- !is.na(S[, pw])
    X <- ws_design(meta[keep, ]); r <- ols_cluster(S[keep, pw], X, meta$osd_id[keep])
    rk <- qr(X)$rank
    rbindlist(lapply(c("flight", "scale", "flight_x_scale"), function(term) {
      i <- match(term, colnames(X))
      # Not testable if the term is in the span of the other columns, or if it varies
      # within fewer than 2 studies: with study fixed effects it is then identified by one
      # cluster, and CR1 cluster-robust errors collapse toward zero (SE ~ 1e-16), the
      # artefact §7 guards against for CARA. Report NA instead.
      g <- meta$osd_id[keep]
      n_var <- sum(tapply(X[, i], droplevels(factor(g)), function(x) length(unique(x)) > 1), na.rm = TRUE)
      est <- qr(X[, -i, drop = FALSE])$rank < rk && n_var >= 2
      data.table(pathway = pw, term = term, estimable = est,
                 estimate = if (est) r$beta[i] else NA_real_, se = if (est) r$se[i] else NA_real_,
                 t = if (est) r$t[i] else NA_real_, p_value = if (est) r$p[i] else NA_real_,
                 n = sum(keep), studies = length(unique(meta$osd_id[keep])))
    }))
  }))[, fdr := p.adjust(p_value, "BH"), by = term][]
}
ws_stored <- fread("tables/reanalysis_2026-10/pathway_within_study_model.tsv")
ws_port <- within_model(stored, m2)
mm <- merge(ws_port, ws_stored, by = c("pathway", "term"), suffixes = c("", ".py"))
port_dmax <- max(abs(mm$estimate - mm$estimate.py)); port_pmax <- max(abs(mm$p_value - mm$p_value.py))
say(sprintf("  R port vs cfd_reanalysis.py: max |d estimate| = %.2e, max |d p| = %.2e over %d rows", port_dmax, port_pmax, nrow(mm)))

# corrected scores: per study, proper log scale, within-study z-scores, no zero-filling
raw <- mapply(load_rnaseq, names(rnaseq_raw_files), rnaseq_raw_files, SIMPLIFY = FALSE)
per_study <- lapply(split(m2$full_sample_id, m2$osd_id), function(ids) {
  osd <- as.character(m2[ids[1], "osd_id"])
  if (m2[ids[1], "omics_type"] == "RNAseq") {
    C <- raw[[osd]][, ids]; C <- round(C[!is.na(rowSums(C)), ])
    d <- DGEList(counts = C)
    d <- calcNormFactors(d[filterByExpr(d, group = m2[ids, "flight"]), , keep.lib.sizes = FALSE])
    E <- cpm(d, log = TRUE)
  } else {
    E <- H[, ids]; E <- E[rowSums(is.na(E)) == 0, ]
  }
  S <- score_matrix(E); rownames(S) <- ids
  attr(S, "n_genes") <- sapply(GS_fixed, function(g) length(intersect(g, rownames(E))))   # unique genes
  S
})
S_new <- do.call(rbind, per_study)[m2$full_sample_id, ]
cover <- sapply(per_study, attr, "n_genes")
# same layout as tables/pathway_scores_gsva.tsv (unnamed index column) so cfd_reanalysis.py can read it
write.table(S_new, file.path(OUT, "pathway_scores_within_study.tsv"), sep = "\t", quote = FALSE, col.names = NA)
fwrite(data.table(pathway = rownames(cover), cover), file.path(OUT, "pathway_gene_coverage_by_study.tsv"), sep = "\t")
n_missing <- sum(cover < 3)
say(sprintf("  per-study pathway coverage: %d of %d pathway x study cells have < 3 genes (score left NA)", n_missing, length(cover)))

ws_new <- within_model(S_new, m2)
rna <- m2$omics_type == "RNAseq"
ws_rna <- within_model(S_new[rna, ], m2[rna, ])
fwrite(ws_new, file.path(OUT, "pathway_within_study_model_rescored.tsv"), sep = "\t")
fwrite(ws_rna, file.path(OUT, "pathway_within_study_model_rescored_rnaseq.tsv"), sep = "\t")
# NB: the model table has a column named `t`, so never filter with `term == t` inside [ ]
cnt <- function(w) sapply(c("flight", "scale", "flight_x_scale"), function(tm) {
  f <- w$fdr[w$term == tm]; if (all(is.na(f))) NA_integer_ else sum(f < 0.05, na.rm = TRUE) })
counts <- rbind(stored_scores = cnt(ws_port), rescored = cnt(ws_new), rescored_rnaseq_only = cnt(ws_rna))
print(counts)
fl <- merge(ws_port[term == "flight", .(pathway, t_orig = t)], ws_new[term == "flight", .(pathway, t_new = t, est_new = estimate, fdr_new = fdr)], by = "pathway")
t_cor <- cor(fl$t_orig, fl$t_new)
say(sprintf("  flight t, stored vs rescored scores: r = %.3f", t_cor))
top <- ws_new[term == "flight"][order(p_value)][1:5]
print(top)

# positive control, defined as in cfd_reanalysis.py §5: plant k x (median over pathways
# of the score SD across all samples) as a flight effect in every pathway. Run on the
# stored scores too, to check the port reproduces §5's 28/28.
X_all <- ws_design(m2)
flt <- ifelse(m2$flight == "FLT", 1, 0)
pos_control <- function(S) {
  sd_med <- median(apply(S, 2, sd, na.rm = TRUE))
  sapply(c(0.25, 0.5), function(k) sum(sapply(colnames(S), function(pw) {
    keep <- !is.na(S[, pw])
    ols_cluster(S[keep, pw] + k * sd_med * flt[keep], X_all[keep, ], m2$osd_id[keep])$p[2] < 0.05
  })))
}
detect_stored <- pos_control(stored); detect <- pos_control(S_new)
sd_ratio <- median(apply(S_new, 2, sd)) / median(apply(stored, 2, sd))
say(sprintf("  positive control: stored %d/28, %d/28; rescored %d/28, %d/28 (0.25 / 0.5 SD); rescored SD is %.2f x the stored SD",
            detect_stored[1], detect_stored[2], detect[1], detect[2], sd_ratio))
fin <- ws_new[term %in% c("flight", "scale", "flight_x_scale"), all(is.finite(se))] && ws_rna[, all(is.finite(se))]
say("  all flight / scale / interaction SEs finite: ", fin)

writeLines(c(log_lines, "", capture.output(sessionInfo())), file.path(OUT, "followups_log.txt"))

# =============================================================================
# REANALYSIS.md section 10
# =============================================================================
fmt <- function(x) format(x, big.mark = ",")
row_cmp <- function(n, label) { r <- cmp[contrast == n]
  sprintf("| %s | %s | %s | %.3f | %s |", label, fmt(r$zero_filled_fdr05), fmt(r$complete_fdr05), r$logFC_cor_common, fmt(r$sig_both)) }
cn <- function(row, tm) { v <- counts[row, tm]; if (is.na(v)) "not estimable" else as.character(v) }
sec <- c(
  "## 10. Follow-ups: no zero-filling at gene level; pathway scores recomputed",
  "",
  "`code/factorial_model/followups_2026-10.R` (shared loaders in `osdr_helpers_2026-10.R`). Tables: `followups/`.",
  "",
  sprintf("**Pooled gene-level fit on genes measured in all 11 studies.** Same pipeline and v2 covariates as §9, but only the %s genes present in every study (%s after the CPM filter), so nothing is zero-filled. Inter-study ρ %.4f (zero-filled 0.7327).",
          fmt(nrow(CC)), fmt(nrow(dge)), r1$rho),
  "",
  "| contrast | zero-filled FDR < 0.05 (of 25,453) | complete-case FDR < 0.05 | logFC r on common genes | significant in both |",
  "|---|---|---|---|---|",
  row_cmp("Flight_within_BRIC", "flight within BRIC"),
  row_cmp("Flight_avg_over_HW", "flight averaged over hardware"),
  row_cmp("Flight_x_HW_CARA", "flight × CARA"),
  row_cmp("Flight_x_HW_VEGGIE", "flight × VEGGIE"),
  row_cmp("CFDv2_g_bl", "Model 2: g_bl"),
  row_cmp("CFDv2_carbon_12h", "Model 2: carbon_12h"),
  "",
  sprintf("- The averaged flight effect from the complete-case fit agrees with the RNA-seq-only fit at logFC r = %.3f.", cc_vs_rna),
  sprintf("- carbon_12h loses most of its genes (%s -> %s): carbon_12h is 1 %% for every BRIC sample and all microarray samples are BRIC, so genes that were zero-filled in the microarray studies looked like a carbon effect. g_bl goes from %s to %s.",
          fmt(cmp[contrast == "CFDv2_carbon_12h"]$zero_filled_fdr05), fmt(cmp[contrast == "CFDv2_carbon_12h"]$complete_fdr05),
          fmt(cmp[contrast == "CFDv2_g_bl"]$zero_filled_fdr05), fmt(cmp[contrast == "CFDv2_g_bl"]$complete_fdr05)),
  "- Removing the zero-filled genes does not rescue either CFD covariate: they stay confounded with study and hardware by design (§4, §9). Microarray and RNA-seq values are still pooled in one voom model here, which is why §9's platform-separated fit remains the preferred analysis.",
  "",
  "**Pathway scores recomputed from the rebuilt matrix.**",
  "",
  sprintf("- Reproduction: the original recipe (log2(x + 1) of the zero-filled matrix, genes z-scored across all 307 samples) applied to the rebuilt matrix reproduces `pathway_scores_gsva.tsv` (max |difference| %.1e; minimum per-pathway r %.6f). That recipe log-transforms the microarray values a second time, and its z-scores mix platforms and studies, so pathway scores carry the same zero-filling artefact as the gene matrix.", repro_max, repro_cor),
  sprintf("- An R port of the §5 within-study model matches `cfd_reanalysis.py` on the stored scores (max |Δ estimate| %.1e, max |Δ p| %.1e).", port_dmax, port_pmax),
  sprintf("- The gene sets in `hormone_carbon_model.py` repeat %d gene entries within their own set (in %d of 28 pathways); the original scores count a repeated gene more than once. The worst case is %s in %s: %d of its %d entries, so the stored %s score is mostly that one gene. The rescored version counts each gene once (`pathway_genesets.tsv`; repeats in `pathway_geneset_repeats.tsv`).", n_dup, n_dup_pw, rep_tab$gene_id[1], rep_tab$pathway[1], rep_tab$times[1], rep_tab$entries[1], rep_tab$pathway[1]),
  sprintf("- Curation check against UniProt: %s. The rescored scores apply this (`pathway_geneset_corrections.tsv`); heat_shock is left with %d genes, all heat shock chaperones, and still lacks the heat shock factors and HSP101.",
          paste(sprintf("%s is removed from %s (%s; %s)", corr$gene_id, corr$pathway, corr$reason, corr$source), collapse = "; "),
          length(unique(GS_fixed[["heat_shock"]]))),
  sprintf("- Rescored per study: RNA-seq as log2-CPM from RSEM expected counts (TMM, `filterByExpr`), microarray on its published log2 scale, each gene z-scored within its study, no zero-filling. %d of %d pathway × study cells have fewer than 3 measured genes and are left missing (`pathway_gene_coverage_by_study.tsv`).", n_missing, length(cover)),
  "",
  "Within-study model (`score ~ study + flight × scale`, CR1 by study, BH per term), pathways with FDR < 0.05 out of 28:",
  "",
  "| scores | flight | scale | flight × scale |",
  "|---|---|---|---|",
  sprintf("| stored (zero-filled, pooled z) — §5 | %s | %s | %s |", cn("stored_scores", "flight"), cn("stored_scores", "scale"), cn("stored_scores", "flight_x_scale")),
  sprintf("| rescored, all 11 studies | %s | %s | %s |", cn("rescored", "flight"), cn("rescored", "scale"), cn("rescored", "flight_x_scale")),
  sprintf("| rescored, RNA-seq studies only (7) | %s | %s | %s |", cn("rescored_rnaseq_only", "flight"), cn("rescored_rnaseq_only", "scale"), cn("rescored_rnaseq_only", "flight_x_scale")),
  "",
  sprintf("- Flight t-statistics from the stored and rescored scores correlate at r = %.3f. Smallest rescored flight p: %s (p = %.2g, FDR = %.2g).",
          t_cor, top$pathway[1], top$p_value[1], top$fdr[1]),
  sprintf("- Positive control (as in §5: a flight effect of 0.25 / 0.5 × the median score SD planted in every pathway): stored scores %d/28 and %d/28 (reproducing §5), rescored %d/28 and %d/28. The rescored SD is %.2f × the stored SD because the stored scores' SD includes between-study and between-platform differences; relative to within-study variation, the §5 control planted larger effects than its label suggests.",
          detect_stored[1], detect_stored[2], detect[1], detect[2], sd_ratio),
  sprintf("- A term is tested only if it varies within at least two studies; otherwise, with study fixed effects, it is identified by a single cluster and cluster-robust errors collapse (as for CARA in §7). Among the RNA-seq studies only OSD-678 has both leaf and rosette samples, so scale is not testable there (%s); flight × scale still varies within every study and is tested.",
          cn("rescored_rnaseq_only", "scale")),
  paste0("- **Reading:** with corrected scores no pathway reaches FDR < 0.05 for flight, as in §5",
         if (counts["rescored", "scale"] > 0) sprintf("; %d pathway(s) now show a scale (rosette vs leaf) difference (%s)", counts["rescored", "scale"],
                                                      paste(ws_new[term == "scale" & fdr < 0.05]$pathway, collapse = ", ")) else "",
         ". But the corrected positive control shows the within-study test detects a 0.25 SD flight effect in only about a quarter of pathways and a 0.5 SD effect in about half, so the pathway-level null now reads as 'no large, consistent flight effect' rather than 'no effect'; effects below about half a within-study SD cannot be excluded at this sample of 11 studies. §5's statement that the null 'is not a powerless test' should be read with this correction. §11 reruns the §7 in-dish analysis and the §8 meta-regression on these scores."),
  "")
rpath <- "tables/reanalysis_2026-10/REANALYSIS.md"
txt <- paste(readLines(rpath), collapse = "\n")
# Insert section 10 after any lower-numbered and before any higher-numbered section,
# removing only the previous copy of section 10 (up to the next "## <number>." heading).
txt <- sub("(?s)\n## 10\\. .*?(?=\n## [0-9]+\\. |\n## Not done here|$)", "", txt, perl = TRUE)
nums <- as.integer(sub("^\n## ([0-9]+)\\. .*", "\\1", regmatches(txt, gregexpr("\n## [0-9]+\\. [^\n]*", txt))[[1]]))
higher <- nums[nums > 10]
marker <- if (length(higher)) paste0("\n## ", min(higher), ". ") else "\n## Not done here"
block <- paste(sec, collapse = "\n")
txt <- if (grepl(marker, txt, fixed = TRUE)) sub(marker, paste0("\n", block, marker), txt, fixed = TRUE) else paste0(txt, "\n", block)
writeLines(txt, rpath)
say("updated REANALYSIS.md section 10")
