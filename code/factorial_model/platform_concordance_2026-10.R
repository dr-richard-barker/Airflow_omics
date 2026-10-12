#!/usr/bin/env Rscript
# =============================================================================
# Is the RNA-seq vs microarray disagreement a platform effect? (REANALYSIS §13)
# =============================================================================
# §9 found that, within BRIC, the RNA-seq and microarray flight effects are
# uncorrelated (logFC r = -0.05). Platform is confounded with study and organ:
# the BRIC RNA-seq studies (OSD-321, -37, -38) are whole seedlings, the BRIC
# microarray studies (OSD-7, -17, -147, -205) are mostly hypocotyl, root, shoot or
# cell culture. So the flight effect is estimated separately in every
# study x organ unit, and the units are compared within and across platforms:
#   - if units agree within a platform but not across, the platform matters;
#   - if they disagree everywhere, the flight response is study-specific (§12),
#     and the platform split is one view of that.
# A matched comparison uses only whole-seedling units (present on both platforms).
#
# Per unit: RNA-seq, voom on RSEM expected counts (TMM, filterByExpr); microarray,
# limma-trend on the published log2 values. Model ~ flight (+ ecotype, light where
# they vary within the unit). Units need >= 2 FLT and >= 2 GC samples.
# Usage (repo root): Rscript code/factorial_model/platform_concordance_2026-10.R
# =============================================================================

suppressPackageStartupMessages({
  library(limma)
  library(edgeR)
  library(data.table)
})

RAW <- path.expand(Sys.getenv("OSDR_RAW", "~/Documents/Airflow_omics_osdr_data/raw_osdr"))
OUT <- "tables/reanalysis_2026-10/followups"
log_lines <- character()
say <- function(...) { s <- paste0(...); cat(s, "\n"); log_lines <<- c(log_lines, s) }
source("code/factorial_model/osdr_helpers_2026-10.R")

meta2 <- fread("tables/reanalysis_2026-10/expression_metadata_with_cfd_v2.tsv", data.table = FALSE)
m2 <- prep_meta(meta2)
m2$organ_simple <- as.character(m2$organ_simple)
H <- outer_merge(c(mapply(load_rnaseq, names(rnaseq_files), rnaseq_files, SIMPLIFY = FALSE),
                   mapply(load_array, names(array_files), array_files, SIMPLIFY = FALSE)))[, m2$full_sample_id]
raw <- mapply(load_rnaseq, names(rnaseq_raw_files), rnaseq_raw_files, SIMPLIFY = FALSE)

# ---- per-unit flight effects ----------------------------------------------------
m2$unit <- paste(m2$osd_id, m2$organ_simple, sep = "|")
units <- split(rownames(m2), m2$unit)
fit_unit <- function(ids) {
  mm <- m2[ids, ]
  if (min(table(factor(mm$flight, levels = c("GC", "FLT")))) < 2) return(NULL)
  vars <- c("ecotype_simple", "light")[sapply(c("ecotype_simple", "light"), function(v) length(unique(mm[[v]])) > 1)]
  X <- drop_aliased(model.matrix(reformulate(c("flight", vars)), droplevels(mm)))
  if (!"flightFLT" %in% colnames(X) || nrow(X) - ncol(X) < 2) return(NULL)
  osd <- as.character(mm$osd_id[1])
  if (mm$omics_type[1] == "RNAseq") {
    C <- raw[[osd]][, ids]; C <- round(C[!is.na(rowSums(C)), ])
    d <- DGEList(counts = C); d <- calcNormFactors(d[filterByExpr(d, group = mm$flight), , keep.lib.sizes = FALSE])
    f <- eBayes(lmFit(voom(d, X), X))
  } else {
    E <- H[, ids]; E <- E[rowSums(is.na(E)) == 0, ]
    f <- eBayes(lmFit(E, X), trend = TRUE)
  }
  r <- topTable(f, coef = "flightFLT", number = Inf, sort.by = "none")
  list(unit = mm$unit[1], osd = osd, organ = mm$organ_simple[1], platform = mm$omics_type[1],
       hardware = as.character(mm$cfd_hardware[1]), n_flt = sum(mm$flight == "FLT"), n_gc = sum(mm$flight == "GC"),
       fdr05 = sum(r$adj.P.Val < 0.05), lfc = setNames(r$logFC, rownames(r)))
}
fits <- Filter(Negate(is.null), lapply(units, fit_unit))
U <- rbindlist(lapply(fits, function(x) x[c("unit", "osd", "organ", "platform", "hardware", "n_flt", "n_gc", "fdr05")]))
setorder(U, hardware, platform, unit)
say("== per-unit flight effects =="); print(U)
fwrite(U, file.path(OUT, "platform_concordance_units.tsv"), sep = "\t")

# ---- pairwise concordance -------------------------------------------------------
L <- setNames(lapply(fits, `[[`, "lfc"), sapply(fits, `[[`, "unit"))
pairs <- combn(names(L), 2, simplify = FALSE)
P <- rbindlist(lapply(pairs, function(p) {
  g <- intersect(names(L[[p[1]]]), names(L[[p[2]]]))
  a <- U[unit == p[1]]; b <- U[unit == p[2]]
  data.table(unit_a = p[1], unit_b = p[2], genes = length(g),
             r = cor(L[[p[1]]][g], L[[p[2]]][g]), rho = cor(L[[p[1]]][g], L[[p[2]]][g], method = "spearman"),
             same_study = a$osd == b$osd, platform = if (a$platform == b$platform) paste("within", a$platform) else "across",
             bric_both = a$hardware == "BRIC" & b$hardware == "BRIC",
             seedling_both = a$organ == "whole_seedling" & b$organ == "whole_seedling")
}))
fwrite(P, file.path(OUT, "platform_concordance_pairs.tsv"), sep = "\t")
summ <- function(sub, label) {
  if (!nrow(sub)) return(NULL)
  data.table(comparison = label, pairs = nrow(sub), median_r = median(sub$r), min_r = min(sub$r), max_r = max(sub$r),
             median_rho = median(sub$rho))
}
dp <- P[same_study == FALSE]   # different studies only: same-study pairs share controls and handling
S <- rbindlist(list(
  summ(dp[bric_both & platform == "within RNAseq"], "BRIC, RNA-seq vs RNA-seq"),
  summ(dp[bric_both & platform == "within microarray"], "BRIC, microarray vs microarray"),
  summ(dp[bric_both & platform == "across"], "BRIC, RNA-seq vs microarray"),
  summ(dp[bric_both & seedling_both & platform == "within RNAseq"], "BRIC whole seedlings, RNA-seq vs RNA-seq"),
  summ(dp[bric_both & seedling_both & platform == "across"], "BRIC whole seedlings, RNA-seq vs microarray"),
  summ(dp[platform == "within RNAseq"], "all hardware, RNA-seq vs RNA-seq"),
  summ(dp[platform == "across"], "all hardware, RNA-seq vs microarray"),
  summ(P[same_study == TRUE], "same study, different organs")))
print(S); fwrite(S, file.path(OUT, "platform_concordance_summary.tsv"), sep = "\t")
writeLines(c(log_lines, "", capture.output(print(S)), "", capture.output(sessionInfo())), file.path(OUT, "platform_concordance_log.txt"))

# ---- REANALYSIS.md section 13 ---------------------------------------------------
row <- function(lbl) { x <- S[comparison == lbl]
  if (!nrow(x)) return(sprintf("| %s | 0 | — | — |", lbl))
  sprintf("| %s | %d | %.2f | %.2f to %.2f |", lbl, x$pairs, x$median_r, x$min_r, x$max_r) }
med <- function(lbl) { x <- S[comparison == lbl]; if (nrow(x)) x$median_r else NA_real_ }
w_rna <- med("BRIC, RNA-seq vs RNA-seq"); w_arr <- med("BRIC, microarray vs microarray"); acr <- med("BRIC, RNA-seq vs microarray")
seed_w <- med("BRIC whole seedlings, RNA-seq vs RNA-seq"); seed_a <- med("BRIC whole seedlings, RNA-seq vs microarray")
unit_tab <- c("| unit | platform | hardware | FLT / GC | genes FDR < 0.05 |", "|---|---|---|---|---|",
              sprintf("| %s | %s | %s | %d / %d | %s |", U$unit, U$platform, U$hardware, U$n_flt, U$n_gc, format(U$fdr05, big.mark = ",")))
sec <- c(
  "## 13. RNA-seq vs microarray: platform or study?",
  "",
  "`code/factorial_model/platform_concordance_2026-10.R`. Tables: `followups/platform_concordance_*`.",
  "",
  "§9 found the RNA-seq and microarray flight effects within BRIC uncorrelated (logFC r = −0.05). Platform is confounded with study and organ: the BRIC RNA-seq studies are whole seedlings, the BRIC microarray studies mostly hypocotyl, root, shoot or cell culture. Here the flight effect is estimated separately in every study × organ unit with at least 2 flight and 2 ground samples (RNA-seq: voom on expected counts; microarray: limma-trend on log2 values), and units from different studies are compared.",
  "",
  unit_tab,
  "",
  "Agreement of flight logFC between units from different studies (Pearson r across shared genes):",
  "",
  "| comparison | pairs | median r | range |",
  "|---|---|---|---|",
  row("BRIC, RNA-seq vs RNA-seq"), row("BRIC, microarray vs microarray"), row("BRIC, RNA-seq vs microarray"),
  row("BRIC whole seedlings, RNA-seq vs RNA-seq"), row("BRIC whole seedlings, RNA-seq vs microarray"),
  row("all hardware, RNA-seq vs RNA-seq"), row("all hardware, RNA-seq vs microarray"),
  "",
  sprintf("For reference, different organs within the same study agree at median r = %.2f. Spearman medians for the first three rows: %.2f, %.2f and %.2f (`platform_concordance_summary.tsv`).",
          med("same study, different organs"), S[comparison == "BRIC, RNA-seq vs RNA-seq"]$median_rho,
          S[comparison == "BRIC, microarray vs microarray"]$median_rho, S[comparison == "BRIC, RNA-seq vs microarray"]$median_rho),
  "",
  sprintf("- **Finding:** within BRIC, the RNA-seq studies share a reproducible flight response (median r %.2f between studies), but the microarray studies agree neither with each other (%.2f) nor with the RNA-seq studies (%.2f; whole seedlings only: %s). The microarray units also detect few flight genes (median %s at FDR < 0.05, against %s in the BRIC RNA-seq units). The platform split in §9 is therefore weak or absent flight signal in the microarray studies, not an opposite response; it can come from the platform's lower sensitivity and smaller groups, or from the earlier flights those studies represent, and this design cannot separate the two. Across hardware, even RNA-seq studies agree poorly (median r %.2f), in line with the study-specific flight response in §12.",
          w_rna, w_arr, acr, if (is.na(seed_a)) "no matched pairs" else sprintf("%.2f", seed_a),
          format(median(U[hardware == "BRIC" & platform == "microarray"]$fdr05), big.mark = ","),
          format(median(U[hardware == "BRIC" & platform == "RNAseq"]$fdr05), big.mark = ","),
          med("all hardware, RNA-seq vs RNA-seq")),
  "- Caveats: one unit per study × organ with small groups (some 3 vs 3 or 3 vs 9), so single-pair correlations are noisy; OSD-7 is APEX-TAGES hardware filed under BRIC (§2); OSD-678's samples with organ 'unknown' are treated as whole seedlings, as in the original pipeline; microarray logFCs are compressed relative to RNA-seq, which changes slopes but not correlations.",
  "")
rpath <- "tables/reanalysis_2026-10/REANALYSIS.md"
txt <- paste(readLines(rpath), collapse = "\n")
# Insert section 13 after any lower-numbered and before any higher-numbered section,
# removing only the previous copy of section 13 (up to the next "## <number>." heading).
txt <- sub("(?s)\n## 13\\. .*?(?=\n## [0-9]+\\. |\n## Not done here|$)", "", txt, perl = TRUE)
nums <- as.integer(sub("^\n## ([0-9]+)\\. .*", "\\1", regmatches(txt, gregexpr("\n## [0-9]+\\. [^\n]*", txt))[[1]]))
higher <- nums[nums > 13]
marker <- if (length(higher)) paste0("\n## ", min(higher), ". ") else "\n## Not done here"
block <- paste(sec, collapse = "\n")
txt <- if (grepl(marker, txt, fixed = TRUE)) sub(marker, paste0("\n", block, marker), txt, fixed = TRUE) else paste0(txt, "\n", block)
writeLines(txt, rpath)
say("updated REANALYSIS.md section 13")
