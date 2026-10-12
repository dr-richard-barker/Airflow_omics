#!/usr/bin/env Rscript
# =============================================================================
# Variance partitioning, redone (REANALYSIS §12)
# =============================================================================
# The stored tables/variance_partition_summary.tsv (fit_factorial_model.R) entered
# categorical factors as fixed effects, which variancePartition does not partition
# correctly, and the summary is internally inconsistent (cfd_hardware median 66 %
# with mean ~0 %; no residual, g_bl or carbon_12h rows). Here every categorical
# factor is a random effect, as the variancePartition authors recommend, and the
# continuous CFD covariates are left out because they are functions of
# hardware x flight x scale (§4). Study is nested in hardware (each study used one
# hardware), so both enter as nested components.
#
#   ~ (1|cfd_hardware) + (1|osd_id) + (1|flight) + (1|flight:cfd_hardware)
#     + (1|flight:osd_id) + (1|organ_simple) + (1|light) + (1|ecotype_simple)
#
# Data: (a) pooled voom on the genes measured in all 11 studies (as §10);
#       (b) RNA-seq studies only, voom on RSEM expected counts (as §9).
# Usage (from the repo root): Rscript code/factorial_model/variance_partition_2026-10.R
# =============================================================================

suppressPackageStartupMessages({
  library(limma)
  library(edgeR)
  library(data.table)
  library(variancePartition)
  library(BiocParallel)
})

RAW <- path.expand(Sys.getenv("OSDR_RAW", "~/Documents/Airflow_omics_osdr_data/raw_osdr"))
OUT <- "tables/reanalysis_2026-10/followups"
NCORES <- as.integer(Sys.getenv("VP_CORES", max(1, parallel::detectCores() - 2)))
log_lines <- character()
say <- function(...) { s <- paste0(...); cat(s, "\n"); log_lines <<- c(log_lines, s) }
source("code/factorial_model/osdr_helpers_2026-10.R")
param <- MulticoreParam(NCORES, progressbar = FALSE)

meta2 <- fread("tables/reanalysis_2026-10/expression_metadata_with_cfd_v2.tsv", data.table = FALSE)
m2 <- prep_meta(meta2)
m2$organ_simple <- factor(m2$organ_simple); m2$ecotype_simple <- factor(m2$ecotype_simple)

# collinearity of the variables in the original formula (canonical correlations)
cc_form <- ~ flight + cfd_hardware + light + organ_simple + ecotype_simple + g_bl + carbon_12h + osd_id
C <- canCorPairs(cc_form, m2)
fwrite(data.table(variable = rownames(C), round(C, 3)), file.path(OUT, "variance_partition_cancor.tsv"), sep = "\t")
say("== canonical correlations, original formula variables ==")
say(sprintf("  cfd_hardware~osd_id %.2f; carbon_12h~osd_id %.2f; g_bl~flight %.2f; carbon_12h~cfd_hardware %.2f",
            C["cfd_hardware", "osd_id"], C["carbon_12h", "osd_id"], C["g_bl", "flight"], C["carbon_12h", "cfd_hardware"]))

vp_form <- ~ (1 | cfd_hardware) + (1 | osd_id) + (1 | flight) + (1 | flight:cfd_hardware) +
  (1 | flight:osd_id) + (1 | organ_simple) + (1 | light) + (1 | ecotype_simple)
voom_design <- function(meta) drop_aliased(model.matrix(F1, meta))

MAXG <- as.integer(Sys.getenv("VP_MAX_GENES", "0"))   # > 0 only for smoke tests
run_vp <- function(label, dge, meta) {
  if (MAXG > 0) dge <- dge[seq_len(min(MAXG, nrow(dge))), ]
  v <- voom(dge, voom_design(meta))
  t0 <- Sys.time()
  vp <- fitExtractVarPartModel(v, vp_form, meta, BPPARAM = param)
  say(sprintf("  %s: %d genes x %d samples, %.1f min", label, nrow(v), ncol(v), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  vp <- as.data.frame(vp) * 100
  fwrite(data.table(gene_id = rownames(vp), vp), file.path(OUT, paste0("variance_partition_", label, ".tsv.gz")), sep = "\t")
  flight_total <- vp[["flight"]] + vp[["flight:cfd_hardware"]] + vp[["flight:osd_id"]]
  s <- data.table(component = colnames(vp), mean_pct = colMeans(vp), median_pct = apply(vp, 2, median),
                  genes_over_10pct = colSums(vp > 10))
  s <- rbind(s, data.table(component = "all flight terms", mean_pct = mean(flight_total),
                           median_pct = median(flight_total), genes_over_10pct = sum(flight_total > 10)))
  s[, data := label][]
}

say("\n== variance partition ==")
H <- outer_merge(c(mapply(load_rnaseq, names(rnaseq_files), rnaseq_files, SIMPLIFY = FALSE),
                   mapply(load_array, names(array_files), array_files, SIMPLIFY = FALSE)))[, m2$full_sample_id]
CC <- H[rowSums(is.na(H)) == 0, ]
d_all <- calcNormFactors(DGEList(counts = CC)[cpm_keep(CC), , keep.lib.sizes = FALSE], method = "TMM")
res_a <- run_vp("pooled_complete", d_all, m2)

rna <- m2$omics_type == "RNAseq"; mr <- droplevels(m2[rna, ])
raw <- outer_merge(mapply(load_rnaseq, names(rnaseq_raw_files), rnaseq_raw_files, SIMPLIFY = FALSE))
R <- round(raw[rowSums(is.na(raw)) == 0, rownames(mr)]); R <- R[grepl("^AT[1-5CM]G[0-9]{5}", rownames(R)), ]
d_rna <- DGEList(counts = R)
d_rna <- calcNormFactors(d_rna[filterByExpr(d_rna, group = interaction(mr$cfd_hardware, mr$flight)), , keep.lib.sizes = FALSE])
res_b <- run_vp("rnaseq", d_rna, mr)

summ <- rbind(res_a, res_b)
fwrite(summ, file.path(OUT, "variance_partition_summary_v2.tsv"), sep = "\t")
print(summ)
writeLines(c(log_lines, "", capture.output(sessionInfo())), file.path(OUT, "variance_partition_log.txt"))

# ---- REANALYSIS.md section 12 ------------------------------------------------
order_c <- c("osd_id", "cfd_hardware", "organ_simple", "light", "ecotype_simple", "flight",
             "flight:cfd_hardware", "flight:osd_id", "all flight terms", "Residuals")
lab <- c(osd_id = "study (within hardware)", cfd_hardware = "hardware", organ_simple = "organ", light = "light",
         ecotype_simple = "ecotype", flight = "flight (common to all studies)", `flight:cfd_hardware` = "flight × hardware",
         `flight:osd_id` = "flight × study", `all flight terms` = "**all flight terms**", Residuals = "residual")
cell <- function(d, comp) { r <- summ[data == d & component == comp]; sprintf("%.1f / %.1f", r$mean_pct, r$median_pct) }
pct <- function(d, comp) summ[data == d & component == comp]$mean_pct
sec <- c(
  "## 12. Variance partitioning, redone",
  "",
  "`code/factorial_model/variance_partition_2026-10.R`. Tables: `followups/variance_partition_*`.",
  "",
  "The stored `tables/variance_partition_summary.tsv` should not be used. `fit_factorial_model.R` entered the categorical factors as fixed effects, which variancePartition does not partition correctly, and the stored summary is internally inconsistent: `cfd_hardware` has a median of 66 % but a mean of ~0 %, which is impossible for non-negative fractions, and there are no residual, g_bl or carbon_12h rows.",
  "",
  sprintf("Here every categorical factor is a random effect, and hardware and study enter as nested components (each study used one hardware; canonical correlation hardware ~ study = %.2f). g_bl and carbon_12h are left out: they are functions of hardware × flight × scale (canonical correlation carbon_12h ~ study = %.2f).",
          C["cfd_hardware", "osd_id"], C["carbon_12h", "osd_id"]),
  "",
  "Percent of each gene's variance, mean / median across genes:",
  "",
  sprintf("| component | pooled, genes in all studies (%s genes, 307 samples) | RNA-seq only (%s genes, %d samples) |",
          format(nrow(d_all), big.mark = ","), format(nrow(d_rna), big.mark = ","), ncol(d_rna)),
  "|---|---|---|",
  sapply(order_c, function(cmp) sprintf("| %s | %s | %s |", lab[[cmp]], cell("pooled_complete", cmp), cell("rnaseq", cmp))),
  "",
  sprintf("- Study and hardware together account for %.0f %% (pooled) and %.0f %% (RNA-seq) of the average gene's variance; all flight terms together for %.0f %% and %.0f %%.",
          pct("pooled_complete", "osd_id") + pct("pooled_complete", "cfd_hardware"), pct("rnaseq", "osd_id") + pct("rnaseq", "cfd_hardware"),
          pct("pooled_complete", "all flight terms"), pct("rnaseq", "all flight terms")),
  sprintf("- Of the flight variance, the part shared by all studies is %.1f %% (pooled) and %.1f %% (RNA-seq); the study-specific part (flight × study) is %.1f %% and %.1f %%. %s",
          pct("pooled_complete", "flight"), pct("rnaseq", "flight"), pct("pooled_complete", "flight:osd_id"), pct("rnaseq", "flight:osd_id"),
          if (pct("pooled_complete", "flight:osd_id") > pct("pooled_complete", "flight") && pct("rnaseq", "flight:osd_id") > pct("rnaseq", "flight"))
            "The flight response is mostly study-specific, which is why pooled flight contrasts are sensitive to which studies are included (§9)."
          else "The common flight component is at least as large as the study-specific one in at least one data set."),
  sprintf("- Genes with more than 10 %% of their variance in the flight effect common to all studies: %d of %s (pooled) and %d of %s (RNA-seq); in flight × study: %s and %s.",
          summ[data == "pooled_complete" & component == "flight"]$genes_over_10pct, format(nrow(d_all), big.mark = ","),
          summ[data == "rnaseq" & component == "flight"]$genes_over_10pct, format(nrow(d_rna), big.mark = ","),
          format(summ[data == "pooled_complete" & component == "flight:osd_id"]$genes_over_10pct, big.mark = ","),
          format(summ[data == "rnaseq" & component == "flight:osd_id"]$genes_over_10pct, big.mark = ",")),
  sprintf("- Hardware and flight × hardware are estimated from only %d hardware groups, so their split from the study components is imprecise; read study + hardware, and the three flight terms, as totals.", nlevels(droplevels(m2$cfd_hardware))),
  "")
rpath <- "tables/reanalysis_2026-10/REANALYSIS.md"
txt <- paste(readLines(rpath), collapse = "\n")
# Insert section 12 after any lower-numbered and before any higher-numbered section,
# removing only the previous copy of section 12 (up to the next "## <number>." heading).
txt <- sub("(?s)\n## 12\\. .*?(?=\n## [0-9]+\\. |\n## Not done here|$)", "", txt, perl = TRUE)
nums <- as.integer(sub("^\n## ([0-9]+)\\. .*", "\\1", regmatches(txt, gregexpr("\n## [0-9]+\\. [^\n]*", txt))[[1]]))
higher <- nums[nums > 12]
marker <- if (length(higher)) paste0("\n## ", min(higher), ". ") else "\n## Not done here"
block <- paste(sec, collapse = "\n")
txt <- if (grepl(marker, txt, fixed = TRUE)) sub(marker, paste0("\n", block, marker), txt, fixed = TRUE) else paste0(txt, "\n", block)
if (MAXG == 0) { writeLines(txt, rpath); say("updated REANALYSIS.md section 12") } else cat(block, "\n")
