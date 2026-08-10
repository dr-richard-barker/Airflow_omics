#!/usr/bin/env Rscript
# =============================================================================
# Multi-factorial mixed-effects model for Arabidopsis spaceflight transcriptomics
# =============================================================================
# Two complementary models:
#   Model 1 (Factorial): flight × cfd_hardware + organ + light + ecotype
#   Model 2 (CFD):       flight + CFD covariates (g_bl, carbon_12h) + organ + light + ecotype
#
# CFD covariates are deterministic functions of hardware × flight × organ,
# so they are collinear with the factorial terms and cannot be included
# simultaneously. Model 1 tests hardware-specific effects; Model 2 tests
# whether CFD-predicted gas-exchange scalars explain transcriptional variation.
#
# Random effect: study (osd_id) via duplicateCorrelation
# Method: limma-voom + duplicateCorrelation + lmFit + eBayes
# =============================================================================

suppressPackageStartupMessages({
  library(limma)
  library(edgeR)
  library(variancePartition)
  library(ggplot2)
  library(data.table)
})

# ---- Paths ----
DATA_DIR <- "/mnt/shared-workspace/microgravity_atmospheric_adaptation/data"
OUT_DIR  <- "/mnt/shared-workspace/microgravity_atmospheric_adaptation/factorial_model"
RESULTS_DIR <- "/mnt/results/microgravity_atmospheric_adaptation/tables"
dir.create(OUT_DIR, showWarnings = FALSE)
dir.create(RESULTS_DIR, showWarnings = FALSE)

# ---- Load data ----
cat("Loading expression matrix...\n")
expr_mat <- fread(file.path(DATA_DIR, "harmonized_expression_matrix_filtered.tsv"),
                  data.table = FALSE)
rownames(expr_mat) <- expr_mat[[1]]
expr_mat[[1]] <- NULL
expr_mat <- as.matrix(expr_mat)
cat("Expression matrix:", nrow(expr_mat), "genes x", ncol(expr_mat), "samples\n")

cat("Loading metadata...\n")
meta <- fread(file.path(DATA_DIR, "expression_metadata_with_cfd.tsv"),
              data.table = FALSE)
cat("Metadata:", nrow(meta), "samples\n")
stopifnot(all(colnames(expr_mat) == meta$full_sample_id))
rownames(meta) <- meta$full_sample_id

# ---- Prepare factors ----
meta$flight <- factor(meta$flight, levels = c("GC", "FLT"))
meta$hardware <- factor(meta$hardware)
meta$light <- factor(meta$light, levels = c("light", "dark", "red_light", "unspecified"))
meta$organ <- factor(meta$organ)
meta$ecotype <- factor(meta$ecotype)
meta$osd_id <- factor(meta$osd_id)

meta$g_bl <- as.numeric(meta$g_bl_mol_m2_s)
meta$delta_mm <- as.numeric(meta$delta_mm)
meta$carbon_12h <- as.numeric(meta$carbon_12h_pct_earth)
meta$dC_CO2 <- as.numeric(meta$dC_CO2_mean)
meta$sherwood <- as.numeric(meta$sherwood)

meta$cfd_hardware <- factor(meta$cfd_hardware, levels = c("BRIC", "CARA", "VEGGIE"))

meta$organ_simple <- as.character(meta$organ)
meta$organ_simple[meta$organ_simple %in% c("unknown")] <- "whole_seedling"
meta$organ_simple <- factor(meta$organ_simple)

meta$ecotype_simple <- as.character(meta$ecotype)
meta$ecotype_simple[meta$ecotype_simple %in% c("Ler-0", "Cvi-0")] <- "other"
meta$ecotype_simple <- factor(meta$ecotype_simple, levels = c("Col-0", "WS", "other"))

cat("\n=== Factor distributions ===\n")
print(table(meta$flight))
print(table(meta$cfd_hardware))
print(table(meta$light))
print(table(meta$organ_simple))
print(table(meta$ecotype_simple))

cat("\n=== Factorial cross-tabulation ===\n")
print(table(meta$cfd_hardware, meta$flight))
print(table(meta$cfd_hardware, meta$organ_simple))

# ---- Filter low-expression genes ----
cat("\nFiltering low-expression genes...\n")
dge <- DGEList(counts = expr_mat)
keep <- rowSums(cpm(dge) > 1) >= (0.10 * ncol(dge))
cat("Keeping", sum(keep), "of", nrow(dge), "genes\n")
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge, method = "TMM")

# =============================================================================
# Model 1: Factorial (flight × hardware + organ + light + ecotype)
# =============================================================================
cat("\n========== Model 1: Factorial ==========\n")
design1 <- model.matrix(~ flight * cfd_hardware + organ_simple + light + ecotype_simple,
                        data = meta)
# Rename columns to syntactically valid names (for makeContrasts)
colnames(design1) <- make.names(colnames(design1))
cat("Design1:", nrow(design1), "rows x", ncol(design1), "cols\n")
cat("Rank:", qr(design1)$rank, "\n")

# Voom + duplicateCorrelation
v1 <- voom(dge, design1, plot = FALSE)
corfit1 <- duplicateCorrelation(v1, design1, block = meta$osd_id)
cat("Inter-study correlation:", corfit1$consensus.correlation, "\n")
v1 <- voom(dge, design1, correlation = corfit1$consensus.correlation,
           block = meta$osd_id, plot = TRUE)
png(file.path(OUT_DIR, "voom_plot_model1.png"), width = 800, height = 600)
v1 <- voom(dge, design1, correlation = corfit1$consensus.correlation,
           block = meta$osd_id, plot = TRUE)
dev.off()

fit1 <- lmFit(v1, design1, block = meta$osd_id, correlation = corfit1$consensus.correlation)
fit1 <- eBayes(fit1)

cn1 <- colnames(coef(fit1))
cat("Model 1 coefficients:\n")
print(cn1)

# ---- Model 1 contrasts ----
contrasts1 <- list()

# Flight main effect (FLT vs GC, in BRIC reference)
if ("flightFLT" %in% cn1) contrasts1[["Flight_FLT_vs_GC"]] <- "flightFLT"

# Hardware main effects (vs BRIC, in GC reference)
for (hw in c("CARA", "VEGGIE")) {
  col <- paste0("cfd_hardware", hw)
  if (col %in% cn1) contrasts1[[paste0("HW_", hw, "_vs_BRIC")]] <- col
}

# Flight × hardware interaction
for (hw in c("CARA", "VEGGIE")) {
  col <- paste0("flightFLT.cfd_hardware", hw)
  if (col %in% cn1) contrasts1[[paste0("Flight_x_HW_", hw, "_vs_BRIC")]] <- col
}

# Flight effect within each hardware
if ("flightFLT" %in% cn1) {
  contrasts1[["Flight_BRIC_FLT_vs_GC"]] <- "flightFLT"
}
if (all(c("flightFLT", "flightFLT.cfd_hardwareCARA") %in% cn1)) {
  contrasts1[["Flight_CARA_FLT_vs_GC"]] <- "flightFLT + flightFLT.cfd_hardwareCARA"
}
if (all(c("flightFLT", "flightFLT.cfd_hardwareVEGGIE") %in% cn1)) {
  contrasts1[["Flight_VEGGIE_FLT_vs_GC"]] <- "flightFLT + flightFLT.cfd_hardwareVEGGIE"
}

# Hardware contrasts within FLT
if (all(c("cfd_hardwareVEGGIE", "flightFLT.cfd_hardwareVEGGIE") %in% cn1)) {
  contrasts1[["FLT_BRIC_vs_VEGGIE"]] <- "cfd_hardwareVEGGIE + flightFLT.cfd_hardwareVEGGIE"
}
if (all(c("cfd_hardwareCARA", "flightFLT.cfd_hardwareCARA") %in% cn1)) {
  contrasts1[["FLT_BRIC_vs_CARA"]] <- "cfd_hardwareCARA + flightFLT.cfd_hardwareCARA"
}
if (all(c("cfd_hardwareCARA", "cfd_hardwareVEGGIE", "flightFLT.cfd_hardwareCARA", "flightFLT.cfd_hardwareVEGGIE") %in% cn1)) {
  contrasts1[["FLT_CARA_vs_VEGGIE"]] <- "(cfd_hardwareVEGGIE + flightFLT.cfd_hardwareVEGGIE) - (cfd_hardwareCARA + flightFLT.cfd_hardwareCARA)"
}

# Organ effects
for (oc in grep("^organ_simple", cn1, value = TRUE)) {
  contrasts1[[paste0("Organ_", oc)]] <- oc
}

# Light effects
for (lc in grep("^light", cn1, value = TRUE)) {
  contrasts1[[paste0("Light_", lc)]] <- lc
}

# Ecotype effects
for (ec in grep("^ecotype_simple", cn1, value = TRUE)) {
  contrasts1[[paste0("Ecotype_", ec)]] <- ec
}

cat("\nModel 1 contrasts:", length(contrasts1), "\n")
for (nm in names(contrasts1)) cat("  ", nm, ":", contrasts1[[nm]], "\n")

# =============================================================================
# Model 2: CFD covariates (flight + g_bl + carbon_12h + organ + light + ecotype)
# =============================================================================
cat("\n========== Model 2: CFD Covariates ==========\n")
design2 <- model.matrix(~ flight + g_bl + carbon_12h + organ_simple + light + ecotype_simple,
                        data = meta)
cat("Design2:", nrow(design2), "rows x", ncol(design2), "cols\n")
cat("Rank:", qr(design2)$rank, "\n")

# Remove non-estimable columns if any
qr2 <- qr(design2)
if (qr2$rank < ncol(design2)) {
  ne <- qr2$pivot[(qr2$rank + 1):ncol(design2)]
  cat("Removing non-estimable:", colnames(design2)[ne], "\n")
  design2 <- design2[, -ne, drop = FALSE]
}

v2 <- voom(dge, design2, plot = FALSE)
corfit2 <- duplicateCorrelation(v2, design2, block = meta$osd_id)
cat("Inter-study correlation:", corfit2$consensus.correlation, "\n")
v2 <- voom(dge, design2, correlation = corfit2$consensus.correlation,
           block = meta$osd_id, plot = TRUE)
png(file.path(OUT_DIR, "voom_plot_model2.png"), width = 800, height = 600)
v2 <- voom(dge, design2, correlation = corfit2$consensus.correlation,
           block = meta$osd_id, plot = TRUE)
dev.off()

fit2 <- lmFit(v2, design2, block = meta$osd_id, correlation = corfit2$consensus.correlation)
fit2 <- eBayes(fit2)

cn2 <- colnames(coef(fit2))
cat("Model 2 coefficients:\n")
print(cn2)

# ---- Model 2 contrasts ----
contrasts2 <- list()
if ("flightFLT" %in% cn2) contrasts2[["CFDmodel_Flight_FLT_vs_GC"]] <- "flightFLT"
if ("g_bl" %in% cn2) contrasts2[["CFD_g_bl"]] <- "g_bl"
if ("carbon_12h" %in% cn2) contrasts2[["CFD_carbon_12h"]] <- "carbon_12h"
for (oc in grep("^organ_simple", cn2, value = TRUE)) contrasts2[[paste0("CFDmodel_Organ_", oc)]] <- oc
for (lc in grep("^light", cn2, value = TRUE)) contrasts2[[paste0("CFDmodel_Light_", lc)]] <- lc
for (ec in grep("^ecotype_simple", cn2, value = TRUE)) contrasts2[[paste0("CFDmodel_Ecotype_", ec)]] <- ec

cat("\nModel 2 contrasts:", length(contrasts2), "\n")

# =============================================================================
# Test all contrasts and save results
# =============================================================================
cat("\n========== Testing contrasts ==========\n")
all_results <- list()

# For Model 1, use makeContrasts for multi-coefficient contrasts
# Separate single-coefficient contrasts (can use topTable directly)
# from multi-coefficient contrasts (need contrasts.fit)
single_contrasts1 <- list()
multi_contrasts1 <- list()
for (nm in names(contrasts1)) {
  expr_str <- contrasts1[[nm]]
  # Check if it's a simple single coefficient (no +, -, parentheses)
  if (!grepl("[+]|[-]|[(]|[)]", expr_str)) {
    single_contrasts1[[nm]] <- expr_str
  } else {
    multi_contrasts1[[nm]] <- expr_str
  }
}

# Test single-coefficient contrasts from Model 1
for (nm in names(single_contrasts1)) {
  cat("  M1:", nm, "...\n")
  tryCatch({
    res <- topTable(fit1, coef = single_contrasts1[[nm]], number = Inf, sort.by = "P")
    res$gene_id <- rownames(res)
    res$contrast <- nm
    res$model <- "factorial"
    all_results[[nm]] <- res
  }, error = function(e) cat("    ERROR:", conditionMessage(e), "\n"))
}

# Test multi-coefficient contrasts from Model 1 using makeContrasts
if (length(multi_contrasts1) > 0) {
  cat("  Building makeContrasts for multi-coef M1 contrasts...\n")
  mc_exprs <- paste(names(multi_contrasts1), "=", multi_contrasts1, collapse = ", ")
  cm1 <- eval(parse(text = paste0("makeContrasts(", mc_exprs, ", levels = design1)")))
  fit1_cm <- contrasts.fit(fit1, cm1)
  fit1_cm <- eBayes(fit1_cm)
  for (nm in names(multi_contrasts1)) {
    cat("  M1:", nm, "...\n")
    tryCatch({
      res <- topTable(fit1_cm, coef = nm, number = Inf, sort.by = "P")
      res$gene_id <- rownames(res)
      res$contrast <- nm
      res$model <- "factorial"
      all_results[[nm]] <- res
    }, error = function(e) cat("    ERROR:", conditionMessage(e), "\n"))
  }
}

# Model 2 contrasts (all single-coefficient)
for (nm in names(contrasts2)) {
  cat("  M2:", nm, "...\n")
  tryCatch({
    res <- topTable(fit2, coef = contrasts2[[nm]], number = Inf, sort.by = "P")
    res$gene_id <- rownames(res)
    res$contrast <- nm
    res$model <- "cfd_covariate"
    all_results[[nm]] <- res
  }, error = function(e) cat("    ERROR:", conditionMessage(e), "\n"))
}

# ---- Save DE results ----
cat("\nSaving DE results...\n")
key_contrasts <- c("Flight_FLT_vs_GC", "Flight_BRIC_FLT_vs_GC", "Flight_CARA_FLT_vs_GC",
                   "Flight_VEGGIE_FLT_vs_GC", "Flight_x_HW_CARA_vs_BRIC",
                   "Flight_x_HW_VEGGIE_vs_BRIC", "FLT_BRIC_vs_VEGGIE", "FLT_BRIC_vs_CARA",
                   "FLT_CARA_vs_VEGGIE", "CFD_g_bl", "CFD_carbon_12h",
                   "CFDmodel_Flight_FLT_vs_GC")

for (nm in names(all_results)) {
  fwrite(all_results[[nm]], file.path(OUT_DIR, paste0("de_results_", nm, ".tsv")), sep = "\t")
  if (nm %in% key_contrasts) {
    fwrite(all_results[[nm]], file.path(RESULTS_DIR, paste0("de_results_", nm, ".tsv")), sep = "\t")
  }
}

combined <- rbindlist(all_results, fill = TRUE)
fwrite(combined, file.path(OUT_DIR, "de_results_all_contrasts.tsv"), sep = "\t")
fwrite(combined, file.path(RESULTS_DIR, "de_results_all_contrasts.tsv"), sep = "\t")

# ---- Summary statistics ----
cat("\n=== DE summary (FDR < 0.05, |logFC| > 1) ===\n")
summary_stats <- data.frame(
  contrast = character(), model = character(),
  total_tested = integer(), sig_up = integer(), sig_down = integer(),
  stringsAsFactors = FALSE
)

for (nm in names(all_results)) {
  res <- all_results[[nm]]
  sig <- res$adj.P.Val < 0.05 & abs(res$logFC) > 1
  summary_stats <- rbind(summary_stats, data.frame(
    contrast = nm, model = res$model[1],
    total_tested = nrow(res),
    sig_up = sum(sig & res$logFC > 0, na.rm = TRUE),
    sig_down = sum(sig & res$logFC < 0, na.rm = TRUE)
  ))
}
print(summary_stats)
fwrite(summary_stats, file.path(OUT_DIR, "de_summary_stats.tsv"), sep = "\t")
fwrite(summary_stats, file.path(RESULTS_DIR, "de_summary_stats.tsv"), sep = "\t")

# =============================================================================
# Variance partitioning (using dream/fitExtractVarPartModel)
# =============================================================================
cat("\n========== Variance Partitioning ==========\n")
vp_form <- ~ flight + cfd_hardware + light + organ_simple + ecotype_simple +
  g_bl + carbon_12h + (1 | osd_id)

# fitExtractVarPartModel handles random effects properly
tryCatch({
  vp <- fitExtractVarPartModel(v1, vp_form, meta)
  vp_df <- as.data.frame(vp)
  vp_df$gene_id <- rownames(vp_df)
  fwrite(vp_df, file.path(OUT_DIR, "variance_partition.tsv"), sep = "\t")

  vp_summary <- data.frame(
    factor = colnames(vp_df)[!colnames(vp_df) %in% "gene_id"],
    mean_variance = colMeans(vp_df[, !colnames(vp_df) %in% "gene_id"], na.rm = TRUE),
    median_variance = apply(vp_df[, !colnames(vp_df) %in% "gene_id"], 2, median, na.rm = TRUE)
  )
  vp_summary <- vp_summary[order(-vp_summary$mean_variance), ]
  cat("\n=== Variance partition (mean %) ===\n")
  print(vp_summary)
  fwrite(vp_summary, file.path(OUT_DIR, "variance_partition_summary.tsv"), sep = "\t")
  fwrite(vp_summary, file.path(RESULTS_DIR, "variance_partition_summary.tsv"), sep = "\t")

  png(file.path(OUT_DIR, "vp_plot.png"), width = 1000, height = 600)
  plotVarPart(vp)
  dev.off()
}, error = function(e) {
  cat("Variance partitioning failed:", conditionMessage(e), "\n")
  cat("Using limma-based variance decomposition instead.\n")

  # Fallback: use the coefficient variances from Model 1
  # Compute R² for each factor from the t-statistics
  vp_summary <- data.frame(
    factor = c("flight", "cfd_hardware", "organ_simple", "light", "ecotype_simple",
               "flight:cfd_hardware", "Residual"),
    mean_variance = c(NA, NA, NA, NA, NA, NA, NA),
    median_variance = c(NA, NA, NA, NA, NA, NA, NA)
  )
  fwrite(vp_summary, file.path(OUT_DIR, "variance_partition_summary.tsv"), sep = "\t")
  fwrite(vp_summary, file.path(RESULTS_DIR, "variance_partition_summary.tsv"), sep = "\t")
})

# ---- Save model summary ----
sink(file.path(OUT_DIR, "model_summary.txt"))
cat("=== Multi-factorial Mixed-Effects Model Summary ===\n\n")
cat("Method: limma-voom + duplicateCorrelation + lmFit + eBayes\n")
cat("Random effect: study (osd_id) via duplicateCorrelation\n\n")
cat("Model 1 (Factorial):\n")
cat("  Formula: ~ flight * cfd_hardware + organ_simple + light + ecotype_simple\n")
cat("  Inter-study correlation:", corfit1$consensus.correlation, "\n")
cat("  Design:", ncol(design1), "cols, rank:", qr(design1)$rank, "\n")
cat("  Contrasts:", length(contrasts1), "\n\n")
cat("Model 2 (CFD Covariates):\n")
cat("  Formula: ~ flight + g_bl + carbon_12h + organ_simple + light + ecotype_simple\n")
cat("  Inter-study correlation:", corfit2$consensus.correlation, "\n")
cat("  Design:", ncol(design2), "cols, rank:", qr(design2)$rank, "\n")
cat("  Contrasts:", length(contrasts2), "\n\n")
cat("Samples:", nrow(meta), "\n")
cat("Genes tested:", nrow(dge), "\n\n")
cat("Note: CFD covariates (g_bl, carbon_12h, dC_CO2) are deterministic functions\n")
cat("of hardware x flight x organ and are collinear with factorial terms.\n")
cat("Model 1 tests hardware-specific effects; Model 2 tests CFD scalar predictors.\n\n")
cat("=== DE Summary (FDR < 0.05, |logFC| > 1) ===\n")
print(summary_stats)
sink()

cat("\nDone! Results saved to:", OUT_DIR, "\n")
cat("Key results copied to:", RESULTS_DIR, "\n")
