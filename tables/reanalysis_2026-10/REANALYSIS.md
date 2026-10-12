# CFD covariate update and pathway re-analysis

307 samples, 11 studies, 28 pathways.

## 1. Mapping check (old LunarLeaf tables -> stored covariates)

Rebuilt the stored covariates from LunarLeaf `lunarleaf_tables_b3567f6`: 25/30 columns reproduce exactly.
Columns that do not reproduce:

| column | n_mismatch |
|---|---|
| carbon_12h_pct_earth | 106 |
| A_start_umol_m2_s | 106 |
| A_end_umol_m2_s | 106 |
| earth_equiv_fan_cm_s | 136 |
| fan_cm_s | 35 |

All expected, from deliberate changes to extract_cfd.py:

- `carbon_12h_pct_earth`: BRIC ground controls use the same sealed hardware, so they get T11's BRIC value, not 100 % (2026-10-08). The CFD class BRIC includes APEX-TAGES (OSD-7), whose 18 ground controls get the same value; TODO: confirm they were sealed
- `A_start_umol_m2_s`: as carbon_12h_pct_earth: BRIC ground controls now carry T11's BRIC values (2026-10-08)
- `A_end_umol_m2_s`: as carbon_12h_pct_earth: BRIC ground controls now carry T11's BRIC values (2026-10-08)
- `earth_equiv_fan_cm_s`: read from T9 instead of hard-coded 2.8 / 11.0 cm/s (2026-10-07)
- `fan_cm_s`: VEGGIE fan is the modelled 8.3 cm/s, not the Earth-equivalent 2.8 (2026-10-07)

Every other column reproduces, which validates the mapping.

## 2. Covariates that change with the corrected LunarLeaf tables

| column | n_mismatch |
|---|---|
| dC_CO2_peak | 307 |
| earth_equiv_fan_cm_s | 307 |
| dC_H2O_peak | 239 |
| dC_O2_peak | 239 |
| g_bl_mol_m2_s | 221 |
| sherwood | 221 |
| dC_H2O_mean | 153 |
| o2_excess_ppm_bl | 153 |
| dC_CO2_mean | 153 |
| delta_mm | 153 |
| dC_O2_mean | 153 |
| dC_CO2_mean_bl | 153 |
| carbon_12h_pct_earth | 106 |
| A_start_umol_m2_s | 106 |
| A_end_umol_m2_s | 106 |
| u_max | 86 |
| surf_dC_H2O_hw | 35 |
| surf_dC_CO2_hw | 35 |
| fan_cm_s | 35 |

g_bl per cell, old -> new: earth/leaf 1.000 -> 0.997; earth/rosette 0.540 -> 0.540; ug/leaf 0.494 -> 0.443; ug/rosette 0.291 -> 0.272

Not updated because their LunarLeaf source tables were not regenerated: co2_drawdown_* / o2_buildup_* (T5), net_assimilation_pct (T12).

OSDR hardware -> CFD class (`extract_cfd.py`). APEX-TAGES is folded into BRIC and EMCS into VEGGIE, so every CFD covariate for them, including the ground-control carbon value, is the class value:

| hardware | cfd_hardware | studies | n |
|---|---|---|---|
| APEX-TAGES | BRIC | OSD-7 | 36 |
| BRIC | BRIC | OSD-147, OSD-17, OSD-205, OSD-321, OSD-37, OSD-38 | 170 |
| CARA | CARA | OSD-120 | 36 |
| EMCS | VEGGIE | OSD-314 | 17 |
| Veggie | VEGGIE | OSD-217, OSD-678 | 48 |

TODO: confirm that the APEX-TAGES (OSD-7) ground controls were sealed like BRIC's, and whether EMCS (OSD-314) is adequately represented by the vented VEGGIE scenario.

## 3. Original method (Pearson, n = 307, BH)

Reproduces the stored `pathway_cfd_correlations.tsv`: 84/84 rows matched, max |Δr| = 1.11e-16.

- `g_bl_mol_m2_s`: FDR < 0.05 in 28/28 pathways with old covariates, 28/28 with new; max |r| 0.467 -> 0.431.
- `carbon_12h_pct_earth`: FDR < 0.05 in 27/28 pathways with old covariates, 27/28 with new; max |r| 0.294 -> 0.585.
- `delta_mm`: FDR < 0.05 in 28/28 pathways with old covariates, 28/28 with new; max |r| 0.460 -> 0.422.

`carbon_12h_pct_earth` by hardware and flight (new covariates): BRIC FLT 1 / GC 1; CARA FLT 90 / GC 100; VEGGIE FLT 100 / GC 100. Within BRIC it no longer differs between flight and ground, so its correlation with pathway scores now reflects which studies used BRIC (a between-study, hardware difference), not a flight effect.

## 4. Is g_bl identifiable separately from flight and scale?

| gravity_key | cfd_scale | 0.272 | 0.443 | 0.54 | 0.997 |
|---|---|---|---|---|---|
| earth | leaf | 0 | 0 | 0 | 68 |
| earth | rosette | 0 | 0 | 86 | 0 |
| ug | leaf | 0 | 68 | 0 | 0 |
| ug | rosette | 85 | 0 | 0 | 0 |

Design rank: study + flight × scale = 14 columns; adding g_bl gives rank 14 of 15 columns. g_bl is a deterministic relabelling of the flight × scale cells, so no model can attribute an effect to g_bl rather than to flight and organ scale in this design. Its correlations with pathways are flight × scale contrasts under another name.

## 5. Within-study model per pathway

`score ~ study + flight × scale` (effect coding), CR1 cluster-robust SEs by study (11 clusters, t on G − 1 df), BH across pathways per term.

- **flight**: FDR < 0.05 in 0/28 pathways.
- **scale**: FDR < 0.05 in 0/28 pathways.
- **flight_x_scale**: FDR < 0.05 in 0/28 pathways.

For comparison, the same model with naive (independent-sample) SEs calls the flight term significant in 7/28 pathways at Bonferroni 0.05: the original significance comes from treating 307 samples as independent when they cluster in 11 studies.
Positive control: a flight effect planted in every pathway is detected (p < 0.05) in 28/28 at 0.25 SD, 28/28 at 0.5 SD, so the null result is not a powerless test; the pathway-level flight effects are small or inconsistent across studies relative to between-study variation.

Of the 28 pathways 'significant' for g_bl by Pearson, 0 have any FDR < 0.05 term (flight, scale or interaction) once within-study structure and clustering are respected.

## 6. Gene-level checks on the stored DE tables

- `Flight_FLT_vs_GC` vs `Flight_BRIC_FLT_vs_GC`: 25453 genes compared; logFC identical: **True** (max |Δ| = 0.00e+00). With BRIC as the reference level, the 'Flight' coefficient is the flight effect within BRIC, not a main effect across hardware.
- `CFD_g_bl` (logFC per 1 mol m⁻² s⁻¹; the covariate spans only 0.709, SD 0.253): |logFC| > 1 & FDR < 0.05 → 3633 genes. The same genes per 1 SD of g_bl: |logFC·SD| > 1 → 813; per full observed range: 3010. The cutoff counts depend on the unit chosen, so for a continuous covariate report FDR plus an effect size per SD or per observed range.
## 7. In-dish covariate (LunarLeaf T19)

`g_bl_dish` (leaf conductance to the dish's own air, in the real BRIC 60 mm / CARA 120 mm dish) is defined for 112 leaf samples in BRIC or CARA hardware from 5 studies. It is NaN for rosette-scale and VEGGIE samples, so it is kept as a separate column rather than replacing `g_bl`.

| cfd_hardware | gravity_key | n | studies | g_bl_dish | g_bl_open_box | note |
|---|---|---|---|---|---|---|
| BRIC | earth | 38 | 4 | 0.33 | 0.997 | indicative: sealed 1 g not steady, beyond Boussinesq limit |
| BRIC | ug | 38 | 4 | 0.185 | 0.443 | steady |
| CARA | earth | 18 | 1 | 0.263 | 0.997 | steady |
| CARA | ug | 18 | 1 | 0.186 | 0.443 | steady |

Again one value per hardware × gravity cell, so it is identifiable only as hardware × flight. CARA comes from OSD-120 alone, so a CARA-vs-BRIC difference in the flight response is also a difference between studies.

Predicted in-dish penalty (1 g − 0 g of g_bl_dish): BRIC 0.145 (its 1 g value is indicative only), CARA 0.077. In the open box it was 0.554 for every hardware.

Within-study model on this subset (`score ~ study + flight + flight × CARA`, CR1 by study, 5 clusters, t on 4 df): flight effect FDR < 0.05 in 0/28 pathways. CARA-vs-BRIC difference in the flight response: **not testable**: CARA is a single study, and cluster-robust errors for a term that varies within one cluster collapse toward zero (an unguarded run reports it 'significant' in most pathways, which is an artefact).

**Reading:** the dish geometry changes the *size* of the predicted gravity penalty (smaller than in the open box, and different for BRIC and CARA), but it adds no identifiable covariate. In this design it is still a hardware × flight label. Testing it as a dose-response needs study-specific dish scenarios across more studies per hardware.

## 8. Meta-regression: does the predicted boundary-layer penalty explain which studies respond?

14 study x scale units from 11 studies; the predictor takes 5 distinct values (it varies only by hardware x scale).

| unit | hardware | FLT/GC | predicted penalty (|ΔC CO₂| µg − Earth) |
|---|---|---|---|
| OSD-120|leaf | CARA | 18/18 | 0.145 |
| OSD-147|leaf | BRIC | 9/9 | 0.151 |
| OSD-17|leaf | BRIC | 6/6 | 0.151 |
| OSD-17|rosette | BRIC | 6/6 | 0.279 |
| OSD-205|leaf | BRIC | 8/8 | 0.151 |
| OSD-217|leaf | VEGGIE | 6/6 | -0.074 |
| OSD-314|rosette | VEGGIE | 11/6 | -0.004 |
| OSD-321|rosette | BRIC | 22/22 | 0.279 |
| OSD-37|rosette | BRIC | 28/28 | 0.279 |
| OSD-38|rosette | BRIC | 3/9 | 0.279 |
| OSD-678|leaf | VEGGIE | 6/6 | -0.074 |
| OSD-678|rosette | VEGGIE | 12/12 | -0.004 |
| OSD-7|leaf | BRIC | 15/15 | 0.151 |
| OSD-7|rosette | BRIC | 3/3 | 0.279 |

- Knapp–Hartung (primary): flight effect scales with the predicted penalty at FDR < 0.05 in **1/28** pathways (smallest p = 0.000822, ga_biosynthesis). The DerSimonian–Laird z-test, which is over-confident with this few units, would have reported 14/28 (smallest p = 1.1e-12).
- Of the FDR-significant pathways, 1 are also significant with a plain VEGGIE-vs-other indicator, and 0 show the slope within BRIC units alone (where the penalty varies only by organ scale). Leave-one-study-out: 1 stay p < 0.05 whichever study is dropped.
- **Reading:** the signal is a VEGGIE-versus-sealed-hardware difference in the flight response, consistent with the boundary-layer physics but not separable from which studies (labs, years, organs) used VEGGIE. There is no evidence yet of a within-hardware dose-response to the predicted penalty.
- Positive control: a planted effect of 1 SD of the unit-level effects per SD of penalty is detected (p < 0.05) in 24/28 pathways.
- Caveats: few units, the penalty varies only by hardware × scale (so it is also confounded with which studies used which hardware), and T10's BRIC/CARA values are the original 3.8 s values, not steady state.

## 9. Gene-level limma refit with the v2 covariates

`code/factorial_model/refit_gene_level_2026-10.R`, from the 11 OSDR processed files (7 RNA-seq, 4 microarray; download manifest with SHA-256 in `Airflow_omics_osdr_data/raw_osdr/MANIFEST.json`). Gene tables: `gene_level/`.

**Reconstructing the missing filtering step.** No script in the package writes `harmonized_expression_matrix_filtered.tsv`. Of two candidate rules, only one reproduces the stored fit:

| rule | genes after CPM filter | inter-study ρ (Model 1) | t vs stored Flight_FLT_vs_GC |
|---|---|---|---|
| drop genes missing in any study | 17,129 | 0.6412 | 0.9917 |
| set missing values to 0 | 25,453 | 0.7327 | 0.9988 |
| stored fit | 25,453 | 0.7324 | — |

Model 2 reproduces too (ρ 0.7767 vs stored 0.7782; t-statistic correlation with the stored g_bl, carbon_12h and flight tables 0.996–0.998).
So the original pipeline **set unmeasured genes to 0**: of 30,833 genes in the outer join, those measured only by RNA-seq enter every microarray sample as zero expression (and genes absent from one RNA-seq study are zero there). This builds study × platform differences into the matrix. It also passes microarray log2 intensities through `DGEList`/`voom` as if they were counts.

**Pooled refit (same pipeline, v2 covariates).** Genes at FDR < 0.05 (|logFC| > 1 in brackets):

| contrast | v1 (stored) | v2 refit |
|---|---|---|
| flight within BRIC (stored name `Flight_FLT_vs_GC`) | 11,783 | 11,675 (715) |
| **flight averaged over BRIC/CARA/VEGGIE** | not fitted | 4,598 (285) |
| flight × CARA | 6,909 | 6,822 |
| flight × VEGGIE | 10,904 | 10,491 |
| Model 2: g_bl | 4,020 | 8,954 |
| Model 2: carbon_12h | 9,561 | 8,449 |

- The within-BRIC and averaged flight effects agree only weakly (logFC correlation 0.352). The stored headline is the BRIC-specific effect.
- g_bl: VIF 13.6 in Model 2. It takes 4 values and recodes flight × scale (§4), so its gene count is not a boundary-layer effect.
- carbon_12h: study membership explains 99.9 % of its variance (v2 sets every BRIC sample to 1 %). Its genes are a BRIC-studies-versus-the-rest contrast, not a carbon dose response. The |logFC| > 1 count is 0 because the unit is 1 % of Earth carbon.

**Platform-separated sensitivity.** RNA-seq: voom on RSEM expected counts, `filterByExpr`, no zero-filling. Microarray: limma-trend on log2 intensities. All microarray samples are BRIC, so hardware terms and carbon_12h drop out there.

| contrast | RNA-seq | microarray |
|---|---|---|
| samples / studies | 213 / 7 | 94 / 4 |
| flight within BRIC | 10,995 | 115 |
| flight averaged over hardware | 3,258 | — |
| g_bl | 7,460 | 656 |
| carbon_12h | 2,598 | not estimable |

- Within BRIC, the RNA-seq and microarray flight effects are uncorrelated (logFC r = -0.050), and the pooled within-BRIC effect does not track the microarray one (r = -0.022). The pooled averaged effect does track the RNA-seq-only one (r = 0.849).
- **Reading:** the gene level agrees with §5 and §8. Neither CFD covariate is identifiable from study and hardware, so neither gene count supports a gas-exchange interpretation. The flight response differs by hardware and by platform/study era, and pooling with zero-filled genes hides that. Report the averaged flight effect, by platform, as the headline.
- Caveats: the inter-study correlation uses one `duplicateCorrelation` pass, matching the original script. The variance-partition section of `fit_factorial_model.R` was not rerun (variancePartition not installed). Pathway (GSVA) scores were not recomputed from the rebuilt matrix.

## 10. Follow-ups: no zero-filling at gene level; pathway scores recomputed

`code/factorial_model/followups_2026-10.R` (shared loaders in `osdr_helpers_2026-10.R`). Tables: `followups/`.

**Pooled gene-level fit on genes measured in all 11 studies.** Same pipeline and v2 covariates as §9, but only the 17,129 genes present in every study (17,129 after the CPM filter), so nothing is zero-filled. Inter-study ρ 0.6412 (zero-filled 0.7327).

| contrast | zero-filled FDR < 0.05 (of 25,453) | complete-case FDR < 0.05 | logFC r on common genes | significant in both |
|---|---|---|---|---|
| flight within BRIC | 11,675 | 10,831 | 0.968 | 9,655 |
| flight averaged over hardware | 4,598 | 3,637 | 0.958 | 3,352 |
| flight × CARA | 6,822 | 6,171 | 0.959 | 5,467 |
| flight × VEGGIE | 10,491 | 8,762 | 0.965 | 8,148 |
| Model 2: g_bl | 8,954 | 7,111 | 0.952 | 6,204 |
| Model 2: carbon_12h | 8,449 | 2,099 | 0.931 | 1,904 |

- The averaged flight effect from the complete-case fit agrees with the RNA-seq-only fit at logFC r = 0.893.
- carbon_12h loses most of its genes (8,449 -> 2,099): carbon_12h is 1 % for every BRIC sample and all microarray samples are BRIC, so genes that were zero-filled in the microarray studies looked like a carbon effect. g_bl goes from 8,954 to 7,111.
- Removing the zero-filled genes does not rescue either CFD covariate: they stay confounded with study and hardware by design (§4, §9). Microarray and RNA-seq values are still pooled in one voom model here, which is why §9's platform-separated fit remains the preferred analysis.

**Pathway scores recomputed from the rebuilt matrix.**

- Reproduction: the original recipe (log2(x + 1) of the zero-filled matrix, genes z-scored across all 307 samples) applied to the rebuilt matrix reproduces `pathway_scores_gsva.tsv` (max |difference| 5.3e-03; minimum per-pathway r 1.000000). That recipe log-transforms the microarray values a second time, and its z-scores mix platforms and studies, so pathway scores carry the same zero-filling artefact as the gene matrix.
- An R port of the §5 within-study model matches `cfd_reanalysis.py` on the stored scores (max |Δ estimate| 4.8e-15, max |Δ p| 2.2e-13).
- The gene sets in `hormone_carbon_model.py` repeat 63 gene entries within their own set (in 20 of 28 pathways); the original scores count a repeated gene more than once. The worst case is AT5G52640 in heat_shock: 21 of its 32 entries, so the stored heat_shock score is mostly that one gene. The rescored version counts each gene once (`pathway_genesets.tsv`; repeats in `pathway_geneset_repeats.tsv`).
- Gene-set curation (UniProt-checked; `pathway_geneset_corrections.tsv`, applied to the rescored scores only): AT4G21840 removed (MSRB8, peptide methionine sulfoxide reductase B8; not a heat shock protein); 12 genes added to heat_shock (HSFA1a, HSFA1b, HSFA1d, HSFA1e, HSFA2, HSFA3, HSFA7a, HSFA7b, HSFB1, HSFB2a, HSFB2b, HSP101 (CLPB1)). heat_shock now has 21 genes, of which 17 are measured in every study (`pathway_gene_coverage_by_study.tsv`).
- Rescored per study: RNA-seq as log2-CPM from RSEM expected counts (TMM, `filterByExpr`), microarray on its published log2 scale, each gene z-scored within its study, no zero-filling. 0 of 308 pathway × study cells have fewer than 3 measured genes and are left missing (`pathway_gene_coverage_by_study.tsv`).

Within-study model (`score ~ study + flight × scale`, CR1 by study, BH per term), pathways with FDR < 0.05 out of 28:

| scores | flight | scale | flight × scale |
|---|---|---|---|
| stored (zero-filled, pooled z) — §5 | 0 | 0 | 0 |
| rescored, all 11 studies | 0 | 3 | 0 |
| rescored, RNA-seq studies only (7) | 0 | not estimable | 0 |

- Flight t-statistics from the stored and rescored scores correlate at r = 0.634. Smallest rescored flight p: cytokinin_signaling (p = 0.023, FDR = 0.4).
- Positive control (as in §5: a flight effect of 0.25 / 0.5 × the median score SD planted in every pathway): stored scores 28/28 and 28/28 (reproducing §5), rescored 8/28 and 15/28. The rescored SD is 0.43 × the stored SD because the stored scores' SD includes between-study and between-platform differences; relative to within-study variation, the §5 control planted larger effects than its label suggests.
- A term is tested only if it varies within at least two studies; otherwise, with study fixed effects, it is identified by a single cluster and cluster-robust errors collapse (as for CARA in §7). Among the RNA-seq studies only OSD-678 has both leaf and rosette samples, so scale is not testable there (not estimable); flight × scale still varies within every study and is tested.
- **Reading:** with corrected scores no pathway reaches FDR < 0.05 for flight, as in §5; 3 pathway(s) now show a scale (rosette vs leaf) difference (ethylene_biosynthesis, ja_signaling, gravitropism_signaling). But the corrected positive control shows the within-study test detects a 0.25 SD flight effect in only about a quarter of pathways and a 0.5 SD effect in about half, so the pathway-level null now reads as 'no large, consistent flight effect' rather than 'no effect'; effects below about half a within-study SD cannot be excluded at this sample of 11 studies. §5's statement that the null 'is not a powerless test' should be read with this correction. §11 reruns the §7 in-dish analysis and the §8 meta-regression on these scores.

## 11. §5, §7 and §8 on the rescored pathway scores

The same scripts (`cfd_reanalysis.py`, `cfd_meta_regression.py`) rerun with `PATHWAY_SCORES=followups/pathway_scores_within_study.tsv` and `REANALYSIS_OUT=followups/rescored` (full report: `followups/rescored/REANALYSIS.md`). Comparison by `code/cfd_covariates/rescored_comparison.py`.

| analysis | stored scores | rescored scores |
|---|---|---|
| §5 within-study model, flight (FDR < 0.05 of 28) | 0 (none) | 0 (none) |
| §5 within-study model, scale (FDR < 0.05 of 28) | 0 (none) | 3 (ethylene_biosynthesis, gravitropism_signaling, ja_signaling) |
| §5 within-study model, flight_x_scale (FDR < 0.05 of 28) | 0 (none) | 0 (none) |
| §5 positive control, 0.25 / 0.5 SD | 28/28, 28/28 | 8/28, 15/28 |
| §7 in-dish subset, flight (FDR < 0.05 of 28) | 0 (none) | 1 (auxin_biosynthesis) |
| §8 Knapp–Hartung meta-regression (FDR < 0.05 of 28) | 1 (ga_biosynthesis) | 0 (none) |
| §8 positive control (detected of 28) | 24 | 20 |

- Per-pathway flight t-statistics (§5) correlate at r = 0.63 between stored and rescored scores; meta-regression slopes (§8) at r = 0.33.
- §8's one stored hit, ga_biosynthesis (Knapp–Hartung p = 0.00082), has p = 0.0067 (FDR = 0.094) on the rescored scores. Its stored score repeats 4 gene(s) in its set and includes zero-filled genes, and the VEGGIE-vs-sealed-hardware signal reported in §8 does not survive correct scoring.
- The rescored scores show a scale (rosette vs leaf) difference in 3 pathway(s), consistent with organ composition once study and platform scale differences are removed.
- **Reading:** on correctly scored pathways, neither the all-study model (§5) nor the meta-regression (§8) finds a flight effect or a CFD dose-response at FDR < 0.05; the in-dish subset (§7: 5 studies, t on 4 df) flags 1 pathway(s) for flight (auxin_biosynthesis), which is a single hit among 28 on a small subset and is not replicated in §5. The positive controls show the pathway-level tests can only rule out large effects. The pathway results support 'no large, consistent pathway-level flight response across these 11 studies', not a boundary-layer mechanism; the CARA-specific flight response remains untestable (one study).

## 12. Variance partitioning, redone

`code/factorial_model/variance_partition_2026-10.R`. Tables: `followups/variance_partition_*`.

The stored `tables/variance_partition_summary.tsv` should not be used. `fit_factorial_model.R` entered the categorical factors as fixed effects, which variancePartition does not partition correctly, and the stored summary is internally inconsistent: `cfd_hardware` has a median of 66 % but a mean of ~0 %, which is impossible for non-negative fractions, and there are no residual, g_bl or carbon_12h rows.

Here every categorical factor is a random effect, and hardware and study enter as nested components (each study used one hardware; canonical correlation hardware ~ study = 1.00). g_bl and carbon_12h are left out: they are functions of hardware × flight × scale (canonical correlation carbon_12h ~ study = 1.00).

Percent of each gene's variance, mean / median across genes:

| component | pooled, genes in all studies (17,129 genes, 307 samples) | RNA-seq only (23,531 genes, 213 samples) |
|---|---|---|
| study (within hardware) | 53.3 / 56.1 | 44.8 / 45.9 |
| hardware | 0.2 / 0.0 | 0.2 / 0.0 |
| organ | 1.5 / 0.0 | 1.3 / 0.0 |
| light | 9.9 / 6.8 | 10.2 / 7.4 |
| ecotype | 6.7 / 3.4 | 7.4 / 3.3 |
| flight (common to all studies) | 0.3 / 0.0 | 0.6 / 0.0 |
| flight × hardware | 2.5 / 0.0 | 2.2 / 0.0 |
| flight × study | 8.3 / 4.4 | 11.9 / 6.5 |
| **all flight terms** | 11.1 / 5.7 | 14.7 / 8.5 |
| residual | 17.3 / 14.8 | 21.4 / 17.5 |

- Study and hardware together account for 54 % (pooled) and 45 % (RNA-seq) of the average gene's variance; all flight terms together for 11 % and 15 %.
- Of the flight variance, the part shared by all studies is 0.3 % (pooled) and 0.6 % (RNA-seq); the study-specific part (flight × study) is 8.3 % and 11.9 %. The flight response is mostly study-specific, which is why pooled flight contrasts are sensitive to which studies are included (§9).
- Genes with more than 10 % of their variance in the flight effect common to all studies: 10 of 17,129 (pooled) and 71 of 23,531 (RNA-seq); in flight × study: 4,597 and 8,899.
- Hardware and flight × hardware are estimated from only 3 hardware groups, so their split from the study components is imprecise; read study + hardware, and the three flight terms, as totals.

## 13. RNA-seq vs microarray: platform or study?

`code/factorial_model/platform_concordance_2026-10.R`. Tables: `followups/platform_concordance_*`.

§9 found the RNA-seq and microarray flight effects within BRIC uncorrelated (logFC r = −0.05). Platform is confounded with study and organ: the BRIC RNA-seq studies are whole seedlings, the BRIC microarray studies mostly hypocotyl, root, shoot or cell culture. Here the flight effect is estimated separately in every study × organ unit with at least 2 flight and 2 ground samples (RNA-seq: voom on expected counts; microarray: limma-trend on log2 values), and units from different studies are compared.

| unit | platform | hardware | FLT / GC | genes FDR < 0.05 |
|---|---|---|---|---|
| OSD-321|whole_seedling | RNAseq | BRIC | 22 / 22 | 16,476 |
| OSD-37|whole_seedling | RNAseq | BRIC | 28 / 28 |  7,381 |
| OSD-38|whole_seedling | RNAseq | BRIC | 3 / 9 |  3,722 |
| OSD-147|hypocotyl | microarray | BRIC | 9 / 9 |      8 |
| OSD-17|undifferentiated cell culture | microarray | BRIC | 6 / 6 |     12 |
| OSD-17|whole_seedling | microarray | BRIC | 6 / 6 |      3 |
| OSD-205|hypocotyl | microarray | BRIC | 8 / 8 |      0 |
| OSD-7|hypocotyl | microarray | BRIC | 5 / 5 |      0 |
| OSD-7|root | microarray | BRIC | 5 / 5 |    118 |
| OSD-7|shoot | microarray | BRIC | 5 / 5 |     12 |
| OSD-7|whole_seedling | microarray | BRIC | 3 / 3 |      0 |
| OSD-120|root | RNAseq | CARA | 18 / 18 |    328 |
| OSD-217|leaf | RNAseq | VEGGIE | 3 / 3 | 13,827 |
| OSD-217|root | RNAseq | VEGGIE | 3 / 3 | 15,030 |
| OSD-314|whole_seedling | RNAseq | VEGGIE | 11 / 6 |  6,455 |
| OSD-678|leaf | RNAseq | VEGGIE | 6 / 6 |      7 |
| OSD-678|whole_seedling | RNAseq | VEGGIE | 12 / 12 |  4,633 |

Agreement of flight logFC between units from different studies (Pearson r across shared genes):

| comparison | pairs | median r | range |
|---|---|---|---|
| BRIC, RNA-seq vs RNA-seq | 3 | 0.40 | 0.26 to 0.41 |
| BRIC, microarray vs microarray | 21 | 0.01 | -0.21 to 0.78 |
| BRIC, RNA-seq vs microarray | 24 | -0.03 | -0.24 to 0.30 |
| BRIC whole seedlings, RNA-seq vs RNA-seq | 3 | 0.40 | 0.26 to 0.41 |
| BRIC whole seedlings, RNA-seq vs microarray | 6 | -0.04 | -0.24 to 0.06 |
| all hardware, RNA-seq vs RNA-seq | 34 | 0.09 | -0.38 to 0.57 |
| all hardware, RNA-seq vs microarray | 72 | 0.02 | -0.24 to 0.30 |

For reference, different organs within the same study agree at median r = 0.35. Spearman medians for the first three rows: 0.42, 0.01 and -0.01 (`platform_concordance_summary.tsv`).

- **Finding:** within BRIC, the RNA-seq studies share a reproducible flight response (median r 0.40 between studies), but the microarray studies agree neither with each other (0.01) nor with the RNA-seq studies (-0.03; whole seedlings only: -0.04). The microarray units also detect few flight genes (median 5.5 at FDR < 0.05, against 7,381 in the BRIC RNA-seq units). The platform split in §9 is therefore weak or absent flight signal in the microarray studies, not an opposite response; it can come from the platform's lower sensitivity and smaller groups, or from the earlier flights those studies represent, and this design cannot separate the two. Across hardware, even RNA-seq studies agree poorly (median r 0.09), in line with the study-specific flight response in §12.
- Caveats: one unit per study × organ with small groups (some 3 vs 3 or 3 vs 9), so single-pair correlations are noisy; OSD-7 is APEX-TAGES hardware filed under BRIC (§2); microarray logFCs are compressed relative to RNA-seq, which changes slopes but not correlations.

## Not done here

The gene-level refit is §9 (`code/factorial_model/refit_gene_level_2026-10.R`), its follow-ups §10 (`followups_2026-10.R`), the rescored §5/§7/§8 comparison §11 (`code/cfd_covariates/rescored_comparison.py`), the variance partition §12 (`variance_partition_2026-10.R`) and the platform comparison §13 (`platform_concordance_2026-10.R`); run them in that order after `cfd_meta_regression.py`. No analysis items from this reanalysis remain open. Open for the authors: the headline claim and the manuscript, README and site text that follow from it; the RNA-seq vs microarray result (§13) is to be reported as a finding.
