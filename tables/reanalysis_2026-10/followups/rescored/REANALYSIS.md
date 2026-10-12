# CFD covariate update and pathway re-analysis

307 samples, 11 studies, 28 pathways.

## 1. Mapping check (old LunarLeaf tables -> stored covariates)

LUNARLEAF_OLD_TABLES not set: mapping check skipped.

## 2. Covariates that change with the corrected LunarLeaf tables

| column | n_mismatch |
|---|---|
| earth_equiv_fan_cm_s | 307 |
| dC_CO2_peak | 307 |
| dC_H2O_peak | 239 |
| dC_O2_peak | 239 |
| g_bl_mol_m2_s | 227 |
| sherwood | 227 |
| dC_CO2_mean_bl | 159 |
| o2_excess_ppm_bl | 159 |
| dC_H2O_mean | 159 |
| dC_CO2_mean | 159 |
| dC_O2_mean | 159 |
| delta_mm | 159 |
| carbon_12h_pct_earth | 105 |
| A_end_umol_m2_s | 99 |
| A_start_umol_m2_s | 99 |
| u_max | 97 |
| surf_dC_H2O_hw | 41 |
| surf_dC_CO2_hw | 41 |
| fan_cm_s | 41 |
| enclosure_type | 35 |
| co2_drawdown_mean_ppm | 17 |
| co2_drawdown_peak_ppm | 17 |
| o2_buildup_mean_ppm | 17 |
| Ra_H2O | 17 |
| has_forced_airflow | 17 |
| net_assimilation_pct | 11 |
| dishmean_CO2_hw | 6 |

g_bl per cell, old -> new: earth/leaf 1.000 -> 0.997; earth/rosette 0.540 -> 0.540; ug/leaf 0.494 -> 0.443; ug/rosette 0.291 -> 0.272

Not updated because their LunarLeaf source tables were not regenerated: co2_drawdown_* / o2_buildup_* (T5), net_assimilation_pct (T12).

OSDR hardware -> CFD class used for the v2 covariates (`extract_cfd.py`). Ground-control enclosures as confirmed by R. Barker (2026-10-12): BRIC ground controls use the sealed flight hardware; APEX-TAGES (OSD-7) ground controls were not sealed, so they keep 100 % carbon while OSD-7 flight samples take the BRIC class; CARA and VEGGIE ground controls grew in the ISS Environmental Simulator (ISSES) chamber at KSC (100 %); EMCS (OSD-314) is its own class with no CFD scenario, so its CFD covariates are missing and it drops out of the CFD analyses:

| hardware | cfd_hardware | studies | n |
|---|---|---|---|
| APEX-TAGES | BRIC | OSD-7 | 36 |
| BRIC | BRIC | OSD-147, OSD-17, OSD-205, OSD-321, OSD-37, OSD-38 | 170 |
| CARA | CARA | OSD-120 | 36 |
| EMCS | EMCS | OSD-314 | 17 |
| Veggie | VEGGIE | OSD-217, OSD-678 | 48 |

## 3. Original method (Pearson, BH; n = 307 for the stored covariates, 290 for v2, which EMCS lacks)

Differs from the stored `pathway_cfd_correlations.tsv` (these are not the stored scores): 84/84 pathway x covariate rows compared, max |Δr| = 6.33e-01.

- `g_bl_mol_m2_s`: FDR < 0.05 in 8/28 pathways with old covariates, 5/28 with new; max |r| 0.210 -> 0.224.
- `carbon_12h_pct_earth`: FDR < 0.05 in 12/28 pathways with old covariates, 0/28 with new; max |r| 0.354 -> 0.060.
- `delta_mm`: FDR < 0.05 in 13/28 pathways with old covariates, 11/28 with new; max |r| 0.311 -> 0.339.

`carbon_12h_pct_earth` by OSDR hardware and flight (new covariates): APEX-TAGES FLT 1 / GC 100; BRIC FLT 1 / GC 1; CARA FLT 90 / GC 100; EMCS FLT NaN / GC NaN; Veggie FLT 100 / GC 100. Within BRIC it no longer differs between flight and ground, so its correlation with pathway scores reflects which studies used sealed hardware (a between-study difference), not a flight effect.

## 4. Is g_bl identifiable separately from flight and scale?

| gravity_key | cfd_scale | 0.272 | 0.443 | 0.54 | 0.997 |
|---|---|---|---|---|---|
| earth | leaf | 0 | 0 | 0 | 68 |
| earth | rosette | 0 | 0 | 80 | 0 |
| ug | leaf | 0 | 68 | 0 | 0 |
| ug | rosette | 74 | 0 | 0 | 0 |

Design rank, 290 samples with g_bl (17 without, excluded): study + flight × scale = 13 columns; adding g_bl gives rank 13 of 14 columns. g_bl is a deterministic relabelling of the flight × scale cells, so no model can attribute an effect to g_bl rather than to flight and organ scale in this design. Its correlations with pathways are flight × scale contrasts under another name.

## 5. Within-study model per pathway

`score ~ study + flight × scale` (effect coding), CR1 cluster-robust SEs by study (11 clusters, t on G − 1 df), BH across pathways per term.

- **flight**: FDR < 0.05 in 0/28 pathways.
- **scale**: FDR < 0.05 in 3/28 pathways.
- **flight_x_scale**: FDR < 0.05 in 0/28 pathways.

For comparison, the same model with naive (independent-sample) SEs calls the flight term significant in 8/28 pathways at Bonferroni 0.05: the original significance comes from treating 307 samples as independent when they cluster in 11 studies.
Positive control: a flight effect planted in every pathway is detected (p < 0.05) in 8/28 at 0.25 SD, 15/28 at 0.5 SD, so the test has limited power: flight effects below about half an SD cannot be excluded, and the null means no large, consistent flight effect rather than no effect.

Of the 5 pathways 'significant' for g_bl by Pearson, 0 have any FDR < 0.05 term (flight, scale or interaction) once within-study structure and clustering are respected.

## 6. Gene-level checks on the stored DE tables

- `Flight_FLT_vs_GC` vs `Flight_BRIC_FLT_vs_GC`: 25453 genes compared; logFC identical: **True** (max |Δ| = 0.00e+00). With BRIC as the reference level, the 'Flight' coefficient is the flight effect within BRIC, not a main effect across hardware.
- `CFD_g_bl` (logFC per 1 mol m⁻² s⁻¹; the covariate spans only 0.709, SD 0.253): |logFC| > 1 & FDR < 0.05 → 3633 genes. The same genes per 1 SD of g_bl: |logFC·SD| > 1 → 813; per full observed range: 3010. The cutoff counts depend on the unit chosen, so for a continuous covariate report FDR plus an effect size per SD or per observed range.
## 7. In-dish covariate (LunarLeaf T19)

`g_bl_dish` (leaf conductance to the dish's own air, in the real BRIC 60 mm / CARA 120 mm dish) is defined for 97 leaf-scale samples (leaf, hypocotyl, root and cell culture are mapped to leaf scale) in BRIC or CARA hardware from 5 studies. It is NaN for rosette-scale and VEGGIE samples, so it is kept as a separate column rather than replacing `g_bl`.

| cfd_hardware | gravity_key | n | studies | g_bl_dish | g_bl_open_box | note |
|---|---|---|---|---|---|---|
| BRIC | earth | 23 | 3 | 0.33 | 0.997 | indicative: sealed 1 g not steady, beyond Boussinesq limit |
| BRIC | ug | 38 | 4 | 0.185 | 0.443 | steady |
| CARA | earth | 18 | 1 | 0.263 | 0.997 | steady |
| CARA | ug | 18 | 1 | 0.186 | 0.443 | steady |

Again one value per hardware × gravity cell, so it is identifiable only as hardware × flight. CARA comes from OSD-120 alone, so a CARA-vs-BRIC difference in the flight response is also a difference between studies.

Predicted in-dish penalty (1 g − 0 g of g_bl_dish): BRIC 0.145 (its 1 g value is indicative only), CARA 0.077. In the open box it was 0.554 for every hardware.

Within-study model on this subset (`score ~ study + flight + flight × CARA`, CR1 by study, 5 clusters, t on 4 df): flight effect FDR < 0.05 in 1/28 pathways. CARA-vs-BRIC difference in the flight response: **not testable**: CARA is a single study, and cluster-robust errors for a term that varies within one cluster collapse toward zero (an unguarded run reports it 'significant' in most pathways, which is an artefact).

**Reading:** the dish geometry changes the *size* of the predicted gravity penalty (smaller than in the open box, and different for BRIC and CARA), but it adds no identifiable covariate. In this design it is still a hardware × flight label. Testing it as a dose-response needs study-specific dish scenarios across more studies per hardware.

## 8. Meta-regression: does the predicted boundary-layer penalty explain which studies respond?

13 study x scale units from 10 studies; the predictor takes 5 distinct values (it varies only by hardware x scale). Excluded for lack of a CFD penalty: OSD-314|rosette (EMCS has no CFD scenario).

| unit | hardware | FLT/GC | predicted penalty (|ΔC CO₂| µg − Earth) |
|---|---|---|---|
| OSD-120|leaf | CARA | 18/18 | 0.145 |
| OSD-147|leaf | BRIC | 9/9 | 0.151 |
| OSD-17|leaf | BRIC | 6/6 | 0.151 |
| OSD-17|rosette | BRIC | 6/6 | 0.279 |
| OSD-205|leaf | BRIC | 8/8 | 0.151 |
| OSD-217|leaf | VEGGIE | 6/6 | -0.074 |
| OSD-321|rosette | BRIC | 22/22 | 0.279 |
| OSD-37|rosette | BRIC | 28/28 | 0.279 |
| OSD-38|rosette | BRIC | 3/9 | 0.279 |
| OSD-678|leaf | VEGGIE | 6/6 | -0.074 |
| OSD-678|rosette | VEGGIE | 12/12 | -0.004 |
| OSD-7|leaf | BRIC | 15/15 | 0.151 |
| OSD-7|rosette | BRIC | 3/3 | 0.279 |

- Knapp–Hartung (primary): flight effect scales with the predicted penalty at FDR < 0.05 in **0/28** pathways (smallest p = 0.00476, ga_biosynthesis). The DerSimonian–Laird z-test, which is over-confident with this few units, would have reported 2/28 (smallest p = 0.00042).
- Of the FDR-significant pathways, 0 are also significant with a plain VEGGIE-vs-other indicator, and 0 show the slope within BRIC units alone (where the penalty varies only by organ scale). Leave-one-study-out: 0 stay p < 0.05 whichever study is dropped.
- **Reading:** no pathway's flight effect tracks the predicted penalty after FDR correction, so the meta-regression gives no support for a boundary-layer dose-response at this number of study units.
- Positive control: a planted effect of 1 SD of the unit-level effects per SD of penalty is detected (p < 0.05) in 21/28 pathways.
- Caveats: few units, the penalty varies only by hardware × scale (so it is also confounded with which studies used which hardware), and T10's BRIC/CARA values are the original 3.8 s values, not steady state.

## Not done here

The gene-level refit is §9 (`code/factorial_model/refit_gene_level_2026-10.R`), its follow-ups §10 (`followups_2026-10.R`), the rescored §5/§7/§8 comparison §11 (`code/cfd_covariates/rescored_comparison.py`), the variance partition §12 (`variance_partition_2026-10.R`) and the platform comparison §13 (`platform_concordance_2026-10.R`); run them in that order after `cfd_meta_regression.py`. No analysis items remain open. The headline claim is in `reports/HEADLINE_2026-10.md`; the manuscript still needs the author ORCID and an archival DOI.
