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

- `carbon_12h_pct_earth`: BRIC ground controls use the same sealed hardware, so they get T11's BRIC value, not 100 % (2026-10-08)
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

## 3. Original method (Pearson, n = 307, BH)

Reproduces the stored `pathway_cfd_correlations.tsv`: 84/84 rows matched, max |Δr| = 7.77e-16.

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

## Not done here

Gene-level limma models (`fit_factorial_model.R`) need `harmonized_expression_matrix_filtered.tsv` from the Biomni workspace. Re-run them with `expression_metadata_with_cfd_v2.tsv`, treating g_bl as non-identifiable from flight × organ, and with effect (sum-to-zero) coding so 'Flight' is a main effect.
