# Headline claim (2026-10)

*Proposed for the manuscript, README and site rewrite. Every number below is traced to `tables/reanalysis_2026-10/REANALYSIS.md` and its tables. It replaces the CFD-correlation headline of the August draft.*

## The claim

> **Spaceflight changes Arabidopsis gene expression substantially, but differently in each experiment. Across 11 OSDR studies, a flight response shared by all studies explains under 1 % of gene variance. Within BRIC hardware, only the RNA-seq studies share a reproducible flight signature. Neither CFD-predicted boundary-layer conductance nor carbon supply explains which studies respond, and that is consistent with the CFD itself: inside a real flight dish, gravity changes leaf gas exchange far less than open-air models suggest.**

## Supporting findings

| # | Finding | Number | Source |
|---|---|---|---|
| 1 | The flight response is study-specific | A flight effect common to all studies explains 0.3 % of gene variance on average (0.6 % in RNA-seq alone). Study-specific flight effects explain 8.3 % (11.9 %), and study explains 53 % (45 %). Only 10 of 17,129 genes (71 of 23,531 in RNA-seq) have more than 10 % of their variance in the common flight effect. | §12; `followups/variance_partition_summary_v2.tsv` |
| 2 | The pooled "flight main effect" was the BRIC effect | The stored `Flight_FLT_vs_GC` is the flight effect within BRIC: 11,675 genes at FDR < 0.05. Averaged over BRIC, CARA and VEGGIE it is 4,598 genes (285 with \|logFC\| > 1), and the two agree only weakly (r = 0.35). | §6, §9; `gene_level/de_summary_v2.tsv`, `gene_level/flight_effect_concordance.tsv` |
| 3 | Only the BRIC RNA-seq studies agree with each other | Between BRIC studies, flight effects correlate at median r = 0.40 for RNA-seq vs RNA-seq, 0.01 for microarray vs microarray, and −0.03 for RNA-seq vs microarray (−0.04 for whole seedlings). The microarray studies detect a median of 5.5 flight genes at FDR < 0.05, against 7,381 in the BRIC RNA-seq studies. Across hardware, RNA-seq studies agree at r = 0.09. | §13; `followups/platform_concordance_summary.tsv`, `followups/platform_concordance_units.tsv` |
| 4 | The CFD covariates are not identifiable | g_bl takes one value per gravity × organ-scale cell (VIF 13.6). Study membership explains 99.9 % of carbon_12h's variance. Their gene counts are flight × scale and BRIC-vs-other contrasts. | §4, §9 |
| 5 | No pathway-level CFD dose-response | Within-study model on correctly scored pathways: 0/28 for flight. Meta-regression on the predicted boundary-layer penalty: 0/28. The positive control detects a planted 0.5 SD flight effect in 15/28 pathways, so only large effects are ruled out. | §5, §10, §11; `followups/pathway_within_study_model_rescored.tsv`, `followups/rescored/pathway_flight_meta_regression.tsv` |
| 6 | The CFD explains why | In the real BRIC and CARA dishes, the leaf's gravity penalty is 1.25–1.8×, not the open box's 2.25×. At 0 g both dishes give nearly the same leaf conductance (0.185 vs 0.186). A sealed BRIC dish runs out of CO₂ in about 3.5 min at either gravity, so BRIC ground controls are carbon-starved too. In a sealed dish the 1 g and 0 g boundary layers converge within minutes. | LunarLeaf-CFD T19, T20, T22; §7 |

## Claims this replaces

| August draft | Status |
|---|---|
| "Atmospheric microenvironment, not microgravity alone, drives plant spaceflight responses" | Not supported: no CFD covariate is identifiable (4), and no dose-response is found (5) |
| g_bl correlates with photosynthesis, auxin, ABA and JA (r up to −0.467, "all FDR < 10⁻⁹") | Pearson on 307 clustered samples. 0/28 survive the within-study model, and the correlations' sign runs against the stated mechanism |
| "Hardware dominates over flight" / 7.1 % of variance to hardware | Hardware is 0.2 %, nested in study. The 7.1 % was a misread of 7.05 × 10⁻¹⁵ |
| "Flight main effect": 476 up / 264 down | This is the within-BRIC effect (2) |
| BRIC limits carbon fixation to ~1 % "in microgravity" | The sealed dish limits carbon at 1 g and 0 g alike; it is a hardware effect, not a flight effect (6) |

## What the claim does not say

- **That hardware has no effect.** Hardware is nested in study, and CARA is a single study (OSD-120), so hardware effects cannot be separated from study effects with these data.
- **That flight has no effect.** Flight changes thousands of genes within studies; it is the shared, reproducible part that is small.
- **That the microarray studies show no biology.** Their weak signal may reflect the platform (sensitivity, smaller groups) or the earlier flights they represent. This design cannot tell which.
- **Open caveats:**
  - OSD-7 (APEX-TAGES) and OSD-314 (EMCS) are folded into the BRIC and VEGGIE classes (§2, TODO).
  - The CFD humidity saturation level is an estimate.

## Abstract-ready paragraph

Spaceflight transcriptome studies of *Arabidopsis* differ in hardware, tissue, platform and flight, and their results rarely agree. We re-analysed 307 samples from 11 NASA OSDR studies with study-aware models and coupled them to lattice-Boltzmann CFD of the BRIC, CARA and VEGGIE growth hardware.
- A flight response shared by all studies explained under 1 % of gene variance, against 8–12 % for study-specific flight effects.
- Within BRIC, flight signatures agreed between the RNA-seq studies (median r = 0.40) but not with the earlier microarray studies, which detected almost no flight response.
- CFD-predicted boundary-layer conductance and carbon supply were confounded with study and did not explain which studies responded.
- Consistent with this, the CFD shows that inside the flight dishes gravity changes leaf gas exchange only 1.25–1.8-fold, and that sealed BRIC dishes become carbon-limited on the ground as well as in orbit.

Spaceflight responses are real but largely experiment-specific, which argues for replicated flights in common hardware rather than pooled meta-analysis.
