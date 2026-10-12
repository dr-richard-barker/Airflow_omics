# What the microgreen-chamber CFD teaches the LunarLeaf → omics pipeline

*Prepared 2026-10-07 for R. Barker. Compares the `microgreen-chamber-cfd` (GitHub, commit `6ac89c1`) with LunarLeaf-CFD (`b3567f6` plus the changes listed in §5, since merged as PRs #1–#7) and its use in this OSDR meta-analysis package.*

*Follow-up analyses added 2026-10-08 (§3.6): chamber-size sensitivity, surface-distribution metrics, a pathway-level reanalysis of this package's CFD statistics, and the leaf inside the real BRIC and CARA dishes.*

*Updated the same day: the affected LunarLeaf tables (T2, T6, T7, T9, T10, T13) and figures (F3, F4, F6–F9) have been regenerated with both fixes (T14 followed later the same day), and `code/cfd_covariates/extract_cfd.py` now reads the corrected values. §3.5 lists what changed.*

*Updated 2026-10-12 to match `tables/reanalysis_2026-10/REANALYSIS.md` §9–§12: the gene-level limma models were refit from the OSDR processed files, the pathway scores were recomputed without zero-filling or repeated genes, and the variance partition was redone. Where those change a statement below, the later result is given alongside.*

## 1. Bottom line

- **The two projects complement each other more than they overlap.** The OpenFOAM project built a well-engineered 3D model of a whole chamber's airflow. That is the scale LunarLeaf doesn't resolve. LunarLeaf already has the physics the chamber model hasn't run yet: buoyancy, CO₂/O₂/H₂O transport, and 0 g vs 1 g.
- **The chamber model's most transferable asset is its verification habits, not the solver.** Applying three of them to LunarLeaf this afternoon found real problems in the numbers that feed the omics covariates.

  1. **The 30 k-step snapshot was not a steady state for any microgravity case.** The original T13 overstated g_bl by **11 % for the 0 g leaf** (0.494 → 0.443 at steady state) and **7 % for the 0 g rosette** (0.291 → 0.272). The canopy was worse: 0 g g_bl fell 0.109 → **0.058** once run for 900 k steps, and even the 1 g canopy moved 10 %. The Earth leaf and rosette were fine (≤ 0.3 %). CARA is still drifting at 26 s.
  2. **The forced-airflow (fan/VEGGIE) scenes have an outlet boundary-condition defect.**
     - The copy outlet has no pressure reference, so mass accumulates (mean ρ rises 1.00 → 1.045 over 21 s) and the real through-flow decays to 71 % of the set fan speed at the 5.2 s snapshot and to 18 % by 26 s. It never settles.
     - `u_max` hides this, because it reads the pinned inlet column.
     - The age-of-air diagnostic exposed it immediately.
     - With a fixed-pressure outlet the flow is steady by 5 s and passes the outlet-age check to within 3 %.
     - **Consequence:** the Earth-equivalent fan speed (the airflow that restores Earth-level surface gradients in µg) is **2.0 cm/s for a leaf, not 2.8**, and **≈ 6.5 cm/s for a rosette, not 11**. At the VEGGIE setting (8.3 cm/s) the old snapshot understated the leaf's g_bl by 17 % (1.525 vs a steady 1.841).
  3. **The grid is not the limiting error at 0 g.** A three-level resolution ladder (dx 0.288 → 0.144 mm) converges monotonically. The production grid is within 2.0 % (0 g) and 5.5 % (1 g) of the extrapolated g_bl, and the Earth/µg contrast moves from 2.25× to 2.19× (§3.3).
  4. **The fix exposed a solver limit.** With the fan now carrying its full flow, air squeezing between rosette leaves and canopy shoots exceeds the lattice-Boltzmann low-Mach limit (Ma > 0.3) at the faster fan speeds, and three runs blew up. Excluding those runs, **the canopy never reaches its Earth level** within the valid range (≤ 6.6 cm/s). The old "≈ 21 cm/s" canopy figure came from extrapolating past the data. A canopy answer at higher speeds needs a finer lattice.
  5. **Follow-up (§3.6): the 0 g conductance depends on the modelled chamber, and the omics CFD statistics do not survive study structure.** Doubling the chamber lowers 0 g g_bl by 33 % (leaf) and 41 % (rosette) with no plateau, while 1 g moves ≤ 7 %. In this package, g_bl is a relabelling of flight × organ; a within-study model finds no pathway with an FDR-significant flight term (on correctly rescored pathways three show a scale difference; REANALYSIS §10), and a meta-regression finds no dose-response on correctly scored pathways (0/28; REANALYSIS §11).

## 2. How the two models compare

| | OpenFOAM: microgreen-chamber-cfd | Richard: LunarLeaf-CFD |
|---|---|---|
| Question | Does the fan and diffuser ventilate a 2.3 L microgreen chamber evenly, and later, how does that change from 1 g to 0 g? | How does gravity change leaf-surface gas gradients and boundary-layer conductance in spaceflight hardware? |
| Code | OpenFOAM v2606 (ESI): `pimpleFoam`, with `buoyantPimpleFoam` templated | Custom TypeScript D2Q9 lattice Boltzmann, plus D2Q5 advection–diffusion per species |
| Dimension / scale | 3D chamber, 0.26–6 M cells, base cell 3.3 mm | 2D vertical slice, 128 × 96 lattice, dx 0.288 mm |
| Physics actually solved | Isothermal flow, k-ω SST vs laminar, mean age of air | Laminar flow, buoyancy from CO₂/O₂/H₂O concentration, stomatal sources/sinks, CO₂-limited photosynthesis feedback |
| Gravity | Planned Phase 3 (Richardson-number sweep). **Not run yet** | 1, 0.38, 0.17, 0 g. **Run** |
| Plants | Planned porous canopy. **Not implemented** | Solid leaf, rosette and canopy shapes with surface fluxes |
| Time | Transient. Found **no steady state**; averages over a window of τ | Transient. Was reported as a single snapshot at 30 k steps (≈ 5.2 s); now run to steady state with a drift check |
| Verification | Mesh ladder m0/m1/m2; Sandberg outlet-age identity; closed-STL and volume checks; turbulence-model sensitivity | 4 textbook gates (Ghia, Strouhal, de Vahl Davis, erfc). The chamber model has none of these |
| Validation vs measurement | None (no sensor data yet) | Bench chamber gas exchange (T1/F5); no flight-hardware measurements |
| Visual outputs | 6-panel age-of-air maps; around 14 GIF animations with pinned linear and log scales; mesh renders | Static figures F1–F10; no animations until today |
| Reproducibility | Template → `generate_case.sh` → `Allrun` with hard gates; `audit_cases.sh` | esbuild one-liners per script; this package has hard-coded `/mnt` and `/workspace` paths and no environment file |

**What they share:**
- the same core framing: gravity is the variable, and convection versus forced ventilation is the mechanism;
- a validation-first culture;
- the same open gap: no in-flight or in-chamber sensor data to validate against.

## 3. What was borrowed, and what it found

The external repository has no licence, so nothing was copied. The techniques were re-implemented with attribution in each file header.

### 3.1 Check for steady state before trusting a snapshot (`analysis/transient_stats.py`, F11, T15)

Phase 1 testing showed that the "steady" chamber never stops moving. The chamber model therefore runs transiently, throws away a start-up window, and quotes window means with error bars corrected for autocorrelation (N_eff = N·dt / 2T_int).

I ported that machinery: N_eff-corrected SE, a Hann-window PSD, sliding-window RMS, and a two-halves-plus-trend drift test. It has 11 regression tests, including an AR(1) case with known N_eff and a three-grid GCI case with known order. I also added `validation/timeseries.ts`, which records the existing diagnostics every 250 steps for 150 k steps (26 s).

| Scenario | g_bl at 30 k snapshot (current table) | Steady value (13–26 s mean) | Snapshot bias | Verdict |
|---|---|---|---|---|
| leaf-earth | 1.001 | 0.997 | +0.4 % | stationary by about 8 s |
| rosette-earth | 0.541 | 0.540 | +0.2 % | stationary by about 6 s |
| **leaf-ug** | **0.494** | **0.443** | **+11.4 %** | stationary only after about 15 s |
| **rosette-ug** | **0.290** | **0.272** | **+6.6 %** | stationary after about 12 s |
| hw-cara | 0.476 | 0.365 (still falling) | +30 % | **drifting** (−6.5 % across the window) |
| leaf-ug-fan-lo | 1.065 | 0.746 (still falling) | +43 % | **drifting** (−30 %); see 3.2 |

Unlike the microgreen chamber, LunarLeaf's leaf scenes relax monotonically to a true steady state; there is no persistent oscillation. So the error-bar half of that method is mostly moot here. The half that matters is "run long enough, then prove it stopped moving."

The physical reason is simple. In microgravity, transport is diffusive and the diffusion time across the 28 mm chamber, about L²/D, is roughly 15 s, three times the snapshot.

**Fix (applied):** the exports now run 150 k steps (900 k for the 0 g canopy) and record the drift over the last 5 s in every table; see §3.5.

**Effect on the omics covariates:** the 0 g leaf g_bl falls 0.494 → 0.443 and the 0 g rosette 0.291 → 0.272. The Earth/µg contrast widens from 2.0× to 2.25× (leaf) and from 1.86× to 1.99× (rosette).

**A caveat this exposed:** at 0 g the steady halo fills the whole box (see the animation). The µg g_bl is therefore partly set by how far away the ambient-held walls are. The domain-size sensitivity run (T17, §3.6) confirms this.

### 3.2 Age of air exposed a boundary-condition defect (`ageOfAir`, `outletBC`, `validation/age_check.ts`)

The chamber model computes a mean-age-of-air scalar and checks Sandberg's identity: the flux-weighted outlet age must equal the residence time τ = V/Q. A run that fails is either unconverged or not conserving mass. I added the same scalar to LunarLeaf's forced-airflow scenes as an opt-in flag.

Fan U = 0.02 (nominal 3.3 cm/s), 0 g:

| | Original copy outlet | New `outletBC: 'pressure'` |
|---|---|---|
| τ over the run (lattice steps) | 6.6 k → 35 k (inflow decays 5×) | 6,384, constant |
| Mean density in the domain | 1.00 → 1.045, rising | 1.00 |
| Outlet age / τ at 26 s | 0.78, still falling | 0.967 (−3 %; diffusive loss at the inlet) |
| g_bl at 5 s / 26 s | 1.07 / 0.64 | 1.233 / 1.233 |
| ΔC_CO₂ (Earth = 0.161) | −0.150 → −0.251 | −0.131 |

`'pressure'` is now the default. `'copy'` is kept only so the original tables can be reproduced (`LUNARLEAF_LEGACY=1` for the fan sweep).

**Corrected fan sweep (regenerated `T6_forced_airflow.csv`; pressure outlet, 150 k steps, every case steady):**

| Fan (cm/s) | ΔC_H₂O, original T6 (copy, 30 k) | ΔC_H₂O, regenerated T6 |
|---|---|---|
| 0 (0 g, still) | 0.231 | 0.240 |
| 1.7 | 0.153 | 0.135 |
| 3.3 | 0.117 | 0.101 |
| 8.3 (VEGGIE setting) | 0.082 | 0.068 |
| 16.6 | 0.063 | 0.049 |
| Earth 1 g (target) | 0.128 | 0.128 |
| **Earth-equivalent fan speed** | **2.8 cm/s** | **2.0 cm/s** |

A leaf in microgravity needs about 2.0 cm/s of real airflow to recover Earth-level surface gradients, not 2.8. The original sweep needed a higher *nominal* speed because, at the 5.2 s snapshot, its copy outlet was delivering only about 70 % of that speed.

The same defect affects the VEGGIE scenes at 8.3 cm/s (`leaf-ug-fan-mid`):
- **Copy outlet:** g_bl reads 1.525 at the snapshot, then drifts to 1.01 by 26 s.
- **Pressure outlet:** g_bl is steady at 1.841, so the snapshot understates VEGGIE's benefit by 17 %.

**Fan by scale (regenerated T9 / F8).** Every run now records its peak Mach number. Runs above Ma 0.3 are flagged and excluded from the crossing, because the lattice-Boltzmann method is only accurate well below the lattice sound speed:

| Scale | Earth-equivalent fan, old T9 | Regenerated T9 | Fastest valid run |
|---|---|---|---|
| leaf | 2.6 cm/s | **2.0 cm/s** (agrees with T6) | 18.3 cm/s (Ma 0.24) |
| rosette | 11.0 cm/s | **≈ 6.5 cm/s** | 11.6 cm/s (Ma 0.22) |
| canopy | NaN (the README's "≈ 21" was an extrapolation) | **not reached** | 6.6 cm/s (Ma 0.26) |

In the canopy, the surface gap flattens above the Earth level: 0.29 at 6.6 cm/s against Earth's 0.26, and still about 0.28 in the faster, out-of-range runs. Three runs blew up: rosette at 36.5 cm/s and canopy at 26.6 and 36.5 cm/s. The old copy outlet hid this problem, because its flow decayed before it got fast enough to break the solver.

### 3.3 Resolution ladder (`refine`, `validation/resolution_ladder.ts`, `analysis/gci.py`, T16)

The m0/m1/m2 ladder showed the mesh was not the limiting error. I added `refine: r` to the leaf builder. It keeps the physics fixed under diffusive scaling (dx/r, dt/r², source /r, buoyancy /r³, g_bl calibration ×r), and at r = 1 it reproduces T13 exactly. Then I ran r = 1, 1.5 and 2 for 26 s at 1 g and 0 g.

Values are means over 13–26 s. The method is three-grid with non-uniform ratios (Celik et al. 2008, *J. Fluids Eng.* 130:078001). Results are in `results/tables/T16_resolution_gci.csv`.

| Quantity | dx 0.144 / 0.192 / 0.288 mm | Observed order p | Extrapolated | GCI (finest) | Error of production grid (r = 1) |
|---|---|---|---|---|---|
| g_bl, 1 g | 0.980 / 0.986 / 0.997 | 0.75 | 0.955 | 3.2 % | 5.5 % |
| g_bl, 0 g | 0.4394 / 0.4407 / 0.4433 | 1.17 | 0.436 | 0.9 % | 2.0 % |
| ΔC_CO₂, 1 g | −0.164 / −0.163 / −0.161 | 1.00 | −0.167 | 2.3 % | 4.5 % |
| ΔC_CO₂, 0 g | −0.3665 / −0.3654 / −0.3635 | 0.65 | −0.372 | 1.8 % | 2.8 % |

How to read this:

- **Convergence is monotone but low-order**, with p of about 0.7–1.2 against the lattice Boltzmann scheme's formal 2. That is typical when a curved leaf is drawn as a staircase of square cells (bounce-back walls), which is first-order accurate. Because p is low, the extrapolated values are indicative rather than exact.
- **The production grid's discretisation error is 2–5.5 %.** At 0 g it is much smaller than the 11 % snapshot bias, so getting to steady state matters more than refining the grid. At 1 g it is about the same size as the other errors.
- **The Earth/µg g_bl contrast barely moves.** It is 2.25× on the production grid and 2.19× after extrapolation.
- **δ was not analysed separately.** The diagnostics print it to only 2 decimal places, which leaves too few significant figures for a GCI. Since δ = D/g_bl, its error follows from g_bl's.

**What this does not test:** the blade is 2.3 mm thick at every level (about 11× life size). Discretisation error and the thickness inflation are separate questions. At r = 2 a thinner blade becomes resolvable, which is the natural next run.

### 3.4 Animations with pinned colour scales (`analysis/animate_fields.py`)

A key visualization standard is to compute one colour range over all frames and all panels, never auto-rescale, and render a log twin so the near-leaf gradient and the far-field plume are both readable. I applied it to the H₂O excess field. Pillow writes the GIFs without dithering, so no ffmpeg is needed.

- `results/animations/leaf_1g_vs_0g.gif` (+ `_log`): at 1 g a plume rises off the leaf within about 1 s; at 0 g a diffusive halo grows until it fills the chamber.
- `results/animations/rosette_1g_vs_0g.gif` (+ `_log`).
- `results/animations/cara_vs_open_0g.gif` (+ `_log`): CARA taped dish vs ambient walls at 0 g; CARA is still darkening at 26 s.
- `results/animations/fan_copy_vs_pressure.gif` and `veggie_fan_copy_vs_pressure.gif` (+ `_log`): the outlet defect at 3.3 and 8.3 cm/s. With the copy outlet the halo grows back around the leaf; with the pressure outlet the wake stays fixed.
- `results/figures/F11_stationarity.png`: g_bl(t) for six scenarios, with the 30 k snapshot marked.
- `results/figures/F12_outlet_bc.png`: through-flow and g_bl over time, copy outlet vs pressure outlet.

### 3.5 What was regenerated

| Output | Change |
|---|---|
| **T13** (feeds `extract_cfd.py`) | Steady-state g_bl: leaf 1 g 0.997, leaf 0 g **0.443**, rosette 1 g 0.540, rosette 0 g **0.272**, canopy 1 g 0.133, canopy 0 g **0.058**. New column `gbl_drift_pct_last_5s`, which is 0.00 % for all rows. |
| **T2**, F3, F4 | Earth→µg surface gaps now 1.7–2.3×, up from 1.5–1.8×. |
| **T6**, F6 | Earth-equivalent fan speed 2.0 cm/s, down from 2.8. |
| **T9**, F8 | leaf 2.0 / rosette ≈ 6.5 / canopy not reached. New `ma_max` column. |
| **T10**, F9 | VEGGIE rows at steady state (leaf CO₂ gap 0.101 → 0.087). BRIC and CARA unchanged byte-for-byte, because they keep their original 22 k steps (BRIC has no steady state). Canopy VEGGIE is at Ma 0.32 and hatched in the figure. |
| **T14** | Leaf-shape sweep at steady state. µg g_bl fell ~10 % (ellipse 0.494 → 0.443, real outline 0.474 → 0.426), but the shape ratios hardly moved (outline −5.6 % / −3.8 %, curl ×0.814 / ×0.832). The curled blade's µg worst spot is now 2.6× the flat 1 g peak (was 2.4×). |
| **T7**, F7 | VEGGIE time series only; the other cases match the original exactly. F7 had been unbuildable since commit `b3567f6` renamed `CARA tape`; the name is now fixed in `hardware_analysis.py`. |
| `code/cfd_covariates/extract_cfd.py` | `earth_equiv_fan_cm_s` is read from T9 (2.0 / 6.5 / NaN). `fan_cm_s` for VEGGIE is 8.3 cm/s; the old 2.8 was the Earth-equivalent speed mislabelled as the VEGGIE fan. T6 and T9 are read with `comment="#"`. `LUNARLEAF_TABLES` overrides the table path. |

**Not regenerated:**
- **T3 / T5:** need the raw Vernier data (`validation/raw/`), which is not on this machine. The T5 ppm values for µg will rise.
- **T12:** transcribed by hand from console output.
- **The slide deck:** `CFD_Spaceflight_Botany.pptx`/`.pdf` was made by hand and still shows the old numbers. (The paper text — `results/MANUSCRIPT.md`, `results/manuscript/manuscript.tex`, `results/DISCUSSION.md` and the `.docx` — has since been updated, with a TODO on T5; reference 5 was corrected in LunarLeaf PR #7.)
- **The omics outputs downstream of `extract_cfd.py`:** since refit locally from the OSDR processed files (REANALYSIS §9–§12).

### 3.6 Follow-up analyses (2026-10-08)

**Chamber size (LunarLeaf T17, F13).** The 0 g halo reaches the ambient-held walls (§3.1), so the same leaf and rosette were run in chambers 1.5× and 2× the production 36.9 × 27.6 mm, at the same lattice spacing and to steady state.

| g_bl | 1× (37 × 28 mm) | 1.5× | 2× (74 × 55 mm) |
|---|---|---|---|
| leaf · 1 g | 0.997 | 1.015 | 1.037 (+4 %) |
| **leaf · 0 g** | **0.443** | **0.344** | **0.297 (−33 %)** |
| rosette · 1 g | 0.540 | 0.501 | 0.505 (−6 %) |
| **rosette · 0 g** | **0.272** | **0.194** | **0.161 (−41 %)** |

At 1 g, buoyant convection sets the boundary layer locally, so g_bl is a property of the plant. At 0 g it is set by diffusion to the walls and keeps falling, as expected in 2-D, where diffusion to a fixed-value boundary has no size-independent limit. In 3-D it would level off, so the 2-D model overstates this sensitivity.

Consequences:
- The absolute 0 g values, and so the Earth/µg ratio (leaf 2.25× → 3.5×, rosette 2.0× → 3.1×), depend on the chamber modelled.
- Only the ordering (0 g < 1 g, rosette < leaf) is robust.
- The hardware scenes should use each system's real dish dimensions. They were supplied later and are used in T19 (below).

**Surface distribution (LunarLeaf T18).** These are the chamber model's tray-plane CoV and target-band fraction, applied to the leaf surface.

| scenario | mean \|ΔC CO₂\| | CoV | p90/p10 | surface worse than Earth-leaf mean |
|---|---|---|---|---|
| leaf · 1 g / 0 g | 0.161 / 0.364 | 0.23 / 0.13 | 1.8 / 1.4 | 31 % / 100 % |
| rosette · 1 g / 0 g | 0.298 / 0.593 | 0.49 / 0.26 | 3.5 / 2.0 | 83 % / 100 % |
| canopy · 1 g / 0 g | 0.364 / 0.829 | 0.46 / 0.37 | 9.7 / 6.0 | 85 % / 91 % |
| VEGGIE leaf · 0 g | 0.087 | 0.36 | 3.5 | 0 % |
| VEGGIE rosette · 0 g | 0.294 | 0.65 | 5.4 | 71 % |

- **Microgravity degrades the whole surface, not just the worst spots.** The CoV roughly halves, and every part of a leaf or rosette is worse than an Earth leaf's average.
- **VEGGIE restores the rosette's mean (0.294 vs 0.298) but not its distribution.** It leaves the most uneven surface of any case, with p90 0.61 against 0.52 on Earth and a worse crown peak. The manuscript's "brings a rosette to about its Earth level" holds only for the mean.

**Pathway-level reanalysis of this package (`code/cfd_covariates/cfd_reanalysis.py`, `code/cfd_covariates/cfd_meta_regression.py`; full write-up in `tables/reanalysis_2026-10/REANALYSIS.md`).**

- **Covariates.** These were rebuilt with `extract_cfd.py`'s own code. The old LunarLeaf tables reproduce 25/30 stored covariate columns exactly. The other five differ by design: the two fan columns, and carbon_12h, A_start and A_end, because BRIC ground controls now get the sealed-dish value. The new covariates are in `expression_metadata_with_cfd_v2.tsv`.
- **Original method.** Pearson, n = 307 (290 with the v2 covariates, which EMCS lacks), BH. It reproduces the stored correlations to 10⁻¹⁵, and still calls 28/28 pathways significant for g_bl with the new covariates.
- **Identifiability.** g_bl takes exactly one value per gravity × scale cell, so adding it does not raise the design rank (REANALYSIS §4). It cannot be attributed separately from flight and organ.
- **Within-study model.** The model is `score ~ study + flight × scale`, with effect coding, cluster-robust SEs by study (11 studies) and BH. **0/28 pathways** have an FDR-significant flight, scale or interaction term on the stored scores (rescored: 0/28 flight, 3/28 scale; REANALYSIS §10).
  - Positive control: a planted 0.25 SD flight effect is detected in 28/28 on the stored scores. On correctly rescored pathways it is detected in only 8/28 (0.25 SD) and 15/28 (0.5 SD), so the null rules out only large effects (REANALYSIS §10).
  - The same model with naive SEs gives 7/28, so the original significance came from treating 307 clustered samples as independent.
- **Meta-regression.** The question: does the predicted boundary-layer penalty explain which studies respond? It uses 14 study × scale units with Knapp–Hartung inference.
  - **1/28 pathways** (GA biosynthesis, p = 8 × 10⁻⁴).
  - That one signal is equally well explained by a VEGGIE-vs-other indicator, and shows no slope within BRIC alone. It is a hardware-associated difference, not a dose-response.
  - A DerSimonian–Laird z-test would have claimed 14/28; it is over-confident with this few units.
  - Updated 2026-10-12: with EMCS (OSD-314) excluded for lack of a CFD scenario, the stored-score result is 6/28 on 13 units; on correctly rescored pathways it is 0/28, with GA biosynthesis at FDR 0.11 (REANALYSIS §8, §11).
- **Gene-level checks (stored DE tables).**
  - `Flight_FLT_vs_GC` equals the BRIC-only contrast for all 25,453 genes.
  - The g_bl |logFC| > 1 count drops from 3,633 to 813 when expressed per SD.
  - The limma models have since been refit from the 11 OSDR processed files (REANALYSIS §9–§10). Neither CFD covariate is identifiable from study and hardware, and a flight effect common to all studies explains only 0.3 % of gene variance (under 0.01 % in RNA-seq alone); the flight response is mostly study-specific (§12).

**Real dish geometry (LunarLeaf T19, F14).** The chamber-size result meant the hardware scenes needed real dimensions. The dishes were rebuilt as their actual cross-sections, from dimensions supplied by R. Barker:
- BRIC 60 mm round: 53 mm wide, ≈ 7.4 mm of air above ≈ 5.6 mm of agar.
- CARA 120 × 120 mm square: ≈ 9.8 mm of air above ≈ 5.2 mm of agar.
- The leaf sits 1 mm above the agar.

| leaf g_bl | 1 g | 0 g | Earth/0 g |
|---|---|---|---|
| open box (T13) | 0.997 | 0.443 | 2.25× |
| BRIC dish, open to cabin | 0.247 | 0.197 | 1.25× |
| CARA dish, open to cabin | 0.251 | 0.176 | 1.43× |
| CARA dish, taped | 0.247 | 0.170 | 1.45× |
| BRIC dish, sealed (vs dish air) | 0.330* | 0.185 | 1.78×* |

The open and taped rows are conductance to cabin air; the sealed BRIC row is to the dish's own air. The omics covariates and the 2026-10-12 figures use conductance to each dish's own air, which for the taped CARA dish is 0.263 / 0.186 (1.41×).

- **The dish is boundary-layer-limited even on Earth.** With so little headspace, convection barely develops at 1 g. The likely main reason is that the agar under the leaf does not exchange gas (not isolated by a separate run). So the gravity penalty in a dish is 1.25–1.8×, smaller than the open box's 2.25×. That fits the reanalysis finding no robust flight effect tied to g_bl.
- **At 0 g, BRIC and CARA give nearly the same leaf conductance to their own air** (0.185 vs 0.186). The hardware difference is in the dish atmosphere (sealed BRIC keeps depleting) rather than in the boundary layer.
- **Headspace sensitivity** (min–max): −18 % / +5 % for BRIC and ±6 % for CARA.
- **Fed into the omics covariates as separate columns.** `extract_cfd.py` adds `g_bl_dish` for the BRIC/CARA leaf-scale samples (97 since 2026-10-12, when the unsealed OSD-7 ground controls were removed; 5 studies), keeping `g_bl` unchanged. Overwriting `g_bl` would mix dish and open-box geometries and reverse the leaf/rosette order. The predicted in-dish penalty (1 g − 0 g) is 0.145 for BRIC (indicative) and 0.077 for CARA, against 0.554 in the open box. It is still one value per hardware × flight cell, so it adds no identifiable covariate (REANALYSIS.md §7). The flight effect on this subset is FDR-significant in 0/28 pathways on the stored scores (1/28, ethylene_signaling, on rescored scores; §11). A CARA-vs-BRIC difference is not testable, because CARA is a single study and cluster-robust errors collapse for a one-cluster term.
- \*The sealed 1 g case is not steady (3.8 % drift) and its depleted air exceeds the Boussinesq limit.
- **Other caveats:**
  - The 2.3 mm (inflated) blade fills about 30 % of the BRIC headspace, and real seedlings are smaller.
  - The 1 mm wall and base thickness is assumed.
  - Rosettes do not fit under the lid yet.

**Sealed BRIC with photosynthesis feedback (LunarLeaf T20, F15).** The T19 sealed runs used a fixed flux, so the dish's CO₂ fell without bound. T20 makes uptake CO₂-limited (stopping at the compensation point) and measures buoyancy from the dish's own mean, the correct reference in a closed dish.

*Original T20 run, with constant transpiration (kept as `T20_dish_feedback_constant_transpiration.csv`); its 1 g values are corrected below.*

| sealed BRIC + feedback | 1 g | 0 g |
|---|---|---|
| assimilation at 26 s / 100 s / 200 s | 84 / 44 / 10 % | 81 / 43 / 11 % |
| dish reaches compensation | ≈ 204 s | ≈ 219 s |
| g_bl to dish air while assimilating | 0.33 → ≈ 0.40 | 0.187 (constant) |

- **The dish runs out of CO₂ in about 3.5 min at either gravity.** The dish-level mass balance does not depend on gravity. T8's analytic estimate is ≈ 7 min; the model time depends on the photosynthesis calibration.
- **The 0 g conductance is well defined** (constant 0.187, matching T19). The 1 g value rose to ≈ 0.40 while photosynthesis was active, but see the correction below: that plateau needs unbounded humidity.
- **The buoyancy fix alone barely changes the earlier 1 g value** (0.327 vs 0.330 at 26 s).
- **CARA (taped) keeps assimilation at 94–96 %.**
- **Implication for the omics covariates (applied 2026-10-08).** BRIC ground controls use the same sealed hardware (confirmed by R. Barker), so they starve as fast as the flight samples. `extract_cfd.py` now gives BRIC ground controls T11's BRIC value (1 %) instead of 100 %. BRIC's carbon deficit is a hardware effect shared by flight and ground, and cancels in FLT − GC. With the fix, `carbon_12h_pct_earth` is mostly a sealed-vs-other hardware label (it also varies between flight and ground within OSD-7 and OSD-120), so its pathway correlations (REANALYSIS §3) now measure differences between studies, not flight. The within-study results are unchanged (0/28). Confirmed 2026-10-12: CARA and VEGGIE ground controls grew in the KSC ISSES chamber and keep 100 %; APEX-TAGES (OSD-7) ground controls were not sealed and keep 100 % too.

**What drives the 1 g conductance, and a correction to T20 (LunarLeaf T21, T22, F16, F17).** Diagnostic runs separated the two candidate explanations for T20's 1 g rise:
- **The rise at 0.5–1 min is convection onset.** It appears equally with fixed uptake.
- **Transpiration then sustained it.** With CO₂ + O₂ buoyancy only, g_bl falls as photosynthesis fades.

But the model held transpiration constant, so the sealed dish accumulated humidity without bound, and after ≈ 3.7 min the flow ran away toward the Mach limit. **T20's late 1 g collapse was this runaway, not an ill-conditioned ratio, and its ≈ 0.40 plateau depended on unbounded humidity.**

With transpiration driven by the leaf-to-air humidity gap (T22; saturation level estimated, not calibrated):
- The dish saturates in 1–2 min, and transpiration stops.
- Convection dies away. The 1 g conductance decays toward the 0 g value: 0.26 at 60 s, 0.21 at 200 s, 0.19 by 5 min. The 0 g value stays 0.187.
- There is no runaway.
- The result holds across a 4× range of the saturation level.

**In a sealed BRIC dish the 1 g and 0 g boundary layers converge within minutes**, so BRIC ground controls see nearly the same leaf-level gas exchange as flight. Real agar releases water, so a real dish is probably near saturation from the start, which strengthens this.

**Wet agar (LunarLeaf T23), a negative result at 1 g.** Holding the agar surface at saturation leaves the 0 g leaf conductance unchanged (BRIC 0.186–0.189, CARA 0.192, the same as dry agar).

At 1 g the agar-to-air humidity contrast pushes the H₂O buoyancy beyond the model's small-perturbation range:
- Sealed BRIC accelerates from about 87 s and diverges at 129.5 s.
- Taped CARA runs near Mach 0.17, then accelerates.

The option is flagged as unsafe at 1 g. Wet agar at 1 g needs a calibrated humidity buoyancy coefficient (real RH values) before it can be modelled.

## 4. Not yet borrowed (recommended next)

1. **Couple chamber scale to leaf scale.** The tray-plane speed (0.30 m/s at the Ø40 mm operating point, control arm; 45.7 % of tray area within 0.2–0.5 m/s) is exactly the far-field boundary condition LunarLeaf's fan scenes need. It is a measured CFD number that could replace a hand-typed fan speed. A joint project writes itself here.
2. **Template → generator → gated run.** the chamber model's `Allrun` refuses to proceed past a failed `checkMesh` and checks that the solver reached `End`. The LunarLeaf equivalent would be one `make results` that runs the gates, then the drift test, then the exports, and refuses to write tables from an unconverged run.
3. **Continuous covariates.** *Partly done (§3.6):* the distribution metrics are now LunarLeaf readouts (T18). Per-scenario metrics are still constant within each hardware × gravity × scale cell, though, so they cannot be separated from flight × organ at the sample level. What would help is study-specific CFD scenarios (each study's real dish size, plant density and age), so that the predicted penalty varies across studies, plus enough studies per hardware to test a dose-response.
4. **Turbulence-model sensitivity.** This was the chamber model's biggest uncertainty (+76 % tray speed between laminar and k-ω SST). It doesn't apply at LunarLeaf's leaf-scale Reynolds numbers, but it does apply to any VEGGIE-scale claim.

## 5. Cross-model insights and recommendations for the chamber CFD model

- **Run the textbook gates before Phase 2.** De Vahl Davis natural convection (Ra 10⁴, Nu 2.24) should be passed by the chamber model's `buoyantPimpleFoam` set-up before any gravity sweep. Ghia is a 10-minute sanity check of its schemes.
- **Species transport for Phase 4.** LunarLeaf's stomatal source/sink formulation and the calibration of g_bl against measured A (T1/T3/T5) are a ready template for future CO₂ and humidity work.
- **Likely bug.** In `templates/transient/system/functions/transientMonitors`, `traySignal` still uses the old inset-tray bounds `(0.0025 0.0308 0) (0.1175 0.1558 1)`. `trayPlane` moved to the flush-tray bounds on 2026-08-16. This probably explains the two different CoV rows in `round2_results.md`.
- **Housekeeping.**
  - The README predates the 2026-08-16 redesign (fan, ports, tray) and the diffuser rounds.
  - About 1,655 PNG frames and a 23 MB `.pptx` are committed without Git LFS (repo is about 524 MB).
  - Mesh independence has not been redone at the Ø40 mm / 11.8 m³/h operating point.
- **The thermal energy balance is the most important result in that repo.** At the original flow, the chamber can reject only about 1.3 W for a 3 K rise, against a 38.4 W LED placeholder. That is the kind of "check the budget before you mesh" result worth highlighting.

**Changes made in LunarLeaf-CFD on 2026-10-07** (since merged in PRs #1–#7) (the 4 gates and 11 analysis tests still pass):
- `src/scenarios/scenarios.ts`:
  - `leafGeometryRefined`, plus the `refine`, `ageOfAir` and `outletBC` options;
  - **`outletBC` now defaults to `'pressure'`.**
- `validation/jobs.ts`: a parallel job runner with a drift check.
- `export_cfd.ts`, `fan_sweep.ts` and `scales.ts` were rewritten on top of it.
- New scripts: `validation/timeseries.ts`, `resolution_ladder.ts` and `age_check.ts`.
- Small fixes to `validation/analyze_model.py` (skips T3/T5 without raw data) and `hardware_analysis.py` (CARA case names).
- New analysis code: `analysis/transient_stats.py`, `gci.py`, `animate_fields.py`, `plot_stationarity.py`, `plot_outlet_bc.py` and `test_transient_stats.py`.
- Outputs:
  - regenerated T2, T6, T7, T9, T10, T13 and F3, F4, F6–F9;
  - new `results/timeseries/`, `resolution/`, `age/`, `animations/`, T15, T16, F11 and F12.
- `results/README.md` and `README.md` are updated with the new numbers and a revision note.

## 6. Side findings in this omics package (flagged, not changed)

- **g_bl has only 4 distinct values** (leaf/rosette × 1 g/0 g; now 0.997 / 0.443 / 0.540 / 0.272) across 307 samples. Its Pearson correlations (n = 307, `hormone_carbon_model.py`) are really a flight × organ contrast with inflated p-values. Confirmed in §3.6: none survives a within-study model with clustering by study.
- **The README said 16 studies; the data contain 11** (307 samples). Fixed in the 2026-10-12 README.
- **`Flight_FLT_vs_GC` is the flight effect within BRIC**, the reference level, not the flight main effect. That is why it is identical to `Flight_BRIC_FLT_vs_GC` (476/264).
- **|logFC| > 1 for continuous covariates depends on units** (per 1 mol m⁻² s⁻¹ of g_bl).
- **Hardware is confounded with study.** CARA comes from OSD-120 alone. Since 2026-10-12, EMCS (OSD-314) is its own class without CFD covariates, and APEX-TAGES (OSD-7) flight samples share the BRIC scenario while its unsealed ground controls keep 100 % carbon.
- **fig6 B/C show identical bars across hardware**, because g_bl does not depend on hardware in T13.
- **Roots are assigned leaf-scale g_bl.**
- **Reproducibility:** hard-coded `/mnt/shared-workspace` and `/workspace` paths (only the LunarLeaf table path can now be overridden); no `requirements.txt` or `sessionInfo`; LunarLeaf commit not pinned; `CITATION.cff` ORCID is a placeholder. Since 2026-10-12 the README lists the run order and requirements, and the refit logs record `sessionInfo()`.

## 7. Reproduce

```bash
cd LunarLeaf-CFD
for s in export_cfd fan_sweep scales hardware_sim timeseries resolution_ladder age_check gates; do npx esbuild validation/$s.ts --bundle --format=esm --platform=node --loader:.json=json --outfile=validation/$s.mjs; done
node validation/gates.mjs                                   # 4 textbook gates
for s in export_cfd fan_sweep scales; do node validation/$s.mjs list | xargs -P 16 -I{} node validation/$s.mjs job {}; node validation/$s.mjs merge; done   # T2/T13, T6, T9/T10
node validation/hardware_sim.mjs                            # T7
python3 validation/analyze_model.py && python3 validation/scales_analysis.py && python3 validation/hardware_analysis.py   # F3-F9
node validation/timeseries.mjs 150000 250 2500 leaf-ug      # one scenario, 26 s
node validation/age_check.mjs 0.02 26 pressure              # Sandberg check
node validation/resolution_ladder.mjs 0 2 26                # g=0, r=2
python3 -m pytest analysis/ -q
python3 analysis/transient_stats.py && python3 analysis/plot_stationarity.py && python3 analysis/gci.py
python3 analysis/animate_fields.py leaf-earth leaf-ug --labels "Leaf · 1 g" "Leaf · 0 g" --out results/animations/leaf_1g_vs_0g
```
