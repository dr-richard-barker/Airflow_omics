# PLAN: NASA OSDR Multiomics Mining for Microgravity Atmospheric Adaptation Model

## Summary

Mine the NASA OSDR Biological Data API for Arabidopsis thaliana multiomics data (RNAseq, microarray, methylation/WGBS, glycomics, proteomics) from Col-0 and WS ecotypes grown on the ISS, build a multi-factorial mixed-effects model (light × organ/tissue × flight × hardware), integrate CFD gas-exchange predictions from LunarLeaf-CFD as continuous covariates, use single-cell RNA-seq atlases to train a variational autoencoder (VAE) for tissue-specific deconvolution of bulk samples, and produce a tissue/organ-specific hormonal and primary carbon allocation adaptation model visualized with ggPlantmap and ggkegg. Deliver as a new standalone npj Microgravity-style LaTeX manuscript with a Zenodo-ready standard package (code + figures + tables + supplementary tables).

## Confirmed Decisions

| Decision point | User choice |
|---|---|
| microRNA gap | Predict from mRNA (psRNATarget/TargetFinder inference) |
| Study scope | Col-0 + WS separately (all studies with either ecotype, ~30+ studies) |
| Deconvolution | scRNA-seq + variational autoencoder (VAE) |
| CFD integration | CFD as covariate (scalar predictors per hardware scenario) |
| Factorial model | Mixed-effects model (fixed: flight/hardware/light/organ/ecotype; random: study/batch) |
| Manuscript | New standalone manuscript |
| Model scope | Full hormone panel (auxin, cytokinin, ethylene, ABA, GA, JA, SA) + primary carbon |
| Deliverables | Standard package (Zenodo-ready code + LaTeX PDF + SVG figures + main/supplementary tables) |

## OSDR Data Landscape (verified via API probing)

**API base**: `https://visualization.osdr.nasa.gov/biodata/api/v2/`
- REST: `/v2/dataset/{OSD-ID}/`, `/v2/dataset/{OSD-ID}/assays/`, `/v2/dataset/{OSD-ID}/files/`
- Query: `/v2/query/metadata/?{field}={value}&{field}={value}` (CSV/TSV/JSON output)
- File download: `/v2/dataset/{OSD-ID}/files/` returns file URLs

**62 Arabidopsis studies, 2,204 samples**. Omics distribution:
- RNAseq: ~30 studies | microarray: ~20 | methylation (WGBS): 4 (OSD-217, OSD-416, OSD-520, OSD-625) | glycomics: 2 (OSD-121, OSD-615) | proteomics: 3 (OSD-16, OSD-38, OSD-522) | **microRNA: 0** (will infer from mRNA)

**Hardware factor** (from `experiment platform` metadata):
- BRIC/BRIC-PDFU: OSD-17, OSD-37, OSD-38, OSD-44, OSD-121, OSD-147, OSD-205, OSD-321
- BRIC-LED: OSD-522
- Veggie: OSD-193, OSD-217, OSD-218, OSD-281, OSD-416, OSD-427, OSD-615, OSD-625
- APEX/TAGES: OSD-7, OSD-16
- EMCS: OSD-223, OSD-251, OSD-314, OSD-346, OSD-437, OSD-480
- CARA (Petri dishes): OSD-120

**Factorial structure** is encoded in sample names (e.g. `Atha_Col-0-PhyD_root_FLT_Alight_Rep1`), not structured metadata. A regex-based parser will extract: ecotype, organ/tissue, flight condition (FLT/GC), light condition, replicate, age/timepoint.

**Organ coverage** (from sample names + metadata): root, shoot/leaf, hypocotyl, whole seedling, cotyledon — highly unbalanced across hardware types.

**Light factor** explicitly present in: OSD-120 (light/dark), OSD-314 (dark), OSD-678 (light/dark). Light is implicit (all studies grown in light unless dark-specified).

## Single-Cell Reference Datasets

| Dataset | Accession | Tissue | Cells | Use |
|---|---|---|---|---|
| Shahan et al. 2020 root atlas | GSE152766 | Root | 110k+ | Primary root VAE training |
| Denyer et al. 2019 | GSE122687 | Root | 12k | Root augmentation |
| Plant Physiology 2023 root+shoot | (published) | Root + Shoot | 70k | Shoot/leaf VAE training |
| UWA ecotype atlas | (UWA repo) | Root tip | 275k (Col-0, C24, Cvi-0, Ler, Ws-2) | Ws-2 ecotype matching |
| Cell cycle reference | GSE262840 | Root | — | Cell-cycle state annotation |

## Implementation Plan (8 subsystems)

### Subsystem 1: OSDR Data Mining & Harmonization (`osdr_miner/`)
**Goal**: Programmatically retrieve all Col-0/WS Arabidopsis OSDR multiomics data and metadata, parse into a harmonized sample × factor metadata table.

**Steps**:
1. Query `/v2/query/metadata/?study.characteristics.organism=arabidopsis thaliana` for all samples across 62 studies (CSV format).
2. For each study, fetch `/v2/dataset/{OSD}/` for study-level metadata (experiment platform, study title, factor names, project type).
3. Filter to spaceflight studies (project type = "Spaceflight Study") with Col-0 or WS ecotype (parse from sample names + ecotype metadata field).
4. Build a regex-based sample-name parser to extract: ecotype (Col-0/WS/other), organ (root/shoot/leaf/hypocotyl/whole_seedling), flight (FLT/GC), light (light/dark/unspecified), replicate, hardware (from study-level experiment platform field).
5. For each assay in each included study, fetch `/v2/dataset/{OSD}/files/` to identify and download processed data files (normalized counts, differential expression tables, raw counts where available).
6. Harmonize gene identifiers across platforms: map Affymetrix probe IDs → AGI locus (ATH1 array via TAIR/Arabidopsis thaliana annotation), RNAseq gene IDs → AGI, WGBS → AGI + context (CG/CHG/CHH), glycomics → glycan epitope profiles, proteomics → AGI.
7. Output: `sample_metadata.tsv` (sample × factor table), per-omics normalized matrices, and a `study_metadata.tsv` (study × hardware/factor table).

**Edge cases**: OSD-219 has empty metadata (fetch via REST endpoint directly); some studies use GSM/SAMN accessions as sample names (parse factors from ISA-Tab metadata instead); microarray and RNAseq need separate normalization before cross-platform harmonization (ComBat or limma::removeBatchEffect for cross-platform, with platform as batch).

### Subsystem 2: Multi-Factorial Mixed-Effects Model (`factorial_model/`)
**Goal**: Fit a mixed-effects model for each gene (and each omics layer) testing the effects of flight, hardware, light, organ, and ecotype on expression/methylation/abundance.

**Design**:
- **Fixed effects**: Flight (FLT vs GC) × Hardware (BRIC/CARA/Veggie/BRIC-LED/APEX/EMCS) × Light (light/dark) × Organ (root/shoot/leaf/hypocotyl/seedling) × Ecotype (Col-0/WS)
- **Random effects**: Study (intercept) + Batch/sequencing run (intercept where available)
- **Framework**: 
  - RNAseq: dream (variancePartition) or limma-voom with duplicateCorrelation for random effects — handles the count→voom→mixed model pipeline
  - microarray: limma with mixed-effects design
  - methylation: DSS or beta-regression mixed model per DMR
  - proteomics/glycomics: limma on log-intensity
- **Contrasts**: Primary contrasts are (1) FLT vs GC within each hardware, (2) hardware differences within FLT, (3) organ-specific flight responses, (4) ecotype × flight interaction.
- **Multiple testing**: BH-FDR < 0.05, with effect-size reporting (log2FC for expression, delta-beta for methylation).
- **Output**: `de_results_{omics}.tsv` per omics layer, `factorial_contrasts.tsv`, summary of significant genes per factor/contrast.

**Assumption**: The unbalanced design means not all factor combinations exist; the mixed model will estimate main effects and available interactions, and we will explicitly report which interaction terms are estimable vs confounded (e.g., light × hardware is only estimable for OSD-120/OSD-678).

### Subsystem 3: CFD Covariate Extraction (`cfd_covariates/`)
**Goal**: Run LunarLeaf-CFD scenarios for each hardware type and extract scalar gas-exchange predictors to use as continuous covariates in the factorial model.

**Steps**:
1. Clone LunarLeaf-CFD repo. The repo already has validated results for BRIC (sealed), CARA (micropore tape), Veggie (vented), and BRIC-LED scenarios across gravity levels.
2. Extract from existing `results/` tables (or re-run headless sweep if needed) the following per-hardware, per-gravity (1g ground control vs µg flight) scalars:
   - Boundary-layer thickness δ (mm)
   - Boundary-layer conductance g_bl (mol m⁻² s⁻¹)
   - Surface-to-bulk ΔC for CO₂, O₂, H₂O
   - Dimensionless numbers: Gr, Ra, Re, Sh, Nu, Pe
   - Predicted carbon fixation rate (% of Earth equivalent)
3. Build a `cfd_covariates.tsv` table: hardware × gravity → scalar predictors.
4. Merge CFD covariates into the sample metadata table so each sample gets the CFD scalars for its hardware × flight condition.
5. In the factorial model (Subsystem 2), add a secondary model variant where g_bl and ΔC_CO₂ replace the categorical hardware factor, testing whether CFD-predicted gas exchange explains the hardware effect.

**Output**: `cfd_covariates.tsv`, `cfd_model_comparison.tsv` (categorical hardware vs CFD covariate model fit comparison).

### Subsystem 4: Single-Cell VAE Deconvolution (`sc_vae_deconv/`)
**Goal**: Train a variational autoencoder on Arabidopsis scRNA-seq atlases to learn latent tissue/cell-type signatures, then project bulk OSDR samples into the latent space to estimate cell-type proportions and assign tissue-specific functions.

**Steps**:
1. Download and preprocess scRNA-seq reference datasets:
   - Root: GSE152766 (Shahan atlas, 110k cells, Col-0) — primary root reference
   - Shoot: Plant Physiology 2023 root+shoot atlas (70k cells) — shoot/leaf reference
   - Ws-2 augmentation: UWA ecotype atlas (Ws-2 root tip cells) for ecotype-specific validation
2. Integrate root + shoot atlases (Seurat/scanpy harmony integration), annotate cell types using published markers (stele, cortex, endodermis, epidermis, columella, QC, xylem, phloem, leaf mesophyll, vasculature, guard cells, etc.).
3. Train a VAE (PyTorch, scVI-style architecture) on the integrated atlas:
   - Encoder: gene expression → latent space (z, 20-30 dims)
   - Decoder: latent space → gene expression
   - Cell-type labels used for semi-supervised training (scVI | scANVI)
   - Latent space captures tissue/cell-type identity
4. Project bulk OSDR RNAseq/microarray samples into the VAE latent space:
   - Pseudobulk the scRNA-seq by cell type → build a cell-type × gene reference matrix
   - Use the VAE decoder + a deconvolution head (BayesPrism-style or direct latent projection) to estimate cell-type proportions per bulk sample
   - Validate deconvolution on pseudobulk samples with known proportions
5. For each bulk sample, output: cell-type proportion estimates, dominant tissue assignment, and tissue-specific expression scores.
6. Cross-reference deconvolved cell-type proportions with the factorial model: test whether flight/hardware effects on specific genes are mediated by shifts in cell-type composition vs within-cell-type expression changes.

**Output**: `cell_type_proportions.tsv`, `deconvolution_validation.tsv`, `tissue_assignments.tsv`.

**Compute**: VAE training on ~180k cells × ~25k genes — moderate GPU/CPU, ~8-16 GB RAM, ~30-60 min on a right-sized machine. Deconvolution of ~500 bulk samples is lightweight.

### Subsystem 5: miRNA Target Inference (`mirna_inference/`)
**Goal**: Since no Arabidopsis microRNA sequencing data exists in OSDR, infer miRNA-target regulatory interactions from the mRNA-seq data.

**Steps**:
1. Download the Arabidopsis miRNA-target prediction database (psRNATarget or sRNAtoolbox) with known Ath miRNAs (miRBase) and their predicted targets.
2. For each OSDR RNAseq sample, compute anti-correlation between miRNA proxy signatures (using target gene expression as inverse proxy) and pathway-level miRNA activity scores.
3. Use the Gene Ontology / KEGG miRNA-target network to identify miRNAs whose predicted targets are differentially expressed in the flight/hardware contrasts.
4. Output a `mirna_activity_scores.tsv` and `mirna_target_network.tsv` to integrate into the systems biology model.

**Output**: `mirna_activity_scores.tsv`, `mirna_target_network.tsv`.

### Subsystem 6: Hormone & Carbon Allocation Systems Biology Model (`systems_biology/`)
**Goal**: Build a tissue/organ-specific hormonal and primary carbon allocation model that integrates the factorial DE results, CFD gas-exchange predictions, and single-cell deconvolution.

**Steps**:
1. **Hormone pathway scoring**: For each bulk sample (deconvolved into cell-type proportions), compute pathway activity scores for the full hormone panel:
   - Auxin: PIN/AUX1/IAA/SAUR/GH3, DR5-responsive genes
   - Cytokinin: ARR-A/ARR-B, CRE1/AHK
   - Ethylene: EIN3/EIL1, ERFs, ACS/ACO
   - ABA: ABI/ABF/SnRK2/RD29A
   - GA: GA20ox/GA3ox/RGA/SLR
   - JA: JAZ/MYC2/LOX/AOS
   - SA: NPR1/PR1/PR5
   - Use GSVA or AUCell for pathway scoring per sample
2. **Carbon allocation**: Score primary carbon metabolism pathways (starch/sucrose, glycolysis, TCA, Calvin cycle, photorespiration) using KEGG ath pathways (ath00100, ath00620, ath00010, ath00020, ath00710, ath00630).
3. **Tissue-specific mapping**: Map hormone/carbon scores onto cell types using the scRNA-seq atlas cell-type × pathway score matrix, then onto organ anatomy using ggPlantmap:
   - Root: `ggPm.At.roottip.longitudinal` (meristem, QC, elongation, maturation zones)
   - Root cross-section: `ggPm.At.roottip.crosssection` (stele, cortex, endodermis, epidermis)
   - Shoot/leaf: `ggPm.At.leaf.topview`, `ggPm.At.leaf.crosssection`
   - Seedling: `ggPm.At.seedling.saltdrought`
   - Rosette: `ggPm.At.3weekrosette.topview`
4. **virtual-root integration**: Feed the auxin pathway scores and CFD-predicted gas exchange into the virtual-root auxin transport model to simulate how altered auxin distribution (from spaceflight transcriptomics) affects root tip patterning under each hardware scenario.
5. **KEGG pathway visualization**: Use ggkegg to overlay differentially expressed genes onto Arabidopsis KEGG pathway maps (carbon metabolism, hormone signal transduction, photosynthesis) with flight/hardware-specific fold changes.
6. **CFD → metabolism linkage**: Test whether CFD-predicted ΔC_CO₂ / g_bl correlates with carbon metabolism pathway scores and hormone pathway scores across hardware types.

**Output**: `hormone_pathway_scores.tsv`, `carbon_pathway_scores.tsv`, `tissue_pathway_maps/` (ggPlantmap + ggkegg figures), `virtual_root_predictions/`.

### Subsystem 7: Visualization & Figures (`figures/`)
**Goal**: Produce all manuscript figures as SVG (per user preference).

**Figures** (planned):
1. **Fig 1**: Study design overview — factorial design matrix (studies × factors), hardware schematic, sample counts per cell
2. **Fig 2**: CFD gas-exchange predictions per hardware × gravity — δ, g_bl, ΔC CO₂/O₂/H₂O (from LunarLeaf-CFD results)
3. **Fig 3**: Multi-factorial model results — volcano plots per contrast, heatmap of significant genes × factors
4. **Fig 4**: Single-cell VAE deconvolution — UMAP of atlas, cell-type proportion shifts (flight vs GC) per hardware, tissue assignment validation
5. **Fig 5**: Hormone pathway activity — ggPlantmap heatmaps (root tip longitudinal + leaf) showing hormone pathway scores per organ × flight × hardware
6. **Fig 6**: Carbon allocation — ggkegg pathway maps with DE overlay, starch/sucrose/TCA pathway activity per hardware
7. **Fig 7**: CFD-transcriptome integration — correlation of g_bl/ΔC with carbon/hormone pathway scores, virtual-root auxin distribution predictions per hardware
8. **Fig 8**: Multi-omics integration — methylation × expression × proteomics concordance for key hormone/carbon genes
9. **Supplementary figures**: per-study QC, deconvolution validation, full hormone panel per cell type, glycomics cell-wall profiles

**Libraries**: Python (seaborn + matplotlib) for data plots; R (ggplot2 + ggPlantmap + ggkegg + ComplexHeatmap) for pathway/anatomy maps. All saved as SVG with `svg.fonttype = 'none'` (Python) or `svglite` (R).

### Subsystem 8: Manuscript & Zenodo Package (`manuscript/`, repo root)
**Goal**: Produce a new standalone npj Microgravity-style LaTeX manuscript and a Zenodo-ready software package.

**Manuscript structure** (npj Microgravity style):
1. Title: "Atmospheric boundary-layer physics drives tissue-specific hormonal and carbon allocation responses in Arabidopsis grown in microgravity"
2. Abstract
3. Introduction — spaceflight gas-exchange problem, boundary-layer mechanism, multi-study OSDR meta-analysis rationale
4. Results
   - CFD gas-exchange predictions across hardware × gravity
   - Multi-factorial transcriptomic response (flight × hardware × organ × light × ecotype)
   - Single-cell deconvolution reveals tissue-specific flight responses
   - Hormone pathway remodeling is hardware- and tissue-dependent
   - Carbon allocation shifts correlate with CFD-predicted gas exchange
   - Multi-omics concordance (methylation, proteomics, glycomics)
   - virtual-root auxin model predicts root tip phenotypes per hardware
5. Discussion — atmospheric design implications for spaceflight growth hardware
6. Methods — OSDR data retrieval, factorial mixed model, CFD covariates, scRNA-seq VAE deconvolution, miRNA inference, pathway scoring, systems biology model
7. Data & code availability — Zenodo DOI + GitHub

**LaTeX**: npj Microgravity LaTeX template (Nature-style). Compiled to PDF.

**Zenodo package structure**:
```
microgravity_atmospheric_adaptation/
├── README.md
├── LICENSE (MIT for code, CC-BY-4.0 for manuscript/figures)
├── .zenodo.json
├── CITATION.cff
├── osdr_miner/           # OSDR API mining + harmonization scripts
├── factorial_model/      # Mixed-effects model scripts
├── cfd_covariates/       # CFD covariate extraction
├── sc_vae_deconv/        # Single-cell VAE deconvolution
├── mirna_inference/      # miRNA target inference
├── systems_biology/      # Hormone/carbon allocation model
├── figures/              # All figures (SVG)
├── tables/               # Main tables (CSV/TSSV)
├── supplementary/        # Supplementary tables and figures
├── manuscript/           # LaTeX source + compiled PDF
│   ├── manuscript.tex
│   ├── manuscript.pdf
│   ├── references.bib
│   └── npj_microgravity.cls
└── data/                 # Processed/harmonized data (not raw OSDR data)
    ├── sample_metadata.tsv
    ├── study_metadata.tsv
    ├── cfd_covariates.tsv
    ├── de_results_*.tsv
    ├── cell_type_proportions.tsv
    ├── hormone_pathway_scores.tsv
    └── carbon_pathway_scores.tsv
```

**Tables** (in `tables/`):
- Table 1: OSDR studies included (OSD-ID, title, hardware, ecotype, organ, assay, n samples, factors)
- Table 2: Factorial model summary (significant genes per factor/contrast)
- Table 3: CFD covariates per hardware × gravity
- Table 4: Cell-type proportion shifts (flight vs GC per hardware)
- Table 5: Top hormone/carbon pathway genes per contrast

**Supplementary tables** (in `supplementary/`):
- S1: Full sample metadata
- S2: Full DE results per omics layer
- S3: miRNA activity scores
- S4: Cell-type proportions per sample
- S5: Pathway scores per sample
- S6: CFD full dimensionless number table

## Compute & Resource Estimate

| Subsystem | Workload | Est. RAM | Est. runtime | Execution target |
|---|---|---|---|---|
| 1. OSDR mining | API queries + file downloads | 2-4 GB | 30-60 min (network-bound) | worker-0 (foreground) |
| 2. Factorial model | Mixed-effects on ~500 samples × 25k genes × 5 omics | 8-16 GB | 1-2 h | Right-sized machine (16 GB) |
| 3. CFD covariates | Extract from existing LunarLeaf-CFD results | 1 GB | 15 min | worker-0 |
| 4. scRNA-seq VAE | Train VAE on ~180k cells × 25k genes | 16-32 GB | 1-3 h (GPU beneficial) | Right-sized machine (32 GB) |
| 5. miRNA inference | Target prediction + correlation | 4 GB | 30 min | worker-0 |
| 6. Systems biology | Pathway scoring + ggPlantmap/ggkegg + virtual-root | 8 GB | 1-2 h | Right-sized machine (16 GB) |
| 7. Figures | Plotting | 4-8 GB | 30-60 min | worker-0 |
| 8. Manuscript | LaTeX compilation | 1 GB | 10 min | worker-0 |

**Total estimated runtime**: ~6-10 h of compute, chunked across multiple foreground/background sessions.

**Machine plan**: 
- worker-0 for lightweight tasks (OSDR mining, CFD extraction, miRNA, figures, LaTeX)
- Create a 32 GB machine for the VAE training and factorial model fitting
- All intermediate checkpoints saved to `/mnt/shared-workspace/`
- All final deliverables saved to `/mnt/results/`

## Key Assumptions

1. OSDR processed data files (normalized counts, DE tables) are available for download via the files endpoint — if only raw FASTQ is available, we will use GeneLab RCP-processed versions or process a subset.
2. The LunarLeaf-CFD existing `results/` tables contain the per-hardware gas-exchange scalars we need; if not, the headless sweep can be re-run.
3. The Arabidopsis scRNA-seq atlases (GSE152766, Plant Phys 2023) are publicly downloadable from GEO/SRA.
4. Cross-platform harmonization (microarray ↔ RNAseq) via ComBat/limma is acceptable for the meta-analytic mixed model, with platform included as a random effect.
5. The unbalanced factorial design limits estimable interaction terms; we will report which terms are confounded and focus interpretation on main effects + estimable interactions.
6. miRNA inference from mRNA is a proxy and will be clearly labeled as predictive, not measured.
7. The virtual-root model's auxin transport parameters can be modulated by transcriptomics-derived hormone pathway scores as a qualitative-to-quantitative bridge.

## Testing & Acceptance Criteria

1. **OSDR mining**: All included studies have complete sample metadata with parsed factors; spot-check 5 studies against OSDR web interface.
2. **Factorial model**: Model converges for >95% of genes; residual diagnostics pass for top genes; known spaceflight response genes (e.g., auxin-related, cell-wall remodeling) recovered in FLT vs GC contrast.
3. **CFD covariates**: CFD scalars match published LunarLeaf-CFD results tables; covariate model explains a significant portion of hardware effect.
4. **VAE deconvolution**: Cell-type proportions sum to ~1 per sample; pseudobulk validation R² > 0.7; known tissue markers enriched in assigned cell types.
5. **Systems biology**: Hormone pathway scores discriminate flight vs GC; CFD g_bl correlates with carbon pathway scores (Spearman |ρ| > 0.4 across hardware).
6. **Manuscript**: LaTeX compiles to PDF without errors; all figures referenced; all tables present.
7. **Zenodo package**: `.zenodo.json` valid; CITATION.cff valid; README includes run instructions; all scripts execute end-to-end on a fresh clone.
