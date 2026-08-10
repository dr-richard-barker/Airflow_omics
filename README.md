# Microgravity Atmospheric Adaptation Model

Multi-omics meta-analysis of NASA OSDR Arabidopsis spaceflight data integrated with CFD gas-exchange modeling, single-cell VAE deconvolution, and systems biology models of hormone and carbon allocation.

## Overview

This repository mines the NASA Open Science Data Repository (OSDR) for *Arabidopsis thaliana* multiomics data (RNA-seq, microarray, methylation, proteomics, glycomics) from Col-0 and WS ecotypes grown on the ISS. It builds a multi-factorial mixed-effects model (flight × hardware × organ × light × ecotype), integrates CFD gas-exchange predictions from [LunarLeaf-CFD](https://github.com/dr-richard-barker/LunarLeaf-CFD) as continuous covariates, uses a scVI variational autoencoder trained on the Shahan root atlas (GSE152766) for tissue-specific deconvolution, infers miRNA activity from mRNA targets, and produces a hormone and primary carbon allocation adaptation model.

## Key Results

- **451 spaceflight samples** across 16 NASA OSDR studies, harmonized to 307 expression-profiled samples (30,831 genes)
- **Hardware dominates over flight**: CARA vs BRIC shows 7,675 downregulated genes vs 264 for flight itself
- **CFD gas-exchange modeling**: BRIC limits carbon fixation to ~1% Earth rates in µg; Veggie maintains ~100%
- **Boundary layer conductance** negatively correlates with photosynthesis (r=-0.467), auxin (r=-0.372), ABA (r=-0.352), JA (r=-0.358) — all FDR < 10⁻⁹
- **scVI VAE deconvolution**: 13 root cell types identified from 6,433 single cells; 307 bulk samples deconvolved
- **miRNA inference**: miR164/167/159 targets suppressed in CARA vs BRIC, implicating auxin/ABA regulation
- **28 pathway gene sets** scored across all samples; 154 significant pathway changes (FDR<0.05)

## Structure

```
Airflow_omics/
├── code/
│   ├── osdr_miner/           # OSDR API mining + harmonization scripts
│   │   ├── mine_osdr.py
│   │   ├── enrich_metadata.py
│   │   ├── finalize_metadata.py
│   │   ├── download_data.py
│   │   └── build_harmonized_matrix.py
│   ├── factorial_model/      # limma-voom mixed-effects model
│   │   └── fit_factorial_model.R
│   ├── cfd_covariates/       # CFD covariate extraction from LunarLeaf-CFD
│   │   └── extract_cfd.py
│   ├── sc_vae_deconv/        # scVI VAE training + NNLS deconvolution
│   │   ├── train_vae_deconvolve.py
│   │   └── run_vae_deconv.py
│   ├── mirna_inference/      # miRNA target inference from mRNA DE
│   │   └── infer_mirna_targets.py
│   └── systems_biology/      # Hormone/carbon pathway scoring (GSVA-like)
│       └── hormone_carbon_model.py
├── figures/                  # SVG and PNG figures + generation script
│   ├── generate_figures.py   # Figure generation script
│   └── fig1_study_overview_variance.png (etc.)
├── manuscript/               # LaTeX manuscript sources + compiled PDF
│   ├── manuscript.tex
│   └── manuscript.pdf
├── docs/                     # GitHub Pages landing page
│   ├── index.html
│   └── assets/
├── CITATION.cff              # Citation info
├── LICENSE                   # MIT license
├── PLAN.md                   # Project design and objectives
└── README.md                 # This overview
```

## Requirements

- Python 3.10+ (pandas, numpy, scipy, scvi-tools, scanpy, torch, matplotlib, seaborn, statsmodels)
- R 4.3+ (limma, edgeR, variancePartition, Seurat, anndata)
- LaTeX (texlive-latex-base, texlive-latex-extra)

## Usage

1. `code/osdr_miner/mine_osdr.py` — Mine OSDR API and build harmonized metadata
2. `code/cfd_covariates/extract_cfd.py` — Extract CFD gas-exchange covariates from LunarLeaf-CFD
3. `code/factorial_model/fit_factorial_model.R` — Fit limma-voom mixed-effects model (2 models, 34 contrasts)
4. `code/sc_vae_deconv/run_vae_deconv.py` — Train scVI VAE and deconvolve bulk samples
5. `code/mirna_inference/infer_mirna_targets.py` — Infer miRNA activity from mRNA DE results
6. `code/systems_biology/hormone_carbon_model.py` — Score 28 hormone/carbon pathways, test differential activity
7. `figures/generate_figures.py` — Generate all 8 main figures (SVG + PNG)
8. `manuscript/manuscript.tex` — Compile with pdflatex

## Data Sources

- NASA OSDR: https://osdr.nasa.gov/bio/repo/
- LunarLeaf-CFD: https://github.com/dr-richard-barker/LunarLeaf-CFD
- virtual-root: https://github.com/dr-richard-barker/virtual-root
- Shahan root atlas: GSE152766 (GEO)

## License

- Code: MIT
- Manuscript and figures: CC-BY-4.0

## Citation

```bibtex
@software{barker2026microgravity,
  title={Computational fluid dynamics-guided multi-omics meta-analysis of Arabidopsis spaceflight adaptation},
  author={Barker, Richard},
  year={2026},
  url={https://github.com/dr-richard-barker/microgravity_atmospheric_adaptation}
}
```
