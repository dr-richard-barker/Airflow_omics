# Airflow_omics: spaceflight transcriptomes and growth-hardware airflow

A study-aware re-analysis of NASA OSDR *Arabidopsis thaliana* spaceflight transcriptomes, coupled to lattice-Boltzmann CFD of the BRIC, CARA and VEGGIE growth hardware ([LunarLeaf-CFD](https://github.com/dr-richard-barker/LunarLeaf-CFD)).

## Headline

**Spaceflight changes Arabidopsis gene expression substantially, but differently in each experiment.**
- Across 11 OSDR studies, a flight response shared by all studies explains under 1 % of gene variance.
- Within BRIC hardware, only the RNA-seq studies share a reproducible flight signature.
- Neither CFD-predicted boundary-layer conductance nor carbon supply explains which studies respond.
- That is consistent with the CFD itself: inside a real flight dish, gravity changes leaf gas exchange far less than open-air models suggest.

The claim, its evidence and what it does not say are in [`reports/HEADLINE_2026-10.md`](reports/HEADLINE_2026-10.md). The full analysis is in [`tables/reanalysis_2026-10/REANALYSIS.md`](tables/reanalysis_2026-10/REANALYSIS.md) (§1–§13).

## Key results

| | Result | Source |
|---|---|---|
| Data | 307 expression samples (153 flight, 154 ground) from 11 OSDR studies: 213 RNA-seq (7 studies) and 94 microarray (4 studies), in four hardware classes (BRIC, CARA, VEGGIE, EMCS). | `tables/reanalysis_2026-10/expression_metadata_with_cfd_v2.tsv` |
| Study-specific response | A flight effect common to all studies explains 0.3 % of gene variance (under 0.01 % in RNA-seq alone); study-specific flight effects 7–11 %; study 44–53 %. | REANALYSIS §12 |
| Flight contrasts | The stored "flight main effect" is the flight effect within BRIC (11,797 genes at FDR < 0.05). Averaged over the four hardware classes it is 2,981 genes. | §6, §9 |
| Platform | BRIC RNA-seq studies agree with each other (median r = 0.40). The BRIC microarray studies agree with nothing (median r ≈ 0) and detect almost no flight genes. | §13 |
| CFD covariates | g_bl takes one value per gravity × organ-scale cell. Study membership explains 86 % of carbon_12h's variance, and its only within-study variation is OSD-7's flight contrast. Neither is identifiable from study and hardware. | §2, §4, §9 |
| Pathways | 0/28 pathways show a flight effect or a CFD dose-response in study-aware models, on correctly computed scores. The tests rule out only large effects. | §10, §11 |
| CFD | In the real BRIC and CARA dishes, gravity changes leaf conductance 1.25–1.8-fold (2.25-fold in an open box). Sealed BRIC dishes become carbon-limited at 1 g as well as in orbit. | LunarLeaf-CFD T19–T22, §7 |

**Superseded results.** The August 2026 draft reported strong CFD–pathway correlations, "hardware dominates over flight" and a 7.1 % hardware variance share. These do not survive the re-analysis; `reports/HEADLINE_2026-10.md` lists each claim and its status.

**Kept for the record, but don't cite:**
- the stored outputs from that pipeline in `tables/`: `de_results_*.tsv`, `pathway_*.tsv`, `variance_partition_summary.tsv`;
- `figures/fig*.png`.

## Structure

```
Airflow_omics/
├── code/
│   ├── osdr_miner/            # OSDR API mining and metadata harmonisation (original pipeline)
│   ├── cfd_covariates/        # CFD covariates from LunarLeaf-CFD, pathway re-analysis (§1–§8, §11)
│   │   ├── extract_cfd.py
│   │   ├── cfd_reanalysis.py
│   │   ├── cfd_meta_regression.py
│   │   └── rescored_comparison.py
│   ├── factorial_model/       # limma-voom models
│   │   ├── fit_factorial_model.R           # original pipeline
│   │   ├── osdr_helpers_2026-10.R          # shared loaders for the 2026-10 refits
│   │   ├── refit_gene_level_2026-10.R      # §9
│   │   ├── followups_2026-10.R             # §10
│   │   ├── variance_partition_2026-10.R    # §12
│   │   └── platform_concordance_2026-10.R  # §13
│   ├── figures/figures_2026-10.py          # figures/2026-10/ and the site figures
│   ├── sc_vae_deconv/         # scVI training and NNLS deconvolution
│   ├── mirna_inference/       # miRNA target inference from mRNA DE
│   └── systems_biology/       # pathway gene sets and scoring
├── tables/reanalysis_2026-10/ # REANALYSIS.md and every table it cites
├── reports/                   # HEADLINE_2026-10.md; CFD_comparison_chamber_model.md
├── figures/2026-10/           # current figures; figures/fig*.png are the superseded August set
├── manuscript/                # LaTeX manuscript and compiled PDF
└── docs/                      # GitHub Pages landing page
```

## Reproduce

The OSDR processed files are downloaded outside the repo, to `~/Documents/Airflow_omics_osdr_data/raw_osdr/`. A `MANIFEST.json` there records the URL and SHA-256 of each file. Run from the repo root, in order:

1. `python3 code/cfd_covariates/cfd_reanalysis.py`. Set `LUNARLEAF_TABLES` (and optionally `LUNARLEAF_OLD_TABLES`) to point at the LunarLeaf-CFD tables.
2. `python3 code/cfd_covariates/cfd_meta_regression.py`
3. `Rscript code/factorial_model/refit_gene_level_2026-10.R`. Set `OSDR_RAW` to point at the downloaded files.
4. `Rscript code/factorial_model/followups_2026-10.R`
5. Rerun steps 1 and 2 with `REANALYSIS_OUT=tables/reanalysis_2026-10/followups/rescored` and `PATHWAY_SCORES=tables/reanalysis_2026-10/followups/pathway_scores_within_study.tsv`, then `python3 code/cfd_covariates/rescored_comparison.py`.
6. `Rscript code/factorial_model/variance_partition_2026-10.R`
7. `Rscript code/factorial_model/platform_concordance_2026-10.R`
8. `python3 code/figures/figures_2026-10.py`

Requirements:
- **Python 3.10+:** pandas, numpy, scipy, matplotlib.
- **R 4.3+:** limma, edgeR, data.table, variancePartition.
- **TeX Live,** to build the manuscript. Build from `manuscript/` with `TEXINPUTS=..:` so that `figures/` resolves.

## Data sources

- NASA OSDR: https://osdr.nasa.gov/bio/repo/
- LunarLeaf-CFD: https://github.com/dr-richard-barker/LunarLeaf-CFD
- Shahan root atlas: GEO GSE152766, used for the deconvolution
- virtual-root: https://github.com/dr-richard-barker/virtual-root

## License

- Code: MIT
- Manuscript and figures: CC-BY-4.0

## Citation

```bibtex
@software{barker2026airflowomics,
  title={Airflow_omics: spaceflight transcriptomes and growth-hardware airflow},
  author={Barker, Richard},
  year={2026},
  url={https://github.com/dr-richard-barker/Airflow_omics}
}
```
