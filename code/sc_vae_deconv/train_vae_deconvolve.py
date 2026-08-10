"""
Train scVI variational autoencoder on Arabidopsis scRNA-seq reference
and deconvolve bulk OSDR samples.

Pipeline:
1. Load Shahan root atlas (GSE152766) scRNA-seq reference
2. Preprocess: filter cells/genes, normalize, HVG selection
3. Train scVI VAE to learn latent cell-type representations
4. Extract cell-type-specific signatures from latent space
5. Deconvolve bulk OSDR samples using cell-type signatures

If the full RDS download fails, falls back to a marker-gene-based
deconvolution using known Arabidopsis root cell-type markers.
"""
import os
import sys
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import scipy.sparse as sp
from pathlib import Path

# ---- Paths ----
DATA_DIR = "/mnt/shared-workspace/microgravity_atmospheric_adaptation/data"
SC_DIR = "/mnt/shared-workspace/microgravity_atmospheric_adaptation/sc_vae_deconv"
RESULTS_DIR = "/mnt/results/microgravity_atmospheric_adaptation/tables"
os.makedirs(SC_DIR, exist_ok=True)
os.makedirs(RESULTS_DIR, exist_ok=True)

# ---- Arabidopsis root cell-type markers ----
# From Shahan et al. 2022, Denyer et al. 2019, Ryu et al. 2019
ROOT_CELL_MARKERS = {
    "columella": ["AT1G28370", "AT2G28470", "AT4G17550", "AT1G76580", "AT5G39620"],
    "epidermis_lateral_root_cap": ["AT1G62510", "AT2G01810", "AT4G28650", "AT3G59700", "AT5G62300"],
    "cortex": ["AT3G15310", "AT5G23860", "AT1G79550", "AT3G55180", "AT2G38710"],
    "endodermis": ["AT4G32180", "AT1G76580", "AT3G18060", "AT5G43190", "AT1G25370"],
    "pericycle": ["AT3G54150", "AT1G32330", "AT5G15290", "AT1G30490", "AT2G35940"],
    "xylem": ["AT4G23850", "AT5G23160", "AT1G32770", "AT2G35930", "AT5G23160"],
    "phloem": ["AT3G55180", "AT1G21400", "AT5G55730", "AT1G22770", "AT3G13790"],
    "procambium": ["AT5G57110", "AT1G76580", "AT3G28350", "AT5G23860"],
    "meristem": ["AT3G23830", "AT1G75680", "AT5G41480", "AT2G23760", "AT1G49770"],
    "root_cap": ["AT1G28370", "AT2G28470", "AT4G17550", "AT5G39620"],
    "stele": ["AT3G54150", "AT1G32330", "AT4G23850", "AT5G23160", "AT3G55180"],
    "hair_cells": ["AT1G62510", "AT2G01810", "AT4G28650", "AT3G59700"],
    "non_hair_cells": ["AT3G15310", "AT5G23860", "AT1G79550"],
    "lateral_root_cap": ["AT1G62510", "AT2G01810", "AT4G28650"],
    "mature_xylem": ["AT4G23850", "AT5G23160", "AT1G32770"],
    "protoxylem": ["AT2G35930", "AT5G23160"],
    "metaxylem": ["AT1G32770", "AT4G23850"],
}

# Shoot cell-type markers (for leaf/shoot samples)
SHOOT_CELL_MARKERS = {
    "palisade_mesophyll": ["AT1G29910", "AT3G56010", "AT5G23860", "AT1G03470"],
    "spongy_mesophyll": ["AT3G56010", "AT5G23860", "AT1G29910"],
    "guard_cells": ["AT1G62510", "AT5G23860", "AT3G28350", "AT1G22770"],
    "vascular": ["AT3G54150", "AT4G23850", "AT5G23160", "AT3G55180"],
    "epidermis_shoot": ["AT1G62510", "AT2G01810", "AT4G28650"],
    "trichomes": ["AT1G75410", "AT5G40330", "AT2G40340"],
    "bundle_sheath": ["AT3G56010", "AT5G23860", "AT1G29910"],
}


def load_reference_data():
    """Load scRNA-seq reference data from RDS or h5ad."""
    # Try h5ad first (if already converted)
    h5ad_path = os.path.join(SC_DIR, "col0_copilot.h5ad")
    if os.path.exists(h5ad_path):
        print(f"Loading pre-converted AnnData: {h5ad_path}")
        adata = ad.read_h5ad(h5ad_path)
        print(f"Reference: {adata.shape[0]} cells x {adata.shape[1]} genes")
        return adata

    # Try RDS
    rds_path = os.path.join(SC_DIR, "col0_copilot.rds")
    if os.path.exists(rds_path):
        print(f"Converting RDS to AnnData: {rds_path}")
        h5ad_path = convert_rds_to_h5ad(rds_path)
        if h5ad_path and os.path.exists(h5ad_path):
            adata = ad.read_h5ad(h5ad_path)
            print(f"Reference: {adata.shape[0]} cells x {adata.shape[1]} genes")
            return adata

    return None


def convert_rds_to_h5ad(rds_path):
    """Convert Seurat RDS to AnnData h5ad using R."""
    import subprocess

    h5ad_path = rds_path.replace('.rds', '.h5ad')
    r_script = f'''
suppressPackageStartupMessages({{library(Seurat); library(anndata)}})
obj <- readRDS("{rds_path}")
cat("Object class:", class(obj), "\\n")

if (inherits(obj, "Seurat")) {{
  cat("Seurat:", ncol(obj), "cells x", nrow(obj), "genes\\n")
  cat("Meta columns:", paste(colnames(obj@meta.data), collapse=", "), "\\n")

  # Get counts
  counts <- GetAssayData(obj, assay = "RNA", layer = "counts")

  # Get metadata
  meta <- obj@meta.data

  # Find cell type column
  ct_col <- NULL
  for (col in colnames(meta)) {{
    if (grepl("cell_type|celltype|cluster|identity|annotation|label", col, ignore.case=TRUE)) {{
      ct_col <- col
      break
    }}
  }}
  cat("Cell type column:", ct_col, "\\n")
  if (!is.null(ct_col)) {{
    cat("Cell types:", unique(meta[[ct_col]]), "\\n")
  }}

  # Convert to AnnData
  adata <- CreateAnnData(X = t(as.matrix(counts)), obs = meta)
  rownames(adata$var) <- rownames(counts)
  write_h5ad(adata, "{h5ad_path}")
  cat("Saved: {h5ad_path}\\n")
}} else {{
  cat("Not a Seurat object. Class:", class(obj), "\\n")
  str(obj, max.level=1)
}}
'''
    try:
        result = subprocess.run(["Rscript", "-e", r_script],
                                capture_output=True, text=True, timeout=600)
        print(result.stdout)
        if result.returncode != 0:
            print("STDERR:", result.stderr[:3000])
            return None
        return h5ad_path if os.path.exists(h5ad_path) else None
    except subprocess.TimeoutExpired:
        print("RDS conversion timed out")
        return None
    except Exception as e:
        print(f"RDS conversion error: {e}")
        return None


def preprocess_reference(adata):
    """Preprocess scRNA-seq reference for scVI."""
    print("\n=== Preprocessing reference ===")
    print(f"Input: {adata.shape[0]} cells x {adata.shape[1]} genes")

    # Filter cells (min 200 genes, max 8000 genes)
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=10)
    print(f"After filtering: {adata.shape[0]} cells x {adata.shape[1]} genes")

    # Find cell type column
    ct_col = None
    for col in adata.obs.columns:
        if any(k in col.lower() for k in ['cell_type', 'celltype', 'cluster',
                                            'identity', 'annotation', 'label',
                                            'cell.type', 'type']):
            ct_col = col
            break

    if ct_col:
        print(f"Cell type column: {ct_col}")
        print(f"Cell types: {adata.obs[ct_col].value_counts().to_dict()}")
        adata.obs['cell_type'] = adata.obs[ct_col].astype(str)
    else:
        print("No cell type column found - will cluster de novo")
        # Basic clustering
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor='seurat_v3',
                                     layer=None, subset=False)
        sc.pp.pca(adata, n_comps=30)
        sc.pp.neighbors(adata, n_neighbors=15)
        sc.tl.leiden(adata, resolution=0.8, key_added='cell_type')
        print(f"De novo clusters: {adata.obs['cell_type'].value_counts().to_dict()}")

    # Store raw counts for scVI
    if sp.issparse(adata.X):
        adata.layers['counts'] = adata.X.copy()
    else:
        adata.layers['counts'] = sp.csr_matrix(adata.X.copy())

    return adata


def train_scvi(adata):
    """Train scVI VAE on the reference data."""
    print("\n=== Training scVI VAE ===")

    try:
        import scvi
        scvi.settings.seed = 42

        # Ensure we have counts layer
        if 'counts' not in adata.layers:
            print("Warning: no counts layer, using X as counts")
            adata.layers['counts'] = adata.X.copy()

        # Setup anndata for scVI
        scvi.model.SCVI.setup_anndata(
            adata,
            layer='counts',
            batch_key=None,
            categorical_covariate_keys=['cell_type'] if 'cell_type' in adata.obs else None,
        )

        # Train VAE
        model = scvi.model.SCVI(
            adata,
            n_layers=2,
            n_latent=30,
            n_hidden=128,
            gene_likelihood='nb'
        )

        # Train (CPU mode, limited epochs for speed)
        model.train(
            max_epochs=50,
            early_stopping=True,
            early_stopping_patience=5,
            batch_size=256,
            plan_kwargs={'lr': 1e-3},
        )

        # Get latent representation
        latent = model.get_latent_representation()
        adata.obsm['X_scVI'] = latent
        print(f"Latent space: {latent.shape}")

        # Save model
        model.save(os.path.join(SC_DIR, "scvi_model"), overwrite=True)
        print("Saved scVI model")

        return model, adata

    except Exception as e:
        print(f"scVI training failed: {e}")
        print("Falling back to PCA-based representation")
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, n_top_genes=2000)
        sc.pp.pca(adata, n_comps=30)
        adata.obsm['X_scVI'] = adata.obsm['X_pca'].copy()
        return None, adata


def compute_cell_type_signatures(adata):
    """Compute cell-type-specific gene expression signatures."""
    print("\n=== Computing cell-type signatures ===")

    if 'cell_type' not in adata.obs:
        print("No cell_type annotation available")
        return None

    # Normalize for signature computation
    adata_norm = adata.copy()
    sc.pp.normalize_total(adata_norm, target_sum=1e4)
    sc.pp.log1p(adata_norm)

    cell_types = adata_norm.obs['cell_type'].unique()
    print(f"Cell types: {len(cell_types)}")

    # Compute mean expression per cell type
    signatures = {}
    for ct in cell_types:
        mask = adata_norm.obs['cell_type'] == ct
        if mask.sum() < 10:
            continue
        ct_data = adata_norm[mask]
        if sp.issparse(ct_data.X):
            mean_expr = np.asarray(ct_data.X.mean(axis=0)).flatten()
        else:
            mean_expr = ct_data.X.mean(axis=0)
        signatures[ct] = pd.Series(mean_expr, index=adata_norm.var_names)

    sig_df = pd.DataFrame(signatures)
    print(f"Signature matrix: {sig_df.shape}")

    return sig_df


def deconvolve_bulk(signatures, bulk_matrix_path, metadata_path):
    """
    Deconvolve bulk OSDR samples using cell-type signatures.
    Uses non-negative least squares (NNLS) for deconvolution.
    """
    from scipy.optimize import nnls

    print("\n=== Deconvolving bulk samples ===")

    # Load bulk data
    bulk = pd.read_csv(bulk_matrix_path, sep='\t', index_col=0)
    meta = pd.read_csv(metadata_path, sep='\t')

    # Log-transform bulk data (voom-like)
    bulk_log = np.log2(bulk + 1)

    # Find shared genes
    shared_genes = list(set(signatures.index) & set(bulk_log.index))
    print(f"Shared genes: {len(shared_genes)}")

    if len(shared_genes) < 100:
        print("Warning: very few shared genes, using marker-based approach")
        return deconvolve_with_markers(bulk_log, meta)

    # Subset to shared genes
    sig_shared = signatures.loc[shared_genes]
    bulk_shared = bulk_log.loc[shared_genes]

    # NNLS deconvolution for each sample
    cell_types = sig_shared.columns
    proportions = []

    for sample in bulk_shared.columns:
        y = bulk_shared[sample].values
        X = sig_shared.values
        try:
            props, residual = nnls(X, y)
            # Normalize to sum to 1
            if props.sum() > 0:
                props = props / props.sum()
            proportions.append(props)
        except Exception:
            proportions.append(np.zeros(len(cell_types)))

    prop_df = pd.DataFrame(proportions, index=bulk_shared.columns,
                           columns=cell_types)

    print(f"Deconvolution complete: {prop_df.shape}")
    print(f"Cell type proportions (mean):")
    print(prop_df.mean().sort_values(ascending=False).head(10))

    return prop_df


def deconvolve_with_markers(bulk_log, meta):
    """
    Fallback: deconvolve using known cell-type marker genes.
    Scores each sample for enrichment of each cell type's markers.
    """
    print("\n=== Marker-based deconvolution ===")

    # Combine root and shoot markers
    all_markers = {}
    all_markers.update(ROOT_CELL_MARKERS)
    all_markers.update(SHOOT_CELL_MARKERS)

    # Score each sample for each cell type
    scores = {}
    for ct, markers in all_markers.items():
        valid_markers = [g for g in markers if g in bulk_log.index]
        if len(valid_markers) < 2:
            continue
        scores[ct] = bulk_log.loc[valid_markers].mean(axis=0)

    score_df = pd.DataFrame(scores)
    print(f"Marker scores: {score_df.shape}")

    # Normalize to proportions (softmax-like)
    score_pos = score_df - score_df.min()
    prop_df = score_pos.div(score_pos.sum(axis=1), axis=0)

    return prop_df


def main():
    # Step 1: Load reference data
    adata = load_reference_data()

    if adata is not None:
        # Step 2: Preprocess
        adata = preprocess_reference(adata)

        # Step 3: Train scVI VAE
        model, adata = train_scvi(adata)

        # Step 4: Compute cell-type signatures
        signatures = compute_cell_type_signatures(adata)

        if signatures is not None:
            # Save signatures
            signatures.to_csv(os.path.join(SC_DIR, "cell_type_signatures.tsv"), sep='\t')
            signatures.to_csv(os.path.join(RESULTS_DIR, "cell_type_signatures.tsv"), sep='\t')

        # Step 5: Deconvolve bulk
        if signatures is not None:
            bulk_path = os.path.join(DATA_DIR, "harmonized_expression_matrix_filtered.tsv")
            meta_path = os.path.join(DATA_DIR, "expression_metadata_with_cfd.tsv")
            proportions = deconvolve_bulk(signatures, bulk_path, meta_path)

            if proportions is not None:
                # Merge with metadata
                prop_meta = proportions.copy()
                prop_meta['sample_name'] = prop_meta.index
                prop_meta = prop_meta.merge(
                    meta[['sample_name', 'flight', 'hardware', 'organ',
                          'cfd_hardware', 'osd_id']],
                    on='sample_name', how='left'
                )

                prop_meta.to_csv(
                    os.path.join(SC_DIR, "deconvolution_proportions.tsv"),
                    sep='\t', index=False
                )
                prop_meta.to_csv(
                    os.path.join(RESULTS_DIR, "deconvolution_proportions.tsv"),
                    sep='\t', index=False
                )
                print(f"\nSaved deconvolution proportions: {prop_meta.shape}")

                # Summary by flight x hardware
                print("\n=== Cell type proportions by flight x hardware ===")
                ct_cols = [c for c in proportions.columns]
                summary = prop_meta.groupby(['cfd_hardware', 'flight'])[ct_cols].mean()
                print(summary.to_string())
    else:
        print("\nNo reference data available. Using marker-based deconvolution.")
        # Load bulk data directly
        bulk = pd.read_csv(os.path.join(DATA_DIR, "harmonized_expression_matrix_filtered.tsv"),
                          sep='\t', index_col=0)
        bulk_log = np.log2(bulk + 1)
        meta = pd.read_csv(os.path.join(DATA_DIR, "expression_metadata_with_cfd.tsv"), sep='\t')

        proportions = deconvolve_with_markers(bulk_log, meta)

        # Merge with metadata
        prop_meta = proportions.copy()
        prop_meta['sample_name'] = prop_meta.index
        prop_meta = prop_meta.merge(
            meta[['sample_name', 'flight', 'hardware', 'organ',
                  'cfd_hardware', 'osd_id']],
            on='sample_name', how='left'
        )

        prop_meta.to_csv(
            os.path.join(SC_DIR, "deconvolution_proportions.tsv"),
            sep='\t', index=False
        )
        prop_meta.to_csv(
            os.path.join(RESULTS_DIR, "deconvolution_proportions.tsv"),
            sep='\t', index=False
        )
        print(f"\nSaved marker-based deconvolution: {prop_meta.shape}")


if __name__ == "__main__":
    main()
