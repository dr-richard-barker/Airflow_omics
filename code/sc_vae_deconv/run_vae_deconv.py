"""
Train scVI VAE on Arabidopsis scRNA-seq reference and deconvolve bulk OSDR samples.
Uses the Shahan root atlas (GSE152766) Col-0 COPILOT sample.
"""
import os
import sys
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import scipy.sparse as sp
from scipy.optimize import nnls

DATA_DIR = "/mnt/shared-workspace/microgravity_atmospheric_adaptation/data"
SC_DIR = "/mnt/shared-workspace/microgravity_atmospheric_adaptation/sc_vae_deconv"
RESULTS_DIR = "/mnt/results/microgravity_atmospheric_adaptation/tables"
os.makedirs(RESULTS_DIR, exist_ok=True)

def main():
    # Step 1: Load reference data
    print("=== Loading scRNA-seq reference ===")
    adata = ad.read_h5ad("/workspace/col0_copilot.h5ad")
    print(f"Reference: {adata.shape[0]} cells x {adata.shape[1]} genes")
    print(f"Cell types: {adata.obs['celltype.anno'].value_counts().to_dict()}")

    # Step 2: Preprocess
    print("\n=== Preprocessing ===")
    # Filter cells
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=10)
    print(f"After filtering: {adata.shape}")

    # Store raw counts
    if sp.issparse(adata.X):
        adata.layers['counts'] = adata.X.copy()
    else:
        adata.layers['counts'] = sp.csr_matrix(adata.X.copy())

    # Use celltype.anno as cell type
    adata.obs['cell_type'] = adata.obs['celltype.anno'].astype(str)
    # Handle NA cells
    adata.obs['cell_type'] = adata.obs['cell_type'].fillna('Unknown')
    print(f"Cell types: {adata.obs['cell_type'].value_counts().to_dict()}")

    # Step 3: Train scVI VAE
    print("\n=== Training scVI VAE ===")
    model = None
    try:
        import scvi
        scvi.settings.seed = 42

        scvi.model.SCVI.setup_anndata(
            adata,
            layer='counts',
            categorical_covariate_keys=['cell_type'],
        )

        model = scvi.model.SCVI(
            adata,
            n_layers=2,
            n_latent=30,
            n_hidden=128,
            gene_likelihood='nb'
        )

        model.train(
            max_epochs=50,
            early_stopping=True,
            early_stopping_patience=5,
            batch_size=256,
            plan_kwargs={'lr': 1e-3},
        )

        latent = model.get_latent_representation()
        adata.obsm['X_scVI'] = latent
        print(f"Latent space: {latent.shape}")
        model.save(os.path.join(SC_DIR, "scvi_model"), overwrite=True)
        print("Saved scVI model")
    except Exception as e:
        print(f"scVI training failed: {e}")
        print("Falling back to PCA representation")
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, n_top_genes=2000)
        sc.pp.pca(adata, n_comps=30)
        adata.obsm['X_scVI'] = adata.obsm['X_pca'].copy()

    # Step 4: Compute cell-type signatures
    print("\n=== Computing cell-type signatures ===")
    adata_norm = adata.copy()
    sc.pp.normalize_total(adata_norm, target_sum=1e4)
    sc.pp.log1p(adata_norm)

    cell_types = adata_norm.obs['cell_type'].unique()
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
    print(f"Cell types in signature: {list(sig_df.columns)}")

    # Save signatures
    sig_df.to_csv(os.path.join(SC_DIR, "cell_type_signatures.tsv"), sep='\t')
    sig_df.to_csv(os.path.join(RESULTS_DIR, "cell_type_signatures.tsv"), sep='\t')

    # Step 5: Deconvolve bulk OSDR samples
    print("\n=== Deconvolving bulk OSDR samples ===")
    bulk = pd.read_csv(os.path.join(DATA_DIR, "harmonized_expression_matrix_filtered.tsv"),
                       sep='\t', index_col=0)
    meta = pd.read_csv(os.path.join(DATA_DIR, "expression_metadata_with_cfd.tsv"), sep='\t')

    # Log-transform bulk
    bulk_log = np.log2(bulk + 1)

    # Find shared genes
    shared_genes = list(set(sig_df.index) & set(bulk_log.index))
    print(f"Shared genes: {len(shared_genes)}")

    if len(shared_genes) < 100:
        print("WARNING: Very few shared genes, using marker-based approach")
        return deconvolve_with_markers(bulk_log, meta)

    # Subset to shared genes
    sig_shared = sig_df.loc[shared_genes]
    bulk_shared = bulk_log.loc[shared_genes]

    # NNLS deconvolution
    cell_types_list = sig_shared.columns
    proportions = []

    for i, sample in enumerate(bulk_shared.columns):
        y = bulk_shared[sample].values
        X = sig_shared.values
        try:
            props, residual = nnls(X, y)
            if props.sum() > 0:
                props = props / props.sum()
            proportions.append(props)
        except Exception:
            proportions.append(np.zeros(len(cell_types_list)))
        if (i + 1) % 50 == 0:
            print(f"  Deconvolved {i+1}/{len(bulk_shared.columns)} samples")

    prop_df = pd.DataFrame(proportions, index=bulk_shared.columns,
                           columns=cell_types_list)
    print(f"\nDeconvolution complete: {prop_df.shape}")
    print(f"Cell type proportions (mean):")
    print(prop_df.mean().sort_values(ascending=False))

    # Merge with metadata
    prop_meta = prop_df.copy()
    prop_meta['full_sample_id'] = prop_meta.index
    prop_meta = prop_meta.merge(
        meta[['full_sample_id', 'flight', 'hardware', 'cfd_hardware',
              'organ', 'osd_id', 'ecotype']],
        on='full_sample_id', how='left'
    )

    # Save
    prop_meta.to_csv(os.path.join(SC_DIR, "deconvolution_proportions.tsv"),
                     sep='\t', index=False)
    prop_meta.to_csv(os.path.join(RESULTS_DIR, "deconvolution_proportions.tsv"),
                     sep='\t', index=False)
    print(f"\nSaved deconvolution proportions: {prop_meta.shape}")

    # Summary by flight x hardware
    print("\n=== Cell type proportions by flight x hardware ===")
    ct_cols = [c for c in prop_df.columns]
    summary = prop_meta.groupby(['cfd_hardware', 'flight'])[ct_cols].mean()
    print(summary.round(3).to_string())

    summary.to_csv(os.path.join(SC_DIR, "deconv_summary_hw_flight.tsv"), sep='\t')
    summary.to_csv(os.path.join(RESULTS_DIR, "deconv_summary_hw_flight.tsv"), sep='\t')

    # Summary by organ
    print("\n=== Cell type proportions by organ ===")
    summary_organ = prop_meta.groupby('organ')[ct_cols].mean()
    print(summary_organ.round(3).to_string())

    summary_organ.to_csv(os.path.join(RESULTS_DIR, "deconv_summary_organ.tsv"), sep='\t')

    return prop_meta


def deconvolve_with_markers(bulk_log, meta):
    """Fallback: deconvolve using known cell-type marker genes."""
    print("Using marker-based deconvolution")
    ROOT_MARKERS = {
        "columella": ["AT1G28370", "AT2G28470", "AT4G17550"],
        "epidermis_lrc": ["AT1G62510", "AT2G01810", "AT4G28650"],
        "cortex": ["AT3G15310", "AT5G23860", "AT1G79550"],
        "endodermis": ["AT4G32180", "AT3G18060", "AT5G43190"],
        "pericycle": ["AT3G54150", "AT1G32330", "AT5G15290"],
        "xylem": ["AT4G23850", "AT5G23160", "AT1G32770"],
        "phloem": ["AT3G55180", "AT1G21400", "AT5G55730"],
        "procambium": ["AT5G57110", "AT3G28350"],
        "meristem": ["AT3G23830", "AT1G75680", "AT5G41480"],
    }

    scores = {}
    for ct, markers in ROOT_MARKERS.items():
        valid = [g for g in markers if g in bulk_log.index]
        if len(valid) < 2:
            continue
        scores[ct] = bulk_log.loc[valid].mean(axis=0)

    score_df = pd.DataFrame(scores)
    score_pos = score_df - score_df.min()
    prop_df = score_pos.div(score_pos.sum(axis=1), axis=0)

    prop_meta = prop_df.copy()
    prop_meta['full_sample_id'] = prop_meta.index
    prop_meta = prop_meta.merge(
        meta[['full_sample_id', 'flight', 'cfd_hardware', 'organ', 'osd_id']],
        on='full_sample_id', how='left'
    )
    prop_meta.to_csv(os.path.join(RESULTS_DIR, "deconvolution_proportions.tsv"),
                     sep='\t', index=False)
    print(f"Saved marker-based deconvolution: {prop_meta.shape}")
    return prop_meta


if __name__ == "__main__":
    main()
