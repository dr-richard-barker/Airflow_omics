#!/usr/bin/env python3
"""
Build harmonized expression matrix from downloaded OSDR processed data.
Combines RNAseq normalized counts and microarray expression across studies,
maps all gene IDs to AGI locus identifiers, and merges with the factorial
metadata.
"""
import pandas as pd, numpy as np, os, re

DATA_DIR = "/mnt/shared-workspace/microgravity_atmospheric_adaptation/data"
RAW_DIR = os.path.join(DATA_DIR, "raw_osdr")
OUTPUT_DIR = DATA_DIR


def load_rnaseq_norm(filepath):
    """Load RNAseq normalized counts, return gene x sample DataFrame."""
    df = pd.read_csv(filepath, index_col=0)
    df.index.name = "gene_id"
    # Ensure AGI format
    df.index = df.index.astype(str).str.strip()
    return df


def load_microarray_norm(filepath):
    """Load microarray normalized expression, extract AGI-mapped probeset values."""
    df = pd.read_csv(filepath)
    # OSD-7 format: TAIR column has AGI IDs, ProbesetID has probe IDs
    # Use TAIR as index, drop annotation columns
    if "TAIR" in df.columns:
        df = df.dropna(subset=["TAIR"])
        df["TAIR"] = df["TAIR"].astype(str).str.strip()
        # Remove rows without valid AGI
        df = df[df["TAIR"].str.match(r"AT[1-5]G\d{5}", na=False)]
        df = df.set_index("TAIR")
        # Drop annotation columns (keep only sample columns)
        annot_cols = [c for c in df.columns if c in
                      ["SYMBOL","GENENAME","REFSEQ","ENTREZID","STRING_id",
                       "GOSLIM_IDS","ProbesetID","TAIR"]]
        sample_cols = [c for c in df.columns if c not in annot_cols]
        df = df[sample_cols]
        # If multiple probesets per gene, take the mean
        df = df.groupby(level=0).mean()
    df.index.name = "gene_id"
    return df


def parse_sample_factors_from_column(colname):
    """Parse sample column name to extract factors."""
    s = str(colname)
    sl = s.lower()
    factors = {"sample_name": s}

    # Ecotype
    if "ws-2" in sl or "ws_2" in sl or "ws-0" in sl or "ws0" in sl or "ws_" in sl or "_ws" in sl or "wassilewskij" in sl:
        factors["ecotype"] = "WS"
    elif "col-0" in sl or "col0" in sl or "col_0" in sl or "wt-col" in sl or "wild type" in sl or "wild_type" in sl:
        factors["ecotype"] = "Col-0"
    elif "ler" in sl:
        factors["ecotype"] = "Ler-0"
    elif "cvi" in sl:
        factors["ecotype"] = "Cvi-0"
    elif "phyd" in sl:
        factors["ecotype"] = "Col-0"  # phyD is a Col-0 mutant
    else:
        factors["ecotype"] = "unknown"

    # Organ
    if "root" in sl:
        factors["organ"] = "root"
    elif "leaf" in sl or "leaves" in sl:
        factors["organ"] = "leaf"
    elif "hypocotyl" in sl:
        factors["organ"] = "hypocotyl"
    elif "shoot" in sl:
        factors["organ"] = "shoot"
    elif "sl-pool" in sl or "sl_pool" in sl or "_sl_" in sl or "seedling" in sl:
        factors["organ"] = "whole_seedling"
    else:
        factors["organ"] = "unknown"

    # Flight
    if "_flt" in sl or "flt_" in sl or "flt-" in sl or "flight" in sl:
        factors["flight"] = "FLT"
    elif "_gc" in sl or "gc_" in sl or "gc-" in sl or "ground" in sl or "ground_control" in sl:
        factors["flight"] = "GC"
    else:
        factors["flight"] = "unknown"

    # Light
    if "dark" in sl:
        factors["light"] = "dark"
    elif "alight" in sl or "a_light" in sl or "_light" in sl or "light_" in sl or "light" in sl:
        factors["light"] = "light"
    elif "red" in sl:
        factors["light"] = "red_light"
    else:
        factors["light"] = "unspecified"

    return factors


def main():
    print("=" * 70)
    print("Building harmonized expression matrix")
    print("=" * 70)

    # Define RNAseq files to load
    rnaseq_files = {
        "OSD-37": "OSD-37/GLDS-37_rna_seq_Normalized_Counts_GLbulkRNAseq.csv",
        "OSD-120": "OSD-120/GLDS-120_rna_seq_Normalized_Counts_GLbulkRNAseq.csv",
        "OSD-217": "OSD-217/GLDS-217_rna_seq_Normalized_Counts.csv",
        "OSD-321": "OSD-321/GLDS-321_rna_seq_Normalized_Counts.csv",
        "OSD-314": "OSD-314/GLDS-314_rna_seq_Normalized_Counts.csv",
        "OSD-678": "OSD-678/GLDS-612_rna_seq_Normalized_Counts_GLbulkRNAseq.csv",
        "OSD-38": "OSD-38/GLDS-38_rna_seq_Normalized_Counts.csv",
    }

    # Microarray files
    array_files = {
        "OSD-7": "OSD-7/GLDS-7_array_normalized_expression_probeset_GLmicroarray.csv",
        "OSD-17": "OSD-17/GLDS-17_array_normalized_expression_probeset_GLmicroarray.csv",
        "OSD-147": "OSD-147/GLDS-147_array_normalized_expression_probeset.csv",
        "OSD-205": "OSD-205/GLDS-205_array_normalized_expression_probeset.csv",
    }

    all_matrices = []
    sample_meta_rows = []

    # Load RNAseq
    print("\n[1/3] Loading RNAseq normalized counts...")
    for osd, relpath in rnaseq_files.items():
        fpath = os.path.join(RAW_DIR, relpath)
        if not os.path.exists(fpath):
            print(f"  {osd}: FILE NOT FOUND ({relpath})")
            continue
        df = load_rnaseq_norm(fpath)
        print(f"  {osd}: {df.shape[0]} genes x {df.shape[1]} samples")

        # Parse factors from column names
        for col in df.columns:
            factors = parse_sample_factors_from_column(col)
            factors["osd_id"] = osd
            factors["omics_type"] = "RNAseq"
            sample_meta_rows.append(factors)

        # Prefix columns with OSD ID to avoid collisions
        df.columns = [f"{osd}__{c}" for c in df.columns]
        all_matrices.append(df)

    # Load microarray
    print("\n[2/3] Loading microarray normalized expression...")
    for osd, relpath in array_files.items():
        fpath = os.path.join(RAW_DIR, relpath)
        if not os.path.exists(fpath):
            print(f"  {osd}: FILE NOT FOUND ({relpath})")
            continue
        df = load_microarray_norm(fpath)
        print(f"  {osd}: {df.shape[0]} genes x {df.shape[1]} samples")

        for col in df.columns:
            factors = parse_sample_factors_from_column(col)
            factors["osd_id"] = osd
            factors["omics_type"] = "microarray"
            sample_meta_rows.append(factors)

        df.columns = [f"{osd}__{c}" for c in df.columns]
        all_matrices.append(df)

    # Merge all matrices on gene ID (outer join)
    print(f"\n[3/3] Merging {len(all_matrices)} matrices...")
    from functools import reduce
    merged = reduce(lambda a, b: pd.merge(a, b, left_index=True, right_index=True, how="outer"), all_matrices)
    print(f"  Merged matrix: {merged.shape[0]} genes x {merged.shape[1]} samples")

    # Fill NaN with 0 for genes not measured in a platform (will handle in model)
    # Actually keep NaN - the model will handle missingness
    merged.index.name = "gene_id"

    # Save harmonized expression matrix
    outpath = os.path.join(OUTPUT_DIR, "harmonized_expression_matrix.tsv")
    # Save as parquet for efficiency (but TSV for Zenodo compatibility)
    merged.to_csv(outpath, sep="\t")
    print(f"\n  Saved harmonized_expression_matrix.tsv ({merged.shape[0]} genes x {merged.shape[1]} samples)")

    # Save sample metadata from expression columns
    sample_meta_df = pd.DataFrame(sample_meta_rows)
    sample_meta_df["full_sample_id"] = sample_meta_df["osd_id"] + "__" + sample_meta_df["sample_name"]
    sample_meta_df.to_csv(os.path.join(OUTPUT_DIR, "expression_sample_metadata.tsv"), sep="\t", index=False)
    print(f"  Saved expression_sample_metadata.tsv ({len(sample_meta_df)} samples)")

    # Summary
    print(f"\n=== Summary ===")
    print(f"  Total genes: {merged.shape[0]}")
    print(f"  Total samples: {merged.shape[1]}")
    print(f"  Studies: {sample_meta_df['osd_id'].nunique()}")
    print(f"  Omics types: {sample_meta_df['omics_type'].value_counts().to_dict()}")
    print(f"  Ecotypes: {sample_meta_df['ecotype'].value_counts().to_dict()}")
    print(f"  Organs: {sample_meta_df['organ'].value_counts().to_dict()}")
    print(f"  Flight: {sample_meta_df['flight'].value_counts().to_dict()}")
    print(f"  Light: {sample_meta_df['light'].value_counts().to_dict()}")


if __name__ == "__main__":
    main()
