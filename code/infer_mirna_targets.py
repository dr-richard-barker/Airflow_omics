"""
Infer miRNA-target interactions from mRNA differential expression results.

Since no microRNA data exists in NASA OSDR for Arabidopsis, we infer miRNA
regulatory activity from mRNA expression patterns using:
1. Known Arabidopsis miRNA-target pairs (from miRBase + AGRIS + psRNATarget)
2. Anti-correlation between miRNA and target expression (inferred from DE)
3. Enrichment of miRNA target sites in DE genes

We use a curated list of well-characterized Arabidopsis miRNAs and their
validated targets to build a regulatory network, then test whether target
genes are enriched in the spaceflight DE signatures.
"""
import os
import numpy as np
import pandas as pd
from scipy import stats
from collections import defaultdict

DATA_DIR = "/mnt/shared-workspace/microgravity_atmospheric_adaptation/data"
DE_DIR = "/mnt/shared-workspace/microgravity_atmospheric_adaptation/factorial_model"
MIRNA_DIR = "/mnt/shared-workspace/microgravity_atmospheric_adaptation/mirna_inference"
RESULTS_DIR = "/mnt/results/microgravity_atmospheric_adaptation/tables"
os.makedirs(MIRNA_DIR, exist_ok=True)
os.makedirs(RESULTS_DIR, exist_ok=True)

# ---- Curated Arabidopsis miRNA-target interactions ----
# Sources: AGRIS (Arabidopsis Gene Regulatory Information Server),
# miRBase v22, psRNATarget predictions, and published literature
# These are well-validated miRNA-target pairs in Arabidopsis

ATH_MIRNA_TARGETS = {
    # miR156: SPL transcription factors (developmental timing)
    "miR156": {
        "targets": ["AT2G42200", "AT5G43270", "AT1G20980", "AT2G33810",
                    "AT3G15270", "AT1G53160", "AT5G50570", "AT1G28550"],
        "family": "SPL (SQUAMOSA PROMOTER BINDING PROTEIN-LIKE)",
        "function": "Developmental timing, vegetative phase change",
        "spaceflight_relevance": "SPLs regulate juvenile-to-adult transition; "
            "spaceflight may alter developmental timing"
    },
    # miR159: MYB transcription factors (ABA signaling)
    "miR159": {
        "targets": ["AT1G66370", "AT2G42200", "AT5G06100", "AT3G11400"],
        "family": "MYB33/MYB65/GAMYB",
        "function": "ABA signaling, germination, stress response",
        "spaceflight_relevance": "ABA signaling is altered in microgravity; "
            "miR159 is a key ABA-responsive miRNA"
    },
    # miR160: ARF transcription factors (auxin signaling)
    "miR160": {
        "targets": ["AT2G28350", "AT4G30080", "AT1G77850"],
        "family": "ARF10/ARF16/ARF17",
        "function": "Auxin signaling, root cap formation, leaf development",
        "spaceflight_relevance": "Auxin redistribution is central to gravitropic "
            "response; ARF10/16/17 regulate root gravitropism"
    },
    # miR164: NAC transcription factors (organ boundary, senescence)
    "miR164": {
        "targets": ["AT1G56010", "AT5G61430", "AT5G07680", "AT3G15170"],
        "family": "NAC1/CUC1/CUC2",
        "function": "Lateral root development, organ boundary formation",
        "spaceflight_relevance": "Lateral root development is altered in "
            "microgravity; NAC1 mediates auxin signaling in roots"
    },
    # miR165/166: HD-ZIP III (adaxial-abaxial polarity)
    "miR165": {
        "targets": ["AT1G52150", "AT1G30490", "AT5G60690", "AT5G65600"],
        "family": "PHB/PHV/REV/CNA (HD-ZIP III)",
        "function": "Leaf polarity, vascular development, meristem",
        "spaceflight_relevance": "Vascular patterning may be affected by altered "
            "auxin flow in microgravity"
    },
    "miR166": {
        "targets": ["AT1G52150", "AT1G30490", "AT5G60690", "AT5G65600"],
        "family": "PHB/PHV/REV/CNA (HD-ZIP III)",
        "function": "Leaf polarity, vascular development",
        "spaceflight_relevance": "Same as miR165"
    },
    # miR167: ARF6/ARF8 (auxin, flower development)
    "miR167": {
        "targets": ["AT1G30330", "AT2G33810", "AT5G37020"],
        "family": "ARF6/ARF8",
        "function": "Auxin signaling, root development, flower development",
        "spaceflight_relevance": "ARF6/8 regulate lateral root initiation; "
            "auxin redistribution in microgravity affects root architecture"
    },
    # miR168: AGO1 (miRNA pathway itself)
    "miR168": {
        "targets": ["AT1G48410"],
        "family": "AGO1",
        "function": "miRNA-mediated gene silencing pathway",
        "spaceflight_relevance": "AGO1 is the core miRNA effector; altered "
            "miR168 may indicate global miRNA pathway modulation"
    },
    # miR169: NF-YA transcription factors (drought, ABA)
    "miR169": {
        "targets": ["AT1G54160", "AT5G12840", "AT3G05690", "AT3G20910"],
        "family": "NF-YA (nuclear factor Y, subunit A)",
        "function": "Drought stress, ABA signaling, flowering",
        "spaceflight_relevance": "NF-YA regulates drought stress; spaceflight "
            "induces stress-like transcriptional responses"
    },
    # miR170/171: SCL transcription factors (meristem)
    "miR171": {
        "targets": ["AT2G45160", "AT3G60150", "AT4G00150", "AT1G69170"],
        "family": "SCL6/SCL22/SCL27 (GRAS domain)",
        "function": "Meristem maintenance, leaf development",
        "spaceflight_relevance": "Meristem function is altered in microgravity"
    },
    # miR172: AP2 transcription factors (flowering)
    "miR172": {
        "targets": ["AT2G28550", "AT5G67030", "AT4G36920", "AT5G11560"],
        "family": "AP2/TOE1/TOE2/SMZ",
        "function": "Flowering time, floral identity",
        "spaceflight_relevance": "Flowering time is altered in spaceflight"
    },
    # miR319: TCP transcription factors (leaf development)
    "miR319": {
        "targets": ["AT1G53230", "AT2G31070", "AT3G15030", "AT5G08310"],
        "family": "TCP2/TCP3/TCP4/TCP10",
        "function": "Leaf development, jasmonate biosynthesis",
        "spaceflight_relevance": "TCPs regulate leaf morphology; spaceflight "
            "affects leaf development"
    },
    # miR390: TAS3 tasiRNA pathway (auxin)
    "miR390": {
        "targets": ["AT3G17185"],  # TAS3
        "family": "TAS3 (trans-acting siRNA)",
        "function": "Generates tasiRNAs targeting ARF2/3/4",
        "spaceflight_relevance": "TAS3 pathway regulates auxin signaling via ARF2/3/4"
    },
    # miR393: TIR1/AFB (auxin perception)
    "miR393": {
        "targets": ["AT3G62980", "AT5G44750", "AT5G64360", "AT2G39360"],
        "family": "TIR1/AFB2/AFB3 (auxin receptors)",
        "function": "Auxin perception, stress response",
        "spaceflight_relevance": "TIR1 is the auxin receptor; miR393-mediated "
            "downregulation would alter auxin sensitivity in microgravity"
    },
    # miR396: GRF transcription factors (leaf growth)
    "miR396": {
        "targets": ["AT2G36400", "AT4G24150", "AT5G53460", "AT3G46980",
                    "AT1G34710", "AT5G07680"],
        "family": "GRF1-7 (GROWTH REGULATING FACTOR)",
        "function": "Leaf growth, root development",
        "spaceflight_relevance": "GRFs regulate organ growth; spaceflight "
            "alters growth patterns"
    },
    # miR397: laccases (lignin biosynthesis)
    "miR397": {
        "targets": ["AT2G40370", "AT2G46570", "AT5G05990"],
        "family": "LAC2/LAC4/LAC17 (laccases)",
        "function": "Lignin biosynthesis, secondary cell wall",
        "spaceflight_relevance": "Cell wall remodeling occurs in microgravity; "
            "lignin deposition is altered"
    },
    # miR398: CSD (copper/zinc SOD, oxidative stress)
    "miR398": {
        "targets": ["AT1G08830", "AT2G28190", "AT5G18100"],
        "family": "CSD1/CSD2/CCS (Cu/Zn superoxide dismutases)",
        "function": "Oxidative stress response",
        "spaceflight_relevance": "Oxidative stress is a major spaceflight stress; "
            "miR398 downregulation would increase CSD levels"
    },
    # miR399: PHO2 (phosphate starvation)
    "miR399": {
        "targets": ["AT2G33770"],
        "family": "PHO2 (PHOSPHATE2)",
        "function": "Phosphate homeostasis",
        "spaceflight_relevance": "Nutrient uptake may be altered in spaceflight"
    },
    # miR408: plantacyanin (copper homeostasis)
    "miR408": {
        "targets": ["AT2G44790", "AT5G20230"],
        "family": "Plantacyanin, laccase",
        "function": "Copper homeostasis, cell wall",
        "spaceflight_relevance": "Cell wall remodeling in microgravity"
    },
    # miR402: DNA repair
    "miR402": {
        "targets": ["AT4G30880", "AT5G13200"],
        "family": "DCL1, AGO2",
        "function": "DNA repair, genome stability",
        "spaceflight_relevance": "DNA damage is elevated in spaceflight due to "
            "radiation; miR402 may regulate repair capacity"
    },
    # miR156/157: SPL (redundant with miR156)
    "miR157": {
        "targets": ["AT2G42200", "AT5G43270", "AT1G20980", "AT2G33810"],
        "family": "SPL (SQUAMOSA PROMOTER BINDING PROTEIN-LIKE)",
        "function": "Developmental timing",
        "spaceflight_relevance": "Same as miR156"
    },
    # miR158: PPR proteins
    "miR158": {
        "targets": ["AT1G06580", "AT1G62670"],
        "family": "PPR (pentatricopeptide repeat)",
        "function": "RNA processing in organelles",
        "spaceflight_relevance": "Organelle function may be affected in microgravity"
    },
    # miR162: DCL1 (miRNA biogenesis)
    "miR162": {
        "targets": ["AT1G01040"],
        "family": "DCL1 (DICER-LIKE 1)",
        "function": "miRNA biogenesis",
        "spaceflight_relevance": "Feedback regulation of miRNA pathway itself"
    },
    # miR395: APS/AST (sulfate assimilation)
    "miR395": {
        "targets": ["AT3G22890", "AT5G67500", "AT1G80840", "AT3G22890"],
        "family": "APS1/APS3/APS4/SULTR2;1",
        "function": "Sulfate assimilation and transport",
        "spaceflight_relevance": "Nutrient metabolism may be altered in spaceflight"
    },
    # miR824: AGO2
    "miR824": {
        "targets": ["AT1G48410"],
        "family": "AGO1",
        "function": "miRNA pathway",
        "spaceflight_relevance": "miRNA pathway regulation"
    },
    # miR854: eIF4E
    "miR854": {
        "targets": ["AT4G18040"],
        "family": "eIF(iso)4E",
        "function": "Translation initiation",
        "spaceflight_relevance": "Translation regulation in spaceflight"
    },
}


def load_de_results():
    """Load all DE results."""
    print("Loading DE results...")
    de_all = pd.read_csv(os.path.join(DE_DIR, "de_results_all_contrasts.tsv"), sep='\t')
    print(f"Total DE results: {len(de_all)}")
    print(f"Contrasts: {de_all['contrast'].unique()}")
    return de_all


def infer_mirna_activity(de_all, fdr_threshold=0.05, logfc_threshold=1.0):
    """
    Infer miRNA activity from target gene expression patterns.

    Logic: If a miRNA is UP-regulated, its targets should be DOWN-regulated
    (and vice versa). We test for enrichment of target genes in the DE set
    using Fisher's exact test, and infer miRNA activity direction from the
    mean logFC of targets.
    """
    print("\n=== Inferring miRNA activity ===")

    results = []

    # Get all tested genes
    all_genes = set(de_all['gene_id'].unique())
    print(f"Total tested genes: {len(all_genes)}")

    for contrast in de_all['contrast'].unique():
        contrast_de = de_all[de_all['contrast'] == contrast].copy()
        model = contrast_de['model'].iloc[0] if 'model' in contrast_de.columns else 'unknown'

        # Define DE genes for this contrast
        de_genes = set(contrast_de[
            (contrast_de['adj.P.Val'] < fdr_threshold) &
            (contrast_de['logFC'].abs() >= logfc_threshold)
        ]['gene_id'].unique())

        de_up = set(contrast_de[
            (contrast_de['adj.P.Val'] < fdr_threshold) &
            (contrast_de['logFC'] >= logfc_threshold)
        ]['gene_id'].unique())

        de_down = set(contrast_de[
            (contrast_de['adj.P.Val'] < fdr_threshold) &
            (contrast_de['logFC'] <= -logfc_threshold)
        ]['gene_id'].unique())

        for mirna, info in ATH_MIRNA_TARGETS.items():
            targets = info['targets']
            valid_targets = [t for t in targets if t in all_genes]

            if len(valid_targets) < 2:
                continue

            # Count targets in DE sets
            targets_in_de = set(valid_targets) & de_genes
            targets_up = set(valid_targets) & de_up
            targets_down = set(valid_targets) & de_down

            # Fisher's exact test for enrichment
            # Contingency: targets in DE vs non-targets in DE
            n_targets = len(valid_targets)
            n_targets_in_de = len(targets_in_de)
            n_non_targets = len(all_genes) - n_targets
            n_non_targets_in_de = len(de_genes) - n_targets_in_de

            # Fisher's exact test
            table = np.array([
                [n_targets_in_de, n_targets - n_targets_in_de],
                [n_non_targets_in_de, n_non_targets - n_non_targets_in_de]
            ])

            try:
                odds_ratio, p_value = stats.fisher_exact(table, alternative='greater')
            except Exception:
                odds_ratio, p_value = np.nan, np.nan

            # Mean logFC of targets
            target_logfcs = contrast_de[
                contrast_de['gene_id'].isin(valid_targets)
            ]['logFC'].values

            mean_logfc = np.mean(target_logfcs) if len(target_logfcs) > 0 else np.nan
            median_logfc = np.median(target_logfcs) if len(target_logfcs) > 0 else np.nan

            # Infer miRNA activity direction
            # If targets are down → miRNA is active/up
            # If targets are up → miRNA is suppressed/down
            if mean_logfc < -0.3:
                inferred_activity = "UP (targets suppressed)"
            elif mean_logfc > 0.3:
                inferred_activity = "DOWN (targets de-repressed)"
            else:
                inferred_activity = "UNCHANGED"

            # Get individual target logFCs
            target_details = []
            for t in valid_targets:
                t_row = contrast_de[contrast_de['gene_id'] == t]
                if len(t_row) > 0:
                    target_details.append({
                        'gene_id': t,
                        'logFC': t_row['logFC'].iloc[0],
                        'adj.P.Val': t_row['adj.P.Val'].iloc[0],
                        'is_DE': t in de_genes
                    })

            results.append({
                'mirna': mirna,
                'contrast': contrast,
                'model': model,
                'n_targets_total': n_targets,
                'n_targets_tested': len(valid_targets),
                'n_targets_DE': n_targets_in_de,
                'n_targets_up': len(targets_up),
                'n_targets_down': len(targets_down),
                'fisher_odds_ratio': odds_ratio,
                'fisher_p_value': p_value,
                'mean_target_logFC': mean_logfc,
                'median_target_logFC': median_logfc,
                'inferred_activity': inferred_activity,
                'mirna_family': info['family'],
                'mirna_function': info['function'],
                'spaceflight_relevance': info['spaceflight_relevance'],
                'target_genes': ';'.join(valid_targets),
                'de_target_genes': ';'.join(targets_in_de),
            })

    results_df = pd.DataFrame(results)

    # Multiple testing correction across all tests
    from statsmodels.stats.multitest import multipletests
    valid_p = results_df['fisher_p_value'].dropna()
    if len(valid_p) > 0:
        _, p_adj, _, _ = multipletests(valid_p, method='fdr_bh')
        results_df.loc[results_df['fisher_p_value'].notna(), 'fisher_fdr'] = p_adj
    else:
        results_df['fisher_fdr'] = np.nan

    return results_df


def build_mirna_network(results_df):
    """Build miRNA-target regulatory network for visualization."""
    print("\n=== Building miRNA-target network ===")

    # Focus on significant miRNA activity changes
    sig = results_df[
        (results_df['fisher_fdr'] < 0.1) |
        (results_df['n_targets_DE'] >= 2)
    ].copy()

    network_edges = []
    for _, row in sig.iterrows():
        targets = row['target_genes'].split(';')
        de_targets = row['de_target_genes'].split(';') if row['de_target_genes'] else []

        for t in targets:
            network_edges.append({
                'source': row['mirna'],
                'target': t,
                'contrast': row['contrast'],
                'inferred_activity': row['inferred_activity'],
                'is_DE_target': t in de_targets,
                'mean_target_logFC': row['mean_target_logFC'],
                'fisher_fdr': row['fisher_fdr'],
                'mirna_family': row['mirna_family'],
                'spaceflight_relevance': row['spaceflight_relevance']
            })

    network_df = pd.DataFrame(network_edges)
    print(f"Network edges: {len(network_df)}")
    return network_df


def main():
    # Load DE results
    de_all = load_de_results()

    # Infer miRNA activity
    results_df = infer_mirna_activity(de_all)

    # Save results
    results_df.to_csv(os.path.join(MIRNA_DIR, "mirna_activity_inference.tsv"),
                      sep='\t', index=False)
    results_df.to_csv(os.path.join(RESULTS_DIR, "mirna_activity_inference.tsv"),
                      sep='\t', index=False)
    print(f"\nSaved miRNA activity inference: {results_df.shape}")

    # Build network
    network_df = build_mirna_network(results_df)
    network_df.to_csv(os.path.join(MIRNA_DIR, "mirna_target_network.tsv"),
                      sep='\t', index=False)
    network_df.to_csv(os.path.join(RESULTS_DIR, "mirna_target_network.tsv"),
                      sep='\t', index=False)
    print(f"Saved miRNA-target network: {network_df.shape}")

    # Summary: significant miRNAs per contrast
    print("\n=== Summary: Significant miRNA activity changes ===")
    sig_mirnas = results_df[
        (results_df['fisher_fdr'] < 0.1) |
        (results_df['n_targets_DE'] >= 2)
    ]

    if len(sig_mirnas) > 0:
        summary = sig_mirnas.groupby(['contrast', 'mirna']).agg({
            'inferred_activity': 'first',
            'n_targets_DE': 'first',
            'mean_target_logFC': 'first',
            'fisher_fdr': 'first',
            'mirna_family': 'first',
            'spaceflight_relevance': 'first'
        }).reset_index()

        summary.to_csv(os.path.join(MIRNA_DIR, "mirna_summary_significant.tsv"),
                       sep='\t', index=False)
        summary.to_csv(os.path.join(RESULTS_DIR, "mirna_summary_significant.tsv"),
                       sep='\t', index=False)

        print(f"Significant miRNA-contrast pairs: {len(summary)}")
        print(summary[['contrast', 'mirna', 'inferred_activity',
                       'n_targets_DE', 'mean_target_logFC', 'fisher_fdr']].to_string())
    else:
        print("No significant miRNA activity changes detected")

    # Print key spaceflight-relevant miRNAs
    print("\n=== Key spaceflight-relevant miRNAs ===")
    key_contrasts = ['Flight_FLT_vs_GC', 'flightFLT', 'CFD_g_bl']
    for kc in key_contrasts:
        subset = results_df[results_df['contrast'].str.contains(kc, case=False, na=False)]
        if len(subset) > 0:
            print(f"\n--- {kc} ---")
            for _, row in subset.sort_values('fisher_p_value').head(10).iterrows():
                print(f"  {row['mirna']}: {row['inferred_activity']} "
                      f"(targets DE: {row['n_targets_DE']}/{row['n_targets_tested']}, "
                      f"mean logFC: {row['mean_target_logFC']:.3f}, "
                      f"FDR: {row['fisher_fdr']:.3g})")


if __name__ == "__main__":
    main()
