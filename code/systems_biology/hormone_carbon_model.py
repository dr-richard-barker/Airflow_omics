"""
Build tissue/organ-specific hormonal and primary carbon allocation
adaptation model using pathway scoring (GSVA-like) on the harmonized OSDR
expression matrix.

Hormone panel: auxin, cytokinin, ethylene, ABA, GA, JA, SA
Plus: primary carbon metabolism (photosynthesis, starch, sucrose, glycolysis)
"""
import os, sys
import numpy as np
import pandas as pd
from scipy import stats

DATA_DIR = "/mnt/shared-workspace/microgravity_atmospheric_adaptation/data"
DE_DIR = "/mnt/shared-workspace/microgravity_atmospheric_adaptation/factorial_model"
SB_DIR = "/mnt/shared-workspace/microgravity_atmospheric_adaptation/systems_biology"
RESULTS_DIR = "/mnt/results/microgravity_atmospheric_adaptation/tables"
os.makedirs(SB_DIR, exist_ok=True)
os.makedirs(RESULTS_DIR, exist_ok=True)

# ---- Hormone and carbon pathway gene sets (curated from KEGG/AraPath/TAIR) ----
PATHWAY_GENESETS = {
    "auxin_biosynthesis": ["AT4G28640","AT1G70560","AT1G04610","AT5G11320","AT1G15050",
        "AT2G33230","AT4G28200","AT1G73560","AT1G48910","AT1G21400","AT1G23320","AT4G24670"],
    "auxin_signaling": ["AT3G62980","AT5G44750","AT5G64360","AT2G39360","AT1G19850",
        "AT3G15500","AT1G51950","AT1G04550","AT1G80390","AT1G28550","AT3G23030","AT2G46990",
        "AT2G01200","AT4G28640","AT1G76580","AT2G01200","AT1G30330","AT2G33810","AT5G37020",
        "AT2G28350","AT4G30080","AT1G77850","AT1G19220","AT5G60450","AT1G53160"],
    "auxin_transport": ["AT2G01180","AT1G23060","AT1G70940","AT2G01420","AT5G16530",
        "AT1G77110","AT1G23080","AT5G15100","AT5G54790","AT2G21050","AT1G77690",
        "AT3G28860","AT4G28770","AT2G36910"],
    "auxin_response_genes": ["AT1G29490","AT2G23170","AT1G80670","AT4G03400","AT2G14960",
        "AT5G54510","AT1G75390","AT1G75590","AT1G75610","AT1G75620","AT1G75630","AT1G75640",
        "AT1G75650","AT2G21210","AT3G03850","AT4G38850","AT5G18060","AT5G18080","AT5G10990"],
    "cytokinin_biosynthesis": ["AT3G63110","AT4G24160","AT3G59570","AT4G28520","AT5G50850",
        "AT1G25410","AT1G68360","AT3G19160","AT5G20040","AT1G67110","AT2G02130","AT2G28305",
        "AT2G40230","AT5G26340","AT1G33450"],
    "cytokinin_signaling": ["AT2G01830","AT5G35750","AT1G27320","AT2G23790","AT5G26010",
        "AT3G29560","AT5G39340","AT3G60510","AT1G20700","AT5G62920","AT4G31920","AT3G16857",
        "AT1G10470","AT3G48105","AT5G62920","AT1G71910","AT5G22300","AT1G49190","AT1G49200",
        "AT2G27070","AT2G25180","AT1G25330","AT2G01760","AT1G84030","AT5G08410","AT3G56380",
        "AT3G04280","AT5G07210","AT3G57100","AT1G85900"],
    "ethylene_biosynthesis": ["AT1G01480","AT1G01500","AT2G22810","AT5G65800","AT4G11280",
        "AT1G62960","AT2G43790","AT3G49700","AT3G51510","AT1G08470","AT1G05010",
        "AT5G51490","AT1G77330","AT1G12010","AT2G19590","AT1G04350"],
    "ethylene_signaling": ["AT2G19590","AT3G23150","AT1G22170","AT2G40940","AT5G66330",
        "AT3G20770","AT3G20780","AT2G27050","AT5G21120","AT1G73730","AT5G10020","AT5G25350",
        "AT3G15210","AT3G14230","AT2G40940","AT1G28360","AT5G58250","AT3G09600","AT1G01520",
        "AT5G47220","AT1G28370","AT4G18450","AT2G31230","AT1G50720","AT3G14230"],
    "aba_biosynthesis": ["AT5G67030","AT1G78390","AT3G14430","AT1G30100","AT4G19170",
        "AT1G13700","AT5G24740","AT1G05060"],
    "aba_signaling": ["AT1G07870","AT1G72770","AT2G26040","AT1G73000","AT2G38310",
        "AT5G05440","AT2G40330","AT4G01026","AT5G67630","AT5G45130","AT4G33950","AT5G60970",
        "AT2G38960","AT4G17870","AT5G59220","AT3G24220","AT1G17550","AT5G08350","AT3G11410",
        "AT5G57170","AT3G56850","AT1G78590","AT2G28900","AT3G54370","AT1G45249","AT4G34000",
        "AT3G19290","AT2G36270","AT3G54370","AT1G49720","AT3G61890"],
    "ga_biosynthesis": ["AT1G79460","AT1G60980","AT5G25900","AT1G15550","AT1G43160",
        "AT1G30040","AT4G25420","AT5G51810","AT5G07200","AT1G78370","AT1G80340","AT1G80330",
        "AT4G25420","AT1G02400","AT1G30040","AT2G34550","AT1G47990","AT1G02400","AT1G60980"],
    "ga_signaling": ["AT3G03450","AT1G14920","AT5G27620","AT1G52830","AT3G03450",
        "AT5G17490","AT1G66970","AT2G01570","AT5G27920"],
    "ja_biosynthesis": ["AT1G72520","AT3G45140","AT1G17440","AT1G72550","AT3G22400",
        "AT1G67530","AT2G06050","AT5G42650","AT3G25760","AT3G25770","AT3G25780","AT1G20510",
        "AT2G25220","AT1G04620","AT2G25220","AT4G23600","AT1G72520"],
    "ja_signaling": ["AT1G17380","AT1G74950","AT3G17860","AT1G48500","AT1G17380",
        "AT1G17380","AT2G34600","AT1G17380","AT1G17380","AT5G13220","AT3G43440","AT5G20900",
        "AT1G19180","AT1G19910","AT5G46760","AT1G71030"],
    "sa_biosynthesis": ["AT2G14610","AT3G10340","AT5G04230","AT3G47340","AT3G47340",
        "AT4G31840","AT2G37040","AT1G18870","AT1G18870","AT4G31840","AT1G75000"],
    "sa_signaling": ["AT1G64280","AT4G26110","AT5G19640","AT4G19670","AT1G74060",
        "AT5G06950","AT1G22070","AT5G10030","AT2G41870","AT1G10120","AT5G06950","AT1G22070",
        "AT5G06950"],
    "photosynthesis_light_reactions": ["ATCG00020","ATCG00270","ATCG00280","ATCG00680",
        "ATCG00220","ATCG00510","ATCG00440","ATCG00030","ATCG00730","ATCG00740","ATCG00750",
        "ATCG00760","ATCG00780","ATCG00500","ATCG00670","ATCG00350","ATCG00360","ATCG00420",
        "ATCG00410","ATCG00490","ATCG00480","AT1G30380","AT1G44575","AT1G06680","AT4G03280",
        "ATCG00120","ATCG00770","ATCG00130","ATCG00480","ATCG00690","ATCG00530"],
    "calvin_cycle": ["AT3G55800","AT5G38430","AT5G38420","AT5G38410","AT1G67090",
        "AT2G39730","AT3G23720","AT3G01530","AT1G42970","AT3G54050","AT2G21170","AT3G01450",
        "AT1G56190","AT3G26650","AT2G36460","AT1G13440","AT1G12900","AT3G55800","AT1G32060"],
    "starch_metabolism": ["AT1G32900","AT5G19220","AT1G05610","AT1G27710","AT4G39210",
        "AT2G21590","AT1G10740","AT5G24300","AT5G65700","AT4G18240","AT1G69830","AT4G09020",
        "AT2G39930","AT4G09020","AT5G04960","AT3G23920","AT4G17090","AT5G64860","AT1G69830",
        "AT1G62330","AT4G25010","AT5G26570","AT2G40840","AT5G64860","AT3G46970","AT5G65700",
        "AT1G10740","AT1G03260"],
    "sucrose_metabolism": ["AT1G73260","AT5G49190","AT4G15210","AT3G43190","AT5G37180",
        "AT2G47180","AT1G12240","AT1G22710","AT2G02860","AT1G09960","AT1G71880","AT5G06170",
        "AT1G66570","AT2G35800","AT3G43190","AT4G09510","AT3G06500","AT3G13784","AT1G55580",
        "AT2G36190","AT3G13784","AT1G72000"],
    "glycolysis": ["AT1G12900","AT2G36460","AT1G13440","AT3G01450","AT1G56190",
        "AT3G23720","AT2G21170","AT1G79550","AT3G12560","AT2G19860","AT4G29130","AT1G50460",
        "AT4G26270","AT2G15220","AT5G57830","AT4G26270","AT5G61540","AT3G52930","AT4G38970",
        "AT2G01130","AT2G30950","AT3G54050","AT1G42970","AT5G08670","AT1G55560","AT2G29560",
        "AT1G74030","AT5G12040","AT5G56350","AT3G22960","AT5G56350"],
    "trehalose_metabolism": ["AT1G23870","AT1G70290","AT1G16980","AT4G27550","AT1G78280",
        "AT4G24040","AT4G39750","AT2G18700","AT5G10100","AT5G51460","AT3G02220","AT1G73970",
        "AT5G25140","AT1G60140"],
    "cell_wall_biosynthesis": ["AT4G18780","AT5G64740","AT5G09870","AT2G21770","AT5G17420",
        "AT5G44030","AT1G02480","AT3G03050","AT4G13430","AT2G32610","AT1G61180","AT3G28180",
        "AT2G20750","AT5G22740","AT4G14130","AT5G57550","AT1G10550","AT4G03210","AT1G14720",
        "AT2G06850","AT5G39260","AT1G69530","AT2G37640","AT3G29030","AT1G65680","AT1G12560",
        "AT2G40610","AT5G56980","AT1G26770","AT3G15370","AT3G15379","AT5G56980","AT4G17030",
        "AT1G71240","AT5G07540"],
    "lignin_biosynthesis": ["AT4G36220","AT2G30490","AT5G54160","AT3G21230","AT1G51680",
        "AT1G65060","AT1G80820","AT3G19730","AT5G48930","AT1G15950","AT1G80820","AT4G34230",
        "AT1G72680","AT3G19450","AT4G36220","AT5G19760","AT2G40370","AT2G46570","AT5G05990",
        "AT5G48930","AT1G71695","AT3G21230","AT5G54160"],
    "oxidative_stress": ["AT1G08830","AT2G28190","AT3G10920","AT5G18100","AT4G25100",
        "AT5G23670","AT1G20630","AT4G35090","AT1G20620","AT1G77120","AT3G09640","AT4G35000",
        "AT4G09010","AT5G39580","AT1G32350","AT2G42830","AT3G11630","AT1G75270","AT5G16710",
        "AT1G28480","AT2G30810","AT4G02520","AT2G29490","AT2G29420","AT2G29450","AT1G02920",
        "AT1G78370","AT2G31570","AT2G43350"],
    "heat_shock": ["AT5G12030","AT5G52640","AT3G12580","AT1G56410","AT5G02500",
        "AT3G09440","AT2G32120","AT5G51440","AT1G16030","AT4G21840","AT2G32120","AT5G52640",
        "AT5G52640","AT5G52640","AT5G52640","AT5G52640","AT5G52640","AT5G52640","AT2G32120",
        "AT5G52640","AT5G52640","AT5G52640","AT5G52640","AT5G52640","AT5G52640","AT5G52640",
        "AT5G52640","AT5G52640","AT5G52640","AT5G52640","AT5G52640","AT5G52640"],
    "dna_damage_response": ["AT1G65480","AT4G21070","AT5G40820","AT5G13810","AT3G48190",
        "AT5G45160","AT3G13470","AT5G20850","AT3G48190","AT1G77510","AT2G31320","AT4G35740",
        "AT3G02490","AT5G60530","AT1G75240"],
    "gravitropism_signaling": ["AT2G03840","AT1G23060","AT2G01180","AT5G15100",
        "AT3G62980","AT2G39360","AT1G70560","AT4G28640","AT1G21400","AT4G24670","AT2G01810",
        "AT4G34390","AT3G25850","AT2G20190","AT3G14370","AT5G25850","AT1G22770","AT3G46600",
        "AT5G59290","AT3G46600"],
}


def load_expression_data():
    print("Loading expression data...")
    expr = pd.read_csv(os.path.join(DATA_DIR, "harmonized_expression_matrix_filtered.tsv"),
                       sep='\t', index_col=0)
    meta = pd.read_csv(os.path.join(DATA_DIR, "expression_metadata_with_cfd.tsv"), sep='\t')
    print(f"Expression: {expr.shape[0]} genes x {expr.shape[1]} samples")
    print(f"Metadata: {meta.shape[0]} samples")
    return expr, meta


def compute_gsva_scores(expr, meta):
    print("\n=== Computing pathway scores (GSVA-like) ===")
    expr_log = np.log2(expr + 1)
    pathway_scores = {}
    for pathway_name, gene_list in PATHWAY_GENESETS.items():
        valid_genes = [g for g in gene_list if g in expr_log.index]
        if len(valid_genes) < 3:
            print(f"  {pathway_name}: only {len(valid_genes)} genes found, skipping")
            continue
        pathway_expr = expr_log.loc[valid_genes]
        gene_means = pathway_expr.mean(axis=1)
        gene_stds = pathway_expr.std(axis=1)
        gene_stds[gene_stds == 0] = 1
        z_scored = pathway_expr.subtract(gene_means, axis=0).div(gene_stds, axis=0)
        scores = z_scored.mean(axis=0)
        pathway_scores[pathway_name] = scores
        print(f"  {pathway_name}: {len(valid_genes)} genes, "
              f"score range [{scores.min():.3f}, {scores.max():.3f}]")
    scores_df = pd.DataFrame(pathway_scores)
    print(f"\nPathway score matrix: {scores_df.shape}")
    return scores_df


def test_differential_pathways(scores_df, meta):
    print("\n=== Testing differential pathway activity ===")
    # Work on a copy to avoid in-place modification issues
    sdf = scores_df.copy()
    sdf['full_sample_id'] = sdf.index.astype(str)
    merged = sdf.merge(
        meta[['full_sample_id', 'flight', 'cfd_hardware', 'organ',
              'light', 'ecotype', 'osd_id']],
        on='full_sample_id', how='left')
    pathway_cols = [c for c in sdf.columns if c not in ('full_sample_id',)]
    print(f"  Merged: {merged.shape}, pathway cols: {len(pathway_cols)}")
    print(f"  Flight distribution: {merged['flight'].value_counts().to_dict()}")
    results = []
    for pathway in pathway_cols:
        # Flight effect
        flt = merged[merged['flight'] == 'FLT'][pathway].dropna()
        gc = merged[merged['flight'] == 'GC'][pathway].dropna()
        if len(flt) > 5 and len(gc) > 5:
            t_stat, p_val = stats.ttest_ind(flt, gc)
            results.append({'pathway': pathway, 'contrast': 'Flight_FLT_vs_GC',
                'mean_FLT': flt.mean(), 'mean_GC': gc.mean(),
                'diff': flt.mean() - gc.mean(), 't_stat': t_stat,
                'p_value': p_val, 'n_FLT': len(flt), 'n_GC': len(gc)})
        # Hardware effect
        for hw in ['BRIC', 'CARA', 'VEGGIE']:
            hw_data = merged[merged['cfd_hardware'] == hw][pathway].dropna()
            other = merged[merged['cfd_hardware'] != hw][pathway].dropna()
            if len(hw_data) > 5 and len(other) > 5:
                t_stat, p_val = stats.ttest_ind(hw_data, other)
                results.append({'pathway': pathway, 'contrast': f'HW_{hw}_vs_other',
                    'mean_group1': hw_data.mean(), 'mean_group2': other.mean(),
                    'diff': hw_data.mean() - other.mean(), 't_stat': t_stat,
                    'p_value': p_val, 'n_group1': len(hw_data), 'n_group2': len(other)})
        # Flight x Hardware
        for hw in ['BRIC', 'CARA', 'VEGGIE']:
            flt_hw = merged[(merged['flight'] == 'FLT') & (merged['cfd_hardware'] == hw)][pathway].dropna()
            gc_hw = merged[(merged['flight'] == 'GC') & (merged['cfd_hardware'] == hw)][pathway].dropna()
            if len(flt_hw) > 3 and len(gc_hw) > 3:
                t_stat, p_val = stats.ttest_ind(flt_hw, gc_hw)
                results.append({'pathway': pathway, 'contrast': f'Flight_{hw}_FLT_vs_GC',
                    'mean_FLT': flt_hw.mean(), 'mean_GC': gc_hw.mean(),
                    'diff': flt_hw.mean() - gc_hw.mean(), 't_stat': t_stat,
                    'p_value': p_val, 'n_FLT': len(flt_hw), 'n_GC': len(gc_hw)})
        # Organ effect
        for organ in ['root', 'leaf', 'shoot', 'whole_seedling']:
            organ_data = merged[merged['organ'] == organ][pathway].dropna()
            other = merged[merged['organ'] != organ][pathway].dropna()
            if len(organ_data) > 5 and len(other) > 5:
                t_stat, p_val = stats.ttest_ind(organ_data, other)
                results.append({'pathway': pathway, 'contrast': f'Organ_{organ}_vs_other',
                    'mean_group1': organ_data.mean(), 'mean_group2': other.mean(),
                    'diff': organ_data.mean() - other.mean(), 't_stat': t_stat,
                    'p_value': p_val, 'n_group1': len(organ_data), 'n_group2': len(other)})
    results_df = pd.DataFrame(results)
    print(f"  Results: {len(results_df)} tests")
    if len(results_df) == 0:
        print("  WARNING: No results generated!")
        return results_df, merged
    from statsmodels.stats.multitest import multipletests
    valid_p = results_df['p_value'].dropna()
    if len(valid_p) > 0:
        _, p_adj, _, _ = multipletests(valid_p, method='fdr_bh')
        results_df.loc[results_df['p_value'].notna(), 'fdr'] = p_adj
    return results_df, merged


def build_hormone_carbon_allocation(scores_df, meta):
    print("\n=== Building hormone-carbon allocation model ===")
    sdf = scores_df.copy()
    sdf['full_sample_id'] = sdf.index.astype(str)
    merged = sdf.merge(
        meta[['full_sample_id', 'flight', 'cfd_hardware', 'organ',
              'light', 'ecotype', 'osd_id', 'g_bl_mol_m2_s',
              'carbon_12h_pct_earth', 'delta_mm']],
        on='full_sample_id', how='left')
    pathway_cols = [c for c in sdf.columns if c not in ('full_sample_id',)]
    summary_hw_flight = merged.groupby(['cfd_hardware', 'flight'])[pathway_cols].mean()
    print("Pathway scores by hardware x flight:")
    print(summary_hw_flight.round(3).to_string())
    summary_organ_flight = merged.groupby(['organ', 'flight'])[pathway_cols].mean()
    cfd_cols = ['g_bl_mol_m2_s', 'carbon_12h_pct_earth', 'delta_mm']
    cfd_correlations = []
    for pathway in pathway_cols:
        for cfd_col in cfd_cols:
            valid = merged[[pathway, cfd_col]].dropna()
            if len(valid) > 10:
                r, p = stats.pearsonr(valid[pathway], valid[cfd_col])
                cfd_correlations.append({'pathway': pathway, 'cfd_covariate': cfd_col,
                    'correlation': r, 'p_value': p, 'n': len(valid)})
    cfd_corr_df = pd.DataFrame(cfd_correlations)
    if len(cfd_corr_df) > 0:
        from statsmodels.stats.multitest import multipletests
        _, p_adj, _, _ = multipletests(cfd_corr_df['p_value'], method='fdr_bh')
        cfd_corr_df['fdr'] = p_adj
    return summary_hw_flight, summary_organ_flight, cfd_corr_df


def main():
    expr, meta = load_expression_data()
    scores_df = compute_gsva_scores(expr, meta)
    scores_df.to_csv(os.path.join(SB_DIR, "pathway_scores_gsva.tsv"), sep='\t')
    scores_df.to_csv(os.path.join(RESULTS_DIR, "pathway_scores_gsva.tsv"), sep='\t')
    print(f"\nSaved pathway scores: {scores_df.shape}")
    diff_results, merged_scores = test_differential_pathways(scores_df, meta)
    diff_results.to_csv(os.path.join(SB_DIR, "differential_pathway_activity.tsv"), sep='\t', index=False)
    diff_results.to_csv(os.path.join(RESULTS_DIR, "differential_pathway_activity.tsv"), sep='\t', index=False)
    print(f"Saved differential pathway results: {diff_results.shape}")
    sig = diff_results[diff_results['fdr'] < 0.05]
    print(f"\nSignificant pathway changes (FDR<0.05): {len(sig)}")
    if len(sig) > 0:
        print(sig[['pathway', 'contrast', 'diff', 'fdr']].sort_values('fdr').to_string())
    summary_hw_flight, summary_organ_flight, cfd_corr = build_hormone_carbon_allocation(scores_df, meta)
    summary_hw_flight.to_csv(os.path.join(SB_DIR, "pathway_summary_hw_flight.tsv"), sep='\t')
    summary_hw_flight.to_csv(os.path.join(RESULTS_DIR, "pathway_summary_hw_flight.tsv"), sep='\t')
    summary_organ_flight.to_csv(os.path.join(SB_DIR, "pathway_summary_organ_flight.tsv"), sep='\t')
    summary_organ_flight.to_csv(os.path.join(RESULTS_DIR, "pathway_summary_organ_flight.tsv"), sep='\t')
    cfd_corr.to_csv(os.path.join(SB_DIR, "pathway_cfd_correlations.tsv"), sep='\t', index=False)
    cfd_corr.to_csv(os.path.join(RESULTS_DIR, "pathway_cfd_correlations.tsv"), sep='\t', index=False)
    print(f"\nSaved pathway summaries and CFD correlations")
    print("\n=== Key hormone-carbon allocation patterns ===")
    print("\nFlight effect on hormone pathways (FLT - GC):")
    hormone_pathways = [p for p in scores_df.columns if any(h in p for h in
                        ['auxin', 'cytokinin', 'ethylene', 'aba', 'ga_', 'ja_', 'sa_'])]
    carbon_pathways = [p for p in scores_df.columns if any(c in p for c in
                        ['photosynthesis', 'calvin', 'starch', 'sucrose', 'glycolysis',
                         'trehalose', 'cell_wall', 'lignin'])]
    for pathway in hormone_pathways + carbon_pathways:
        flt = merged_scores[merged_scores['flight'] == 'FLT'][pathway].dropna()
        gc = merged_scores[merged_scores['flight'] == 'GC'][pathway].dropna()
        if len(flt) > 5 and len(gc) > 5:
            diff = flt.mean() - gc.mean()
            _, p = stats.ttest_ind(flt, gc)
            sig_str = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""
            print(f"  {pathway:40s}: {diff:+.3f} {sig_str}")
    print("\nCFD correlations with hormone/carbon pathways:")
    for _, row in cfd_corr[
        (cfd_corr['pathway'].isin(hormone_pathways + carbon_pathways)) &
        (cfd_corr['fdr'] < 0.05)
    ].sort_values('fdr').iterrows():
        print(f"  {row['pathway']:40s} vs {row['cfd_covariate']:25s}: "
              f"r={row['correlation']:+.3f} (FDR={row['fdr']:.3g})")


if __name__ == "__main__":
    main()
