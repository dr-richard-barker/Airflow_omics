"""
Generate all main figures for the npj Microgravity manuscript.
All figures saved as SVG (and PNG for preview).
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from matplotlib.patches import FancyBboxPatch
import matplotlib.patches as mpatches

# ---- Configuration ----
matplotlib.rcParams['font.family'] = ['Liberation Sans', 'Arimo', 'DejaVu Sans']
matplotlib.rcParams['svg.fonttype'] = 'none'  # Keep SVG text editable
matplotlib.rcParams['pdf.fonttype'] = 42

DATA_DIR = "/mnt/shared-workspace/microgravity_atmospheric_adaptation/data"
DE_DIR = "/mnt/shared-workspace/microgravity_atmospheric_adaptation/factorial_model"
SB_DIR = "/mnt/shared-workspace/microgravity_atmospheric_adaptation/systems_biology"
MIRNA_DIR = "/mnt/shared-workspace/microgravity_atmospheric_adaptation/mirna_inference"
CFD_DIR = "/mnt/shared-workspace/microgravity_atmospheric_adaptation/cfd_covariates"
FIG_DIR = "/mnt/results/microgravity_atmospheric_adaptation/figures"
os.makedirs(FIG_DIR, exist_ok=True)

# Phylo color palette
COLORS = {
    'black': '#000000',
    'cream': '#ECE9E2',
    'white': '#FAF9F3',
    'yellow': '#E9ED4C',
    'orange': '#FF9400',
    'green': '#75A025',
    'pink': '#FD9BED',
    'blue': '#0279EE',
}

# Hardware colors
HW_COLORS = {'BRIC': '#FF9400', 'CARA': '#0279EE', 'VEGGIE': '#75A025'}
FLIGHT_COLORS = {'FLT': '#FD9BED', 'GC': '#0279EE'}


def save_fig(fig, name):
    """Save figure as SVG and PNG."""
    svg_path = os.path.join(FIG_DIR, f"{name}.svg")
    png_path = os.path.join(FIG_DIR, f"{name}.png")
    fig.savefig(svg_path, format='svg', bbox_inches='tight', dpi=300)
    fig.savefig(png_path, format='png', bbox_inches='tight', dpi=300)
    plt.close(fig)
    print(f"  Saved: {name}.svg + {name}.png")


# ============================================================
# Figure 1: Study overview + variance partition
# ============================================================
def fig1_study_overview_variance():
    print("Generating Figure 1: Study overview + variance partition")
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    # Panel A: Sample counts by hardware x flight
    meta = pd.read_csv(os.path.join(DATA_DIR, "expression_metadata_with_cfd.tsv"), sep='\t')
    ct = pd.crosstab(meta['cfd_hardware'], meta['flight'])
    ct = ct.reindex(index=['BRIC', 'CARA', 'VEGGIE'], columns=['GC', 'FLT'], fill_value=0)
    ct.plot(kind='bar', ax=axes[0], color=[FLIGHT_COLORS['GC'], FLIGHT_COLORS['FLT']],
            edgecolor='black', linewidth=0.5, width=0.7)
    axes[0].set_title('A. Samples by hardware and flight', fontsize=11, fontweight='bold')
    axes[0].set_xlabel('Hardware')
    axes[0].set_ylabel('Number of samples')
    axes[0].legend(title='Flight', labels=['Ground control', 'Spaceflight'])
    axes[0].set_xticklabels(axes[0].get_xticklabels(), rotation=0)

    # Panel B: Omics type distribution
    omics_counts = meta['omics_type'].value_counts()
    axes[1].barh(range(len(omics_counts)), omics_counts.values,
                 color=COLORS['green'], edgecolor='black', linewidth=0.5)
    axes[1].set_yticks(range(len(omics_counts)))
    axes[1].set_yticklabels(omics_counts.index)
    axes[1].set_title('B. Omics data types', fontsize=11, fontweight='bold')
    axes[1].set_xlabel('Number of samples')

    # Panel C: Variance partition
    vp = pd.read_csv(os.path.join(DE_DIR, "variance_partition_summary.tsv"), sep='\t')
    vp = vp.sort_values('mean_variance', ascending=True)
    colors_vp = []
    for f in vp['factor']:
        if 'osd_id' in f.lower() or 'study' in f.lower():
            colors_vp.append(COLORS['orange'])
        elif 'light' in f.lower():
            colors_vp.append(COLORS['yellow'])
        elif 'ecotype' in f.lower():
            colors_vp.append(COLORS['pink'])
        elif 'organ' in f.lower():
            colors_vp.append(COLORS['green'])
        elif 'flight' in f.lower():
            colors_vp.append(COLORS['blue'])
        else:
            colors_vp.append(COLORS['black'])

    axes[2].barh(range(len(vp)), vp['mean_variance'],
                 color=colors_vp, edgecolor='black', linewidth=0.5)
    axes[2].set_yticks(range(len(vp)))
    axes[2].set_yticklabels(vp['factor'])
    axes[2].set_title('C. Variance partition', fontsize=11, fontweight='bold')
    axes[2].set_xlabel('Mean variance fraction (%)')

    plt.tight_layout()
    save_fig(fig, "fig1_study_overview_variance")


# ============================================================
# Figure 2: DE volcano plots (flight + hardware)
# ============================================================
def fig2_volcano_plots():
    print("Generating Figure 2: DE volcano plots")
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))

    contrasts = [
        ('Flight_FLT_vs_GC', 'A. Flight: FLT vs GC', 'factorial'),
        ('HW_CARA_vs_BRIC', 'B. Hardware: CARA vs BRIC', 'factorial'),
        ('HW_VEGGIE_vs_BRIC', 'C. Hardware: VEGGIE vs BRIC', 'factorial'),
        ('CFD_g_bl', 'D. CFD: boundary layer conductance', 'cfd'),
    ]

    for idx, (contrast, title, model_type) in enumerate(contrasts):
        ax = axes[idx // 2][idx % 2]
        de = pd.read_csv(os.path.join(DE_DIR, f"de_results_all_contrasts.tsv"), sep='\t')
        de = de[(de['contrast'] == contrast)]

        # Filter significant
        de['neg_log10_p'] = -np.log10(de['adj.P.Val'].clip(lower=1e-300))
        sig = de[(de['adj.P.Val'] < 0.05) & (de['logFC'].abs() >= 1)]
        ns = de[~de.index.isin(sig.index)]

        # Plot
        ax.scatter(ns['logFC'], ns['neg_log10_p'], s=1, alpha=0.1, c='grey', rasterized=True)
        up = sig[sig['logFC'] > 0]
        down = sig[sig['logFC'] < 0]
        ax.scatter(up['logFC'], up['neg_log10_p'], s=3, alpha=0.3, c=COLORS['orange'], rasterized=True)
        ax.scatter(down['logFC'], down['neg_log10_p'], s=3, alpha=0.3, c=COLORS['blue'], rasterized=True)

        ax.set_title(title, fontsize=11, fontweight='bold')
        ax.set_xlabel('log2 fold change')
        ax.set_ylabel('-log10(FDR)')
        ax.axhline(-np.log10(0.05), color='grey', linestyle='--', linewidth=0.5)
        ax.axvline(1, color='grey', linestyle='--', linewidth=0.5)
        ax.axvline(-1, color='grey', linestyle='--', linewidth=0.5)

        # Annotate counts
        ax.text(0.02, 0.98, f'Up: {len(up)}\nDown: {len(down)}',
                transform=ax.transAxes, va='top', fontsize=9,
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    plt.tight_layout()
    save_fig(fig, "fig2_volcano_plots")


# ============================================================
# Figure 3: CFD-pathway correlation heatmap
# ============================================================
def fig3_cfd_pathway_heatmap():
    print("Generating Figure 3: CFD-pathway correlation heatmap")
    cfd_corr = pd.read_csv(os.path.join(SB_DIR, "pathway_cfd_correlations.tsv"), sep='\t')

    # Pivot to matrix
    corr_matrix = cfd_corr.pivot(index='pathway', columns='cfd_covariate', values='correlation')

    # Order pathways by category
    pathway_order = [
        'auxin_biosynthesis', 'auxin_signaling', 'auxin_transport', 'auxin_response_genes',
        'cytokinin_biosynthesis', 'cytokinin_signaling',
        'ethylene_biosynthesis', 'ethylene_signaling',
        'aba_biosynthesis', 'aba_signaling',
        'ga_biosynthesis', 'ga_signaling',
        'ja_biosynthesis', 'ja_signaling',
        'sa_biosynthesis', 'sa_signaling',
        'photosynthesis_light_reactions', 'calvin_cycle',
        'starch_metabolism', 'sucrose_metabolism', 'glycolysis', 'trehalose_metabolism',
        'cell_wall_biosynthesis', 'lignin_biosynthesis',
        'oxidative_stress', 'heat_shock', 'dna_damage_response',
        'gravitropism_signaling',
    ]
    corr_matrix = corr_matrix.reindex(index=[p for p in pathway_order if p in corr_matrix.index])

    # Rename columns
    col_names = {
        'g_bl_mol_m2_s': 'Boundary layer\nconductance (g_bl)',
        'carbon_12h_pct_earth': 'Carbon fixation\n(% Earth)',
        'delta_mm': 'Boundary layer\nthickness (delta)',
    }
    corr_matrix = corr_matrix.rename(columns=col_names)

    fig, ax = plt.subplots(figsize=(8, 12))
    sns.heatmap(corr_matrix, annot=True, fmt='.2f', cmap='RdBu_r', center=0,
                vmin=-0.5, vmax=0.5, ax=ax, linewidths=0.5,
                cbar_kws={'label': 'Pearson correlation'})
    ax.set_title('CFD covariate correlations with hormone & carbon pathways',
                 fontsize=12, fontweight='bold')
    ax.set_ylabel('Pathway')
    ax.set_xlabel('')

    plt.tight_layout()
    save_fig(fig, "fig3_cfd_pathway_heatmap")


# ============================================================
# Figure 4: Hormone-carbon allocation by hardware x flight
# ============================================================
def fig4_hormone_carbon_allocation():
    print("Generating Figure 4: Hormone-carbon allocation")
    summary = pd.read_csv(os.path.join(SB_DIR, "pathway_summary_hw_flight.tsv"), sep='\t')

    # Select key hormone and carbon pathways
    key_pathways = [
        'auxin_signaling', 'auxin_transport',
        'cytokinin_signaling',
        'ethylene_signaling',
        'aba_signaling',
        'ga_signaling',
        'ja_signaling',
        'sa_signaling',
        'photosynthesis_light_reactions', 'calvin_cycle',
        'starch_metabolism', 'glycolysis',
        'cell_wall_biosynthesis', 'lignin_biosynthesis',
        'gravitropism_signaling',
    ]

    fig, axes = plt.subplots(3, 5, figsize=(20, 12))
    axes = axes.flatten()

    for idx, pathway in enumerate(key_pathways):
        if idx >= len(axes):
            break
        ax = axes[idx]
        if pathway not in summary.columns:
            ax.set_visible(False)
            continue

        # Create grouped bar plot
        hw_order = ['BRIC', 'CARA', 'VEGGIE']
        flt_vals = []
        gc_vals = []
        for hw in hw_order:
            row_flt = summary[(summary['cfd_hardware'] == hw) & (summary['flight'] == 'FLT')]
            row_gc = summary[(summary['cfd_hardware'] == hw) & (summary['flight'] == 'GC')]
            flt_vals.append(row_flt[pathway].values[0] if len(row_flt) > 0 else 0)
            gc_vals.append(row_gc[pathway].values[0] if len(row_gc) > 0 else 0)

        x = np.arange(len(hw_order))
        width = 0.35
        ax.bar(x - width/2, gc_vals, width, label='GC', color=FLIGHT_COLORS['GC'],
               edgecolor='black', linewidth=0.5)
        ax.bar(x + width/2, flt_vals, width, label='FLT', color=FLIGHT_COLORS['FLT'],
               edgecolor='black', linewidth=0.5)
        ax.set_xticks(x)
        ax.set_xticklabels(hw_order, fontsize=9)
        ax.set_title(pathway.replace('_', ' '), fontsize=9, fontweight='bold')
        ax.legend(fontsize=7)

    plt.suptitle('Pathway activity by hardware x flight condition',
                 fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    save_fig(fig, "fig4_hormone_carbon_allocation")


# ============================================================
# Figure 5: miRNA-target regulatory network
# ============================================================
def fig5_mirna_network():
    print("Generating Figure 5: miRNA-target network")
    network = pd.read_csv(os.path.join(MIRNA_DIR, "mirna_target_network.tsv"), sep='\t')

    # Focus on flight-related contrasts
    flight_contrasts = ['Flight_FLT_vs_GC', 'HW_CARA_vs_BRIC', 'HW_VEGGIE_vs_BRIC',
                        'FLT_BRIC_vs_CARA', 'FLT_BRIC_vs_VEGGIE', 'CFD_g_bl']
    net = network[network['contrast'].isin(flight_contrasts)]

    # Create a summary: for each miRNA, show inferred activity across contrasts
    summary = net.groupby(['source', 'contrast']).agg({
        'inferred_activity': 'first',
        'mean_target_logFC': 'first',
        'is_DE_target': 'mean',
        'fisher_fdr': 'first'
    }).reset_index()

    # Pivot for heatmap
    activity_map = {'UP (targets suppressed)': -1, 'DOWN (targets de-repressed)': 1, 'UNCHANGED': 0}
    summary['activity_score'] = summary['inferred_activity'].map(activity_map)

    pivot = summary.pivot(index='source', columns='contrast', values='activity_score')

    # Order miRNAs by mean activity
    pivot['mean'] = pivot.mean(axis=1)
    pivot = pivot.sort_values('mean', ascending=True)
    pivot = pivot.drop('mean', axis=1)

    fig, ax = plt.subplots(figsize=(10, 8))
    # Custom colormap: blue (UP/suppressed), white (unchanged), orange (DOWN/de-repressed)
    from matplotlib.colors import ListedColormap, BoundaryNorm
    cmap = ListedColormap([COLORS['blue'], '#FFFFFF', COLORS['orange']])
    bounds = [-1.5, -0.5, 0.5, 1.5]
    norm = BoundaryNorm(bounds, cmap.N)

    sns.heatmap(pivot, cmap=cmap, norm=norm, annot=False, ax=ax,
                linewidths=1, linecolor='grey',
                cbar_kws={'label': 'Inferred miRNA activity',
                          'ticks': [-1, 0, 1],
                          'format': '%d'})
    ax.set_title('Inferred miRNA activity across spaceflight contrasts',
                 fontsize=12, fontweight='bold')
    ax.set_ylabel('miRNA')
    ax.set_xlabel('Contrast')

    # Custom colorbar labels
    cbar = ax.collections[0].colorbar
    cbar.ax.set_yticklabels(['UP\n(targets↓)', 'Unchanged', 'DOWN\n(targets↑)'])
    plt.tight_layout()
    save_fig(fig, "fig5_mirna_network")


# ============================================================
# Figure 6: CFD boundary layer schematic + gas exchange
# ============================================================
def fig6_cfd_gas_exchange():
    print("Generating Figure 6: CFD gas exchange")
    cfd = pd.read_csv(os.path.join(CFD_DIR, "cfd_covariates_by_scenario.tsv"), sep='\t')

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    # Panel A: Carbon fixation by hardware and gravity
    if 'carbon_12h_pct_earth' in cfd.columns:
        pivot = cfd.pivot_table(index='cfd_hardware', columns='gravity_key', values='carbon_12h_pct_earth')
        pivot = pivot.reindex(index=['BRIC', 'CARA', 'VEGGIE'])
        pivot.plot(kind='bar', ax=axes[0], color=[FLIGHT_COLORS['GC'], FLIGHT_COLORS['FLT']],
                   edgecolor='black', linewidth=0.5, width=0.7)
        axes[0].set_title('A. Carbon fixation capacity', fontsize=11, fontweight='bold')
        axes[0].set_xlabel('Hardware')
        axes[0].set_ylabel('Carbon fixation (% Earth)')
        axes[0].legend(title='Gravity', labels=['1g (Earth)', 'µg (ISS)'])
        axes[0].set_xticklabels(axes[0].get_xticklabels(), rotation=0)

    # Panel B: Boundary layer conductance by hardware and gravity
    if 'g_bl_mol_m2_s' in cfd.columns:
        pivot2 = cfd.pivot_table(index='cfd_hardware', columns='gravity_key', values='g_bl_mol_m2_s')
        pivot2 = pivot2.reindex(index=['BRIC', 'CARA', 'VEGGIE'])
        pivot2.plot(kind='bar', ax=axes[1], color=[FLIGHT_COLORS['GC'], FLIGHT_COLORS['FLT']],
                    edgecolor='black', linewidth=0.5, width=0.7)
        axes[1].set_title('B. Boundary layer conductance', fontsize=11, fontweight='bold')
        axes[1].set_xlabel('Hardware')
        axes[1].set_ylabel('g_bl (mol m⁻² s⁻¹)')
        axes[1].legend(title='Gravity', labels=['1g (Earth)', 'µg (ISS)'])
        axes[1].set_xticklabels(axes[1].get_xticklabels(), rotation=0)

    # Panel C: Boundary layer thickness by hardware and gravity
    if 'delta_mm' in cfd.columns:
        pivot3 = cfd.pivot_table(index='cfd_hardware', columns='gravity_key', values='delta_mm')
        pivot3 = pivot3.reindex(index=['BRIC', 'CARA', 'VEGGIE'])
        pivot3.plot(kind='bar', ax=axes[2], color=[FLIGHT_COLORS['GC'], FLIGHT_COLORS['FLT']],
                    edgecolor='black', linewidth=0.5, width=0.7)
        axes[2].set_title('C. Boundary layer thickness', fontsize=11, fontweight='bold')
        axes[2].set_xlabel('Hardware')
        axes[2].set_ylabel('δ (mm)')
        axes[2].legend(title='Gravity', labels=['1g (Earth)', 'µg (ISS)'])
        axes[2].set_xticklabels(axes[2].get_xticklabels(), rotation=0)

    plt.tight_layout()
    save_fig(fig, "fig6_cfd_gas_exchange")


# ============================================================
# Figure 7: Pathway score heatmap (samples x pathways)
# ============================================================
def fig7_pathway_heatmap():
    print("Generating Figure 7: Pathway score heatmap")
    scores = pd.read_csv(os.path.join(SB_DIR, "pathway_scores_gsva.tsv"), sep='\t', index_col=0)
    meta = pd.read_csv(os.path.join(DATA_DIR, "expression_metadata_with_cfd.tsv"), sep='\t')

    # Merge with metadata for ordering
    scores['full_sample_id'] = scores.index
    merged = scores.merge(meta[['full_sample_id', 'flight', 'cfd_hardware', 'organ']],
                          on='full_sample_id', how='left')

    # Order samples by hardware, then flight, then organ
    merged['sort_key'] = (merged['cfd_hardware'].fillna('Z') + '_' +
                          merged['flight'].fillna('Z') + '_' +
                          merged['organ'].fillna('Z'))
    merged = merged.sort_values('sort_key')

    pathway_cols = [c for c in scores.columns if c != 'full_sample_id']
    data = merged[pathway_cols].T

    fig, ax = plt.subplots(figsize=(20, 10))
    sns.heatmap(data, cmap='RdBu_r', center=0, ax=ax,
                xticklabels=False, yticklabels=True,
                cbar_kws={'label': 'GSVA-like z-score'},
                linewidths=0)
    ax.set_title('Pathway activity across all samples (ordered by hardware → flight → organ)',
                 fontsize=12, fontweight='bold')
    ax.set_ylabel('Pathway')
    ax.set_xlabel('Samples (n=307)')

    # Add hardware group annotations
    hw_groups = merged.groupby('cfd_hardware').indices
    for hw, indices in hw_groups.items():
        if len(indices) > 0:
            start = indices[0]
            end = indices[-1] + 1
            ax.axvline(start, color='black', linewidth=1)
            ax.text((start + end) / 2, len(pathway_cols) + 1, hw,
                    ha='center', va='bottom', fontsize=10, fontweight='bold')

    plt.tight_layout()
    save_fig(fig, "fig7_pathway_heatmap")


# ============================================================
# Figure 8: DE summary bar chart
# ============================================================
def fig8_de_summary():
    print("Generating Figure 8: DE summary")
    de_all = pd.read_csv(os.path.join(DE_DIR, "de_results_all_contrasts.tsv"), sep='\t')

    # Count DE genes per contrast
    de_summary = de_all.groupby('contrast').apply(
        lambda x: pd.Series({
            'up': ((x['adj.P.Val'] < 0.05) & (x['logFC'] >= 1)).sum(),
            'down': ((x['adj.P.Val'] < 0.05) & (x['logFC'] <= -1)).sum(),
        })
    ).reset_index()

    de_summary = de_summary.sort_values('up', ascending=True)

    fig, ax = plt.subplots(figsize=(10, 12))
    y = np.arange(len(de_summary))
    ax.barh(y, de_summary['up'], color=COLORS['orange'], label='Upregulated',
            edgecolor='black', linewidth=0.5)
    ax.barh(y, -de_summary['down'], color=COLORS['blue'], label='Downregulated',
            edgecolor='black', linewidth=0.5)
    ax.set_yticks(y)
    ax.set_yticklabels(de_summary['contrast'], fontsize=8)
    ax.set_xlabel('Number of DE genes (FDR < 0.05, |logFC| > 1)')
    ax.set_title('Differential expression summary across all contrasts',
                 fontsize=12, fontweight='bold')
    ax.legend(loc='lower right')
    ax.axvline(0, color='black', linewidth=0.5)

    plt.tight_layout()
    save_fig(fig, "fig8_de_summary")


def main():
    print("=" * 60)
    print("Generating all manuscript figures")
    print("=" * 60)

    fig1_study_overview_variance()
    fig2_volcano_plots()
    fig3_cfd_pathway_heatmap()
    fig4_hormone_carbon_allocation()
    fig5_mirna_network()
    fig6_cfd_gas_exchange()
    fig7_pathway_heatmap()
    fig8_de_summary()

    print("\n" + "=" * 60)
    print(f"All figures saved to {FIG_DIR}")
    print("=" * 60)


if __name__ == "__main__":
    main()
