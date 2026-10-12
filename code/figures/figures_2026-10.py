"""
Figures for the 2026-10 rewrite (manuscript, README, site).

Every panel is drawn from a committed table in this repo; nothing is read from outside it.
Run from the repo root:

    python3 code/figures/figures_2026-10.py

Writes figures/2026-10/fig{1..5}_*.png and copies them to docs/assets/figures/.

Palette (validated for colour-vision separation and contrast on white):
accent blue #3B6EA5, teal #2A9D8F, orange #C46A1F.
"""
import os
import shutil

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Patch
from scipy import stats

TAB = "tables/reanalysis_2026-10"
OUT = "figures/2026-10"
SITE = "docs/assets/figures"
os.makedirs(OUT, exist_ok=True)

BLUE, TEAL, ORANGE = "#3B6EA5", "#2A9D8F", "#C46A1F"
INK, MUTED, RULE = "#16212E", "#56657A", "#D9E0E8"
HW_COL = {"BRIC": BLUE, "CARA": TEAL, "VEGGIE": ORANGE, "EMCS": "#7A8794"}  # EMCS: neutral gray, no CFD scenario
HW_ORDER = {"BRIC": 0, "CARA": 1, "VEGGIE": 2, "EMCS": 3}
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
    "xtick.color": MUTED, "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": False, "savefig.dpi": 300, "savefig.bbox": "tight", "savefig.facecolor": "white",
})


def save(fig, name):
    path = os.path.join(OUT, name)
    fig.savefig(path)
    plt.close(fig)
    print("wrote", path)
    return path


meta = pd.read_csv(os.path.join(TAB, "expression_metadata_with_cfd_v2.tsv"), sep="\t")

# ---- Fig 1: study design ------------------------------------------------------
d = (meta.groupby(["osd_id", "hardware", "cfd_hardware", "omics_type"])
     .agg(flt=("flight", lambda x: (x == "FLT").sum()), gc=("flight", lambda x: (x == "GC").sum()))
     .reset_index())
d["order"] = d.cfd_hardware.map(HW_ORDER)
d = d.sort_values(["order", "omics_type", "osd_id"], ascending=[True, False, True]).reset_index(drop=True)
fig, ax = plt.subplots(figsize=(6.4, 3.9))
y = np.arange(len(d))[::-1]
for yi, r in zip(y, d.itertuples()):
    c = HW_COL[r.cfd_hardware]
    ax.barh(yi, r.flt, color=c, height=0.62)
    ax.barh(yi, r.gc, left=r.flt + 0.6, color=c, alpha=0.45, height=0.62)
    plat = "RNA-seq" if r.omics_type == "RNAseq" else "microarray"
    ax.text(r.flt + r.gc + 1.6, yi, f"{r.flt} / {r.gc}   {plat}", va="center", fontsize=8, color=INK)
ax.set_yticks(y)
ax.set_yticklabels([f"{r.osd_id}  ({r.hardware})" for r in d.itertuples()], fontsize=8)
ax.set_xlabel("samples (solid: flight, light: ground control)")
ax.set_xlim(0, d.eval("flt + gc").max() + 26)
ax.legend(handles=[Patch(color=c, label=f"{hw} class") for hw, c in HW_COL.items()], frameon=False, loc="lower right", fontsize=8)
ax.set_title(f"{len(meta)} samples, {meta.osd_id.nunique()} OSDR studies, grouped by hardware class",
             loc="left", fontsize=10, color=INK)
f1 = save(fig, "fig1_study_design.png")

# ---- Fig 2: variance partition ------------------------------------------------
vp = pd.read_csv(os.path.join(TAB, "followups/variance_partition_summary_v2.tsv"), sep="\t")
# gene and sample counts per data set, from the variance-partition run log
import re
VPN = {m.group(1): (int(m.group(2)), int(m.group(3))) for m in
       re.finditer(r"(pooled_complete|rnaseq): (\d+) genes x (\d+) samples", open(os.path.join(TAB, "followups/variance_partition_log.txt")).read())}
comp = [("osd_id", "study (within hardware)"), ("light", "light"), ("ecotype_simple", "ecotype"),
        ("organ_simple", "organ"), ("cfd_hardware", "hardware"), ("flight", "flight, common to all studies"),
        ("flight:cfd_hardware", "flight × hardware"), ("flight:osd_id", "flight × study"), ("Residuals", "residual")]
fig, ax = plt.subplots(figsize=(6.4, 3.6))
y = np.arange(len(comp))[::-1]
for k, (data, col, lab) in enumerate([("pooled_complete", BLUE, "pooled ({:,} genes, {} samples)".format(*VPN["pooled_complete"])),
                                      ("rnaseq", TEAL, "RNA-seq only ({:,} genes, {} samples)".format(*VPN["rnaseq"]))]):
    v = vp[vp.data == data].set_index("component").mean_pct
    vals = [v[c] for c, _ in comp]
    ax.barh(y + (0.19 if k == 0 else -0.19), vals, height=0.36, color=col, label=lab)
    for yi, val in zip(y, vals):
        ax.text(val + 0.8, yi + (0.19 if k == 0 else -0.19), f"{val:.1f}", va="center", fontsize=8, color=INK)
ax.set_yticks(y)
ax.set_yticklabels([l for _, l in comp])
ax.get_yticklabels()[5].set_fontweight("bold")
ax.set_xlabel("mean % of each gene's variance")
ax.set_xlim(0, 62)
ax.legend(frameon=False, fontsize=8, loc="lower right")
ax.set_title("A flight response shared by all studies explains < 1 % of gene variance", loc="left", fontsize=10, color=INK)
f2 = save(fig, "fig2_variance_partition.png")

# ---- Fig 3: platform concordance --------------------------------------------
U = pd.read_csv(os.path.join(TAB, "followups/platform_concordance_units.tsv"), sep="\t")
P = pd.read_csv(os.path.join(TAB, "followups/platform_concordance_pairs.tsv"), sep="\t")
U["order"] = U.hardware.map(HW_ORDER) * 2 + (U.platform == "microarray")
U = U.sort_values(["order", "unit"]).reset_index(drop=True)
names = list(U.unit)
M = pd.DataFrame(np.nan, index=names, columns=names)
for r in P.itertuples():
    M.loc[r.unit_a, r.unit_b] = M.loc[r.unit_b, r.unit_a] = r.r
cmap = LinearSegmentedColormap.from_list("div", [ORANGE, "#EEEEEE", BLUE])
fig, ax = plt.subplots(figsize=(7.8, 7.0))
im = ax.imshow(M.values, cmap=cmap, vmin=-0.8, vmax=0.8)
lab = [f"{u.split('|')[0]} {u.split('|')[1].replace('undifferentiated cell culture', 'cell culture').replace('whole_seedling', 'seedling')}"
       for u in names]
ax.set_xticks(range(len(names)))
ax.set_yticks(range(len(names)))
ax.set_xticklabels(lab, rotation=90, fontsize=8)
ax.set_yticklabels(lab, fontsize=8)
for i, u in enumerate(U.itertuples()):
    c = HW_COL[u.hardware]
    ax.get_yticklabels()[i].set_color(c)
    ax.get_xticklabels()[i].set_color(c)
for i in range(len(names)):
    for j in range(len(names)):
        v = M.values[i, j]
        if np.isfinite(v):
            ax.text(j, i, f"{v:.2f}".replace("0.", ".").replace("-.", "−."), ha="center", va="center", fontsize=8,
                    color="white" if abs(v) > 0.5 else INK)
groups = U.groupby(["hardware", "platform"], sort=False).size()
edge = 0
for (hw, plat), n in groups.items():
    edge += n
    if edge < len(names):
        ax.axhline(edge - 0.5, color="white", lw=2)
        ax.axvline(edge - 0.5, color="white", lw=2)
cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
cb.set_label("Pearson r of flight logFC", fontsize=8)
for spine in ax.spines.values():
    spine.set_visible(False)
ax.set_title("Flight effects agree between BRIC RNA-seq studies, not with the microarray studies\n"
             "(labels coloured by hardware: blue BRIC, teal CARA, orange VEGGIE, gray EMCS; blocks by hardware × platform)",
             loc="left", fontsize=9, color=INK)
f3 = save(fig, "fig3_platform_concordance.png")

# ---- Fig 4: pathway flight effects (rescored, within-study model) ------------
W = pd.read_csv(os.path.join(TAB, "followups/pathway_within_study_model_rescored.tsv"), sep="\t")
W = W[W.term == "flight"].sort_values("estimate").reset_index(drop=True)
tc = stats.t.ppf(0.975, df=W.studies.iloc[0] - 1)
fig, ax = plt.subplots(figsize=(6.0, 5.6))
y = np.arange(len(W))
ax.axvline(0, color=MUTED, lw=0.8)
ax.hlines(y, W.estimate - tc * W.se, W.estimate + tc * W.se, color=BLUE, lw=1.6)
ax.plot(W.estimate, y, "o", color=BLUE, ms=5, mec="white", mew=1)
ax.set_yticks(y)
ax.set_yticklabels(W.pathway.str.replace("_", " "), fontsize=8.5)
ax.set_xlabel("flight effect (FLT − GC), within-study SD units, 95 % CI")
ax.set_title(f"No pathway responds to flight after accounting for study (0/{len(W)} at FDR < 0.05)",
             loc="left", fontsize=10, color=INK)
ax.text(1.0, -0.11, f"score ~ study + flight × scale; cluster-robust SEs, {W.studies.iloc[0]} studies; unadjusted 95 % CIs; smallest FDR {W.fdr.min():.2f}",
        transform=ax.transAxes, ha="right", va="top", fontsize=8, color=MUTED)
f4 = save(fig, "fig4_pathway_flight_effects.png")

# ---- Fig 5: CFD gravity penalty, open box vs real dish ------------------------
leaf = meta[meta.cfd_scale == "leaf"]
def cell(hw, gk, col):
    v = leaf[(leaf.cfd_hardware == hw) & (leaf.gravity_key == gk)][col].dropna().unique()
    assert len(v) == 1, (hw, gk, col, v)
    return float(v[0])
rows = [("open box (leaf)", cell("BRIC", "earth", "g_bl_mol_m2_s"), cell("BRIC", "ug", "g_bl_mol_m2_s"), MUTED, ""),
        ("BRIC 60 mm dish, sealed", cell("BRIC", "earth", "g_bl_dish_mol_m2_s"), cell("BRIC", "ug", "g_bl_dish_mol_m2_s"), BLUE, "*"),
        ("CARA 120 mm dish, taped", cell("CARA", "earth", "g_bl_dish_mol_m2_s"), cell("CARA", "ug", "g_bl_dish_mol_m2_s"), TEAL, "")]
fig, ax = plt.subplots(figsize=(6.2, 2.6))
for i, (lab, g1, g0, c, star) in enumerate(rows[::-1]):
    ax.plot([g0, g1], [i, i], color=c, lw=2.2)
    ax.plot(g0, i, "o", color="white", mec=c, mew=2, ms=7)
    ax.plot(g1, i, "o", color=c, ms=7)
    ax.text(g1 + 0.02, i, f"1 g {g1:.3g}{star}   ×{g1 / g0:.2f}", va="center", fontsize=8, color=INK)
    ax.text(g0 - 0.02, i, f"0 g {g0:.3g}", va="center", ha="right", fontsize=8, color=INK)
ax.set_yticks(range(len(rows)))
ax.set_yticklabels([r[0] for r in rows[::-1]])
ax.set_xlim(0, 1.25)
ax.set_xlabel("leaf boundary-layer conductance g$_{bl}$ (mol m$^{-2}$ s$^{-1}$)")
ax.set_title("Inside a flight dish, gravity changes leaf gas exchange far less than in open air",
             loc="left", fontsize=10, color=INK)
ax.text(0.99, -0.30, "LunarLeaf-CFD T13 (open box) and T19 (dish, to the dish's own air). *Sealed 1 g value indicative only.",
        transform=ax.transAxes, ha="right", fontsize=8, color=MUTED)
f5 = save(fig, "fig5_cfd_dish_penalty.png")

os.makedirs(SITE, exist_ok=True)
for f in (f1, f2, f3, f4, f5):
    shutil.copy(f, os.path.join(SITE, os.path.basename(f)))
print("copied to", SITE)
