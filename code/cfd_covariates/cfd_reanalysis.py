"""
Re-derive the CFD covariates from the corrected LunarLeaf-CFD tables and redo the
pathway-level CFD statistics (2026-10-07).

Why: LunarLeaf-CFD regenerated its tables after two defects were found (0 g scenes
reported before steady state; a fan outlet that let the through-flow decay). The
covariates in tables/expression_metadata_with_cfd.tsv were built from the old tables.
Separately, the original pathway-vs-CFD test (hormone_carbon_model.py) ran Pearson
correlations on n = 307 samples, but g_bl takes only four values, one per
gravity x scale cell, and samples cluster by study.

What runs here (everything local to this package; no expression matrix needed):
  1. Rebuild the per-scenario covariates with extract_cfd.py's own code, once from
     the old LunarLeaf tables and once from the new ones.
  2. VALIDATE: the old tables must reproduce the stored per-sample covariates.
  3. Map the new covariates onto the 307 samples.
  4. VALIDATE: the original Pearson method on the stored covariates must reproduce
     tables/pathway_cfd_correlations.tsv. Then rerun it on the new covariates.
  5. Identifiability: g_bl vs (gravity, scale). Then a within-study model per
     pathway: score ~ study + flight * scale (sum-to-zero coding), with
     cluster-robust (CR1) standard errors by study and BH-FDR.
  6. From the stored DE tables: is 'Flight_FLT_vs_GC' identical to the BRIC-only
     contrast, and how does the |logFC| > 1 call for g_bl depend on its units?

The gene-level limma models (fit_factorial_model.R) are NOT refit: they need
harmonized_expression_matrix_filtered.tsv, which only exists in the original
Biomni workspace.

    LUNARLEAF_TABLES=~/Documents/LunarLeaf-CFD/results/tables \\
    LUNARLEAF_OLD_TABLES=<export of LunarLeaf-CFD@b3567f6 results/tables> \\
    python3 code/cfd_covariates/cfd_reanalysis.py
"""
import os
import re
import sys

import warnings

import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.dirname(os.path.dirname(HERE))  # code/cfd_covariates/ -> repo root
TAB = os.path.join(PKG, "tables")
# REANALYSIS_OUT / PATHWAY_SCORES let the same analysis run on rescored pathway scores
# (followups_2026-10.R) without touching the primary outputs.
OUT = os.path.expanduser(os.environ.get("REANALYSIS_OUT", os.path.join(TAB, "reanalysis_2026-10")))
SCORES = os.path.expanduser(os.environ.get("PATHWAY_SCORES", os.path.join(TAB, "pathway_scores_gsva.tsv")))
os.makedirs(OUT, exist_ok=True)
# numpy 2.0 + macOS Accelerate emits spurious matmul overflow warnings on well-conditioned
# products; coefficients were checked against np.linalg.lstsq (max |diff| ~ 4e-15).
warnings.filterwarnings("ignore", message=".*encountered in matmul")

NEW_TABLES = os.path.expanduser(os.environ.get("LUNARLEAF_TABLES", "~/Documents/LunarLeaf-CFD/results/tables"))
OLD_TABLES = os.path.expanduser(os.environ.get("LUNARLEAF_OLD_TABLES", ""))

KEYS = ["cfd_hardware", "gravity_key", "cfd_scale"]


# ---------------------------------------------------------------- 1. scenarios
def build_scenarios(tables_dir: str) -> pd.DataFrame:
    """Run extract_cfd.py sections 1-3 (scenario table only) against `tables_dir`."""
    src = open(os.path.join(HERE, "extract_cfd.py"), encoding="utf-8").read()
    start = src.index("import os")
    end = src.index("# 4. Map CFD covariates to each OSDR sample")
    block = src[start:end]
    # Keep the original code; only redirect its two output side effects.
    block = block.replace("os.makedirs(OUT, exist_ok=True)", "")
    block = block.replace('scenario_df.to_csv(f"{OUT}/cfd_covariates_by_scenario.tsv", sep="\\t", index=False)', "")
    env_before = os.environ.get("LUNARLEAF_TABLES")
    os.environ["LUNARLEAF_TABLES"] = tables_dir
    ns: dict = {"__name__": "extract_cfd_scenarios"}
    stdout = sys.stdout
    try:
        sys.stdout = open(os.devnull, "w")
        exec(compile(block, "extract_cfd.py[1-3]", "exec"), ns)
    finally:
        sys.stdout.close()
        sys.stdout = stdout
        if env_before is None:
            os.environ.pop("LUNARLEAF_TABLES", None)
        else:
            os.environ["LUNARLEAF_TABLES"] = env_before
    sc = ns["scenario_df"].rename(columns={"scale": "cfd_scale"})
    return sc


def apply(meta: pd.DataFrame, sc: pd.DataFrame) -> pd.DataFrame:
    """Replace the stored covariates with the scenario table's, and append any covariate
    columns the scenario table has that the stored metadata does not (e.g. g_bl_dish)."""
    cols = [c for c in sc.columns if c not in KEYS and c in meta.columns and c != "gravity_g"]
    extra = [c for c in sc.columns if c not in KEYS and c not in meta.columns and c != "gravity_g"]
    out = meta.drop(columns=cols).merge(sc[KEYS + cols + extra], on=KEYS, how="left", validate="many_to_one")
    return out[list(meta.columns) + extra]


def compare(a: pd.DataFrame, b: pd.DataFrame, cols) -> pd.DataFrame:
    rows = []
    for c in cols:
        x, y = a[c], b[c]
        if pd.api.types.is_numeric_dtype(x) and pd.api.types.is_numeric_dtype(y):
            same = np.isclose(x.astype(float), y.astype(float), rtol=1e-6, atol=1e-9, equal_nan=True)
        else:
            same = (x.astype(str) == y.astype(str)).to_numpy()
        rows.append({"column": c, "n_mismatch": int((~same).sum())})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- stats helpers
def md(df: pd.DataFrame, index: bool = False) -> str:
    """Minimal markdown table (avoids the optional `tabulate` dependency)."""
    d = df.reset_index() if index else df
    cols = [str(c) for c in d.columns]
    fmt = lambda v: f"{v:.4g}" if isinstance(v, float) else str(v)  # noqa: E731
    lines = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    lines += ["| " + " | ".join(fmt(v) for v in row) + " |" for row in d.itertuples(index=False)]
    return "\n".join(lines)


def bh(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(q, 0, 1)
    return out


def pearson_table(df: pd.DataFrame, pathways, covs) -> pd.DataFrame:
    """The original method (hormone_carbon_model.py): Pearson per pathway x covariate, BH."""
    rows = []
    for pw in pathways:
        for cv in covs:
            v = df[[pw, cv]].dropna()
            if len(v) > 10:
                r, p = stats.pearsonr(v[pw], v[cv])
                rows.append({"pathway": pw, "cfd_covariate": cv, "correlation": r, "p_value": p, "n": len(v)})
    t = pd.DataFrame(rows)
    t["fdr"] = bh(t["p_value"].to_numpy())
    return t


def ols_cluster(y: np.ndarray, X: np.ndarray, groups: np.ndarray):
    """OLS with CR1 cluster-robust covariance; t-tests on G - 1 df."""
    n, k = X.shape
    XtX_inv = np.linalg.pinv(X.T @ X)
    beta = XtX_inv @ X.T @ y
    resid = y - X @ beta
    meat = np.zeros((k, k))
    gs = np.unique(groups)
    for g in gs:
        idx = groups == g
        s = X[idx].T @ resid[idx]
        meat += np.outer(s, s)
    G = len(gs)
    c = (G / (G - 1)) * ((n - 1) / (n - k))
    V = c * XtX_inv @ meat @ XtX_inv
    se = np.sqrt(np.diag(V))
    t = beta / se
    p = 2 * stats.t.sf(np.abs(t), df=G - 1)
    return beta, se, t, p


def within_study_design(df: pd.DataFrame):
    """score ~ study + flight * scale, sum-to-zero (effect) coding for flight and scale,
    so 'flight' is the main effect averaged over scales, not the effect at a reference."""
    f = np.where(df["flight"] == "FLT", 0.5, -0.5)          # FLT - GC
    s = np.where(df["cfd_scale"] == "rosette", 0.5, -0.5)   # rosette - leaf
    study = pd.get_dummies(df["osd_id"], drop_first=True, dtype=float)
    X = np.column_stack([np.ones(len(df)), f, s, f * s, study.to_numpy()])
    names = ["intercept", "flight", "scale", "flight_x_scale"] + list(study.columns)
    return X, names


# ---------------------------------------------------------------- main
def main():
    report = []
    log = lambda s="": (print(s), report.append(s))  # noqa: E731

    meta = pd.read_csv(os.path.join(TAB, "expression_metadata_with_cfd.tsv"), sep="\t")
    scores = pd.read_csv(SCORES, sep="\t", index_col=0)
    pathways = list(scores.columns)
    log(f"# CFD covariate update and pathway re-analysis\n\n{len(meta)} samples, {meta.osd_id.nunique()} studies, {len(pathways)} pathways.\n")

    sc_new = build_scenarios(NEW_TABLES)
    cfd_cols = [c for c in sc_new.columns if c not in KEYS and c in meta.columns and c != "gravity_g"]

    # 2. validate mapping against the stored covariates
    log("## 1. Mapping check (old LunarLeaf tables -> stored covariates)\n")
    if OLD_TABLES and os.path.isdir(OLD_TABLES):
        sc_old = build_scenarios(OLD_TABLES)
        rebuilt_old = apply(meta, sc_old)
        cmp = compare(meta, rebuilt_old, cfd_cols)
        bad = cmp[cmp.n_mismatch > 0]
        log(f"Rebuilt the stored covariates from LunarLeaf `{os.path.basename(OLD_TABLES)}`: "
            f"{len(cfd_cols) - len(bad)}/{len(cfd_cols)} columns reproduce exactly.")
        if len(bad):
            log("Columns that do not reproduce:\n")
            log(md(bad))
            by_design = {
                "earth_equiv_fan_cm_s": "read from T9 instead of hard-coded 2.8 / 11.0 cm/s (2026-10-07)",
                "fan_cm_s": "VEGGIE fan is the modelled 8.3 cm/s, not the Earth-equivalent 2.8 (2026-10-07)",
                "carbon_12h_pct_earth": "BRIC ground controls use the same sealed hardware, so they get T11's BRIC value, not 100 % (2026-10-08). The CFD class BRIC includes APEX-TAGES (OSD-7), whose 18 ground controls get the same value; TODO: confirm they were sealed",
                "A_start_umol_m2_s": "as carbon_12h_pct_earth: BRIC ground controls now carry T11's BRIC values (2026-10-08)",
                "A_end_umol_m2_s": "as carbon_12h_pct_earth: BRIC ground controls now carry T11's BRIC values (2026-10-08)",
            }
            if set(bad.column) <= set(by_design):
                log("\nAll expected, from deliberate changes to extract_cfd.py:\n")
                for col in bad.column:
                    log(f"- `{col}`: {by_design[col]}")
                log("\nEvery other column reproduces, which validates the mapping.")
            else:
                log("\n**Unexpected mismatches:** " + ", ".join(sorted(set(bad.column) - set(by_design))))
        log("")
    else:
        log("LUNARLEAF_OLD_TABLES not set: mapping check skipped.\n")

    # 3. apply new covariates
    meta_new = apply(meta, sc_new)
    meta_new.to_csv(os.path.join(OUT, "expression_metadata_with_cfd_v2.tsv"), sep="\t", index=False)
    changed = compare(meta, meta_new, cfd_cols)
    changed = changed[changed.n_mismatch > 0].sort_values("n_mismatch", ascending=False)
    changed.to_csv(os.path.join(OUT, "covariate_changes.tsv"), sep="\t", index=False)
    log("## 2. Covariates that change with the corrected LunarLeaf tables\n")
    log(md(changed))
    g = meta_new.groupby(["gravity_key", "cfd_scale"])["g_bl_mol_m2_s"].first()
    g_old = meta.groupby(["gravity_key", "cfd_scale"])["g_bl_mol_m2_s"].first()
    log("\ng_bl per cell, old -> new: " + "; ".join(f"{k[0]}/{k[1]} {g_old[k]:.3f} -> {g[k]:.3f}" for k in g.index))
    log("\nNot updated because their LunarLeaf source tables were not regenerated: "
        "co2_drawdown_* / o2_buildup_* (T5), net_assimilation_pct (T12).\n")
    hwmap = (meta.groupby(["hardware", "cfd_hardware"])
             .agg(studies=("osd_id", lambda x: ", ".join(sorted(x.unique()))), n=("osd_id", "size")).reset_index())
    log("OSDR hardware -> CFD class (`extract_cfd.py`). APEX-TAGES is folded into BRIC and EMCS into VEGGIE, so "
        "every CFD covariate for them, including the ground-control carbon value, is the class value:\n")
    log(md(hwmap))
    log("\nTODO: confirm that the APEX-TAGES (OSD-7) ground controls were sealed like BRIC's, and whether EMCS (OSD-314) "
        "is adequately represented by the vented VEGGIE scenario.\n")

    # 4. original Pearson method: validate, then rerun on new covariates
    covs = ["g_bl_mol_m2_s", "carbon_12h_pct_earth", "delta_mm"]
    s = scores.copy()
    s["full_sample_id"] = s.index.astype(str)
    old_df = s.merge(meta, on="full_sample_id", how="left")
    new_df = s.merge(meta_new, on="full_sample_id", how="left")
    p_old = pearson_table(old_df, pathways, covs)
    stored = pd.read_csv(os.path.join(TAB, "pathway_cfd_correlations.tsv"), sep="\t")
    chk = stored.merge(p_old, on=["pathway", "cfd_covariate"], suffixes=("_stored", "_rebuilt"))
    max_dr = float(np.max(np.abs(chk.correlation_stored - chk.correlation_rebuilt)))
    log("## 3. Original method (Pearson, n = 307, BH)\n")
    if max_dr < 1e-6:
        log(f"Reproduces the stored `pathway_cfd_correlations.tsv`: {len(chk)}/{len(stored)} rows matched, max |Δr| = {max_dr:.2e}.\n")
    else:  # expected when PATHWAY_SCORES points at rescored scores
        log(f"Differs from the stored `pathway_cfd_correlations.tsv` (these are not the stored scores): {len(chk)}/{len(stored)} "
            f"pathway x covariate rows compared, max |Δr| = {max_dr:.2e}.\n")
    p_new = pearson_table(new_df, pathways, covs)
    p_new.to_csv(os.path.join(OUT, "pathway_cfd_correlations_v2_pearson.tsv"), sep="\t", index=False)
    for cv in covs:
        a, b = p_old[p_old.cfd_covariate == cv], p_new[p_new.cfd_covariate == cv]
        log(f"- `{cv}`: FDR < 0.05 in {int((a.fdr < 0.05).sum())}/{len(a)} pathways with old covariates, "
            f"{int((b.fdr < 0.05).sum())}/{len(b)} with new; max |r| {a.correlation.abs().max():.3f} -> {b.correlation.abs().max():.3f}.")
    cv = meta_new.groupby(["cfd_hardware", "flight"])["carbon_12h_pct_earth"].first().unstack()
    log("\n`carbon_12h_pct_earth` by hardware and flight (new covariates): "
        + "; ".join(f"{hw} FLT {cv.loc[hw, 'FLT']:g} / GC {cv.loc[hw, 'GC']:g}" for hw in cv.index)
        + ". Within BRIC it no longer differs between flight and ground, so its correlation with pathway scores "
        "now reflects which studies used BRIC (a between-study, hardware difference), not a flight effect.")
    log("")

    # 5. identifiability + within-study model
    log("## 4. Is g_bl identifiable separately from flight and scale?\n")
    ct = pd.crosstab([meta_new.gravity_key, meta_new.cfd_scale], meta_new.g_bl_mol_m2_s)
    log(md(ct, index=True))
    X, names = within_study_design(new_df)
    Xg = np.column_stack([X, new_df["g_bl_mol_m2_s"].to_numpy()])
    log(f"\nDesign rank: study + flight × scale = {np.linalg.matrix_rank(X)} columns; adding g_bl gives rank "
        f"{np.linalg.matrix_rank(Xg)} of {Xg.shape[1]} columns. g_bl is a deterministic relabelling of the flight × "
        "scale cells, so no model can attribute an effect to g_bl rather than to flight and organ scale in this design. "
        "Its correlations with pathways are flight × scale contrasts under another name.\n")

    log("## 5. Within-study model per pathway\n")
    log("`score ~ study + flight × scale` (effect coding), CR1 cluster-robust SEs by study "
        f"({new_df.osd_id.nunique()} clusters, t on G − 1 df), BH across pathways per term.\n")
    rows = []
    for pw in pathways:
        y = new_df[pw].to_numpy(float)
        beta, se, t, p = ols_cluster(y, X, new_df["osd_id"].to_numpy())
        for term in ["flight", "scale", "flight_x_scale"]:
            i = names.index(term)
            rows.append({"pathway": pw, "term": term, "estimate": beta[i], "se": se[i], "t": t[i], "p_value": p[i]})
    ws = pd.DataFrame(rows)
    ws["fdr"] = np.nan
    for term in ws.term.unique():
        m = ws.term == term
        ws.loc[m, "fdr"] = bh(ws.loc[m, "p_value"].to_numpy())
    ws.to_csv(os.path.join(OUT, "pathway_within_study_model.tsv"), sep="\t", index=False)
    for term in ["flight", "scale", "flight_x_scale"]:
        sub = ws[ws.term == term]
        log(f"- **{term}**: FDR < 0.05 in {int((sub.fdr < 0.05).sum())}/{len(sub)} pathways.")
    # Checks on this model: clustering vs naive SEs, and a positive control (L008:
    # a null result is only informative if the test can detect a real effect).
    naive_sig, planted = 0, {0.25: 0, 0.5: 0}
    sd_med = float(new_df[pathways].std().median())
    flt = (new_df["flight"] == "FLT").to_numpy(float)
    groups = new_df["osd_id"].to_numpy()
    n_obs, k_par = X.shape
    for pw in pathways:
        y = new_df[pw].to_numpy(float)
        b_ = np.linalg.lstsq(X, y, rcond=None)[0]
        r_ = y - X @ b_
        se_ = np.sqrt(np.diag((r_ @ r_ / (n_obs - k_par)) * np.linalg.pinv(X.T @ X)))
        naive_sig += 2 * stats.t.sf(abs(b_[1] / se_[1]), n_obs - k_par) < 0.05 / len(pathways)
        for eff in planted:
            planted[eff] += ols_cluster(y + eff * sd_med * flt, X, groups)[3][1] < 0.05
    log(f"\nFor comparison, the same model with naive (independent-sample) SEs calls the flight term significant "
        f"in {naive_sig}/{len(pathways)} pathways at Bonferroni 0.05: the original significance comes from treating "
        "307 samples as independent when they cluster in 11 studies.")
    log("Positive control: a flight effect planted in every pathway is detected (p < 0.05) in "
        + ", ".join(f"{v}/{len(pathways)} at {k} SD" for k, v in planted.items())
        + (", so the null result is not a powerless test; the pathway-level flight effects are small or inconsistent "
           "across studies relative to between-study variation."
           if min(planted.values()) >= 0.9 * len(pathways) else
           ", so the test has limited power: flight effects below about half an SD cannot be excluded, and the null "
           "means no large, consistent flight effect rather than no effect."))
    g_sig = set(p_new[(p_new.cfd_covariate == "g_bl_mol_m2_s") & (p_new.fdr < 0.05)].pathway)
    any_ws = set(ws[(ws.fdr < 0.05)].pathway)
    log(f"\nOf the {len(g_sig)} pathways 'significant' for g_bl by Pearson, {len(g_sig & any_ws)} have any FDR < 0.05 "
        "term (flight, scale or interaction) once within-study structure and clustering are respected.\n")

    # 7 (reported after 6). In-dish covariate from LunarLeaf T19
    dish_section = []
    dl = dish_section.append
    if "g_bl_dish_mol_m2_s" in meta_new.columns and meta_new["g_bl_dish_mol_m2_s"].notna().any():
        sub = new_df[new_df["g_bl_dish_mol_m2_s"].notna()].copy()
        dl("## 7. In-dish covariate (LunarLeaf T19)\n")
        dl(f"`g_bl_dish` (leaf conductance to the dish's own air, in the real BRIC 60 mm / CARA 120 mm dish) is defined for "
           f"{len(sub)} leaf samples in BRIC or CARA hardware from {sub.osd_id.nunique()} studies. It is NaN for rosette-scale "
           "and VEGGIE samples, so it is kept as a separate column rather than replacing `g_bl`.\n")
        ct = sub.groupby(["cfd_hardware", "gravity_key"]).agg(n=("osd_id", "size"), studies=("osd_id", "nunique"),
                                                            g_bl_dish=("g_bl_dish_mol_m2_s", "first"),
                                                            g_bl_open_box=("g_bl_mol_m2_s", "first"),
                                                            note=("g_bl_dish_note", "first"))
        dl(md(ct, index=True))
        hw_by_study = sub.groupby("osd_id")["cfd_hardware"].first()
        cara = sorted(hw_by_study[hw_by_study == "CARA"].index)
        dl(f"\nAgain one value per hardware × gravity cell, so it is identifiable only as hardware × flight. CARA comes from "
           f"{', '.join(cara)} alone, so a CARA-vs-BRIC difference in the flight response is also a difference between studies.\n")
        f = np.where(sub["flight"] == "FLT", 0.5, -0.5)
        c = (sub["cfd_hardware"] == "CARA").to_numpy(float)
        study = pd.get_dummies(sub["osd_id"], drop_first=True, dtype=float)
        Xs = np.column_stack([np.ones(len(sub)), f, f * c, study.to_numpy()])
        G = sub.osd_id.nunique()
        rows_ = []
        for pw in pathways:
            b_, se_, t_, p_ = ols_cluster(sub[pw].to_numpy(float), Xs, sub["osd_id"].to_numpy())
            rows_.append({"pathway": pw, "flight_BRIC": b_[1], "p_flight": p_[1], "flight_x_CARA": b_[2], "p_x_CARA": p_[2]})
        ds = pd.DataFrame(rows_)
        ds["fdr_flight"] = bh(ds.p_flight.to_numpy())
        # Cluster-robust SEs are invalid for a term that varies within a single cluster:
        # its variance collapses and p-values become spuriously small. CARA is one study.
        n_cara_clusters = sub.loc[sub.cfd_hardware == "CARA", "osd_id"].nunique()
        if n_cara_clusters < 2:
            ds["p_x_CARA"] = np.nan
            ds["fdr_x_CARA"] = np.nan
        else:
            ds["fdr_x_CARA"] = bh(ds.p_x_CARA.to_numpy())
        ds.to_csv(os.path.join(OUT, "pathway_dish_subset_model.tsv"), sep="\t", index=False)
        pen = ct["g_bl_dish"].unstack()
        open_pen = sub.loc[sub.gravity_key == "earth", "g_bl_mol_m2_s"].iloc[0] - sub.loc[sub.gravity_key == "ug", "g_bl_mol_m2_s"].iloc[0]
        dl(f"Predicted in-dish penalty (1 g − 0 g of g_bl_dish): BRIC {pen.loc['BRIC', 'earth'] - pen.loc['BRIC', 'ug']:.3f} "
           f"(its 1 g value is indicative only), CARA {pen.loc['CARA', 'earth'] - pen.loc['CARA', 'ug']:.3f}. In the open box it was "
           f"{open_pen:.3f} for every hardware.\n")
        x_txt = ("**not testable**: CARA is a single study, and cluster-robust errors for a term that varies within one "
                 "cluster collapse toward zero (an unguarded run reports it 'significant' in most pathways, which is an artefact)"
                 if n_cara_clusters < 2 else
                 f"FDR < 0.05 in {int((ds.fdr_x_CARA < 0.05).sum())}/{len(ds)}")
        dl(f"Within-study model on this subset (`score ~ study + flight + flight × CARA`, CR1 by study, {G} clusters, t on {G - 1} df): "
           f"flight effect FDR < 0.05 in {int((ds.fdr_flight < 0.05).sum())}/{len(ds)} pathways. CARA-vs-BRIC difference in the "
           f"flight response: {x_txt}.\n")
        dl("**Reading:** the dish geometry changes the *size* of the predicted gravity penalty (smaller than in the open box, and "
           "different for BRIC and CARA), but it adds no identifiable covariate. In this design it is still a hardware × flight "
           "label. Testing it as a dose-response needs study-specific dish scenarios across more studies per hardware.")

    # 6. stored DE tables: contrast label and logFC units
    log("## 6. Gene-level checks on the stored DE tables\n")
    try:
        a = pd.read_csv(os.path.join(TAB, "de_results_Flight_FLT_vs_GC.tsv"), sep="\t")
        b = pd.read_csv(os.path.join(TAB, "de_results_Flight_BRIC_FLT_vs_GC.tsv"), sep="\t")
        gene = [c for c in a.columns if c.lower() in ("gene", "gene_id", "geneid", "symbol", "id")]
        key = gene[0] if gene else a.columns[0]
        mm = a.merge(b, on=key, suffixes=("_main", "_bric"))
        same = np.allclose(mm["logFC_main"], mm["logFC_bric"], equal_nan=True)
        log(f"- `Flight_FLT_vs_GC` vs `Flight_BRIC_FLT_vs_GC`: {len(mm)} genes compared; logFC identical: **{same}** "
            f"(max |Δ| = {np.nanmax(np.abs(mm.logFC_main - mm.logFC_bric)):.2e}). With BRIC as the reference level, the "
            "'Flight' coefficient is the flight effect within BRIC, not a main effect across hardware.")
    except Exception as e:  # noqa: BLE001
        log(f"- Contrast check skipped: {e}")
    try:
        d = pd.read_csv(os.path.join(TAB, "de_results_CFD_g_bl.tsv"), sep="\t")
        sig = d["adj.P.Val"] < 0.05
        sd_old = meta["g_bl_mol_m2_s"].std()
        span_old = meta["g_bl_mol_m2_s"].max() - meta["g_bl_mol_m2_s"].min()
        log(f"- `CFD_g_bl` (logFC per 1 mol m⁻² s⁻¹; the covariate spans only {span_old:.3f}, SD {sd_old:.3f}): "
            f"|logFC| > 1 & FDR < 0.05 → {int((sig & (d.logFC.abs() > 1)).sum())} genes. "
            f"The same genes per 1 SD of g_bl: |logFC·SD| > 1 → {int((sig & ((d.logFC * sd_old).abs() > 1)).sum())}; "
            f"per full observed range: {int((sig & ((d.logFC * span_old).abs() > 1)).sum())}. The cutoff counts depend on the unit chosen, "
            "so for a continuous covariate report FDR plus an effect size per SD or per observed range.")
    except Exception as e:  # noqa: BLE001
        log(f"- logFC-unit check skipped: {e}")

    for line in dish_section:
        log(line)
    log("\n## Not done here\n\nThe gene-level refit is §9 (`code/factorial_model/refit_gene_level_2026-10.R`), its follow-ups §10 (`followups_2026-10.R`), the rescored §5/§7/§8 comparison §11 (`code/cfd_covariates/rescored_comparison.py`), the variance partition §12 (`variance_partition_2026-10.R`) and the platform comparison §13 (`platform_concordance_2026-10.R`); run them in that order after `cfd_meta_regression.py`. No analysis items from this reanalysis remain open. Open for the authors: the headline claim and the manuscript, README and site text that follow from it; the RNA-seq vs microarray result (§13) is to be reported as a finding.")
    # keep sections written later by other scripts (§9 onward: the R refits and §11)
    rpath = os.path.join(OUT, "REANALYSIS.md")
    text = "\n".join(report) + "\n"
    if os.path.exists(rpath):
        later = re.search(r"\n## (?:9|1[0-9])\. .*?(?=\n## Not done here|\Z)", open(rpath).read(), flags=re.S)
        if later and "\n## Not done here" in text:
            text = text.replace("\n## Not done here", later.group(0) + "\n## Not done here", 1)
    open(rpath, "w").write(text)
    print(f"\nwrote {OUT}/")


if __name__ == "__main__":
    main()
