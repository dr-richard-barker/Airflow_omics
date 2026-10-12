"""
Does the CFD-predicted boundary-layer penalty explain WHICH studies respond most to flight?

g_bl and every other per-scenario CFD covariate is constant within a hardware x gravity x
scale cell, so at the sample level it cannot be separated from flight x organ
(cfd_reanalysis.py, section 4). What IS identifiable is between-study heterogeneity of
the within-study flight effect. Each study x organ-scale unit has a flight effect
(FLT - GC, within the study) and a physics-predicted penalty for its own hardware and
scale. A random-effects meta-regression asks whether units with a larger predicted
penalty show larger flight responses.

  unit      = study x cfd_scale, with >= 2 FLT and >= 2 GC samples
  effect    = mean(FLT) - mean(GC) of each pathway score, Welch SE
  predictor = |surface dC_CO2| in microgravity for the unit's hardware and scale
              minus the Earth value for that scale (surf_dC_CO2_hw, from LunarLeaf T10/T2)
  model     = random-effects meta-regression, tau2 by DerSimonian-Laird, Knapp-Hartung
              (t, k - 2 df) inference because k is small; BH across pathways.
  robustness= the same with a VEGGIE indicator instead of the penalty; BRIC-only units
              (where the penalty still varies by scale); leave-one-study-out range.

Caveats, stated in the output: few units; the predictor varies only by hardware x
scale; T10's BRIC/CARA values come from the original 3.8 s runs (a sealed dish has
no steady state), while VEGGIE is at steady state.

    python3 code/cfd_covariates/cfd_meta_regression.py     # needs tables/reanalysis_2026-10/ from cfd_reanalysis.py
"""
import os
import re

import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
TAB = os.path.join(os.path.dirname(os.path.dirname(HERE)), "tables")  # code/cfd_covariates/ -> repo root
# REANALYSIS_OUT / PATHWAY_SCORES let the same analysis run on rescored pathway scores
# (followups_2026-10.R) without touching the primary outputs.
OUT = os.path.expanduser(os.environ.get("REANALYSIS_OUT", os.path.join(TAB, "reanalysis_2026-10")))
SCORES = os.path.expanduser(os.environ.get("PATHWAY_SCORES", os.path.join(TAB, "pathway_scores_gsva.tsv")))


def bh(p):
    p = np.asarray(p, float)
    n = len(p)
    o = np.argsort(p)
    q = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    out = np.empty(n)
    out[o] = np.clip(q, 0, 1)
    return out


def dl_meta_regression(y, v, x):
    """Random-effects meta-regression y_i = b0 + b1 x_i + u_i + e_i, var(e_i) = v_i,
    u_i ~ N(0, tau2), tau2 by the DerSimonian-Laird moment estimator.
    Returns b1, se(b1), p (normal), tau2, Q_E (residual heterogeneity)."""
    X = np.column_stack([np.ones_like(x), x])
    W = np.diag(1 / v)
    XtWX_inv = np.linalg.inv(X.T @ W @ X)
    b_fe = XtWX_inv @ X.T @ W @ y
    QE = float((y - X @ b_fe) @ W @ (y - X @ b_fe))
    k, p = X.shape
    P = W - W @ X @ XtWX_inv @ X.T @ W
    tau2 = max(0.0, (QE - (k - p)) / np.trace(P))
    Ws = np.diag(1 / (v + tau2))
    Vb = np.linalg.inv(X.T @ Ws @ X)
    b = Vb @ X.T @ Ws @ y
    se = float(np.sqrt(Vb[1, 1]))
    z = b[1] / se
    return float(b[1]), se, float(2 * stats.norm.sf(abs(z))), tau2, QE


def kh_meta_regression(y, v, x):
    """Knapp-Hartung inference on the DL fit: SE scaled by the weighted residual mean
    square (floored at 1), t on k - 2 df. Returns slope, KH p, tau2."""
    _, _, _, tau2, _ = dl_meta_regression(y, v, x)
    X = np.column_stack([np.ones_like(x), x])
    w = 1 / (v + tau2)
    Vb = np.linalg.inv(X.T @ (w[:, None] * X))
    b = Vb @ X.T @ (w * y)
    k = len(y)
    q = float(np.sum(w * (y - X @ b) ** 2) / (k - 2))
    se = np.sqrt(max(q, 1.0) * Vb[1, 1])
    return float(b[1]), float(2 * stats.t.sf(abs(b[1] / se), k - 2)), tau2


def main():
    meta = pd.read_csv(os.path.join(OUT, "expression_metadata_with_cfd_v2.tsv"), sep="\t")
    scores = pd.read_csv(SCORES, sep="\t", index_col=0)
    pathways = list(scores.columns)
    scores["full_sample_id"] = scores.index.astype(str)
    df = scores.merge(meta, on="full_sample_id", how="left")

    units = []
    for (study, scale), g in df.groupby(["osd_id", "cfd_scale"]):
        f, c = g[g.flight == "FLT"], g[g.flight == "GC"]
        if len(f) < 2 or len(c) < 2:
            continue
        hw = f.cfd_hardware.mode().iat[0]
        pen = abs(f.surf_dC_CO2_hw.mean()) - abs(c.surf_dC_CO2_hw.mean())
        d_gbl = f.g_bl_mol_m2_s.mean() - c.g_bl_mol_m2_s.mean()
        row = {"unit": f"{study}|{scale}", "study": study, "scale": scale, "hardware": hw,
               "n_flt": len(f), "n_gc": len(c), "penalty_dCO2": pen, "delta_g_bl": d_gbl}
        for pw in pathways:
            row[f"{pw}__eff"] = f[pw].mean() - c[pw].mean()
            row[f"{pw}__var"] = f[pw].var(ddof=1) / len(f) + c[pw].var(ddof=1) / len(c)
        units.append(row)
    U = pd.DataFrame(units)
    U[["unit", "hardware", "n_flt", "n_gc", "penalty_dCO2", "delta_g_bl"]].to_csv(
        os.path.join(OUT, "meta_regression_units.tsv"), sep="\t", index=False)

    rows = []
    x = U.penalty_dCO2.to_numpy(float)
    veg = (U.hardware == "VEGGIE").to_numpy(float)
    bric = (U.hardware == "BRIC").to_numpy()
    for pw in pathways:
        y = U[f"{pw}__eff"].to_numpy(float)
        v = np.maximum(U[f"{pw}__var"].to_numpy(float), 1e-12)
        b1, p_kh, tau2 = kh_meta_regression(y, v, x)
        p_dl = dl_meta_regression(y, v, x)[2]
        p_veg = kh_meta_regression(y, v, veg)[1]
        p_bric = kh_meta_regression(y[bric], v[bric], x[bric])[1] if bric.sum() >= 4 else np.nan
        loso = [kh_meta_regression(y[U.study != st], v[U.study != st], x[U.study != st])[1] for st in U.study.unique()]
        rows.append({"pathway": pw, "slope_per_unit_penalty": b1, "p_value": p_kh, "p_dl_z": p_dl, "tau2": tau2,
                     "p_veggie_indicator": p_veg, "p_bric_only_units": p_bric, "p_loso_max": max(loso), "k_units": len(y)})
    R = pd.DataFrame(rows)
    R["fdr"] = bh(R.p_value)
    R.sort_values("p_value").to_csv(os.path.join(OUT, "pathway_flight_meta_regression.tsv"), sep="\t", index=False)

    # Positive control: plant a flight effect proportional to the predicted penalty.
    sd_eff = float(np.median([U[f"{pw}__eff"].std() for pw in pathways]))
    detected = 0
    for pw in pathways:
        y = U[f"{pw}__eff"].to_numpy(float) + (x - x.mean()) / (x.std() or 1) * sd_eff
        detected += kh_meta_regression(y, np.maximum(U[f"{pw}__var"].to_numpy(float), 1e-12), x)[1] < 0.05

    lines = [
        "## 8. Meta-regression: does the predicted boundary-layer penalty explain which studies respond?\n",
        f"{len(U)} study x scale units from {U.study.nunique()} studies; the predictor takes "
        f"{U.penalty_dCO2.round(4).nunique()} distinct values (it varies only by hardware x scale).\n",
        "| unit | hardware | FLT/GC | predicted penalty (|ΔC CO₂| µg − Earth) |",
        "|---|---|---|---|",
    ] + [f"| {r.unit} | {r.hardware} | {r.n_flt}/{r.n_gc} | {r.penalty_dCO2:.3f} |" for r in U.itertuples()] + [
        "",
        f"- Knapp–Hartung (primary): flight effect scales with the predicted penalty at FDR < 0.05 in "
        f"**{int((R.fdr < 0.05).sum())}/{len(R)}** pathways (smallest p = {R.p_value.min():.3g}, "
        f"{R.sort_values('p_value').pathway.iat[0]}). The DerSimonian–Laird z-test, which is over-confident "
        f"with this few units, would have reported {int((bh(R.p_dl_z) < 0.05).sum())}/{len(R)} (smallest p = {R.p_dl_z.min():.2g}).",
        f"- Of the FDR-significant pathways, {int(((R.fdr < 0.05) & (R.p_veggie_indicator < 0.05)).sum())} are also "
        f"significant with a plain VEGGIE-vs-other indicator, and {int(((R.fdr < 0.05) & (R.p_bric_only_units < 0.05)).sum())} "
        "show the slope within BRIC units alone (where the penalty varies only by organ scale). "
        f"Leave-one-study-out: {int(((R.fdr < 0.05) & (R.p_loso_max < 0.05)).sum())} stay p < 0.05 whichever study is dropped.",
        ("- **Reading:** the signal is a VEGGIE-versus-sealed-hardware difference in the flight response, consistent "
         "with the boundary-layer physics but not separable from which studies (labs, years, organs) used VEGGIE. "
         "There is no evidence yet of a within-hardware dose-response to the predicted penalty."
         if (R.fdr < 0.05).any() else
         "- **Reading:** no pathway's flight effect tracks the predicted penalty after FDR correction, so the "
         "meta-regression gives no support for a boundary-layer dose-response at this number of study units."),
        f"- Positive control: a planted effect of 1 SD of the unit-level effects per SD of penalty is detected "
        f"(p < 0.05) in {detected}/{len(R)} pathways.",
        "- Caveats: few units, the penalty varies only by hardware × scale (so it is also confounded with which "
        "studies used which hardware), and T10's BRIC/CARA values are the original 3.8 s values, not steady state.",
    ]
    # Insert before the closing "Not done here" note so the sections stay in order.
    rpath = os.path.join(OUT, "REANALYSIS.md")
    text = open(rpath).read()
    block = "\n".join(lines) + "\n"
    higher = [int(k) for k in re.findall(r"\n## (\d+)\. ", text) if int(k) > 8]   # keep §8 ahead of §9+
    marker = f"\n## {min(higher)}. " if higher else "\n## Not done here"
    text = text.replace(marker, "\n" + block + marker, 1) if marker in text else text + "\n" + block
    open(rpath, "w").write(text)
    print("\n".join(lines))
    print(f"\nwrote {OUT}/pathway_flight_meta_regression.tsv, meta_regression_units.tsv; updated REANALYSIS.md")


if __name__ == "__main__":
    main()
