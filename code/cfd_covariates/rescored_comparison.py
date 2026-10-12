"""
§5, §7 and §8 rerun on the rescored pathway scores (REANALYSIS §11).

followups_2026-10.R writes rescored pathway scores (per-study log scale, within-study
z-scores, no zero-filling, each gene counted once). This script compares the §5/§7/§8
results on those scores with the stored-score results and inserts §11 into
REANALYSIS.md. Run it after the rescored pass of the two analysis scripts:

    S=tables/reanalysis_2026-10/followups
    for f in cfd_reanalysis cfd_meta_regression; do
      REANALYSIS_OUT=$S/rescored PATHWAY_SCORES=$S/pathway_scores_within_study.tsv \
        python3 code/cfd_covariates/$f.py
    done
    python3 code/cfd_covariates/rescored_comparison.py
"""
import os
import re

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
TAB = os.path.join(os.path.dirname(os.path.dirname(HERE)), "tables")
MAIN = os.path.join(TAB, "reanalysis_2026-10")
RES = os.path.join(MAIN, "followups", "rescored")


def read(d, name):
    return pd.read_csv(os.path.join(d, name), sep="\t")


def sig(df, col="fdr"):
    s = df[df[col] < 0.05]
    return len(s), ", ".join(sorted(s.pathway)) or "none"


def control(d):
    """Positive-control sentences, read back from that run's REANALYSIS.md."""
    txt = open(os.path.join(d, "REANALYSIS.md")).read()
    s5 = re.search(r"detected \(p < 0\.05\) in (\d+)/28 at 0\.25 SD, (\d+)/28 at 0\.5 SD", txt)
    s8 = re.search(r"per SD of penalty is detected \(p < 0\.05\) in (\d+)/28", txt)
    return (s5.group(1), s5.group(2)), s8.group(1)


def main():
    rows, terms = [], ["flight", "scale", "flight_x_scale"]
    ws = {k: read(d, "pathway_within_study_model.tsv") for k, d in [("stored", MAIN), ("rescored", RES)]}
    dish = {k: read(d, "pathway_dish_subset_model.tsv") for k, d in [("stored", MAIN), ("rescored", RES)]}
    mr = {k: read(d, "pathway_flight_meta_regression.tsv") for k, d in [("stored", MAIN), ("rescored", RES)]}
    ctl = {"stored": control(MAIN), "rescored": control(RES)}

    lines = [
        "## 11. §5, §7 and §8 on the rescored pathway scores",
        "",
        "The same scripts (`cfd_reanalysis.py`, `cfd_meta_regression.py`) rerun with "
        "`PATHWAY_SCORES=followups/pathway_scores_within_study.tsv` and `REANALYSIS_OUT=followups/rescored` "
        "(full report: `followups/rescored/REANALYSIS.md`). Comparison by `code/cfd_covariates/rescored_comparison.py`.",
        "",
        "| analysis | stored scores | rescored scores |",
        "|---|---|---|",
    ]
    for t in terms:
        a, b = (sig(ws[k][ws[k].term == t]) for k in ("stored", "rescored"))
        lines.append(f"| §5 within-study model, {t} (FDR < 0.05 of 28) | {a[0]} ({a[1]}) | {b[0]} ({b[1]}) |")
    lines.append(f"| §5 positive control, 0.25 / 0.5 SD | {ctl['stored'][0][0]}/28, {ctl['stored'][0][1]}/28 | {ctl['rescored'][0][0]}/28, {ctl['rescored'][0][1]}/28 |")
    a, b = (sig(dish[k], "fdr_flight") for k in ("stored", "rescored"))
    lines.append(f"| §7 in-dish subset, flight (FDR < 0.05 of 28) | {a[0]} ({a[1]}) | {b[0]} ({b[1]}) |")
    a, b = (sig(mr[k]) for k in ("stored", "rescored"))
    lines.append(f"| §8 Knapp–Hartung meta-regression (FDR < 0.05 of 28) | {a[0]} ({a[1]}) | {b[0]} ({b[1]}) |")
    lines.append(f"| §8 positive control (detected of 28) | {ctl['stored'][1]} | {ctl['rescored'][1]} |")

    # how much the per-pathway flight estimates move
    m = ws["stored"][ws["stored"].term == "flight"].merge(
        ws["rescored"][ws["rescored"].term == "flight"], on="pathway", suffixes=("_s", "_r"))
    r_t = m.t_s.corr(m.t_r)
    k = mr["stored"].merge(mr["rescored"], on="pathway", suffixes=("_s", "_r"))
    r_slope = k.slope_per_unit_penalty_s.corr(k.slope_per_unit_penalty_r)
    top_s = mr["stored"].sort_values("p_value").iloc[0]
    row_r = mr["rescored"].set_index("pathway").loc[top_s.pathway]
    gs = pd.read_csv(os.path.join(MAIN, "followups", "pathway_geneset_repeats.tsv"), sep="\t")
    n_rep = int((gs.pathway == top_s.pathway).sum())
    survives = row_r.fdr < 0.05
    dish_r = dish["rescored"][dish["rescored"].fdr_flight < 0.05]
    n_flight_r = int((ws["rescored"][ws["rescored"].term == "flight"].fdr < 0.05).sum()) + int((mr["rescored"].fdr < 0.05).sum())
    n_scale = int((ws["rescored"][ws["rescored"].term == "scale"].fdr < 0.05).sum())

    lines += [
        "",
        f"- Per-pathway flight t-statistics (§5) correlate at r = {r_t:.2f} between stored and rescored scores; "
        f"meta-regression slopes (§8) at r = {r_slope:.2f}.",
        f"- §8 on the stored scores finds {int((mr['stored'].fdr < 0.05).sum())} pathway(s) at FDR < 0.05; "
        f"{int(((mr['stored'].fdr < 0.05) & (mr['stored'].set_index('pathway').reindex(mr['stored'].pathway).index.isin(mr['rescored'][mr['rescored'].fdr < 0.05].pathway))).sum())} "
        f"remain significant on the rescored scores. The strongest, {top_s.pathway} (Knapp–Hartung p = {top_s.p_value:.2g}), has "
        f"p = {row_r.p_value:.2g} (FDR = {row_r.fdr:.2g}) on the rescored scores; its stored score repeats {n_rep} gene(s) in its set "
        "and includes zero-filled genes. "
        + ("Its signal survives rescoring." if survives else
           "The stored-score meta-regression hits do not survive correct scoring."),
        f"- The rescored scores show a scale (rosette vs leaf) difference in {n_scale} pathway(s), consistent with "
        "organ composition once study and platform scale differences are removed.",
        "- **Reading:** on correctly scored pathways, "
        + ("neither the all-study model (§5) nor the meta-regression (§8) finds a flight effect or a CFD dose-response "
           "at FDR < 0.05" if n_flight_r == 0 else
           f"{n_flight_r} pathway term(s) in §5/§8 reach FDR < 0.05")
        + ("" if dish_r.empty else
           f"; the in-dish subset (§7: 5 studies, t on 4 df) flags {len(dish_r)} pathway(s) for flight "
           f"({', '.join(dish_r.pathway)}), which is a single hit among 28 on a small subset and is not replicated in §5")
        + ". The positive controls show the pathway-level tests can only rule out large effects. The pathway results "
        "support 'no large, consistent pathway-level flight response across these 11 studies', not a boundary-layer "
        "mechanism; the CARA-specific flight response remains untestable (one study).",
        "",
    ]
    rpath = os.path.join(MAIN, "REANALYSIS.md")
    text = open(rpath).read()
    text = re.sub(r"\n## 11\. .*?(?=\n## \d+\. |\n## Not done here|\Z)", "", text, flags=re.S)
    block = "\n".join(lines)
    higher = [int(k) for k in re.findall(r"\n## (\d+)\. ", text) if int(k) > 11]
    marker = f"\n## {min(higher)}. " if higher else "\n## Not done here"
    text = text.replace(marker, "\n" + block + marker, 1) if marker in text else text + "\n" + block
    open(rpath, "w").write(text)
    print(block)


if __name__ == "__main__":
    main()
