"""Generates supp_app_tables_sc.tex and supp_app_tables_fin.tex (tables of Section S10) from the result files of
sim/applications/single_cell/results and sim/applications/financial_v3/results.

    python gen_supp_app.py <applications dir> supp_app_tables.tex

Numbers are only read and rounded; nothing is recomputed except the
within-regime co-exceedance and innovation standard deviation, which use the
same formulas as financial_v3/report_v3.py (section 4) on prepared/<ds>/X_innov.csv.
"""
import os, sys
import numpy as np, pandas as pd

APP = sys.argv[1] if len(sys.argv) > 1 else "../applications"
OUT = sys.argv[2] if len(sys.argv) > 2 else "supp_app_tables.tex"
SC = os.path.join(APP, "single_cell", "results")
FN = os.path.join(APP, "financial_v3")
L = []


def f3(x):
    return "" if pd.isna(x) else ("$<$0.001" if x < 0.0005 else f"{x:.3f}")


def table(caption, label, colspec, header, rows, size="\\small", sep=None):
    L.extend(["\\begin{table}[ht]", f"\\caption{{{caption}\\label{{{label}}}}}", "\\begin{center}", size + (f"\\setlength{{\\tabcolsep}}{{{sep}}}" if sep else ""),
              f"\\begin{{tabular}}{{{colspec}}}", "\\toprule", header + " \\\\", "\\midrule"])
    L.extend(r + " \\\\" for r in rows)
    L.extend(["\\bottomrule", "\\end{tabular}", "\\end{center}", "\\end{table}", ""])


# ------------------------------------------------------------- single cell
st = pd.read_csv(os.path.join(SC, "animal_structure.csv"), encoding="utf-8-sig")
st = st.sort_values(["weeks", "verified_animal_block"])
stage = {1: "12", 2: "18--20", 3: "30"}
rows = [f"{r.verified_animal_block.replace('_', '\\_')} & {stage[r.group_number]} & {r.weeks} & {r.sex} & {r.tumors} & "
        f"{r.plates} & {r.available_cells} & {'fitting' if r.split == 'fit' else 'test'}" for r in st.itertuples()]
table("Mice of the single-cell application. Stage: group of Section~\\ref{sec:app-sc} (weeks after tumour initiation); "
      "Tumours and Plates: number of tumours and sequencing plates contributing cells; Cells: cells that pass the quality "
      "control of \\citet{marjanovic2020}; Role: role in the main split. Thirty cells are sampled from every mouse.",
      "tab:S_sc_mice", "lcclcccl", "Mouse & Stage & Weeks & Sex & Tumours & Plates & Cells & Role", rows)

cmp_ = pd.read_csv(os.path.join(SC, "comparison.csv"), encoding="utf-8-sig").set_index("configuration")
meth = ["MCA", "MCA_l1", "MCA_quad", "MCA_knn", "MCA_rf", "Sum", "Comb", "DISCO", "MMD", "MST"]
name = {"primary": "Main", "genes_50": "50 genes", "genes_200": "200 genes", "cells_20": "20 cells per mouse",
        "cells_40": "40 cells per mouse", "middle_20w_only": "20-week mice only", "one_tumor": "One tumour per 30-week mouse"}
name.update({f"split_{i:02d}": f"Split {i}" for i in range(1, 11)})
rows = []
for c in cmp_.index:
    rows.append(f"{name[c]} & " + " & ".join(f3(cmp_.loc[c, m]) for m in meth))
table("Exact whole-mouse permutation p-values in the single-cell application for MCA-kDET, its single learners, "
      "the sum and combined statistics, and the three competitors. The sensitivity analyses vary one element of the main "
      "analysis; Splits 1--10 repeat the main analysis with fitting seeds 20261028--20261037. With 8 test mice there are 560 "
      "assignments (210 when only the 20-week mice form the middle group), and the smallest attainable p-value is 1/560.",
      "tab:S_sc_p", "lcccccccccc",
      "Analysis & MCA & L1 & Quad & kNN & RF & Sum & Comb & DISCO & MMD & MST",
      rows, size="\\footnotesize", sep="3.5pt")

auc = pd.read_csv(os.path.join(SC, "animal_auc.csv"), encoding="utf-8-sig")
pv = auc.pivot_table(index=["target", "animal"], columns="learner", values="auc").reset_index()
rows = [f"{stage[r.target]} & {r.animal.replace('_', '\\_')} & {r.l1:.3f} & {r.quad:.3f} & {r.knn:.3f} & {r.rf:.3f}"
        for r in pv.itertuples()]
table("Held-out one-versus-rest AUC of each test mouse in the main single-cell analysis: the score of the mouse's stage, "
      "computed on its 30 cells against the test cells of the other stages. Scores were fitted on other mice. The values "
      "describe the fitted classifiers and are not tests.",
      "tab:S_sc_auc", "llcccc", "Stage & Mouse & L1 & Quad & kNN & RF", rows)

# ------------------------------------------------------------- finance
NSC = len(L)
R = os.path.join(FN, "results")
cal = pd.read_csv(os.path.join(R, "pvalues_jobs_calibration.csv"))
show = ["MCA", "MCA_l1", "MCA_quad", "MCA_knn", "MCA_rf", "DISCO", "MMD", "MST", "Schott"]
rows = []
for kind, kname in (("placebo", "Shifted regimes"), ("bootstrap", "Bootstrap null")):
    for rep in ("innov", "dir"):
        g = cal[(cal.kind == kind) & (cal.representation == rep)]
        rate = g.groupby("method").p_value.apply(lambda p: (p <= 0.05).mean())
        rows.append(f"{kname} & {'Innovations' if rep == 'innov' else 'Directions'} & " + " & ".join(f"{rate[m]:.3f}" for m in show))
table("Rejection rates at $\\alpha = 0.05$ of the week-label permutation of Algorithm~\\ref{M-alg:mca} under two nulls, "
      "$d = 30$, 200 replications each, $B = 199$. Shifted regimes: the VIX regime sequence is shifted circularly by 200 "
      "offsets between 52 and $T - 52$ weeks, with the real data and split. Bootstrap null: weekly innovation vectors are "
      "resampled independently, returns are regenerated from the fitted GARCH(1,1) filters, and the filters are re-estimated. "
      "Schott's test uses its asymptotic p-value. The 95\\% Monte Carlo range for 200 replications is $[0.020, 0.080]$.",
      "tab:S_fin_cal", "llccccccccc",
      "Null & Data & MCA & L1 & Quad & kNN & RF & DISCO & MMD & MST & Schott", rows, size="\\footnotesize", sep="3.5pt")

sp = pd.read_csv(os.path.join(R, "summary_strata_placebo.csv"))
rows = []
for Lw in (26, 52, 104):
    for rep in ("innov", "dir"):
        g = sp[(sp.L == Lw) & (sp.representation == rep)].set_index("method")
        rows.append(f"{Lw} & {'Innovations' if rep == 'innov' else 'Directions'} & " + " & ".join(f"{g.loc[m, 'rate']:.3f}" for m in show))
table("Rejection rates at $\\alpha = 0.05$ of the within-era permutation (labels permuted within calendar windows of $L$ "
      "weeks) under the 200 shifted regime sequences of Table~\\ref{tab:S_fin_cal}, $d = 30$, $B = 199$. Schott's test uses "
      "its asymptotic p-value, which the permutation does not affect.",
      "tab:S_fin_strata_cal", "llccccccccc",
      "$L$ & Data & MCA & L1 & Quad & kNN & RF & DISCO & MMD & MST & Schott", rows, size="\\footnotesize", sep="3.5pt")

sm = pd.read_csv(os.path.join(R, "summary_strata_main.csv")).set_index("setting")
rows = []
for ds, d in (("ind10", 10), ("ind30", 30), ("ind49", 49)):
    for rep in ("innov", "dir"):
        for Lw in ((26, 52, 104) if ds == "ind30" else (52,)):
            s = f"{ds} / {rep} / L = {Lw}"
            rows.append(f"{d} & {'Innovations' if rep == 'innov' else 'Directions'} & {Lw} & " +
                        " & ".join(f3(sm.loc[s, m]) for m in show))
table("P-values of the within-era permutation for the real VIX regimes, $B = 4999$. Schott's test uses its asymptotic "
      "p-value.", "tab:S_fin_strata", "cclccccccccc",
      "$d$ & Data & $L$ & MCA & L1 & Quad & kNN & RF & DISCO & MMD & MST & Schott", rows, size="\\footnotesize", sep="3.5pt")

wl = pd.read_csv(os.path.join(R, "pvalues_jobs_main.csv"))
wl["rep"] = wl.representation + np.where(wl.label_col == "regime", "", "_alt")
order = [("ind30", "raw"), ("ind30", "innov"), ("ind30", "dir"), ("ind30", "innov_pre"), ("ind30", "dir_pre"),
         ("ind30", "innov_alt"), ("ind30", "dir_alt"), ("ind10", "raw"), ("ind10", "innov"), ("ind10", "dir"),
         ("ind49", "raw"), ("ind49", "innov"), ("ind49", "dir")]
rname = {"raw": "Returns", "innov": "Innovations", "dir": "Directions", "innov_pre": "Innov., pre-sample",
         "dir_pre": "Dir., pre-sample", "innov_alt": "Innov., VIX 15/25", "dir_alt": "Dir., VIX 15/25"}
rows = []
for ds, rep in order:
    g = wl[(wl.dataset == ds) & (wl.rep == rep)].set_index("method").p_value
    rows.append(f"{ds[3:]} & {rname[rep]} & " + " & ".join(f3(g[m]) for m in show))
table("P-values of the week-label permutation of Algorithm~\\ref{M-alg:mca} for all settings of the protocol, $B = 4999$. "
      "Table~\\ref{tab:S_fin_cal} shows that this calibration is anti-conservative for MCA-kDET on these data; the values are "
      "reported for completeness and are not used for inference. Pre-sample: GARCH parameters fitted to 1990--1999; VIX 15/25: regimes defined by the thresholds 15 and 25.",
      "tab:S_fin_week", "clccccccccc",
      "$d$ & Data & MCA & L1 & Quad & kNN & RF & DISCO & MMD & MST & Schott", rows, size="\\footnotesize", sep="3.5pt")

rot = pd.read_csv(os.path.join(R, "summary_rotation.csv"))
rows = []
for rep in ("innov", "dir"):
    for m in show:
        x = rot[(rot.representation == rep) & (rot.method == m)].iloc[0]
        nm = {"MCA": "MCA-kDET", "MCA_l1": "L1", "MCA_quad": "Quad", "MCA_knn": "kNN", "MCA_rf": "RF"}.get(m, m)
        if m in ("DISCO",):
            vals = f"{x.T_obs:.1f} & {x.ref_median:.1f} & {x.ref_q95:.1f} & {x.ref_q99:.1f}"
        elif m == "MMD":
            vals = f"{x.T_obs:.4f} & {x.ref_median:.4f} & {x.ref_q95:.4f} & {x.ref_q99:.4f}"
        elif m == "MST":
            vals = f"{x.T_obs:.0f} & {x.ref_median:.0f} & {x.ref_q95:.1f} & {x.ref_q99:.1f}"
        else:
            vals = f"{x.T_obs:.2f} & {x.ref_median:.2f} & {x.ref_q95:.2f} & {x.ref_q99:.2f}"
        rows.append(f"{'Innovations' if rep == 'innov' else 'Directions'} & {nm} & {vals} & {x.p_rotation:.4f}")
table("Rotation test, $d = 30$: observed statistic, median and upper quantiles of its values over the 1292 circular "
      "shifts of the regime sequence, and the rotation p-value. For MST the statistic is the number of edges of the "
      "minimum spanning tree that join two observations of the same regime.",
      "tab:S_fin_rot", "llccccc", "Data & Method & Observed & Median & 95\\% & 99\\% & p-value", rows, size="\\footnotesize", sep="3.5pt")

eco = []
for ds in ("ind10", "ind30", "ind49"):
    wk = pd.read_csv(os.path.join(FN, "prepared", "ind30", "weeks.csv")); reg = wk.regime.values
    Z = pd.read_csv(os.path.join(FN, "prepared", ds, "X_innov.csv")).values; d = Z.shape[1]; iu = np.triu_indices(d, 1)
    R1 = np.corrcoef(Z[reg == 1].T); w = np.linalg.solve(R1, np.ones(d)); w /= w.sum(); pred = w @ R1 @ w
    for g in (1, 2, 3):
        Zg = Z[reg == g]; Rg = np.corrcoef(Zg.T); ev = np.sort(np.linalg.eigvalsh(Rg))[::-1]
        Lw = (Zg < np.quantile(Zg, 0.05, axis=0)).astype(float); cow = (Lw.T @ Lw / len(Lw))[iu].mean()
        eco.append(dict(ds=ds, d=d, g=g, n=len(Zg), sd=Zg.std(0).mean(), corr=Rg[iu].mean(),
                        absr=ev[:int(np.ceil(d / 5))].sum() / d, co=cow / 0.0025, mv=(w @ Rg @ w) / pred))
eco = pd.DataFrame(eco)
eco.to_csv("supp_app_economic.csv", index=False)
gname = {1: "Low", 2: "Medium", 3: "High"}
rows = [f"{r.d} & {gname[r.g]} & {r.n} & {r.sd:.2f} & {r.corr:.2f} & {r.absr:.2f} & {r.co:.1f} & {r.mv:.2f}" for r in eco.itertuples()]
table("Descriptive statistics of the GARCH innovations by VIX regime, all weeks of 2000--2026. Sd: average standard "
      "deviation; Corr: average pairwise correlation; Absorption: share of the $\\lceil d/5 \\rceil$ largest eigenvalues of "
      "the correlation matrix; Co-exceedance: average over pairs of industries of the probability that both fall below their "
      "own 5\\% quantile within the regime, divided by its value $0.0025$ under independence; MV variance: variance, under the "
      "regime's correlation matrix, of the minimum-variance portfolio built from the low-VIX correlation matrix, relative to "
      "its low-VIX variance.",
      "tab:S_fin_eco", "clcccccc", "$d$ & Regime & Weeks & Sd & Corr & Absorption & Co-exceedance & MV variance", rows)

base = OUT[:-4] if OUT.endswith(".tex") else OUT
open(base + "_sc.tex", "w", encoding="utf-8").write("\n".join(L[:NSC]) + "\n")
open(base + "_fin.tex", "w", encoding="utf-8").write("\n".join(L[NSC:]) + "\n")
print(eco.round(3).to_string())
print("wrote", base + "_sc.tex", base + "_fin.tex")
