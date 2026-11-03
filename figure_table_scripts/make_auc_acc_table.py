"""Reviewer item G2: power of the AUC statistic and of the multiclass accuracy computed
from the SAME one-versus-rest fits and the SAME permutations. Sources:
  shape alternatives  sim/motivated_v4_20261102/full_B200/report_tables/all_method_rates_and_intervals.csv
                      (mca_rf vs macc_rf)
  rotated (M3)-(M5)   sim/results_v2_quick/E10_summary.csv (library: mca vs macc_ovr; forest: mca_rf vs macc_rf)
Writes latex/supp_aucacc.tex. No simulation is run."""
import csv
from pathlib import Path

HERE = Path(__file__).resolve().parent
SIM = HERE.parents[2]
R = {}
for r in csv.DictReader(open(SIM / "motivated_v4_20261102/full_B200/report_tables/all_method_rates_and_intervals.csv", encoding="utf-8")):
    R[(r["model"], r["target"], r["d"], r["method"])] = float(r["rate"])
SH = [("N1_bimodal", "(N1)"), ("N2_scale_mix", "(N2)"), ("N4_mixed", "(N4)"), ("N5_skewmix", "(N5)"),
      ("N6_ring", "(N6)"), ("N3_condvar", "(N3)")]
rows = []
for m, lab in SH:
    for t in ["0.7", "0.8"]:
        g = lambda x: R[(m, t, "20", x)]
        rows.append(f"{lab} & oracle AUC {float(t):.2f} & -- & -- & {g('mca_rf'):.2f} & {g('macc_rf'):.2f} \\\\")
rows[-1] += "[2pt]"
e10 = list(csv.DictReader(open(SIM / "results_v2_quick/E10_summary.csv", encoding="utf-8")))
LAB = {"C2": "(M4)", "D": "(M5)", "C": "(M3)"}
for m in ["C2", "D", "C"]:
    for r in sorted((r for r in e10 if r["model"] == m), key=lambda r: float(r["beta"])):
        assert r["rot"] == "1" and r["reps"] == "100"
        b = float(r["beta"])
        setting = f"$s = {-(-round(b * 1e6) // 10000)}$" if m == "D" else f"$\\beta = {b:.2f}$"
        f = lambda x: f"{float(r[x]):.2f}"
        rows.append(f"{LAB[m]}, rotated & {setting} & {f('mca')} & {f('macc_ovr')} & {f('mca_rf')} & {f('macc_rf')} \\\\")
    if m != "C":
        rows[-1] += "[2pt]"

tex = r"""\begin{table}[tbp]
\renewcommand{\baselinestretch}{1}\small
\caption{AUC and accuracy computed from the same one-versus-rest fits with the same permutations, $k = 3$ (4 in (N4)), $n = 300$, $B = 200$, 100 replications per setting. The accuracy assigns a test observation to the population with the largest score of a learner; for the library, the maximum accuracy over the four learners is used. Upper part: shape alternatives with $d = 20$ (the library there contains further components, so only the forest is compared). Lower part: (M3)--(M5) after multiplication of all populations by one fixed random orthogonal matrix.\label{tab:aucacc}}
\begin{center}
\footnotesize\setlength{\tabcolsep}{5pt}
\begin{tabular}{llcccc}
\toprule
 & & \multicolumn{2}{c}{Library} & \multicolumn{2}{c}{Forest} \\
Model & Setting & AUC & Accuracy & AUC & Accuracy \\
\midrule
""" + "\n".join(rows) + r"""
\bottomrule
\end{tabular}
\end{center}
\end{table}
"""
(HERE.parent / "latex" / "supp_aucacc.tex").write_text(tex, encoding="utf-8")
print(len(rows), "rows")
