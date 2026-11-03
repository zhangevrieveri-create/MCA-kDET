"""Supplementary table with every method on the same replications (E2 + E9,
results_v2_quick; E9 recomputes the added competitors on the E2 data).
Writes latex/supp_allmethods.tex. No simulation is run."""
import csv
from pathlib import Path

HERE = Path(__file__).resolve().parent
SIM = HERE.parents[2]
rd = lambda f: list(csv.DictReader(open(SIM / "results_v2_quick" / f, encoding="utf-8")))
key = lambda r: (r["model"], r["beta"], r["d"], r["k"])
e2 = {key(r): r for r in rd("E2_summary.csv")}
e9 = {key(r): r for r in rd("E9_summary.csv")}
ORDER = ["C2", "D", "C", "A", "B"]                         # as in the main text, then (M1), (M2)
LAB = {"A": "(M1)", "B": "(M2)", "C": "(M3)", "C2": "(M4)", "D": "(M5)"}
OLD = ["mca", "macc", "disco", "mmd", "mst", "schott"]
NEW = ["mmd_max", "kmd", "songchen", "songchen_A", "maxmean", "clx_mean", "clx_cov", "schott_perm"]
f2 = lambda x: "%.2f" % float(x)
f3 = lambda x: "%.3f" % float(x)


def xlab(m, b):
    return "$s = %d$" % -(-round(float(b) * 1e6) // 10000) if m == "D" else "$\\beta = %s$" % f2(b)


rows = []
for m in ORDER:
    for k in sorted((k for k in e2 if k[0] == m), key=lambda k: float(k[1])):
        a, b = e2[k], e9[k]
        assert a["reps"] == b["reps"]
        rows.append(" & ".join([LAB[m], xlab(m, k[1])] + [f2(a[c]) for c in OLD] + [f2(b[c]) for c in NEW]) + " \\\\")
    if m != ORDER[-1]:
        rows[-1] += "[2pt]"
nulls = []
for k in sorted((k for k in e9 if k[1] == "0"), key=lambda k: (ORDER.index(k[0]), int(k[2]))):
    r = e9[k]
    nulls.append(" & ".join([LAB[k[0]], "$d = %s$" % k[2]] + ["--"] * 6 + [f3(r[c]) for c in NEW]) + " \\\\")

head = ("Model & Setting & MCA & Acc. & DISCO & MMD & MST & Schott & MMD-max & KMD & SC-S & SC-SA & "
        "Max-mean & CLX-mean & Max-cov & Schott-p \\\\")
tex = r"""\begin{table}[tbp]
\renewcommand{\baselinestretch}{1}\small
\caption{Power at $\alpha = 0.05$ of every method on the same replications in models (M1)--(M5), $k = 3$, $n = 300$, 100 replications per setting (upper part), and sizes of the added competitors under the corresponding null distributions, 200 replications per setting (lower part; the sizes of the first six methods are in Table~\ref{M-tab:size} and Table~\ref{tab:S-size2}). Acc.: multiclass accuracy test; SC-S, SC-SA: the statistics of Song and Chen; Max-mean and Max-cov: coordinatewise max-type mean and covariance statistics; CLX-mean: the max-type mean statistic after a label-free inverse pooled-covariance transform; Schott-p: Schott's statistic calibrated by permutation. All except Schott's test use $B = 200$ permutations.\label{tab:allmethods}}
\begin{center}
\renewcommand{\baselinestretch}{1}\scriptsize\setlength{\tabcolsep}{2.2pt}
\resizebox{\textwidth}{!}{\begin{tabular}{llcccccccccccccc}
\toprule
""" + head + "\n\\midrule\n" + "\n".join(rows) + "\n\\midrule\n" + "\n".join(nulls) + r"""
\bottomrule
\end{tabular}}
\end{center}
\end{table}
"""
(HERE.parent / "latex" / "supp_allmethods.tex").write_text(tex, encoding="utf-8")
print(len(rows), "power rows,", len(nulls), "null rows")
