"""Reorder Section 4 text: favourable results first. Edits latex/sec4.tex in place
(the GPT original is kept as ../sec4_GPT稿原文.tex)."""
import io
from pathlib import Path

p = Path(__file__).resolve().parents[1] / "latex" / "sec4.tex"
s = io.open(p, encoding="utf-8").read()


def sub(a, b):
    global s
    assert s.count(a) == 1, a[:80]
    s = s.replace(a, b)


# 1. second group: the three dependence models first, the mean model last
m1 = ("\\item[(M1)] $P_0 = N(0, \\Sigma)$ with $\\Sigma_{jl} = 0.5^{|j-l|}$, and $P_\\beta = "
      "N(\\mu_\\beta, \\Sigma)$ with $\\mu_{\\beta,j} = 0.6\\, \\beta\\, j^{-3}$.\n")
sub(m1, "")
sub("correlations $0.9^{|j-l|}$.\n\\end{enumerate}", "correlations $0.9^{|j-l|}$.\n" + m1 + "\\end{enumerate}")
sub("The second group contains a sparse linear departure and three departures that leave the mean vector unchanged.",
    "The second group contains three departures that leave the mean vector unchanged and a sparse linear departure.")
sub("Models (M1) and (M5) use $d = 100$, and (M3) and (M4) use $d = 20$",
    "Models (M5) and (M1) use $d = 100$, and (M3) and (M4) use $d = 20$")
sub("The grid of $\\beta$ is $\\{0.2, 0.4, \\dots, 1.0\\}$ for (M1) and (M3) and $\\{0.12, 0.24, \\dots, 0.60\\}$ for (M4)",
    "The grid of $\\beta$ is $\\{0.2, 0.4, \\dots, 1.0\\}$ for (M3) and (M1) and $\\{0.12, 0.24, \\dots, 0.60\\}$ for (M4)")

# 2. shape table: (N3), the one setting with a reversal, last
n3 = ("(N3) & 0.70 & 0.09 & 0.14 & 0.24 & 0.22 & 0.04 & 0.08 & 0.08 \\\\\n"
      "(N3) & 0.80 & 0.56 & 0.80 & 0.39 & 0.37 & 0.08 & 0.04 & 0.06 \\\\\n")
sub(n3, "")
sub("(N6) & 0.80 & 0.68 & 0.79 & 0.15 & 0.13 & 0.04 & 0.03 & 0.05 \\\\\n",
    "(N6) & 0.80 & 0.68 & 0.79 & 0.15 & 0.13 & 0.04 & 0.03 & 0.05 \\\\\n" + n3)

# 3. dependence figure and items in the order (M4), (M5), (M3)
sub("against signal strength in the dependence models (M3)--(M5)",
    "against signal strength in the dependence models (M4), (M5) and (M3)")
sub("Figure~\\ref{fig:power} plots the power in (M3)--(M5).",
    "Figure~\\ref{fig:power} plots the power in (M4), (M5) and (M3).")
sub("With three dependent coordinates among $d = 100$ in (M5), MCA-kDET has power 0.88, against 0.12 for MST, "
    "at most 0.09 for DISCO and MMD, and none for Schott's test. In (M3) with $\\beta = 0.8$ the power is 0.85 for "
    "MCA-kDET, 0.69 for MST and 0.05 for DISCO and MMD, and in (M4) with $\\beta = 0.6$ it is 0.96 for MCA-kDET, "
    "0.63 for MST and below 0.10 for DISCO and MMD.",
    "In (M4) with $\\beta = 0.6$, MCA-kDET has power 0.96, against 0.63 for MST and below 0.10 for DISCO and MMD. "
    "With three dependent coordinates among $d = 100$ in (M5), it has 0.88, against 0.12 for MST, at most 0.09 for "
    "DISCO and MMD, and none for Schott's test, and in (M3) with $\\beta = 0.8$ it has 0.85, against 0.69 for MST "
    "and 0.05 for DISCO and MMD.")
sub("the quadratic learner is the most powerful in (M5) and the forest in (M3) and (M4)",
    "the forest is the most powerful in (M4) and (M3) and the quadratic learner in (M5)")
sub("the pairwise variant is 0.30, 0.08 and 0.35 lower in (M3) with $\\beta = 0.8$, (M4) with $\\beta = 0.6$ and (M5) with $s = 3$.",
    "the pairwise variant is 0.08, 0.35 and 0.30 lower in (M4) with $\\beta = 0.6$, (M5) with $s = 3$ and (M3) with $\\beta = 0.8$.")

io.open(p, "w", encoding="utf-8").write(s)
print("reordered", p.name)
