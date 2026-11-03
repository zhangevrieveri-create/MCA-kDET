"""Second trimming round (main text 39 -> 38 pages). Results unchanged; the (M3)
dimension sentence moves to the supplement."""
import io
from pathlib import Path

L = Path(__file__).resolve().parents[1] / "latex"


def edit(name, pairs):
    p = L / name
    s = io.open(p, encoding="utf-8").read()
    for a, b in pairs:
        assert s.count(a) == 1, (name, a[:70])
        s = s.replace(a, b)
    io.open(p, "w", encoding="utf-8").write(s)


edit("sec4.tex", [
    (r"""Table~\ref{tab:size} reports empirical sizes under the three distinct null distributions, represented by (M1), (M3) and (M5), for $k = 3$, 5 and 10 with equal sample sizes and for $k = 3$ with sample sizes $(150, 300, 600)$. MCA-kDET has rejection rates between 0.043 and 0.054 at $\alpha = 0.05$, with no systematic change with $k$ or with unequal sample sizes.""",
     r"""Table~\ref{tab:size} reports empirical sizes under the null distributions of (M1), (M3) and (M5), for $k = 3$, 5 and 10 and for unequal sample sizes $(150, 300, 600)$. MCA-kDET has rejection rates between 0.043 and 0.054, with no systematic change with $k$ or the sample sizes."""),
    (r"""\textit{Serial dependence.} Theorem~\ref{thm:valid} does not cover serially dependent observations such as the weekly returns of Section~\ref{sec:app-fin}. Standardizing each coordinate by a univariate GARCH(1,1) filter \citep{bollerslev1986} and splitting the data into contiguous blocks separated by a buffer keeps the size of MCA-kDET between 0.042 and 0.058 under a GARCH(1,1) null with heavy-tailed innovations and regimes from a persistent Markov chain independent of the returns, against 0.192--0.861 without filtering (Section~\ref{S-sec:S9}, Table~\ref{S-tab:garch}).""",
     r"""\textit{Serial dependence.} Theorem~\ref{thm:valid} does not cover serially dependent data such as the returns of Section~\ref{sec:app-fin}. Under a GARCH(1,1) null with heavy-tailed innovations and persistent regimes, standardizing each coordinate by a univariate GARCH(1,1) filter \citep{bollerslev1986} and splitting into buffered contiguous blocks keeps the size between 0.042 and 0.058, against 0.192--0.861 without filtering (Table~\ref{S-tab:garch})."""),
    (r"""caption{Power at $\alpha = 0.05$ against the moment-matched shape alternatives, with $d = 20$, $n = 300$ per population and 100 replications per setting. RF is the test that uses the forest alone; SC-S and SC-SA are the two statistics of \citet{song2022}. The max-type mean and covariance tests, MMD and Schott's test have power at most 0.14 in every row (Table~\ref{S-tab:newshape-all}).""",
     r"""caption{Power at $\alpha = 0.05$ against the shape alternatives, $d = 20$, $n = 300$, 100 replications per setting. RF: the forest alone. The max-type tests, MMD and Schott's test have power at most 0.14 in every row (Table~\ref{S-tab:newshape-all})."""),
    (r"""caption{Power at $\alpha = 0.05$ against signal strength in the dependence models (M4), (M5) and (M3), with $k = 3$ and $n = 300$ per population. The horizontal axis is $\beta$ in (M3) and (M4) and the number of dependent coordinates $s$ in (M5). MMD, MMD-max, Schott's test, the max-type mean test and the accuracy test are omitted; all six models are shown in Supplementary Figure~\ref{S-fig:power-all}, and every method in Supplementary Table~\ref{S-tab:allmethods}.""",
     r"""caption{Power at $\alpha = 0.05$ in the dependence models, $k = 3$, $n = 300$; $s$ is the number of dependent coordinates. The other competitors are omitted; every method is in Supplementary Table~\ref{S-tab:allmethods}."""),
    (r"""Figure~\ref{fig:power} plots the power in (M4), (M5) and (M3). The results can be summarized as follows.""",
     r"""Figure~\ref{fig:power} plots the power in (M4), (M5) and (M3)."""),
    (r"""\item MCA-kDET has substantial power in all three models, where distance and kernel tests have almost none. In (M4) with $\beta = 0.6$, MCA-kDET has power 0.96, against 0.63 for MST and below 0.10 for DISCO and MMD. With three dependent coordinates among $d = 100$ in (M5), it has 0.88, against 0.12 for MST, at most 0.09 for DISCO and MMD, and none for Schott's test, and in (M3) with $\beta = 0.8$ it has 0.85, against 0.69 for MST and 0.05 for DISCO and MMD.""",
     r"""\item MCA-kDET has substantial power in all three models, where distance and kernel tests have almost none: 0.96 in (M4) with $\beta = 0.6$, 0.88 in (M5) with three dependent coordinates among $d = 100$, and 0.85 in (M3) with $\beta = 0.8$, against 0.63, 0.12 and 0.69 for MST and at most 0.10 for DISCO and MMD."""),
    (r"""In the sparse mean model (M1), the library pays the price of Remark~\ref{rem:price}: the L1 test is more powerful than the library at $\beta = 0.8$, 0.80 against 0.65, both above DISCO and MMD (0.53), and the two meet at $\beta = 1$ (0.95 and 0.94). The maximum over four learners raises the Bonferroni critical value of the standardized AUC from 2.13 to 2.64. The max-type mean test, which targets exactly this alternative and uses all observations, is more powerful at moderate signals (0.99 at $\beta = 0.8$); the two have the same detection order (Corollary~\ref{cor:sparse}), and the classifier test pays for learning on half of the data.""",
     r"""In the sparse mean model (M1), the library pays the price of Remark~\ref{rem:price}: at $\beta = 0.8$ the L1 test has power 0.80 and the library 0.65, both above DISCO and MMD (0.53), and the two meet at $\beta = 1$. The max-type mean test, which targets this alternative and uses all observations, has 0.99 at $\beta = 0.8$; the two have the same detection order (Corollary~\ref{cor:sparse}), and the classifier test pays for learning on half of the data."""),
    (r""" In (M3) with $n = 300$, no test has power above 0.09 once $d \ge 50$, whereas MCA-kDET has 0.99 for $d = 10$, with five dependent pairs, and 0.70 for $d = 20$, with ten.""", ""),
])

edit("supp_moresim.tex", [
    (r"""\textit{Library breadth.}""",
     r"""\textit{Dimension in (M3).} With $n = 300$, MCA-kDET has power 0.99 for $d = 10$, with five dependent pairs, and 0.70 for $d = 20$, with ten, but no test has power above 0.09 once $d \ge 50$ (Figure~\ref{fig:dimension}).

\textit{Library breadth.}"""),
])
print("round 2 done")
