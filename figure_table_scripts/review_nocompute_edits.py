"""ARS review, round 1: the revision items that need no computation (R2, R3, R6, the text
part of R4, S1, S2, S3, S5 and the acronym check). R5 (CLX sentence) is NOT applied here;
it awaits the author's confirmation. Backup: ../审稿修改前备份_20261103/."""
import io
from pathlib import Path

L = Path(__file__).resolve().parents[1] / "latex"


def edit(name, pairs):
    p = L / name
    s = io.open(p, encoding="utf-8").read()
    for a, b in pairs:
        assert s.count(a) == 1, (name, a[:80])
        s = s.replace(a, b)
    io.open(p, "w", encoding="utf-8").write(s)


# ---------------------------------------------------------------- abstract (R2, S2, S3)
edit("abstract.tex", [
    (r"""In simulations where the populations share their means and covariances but differ in shape, the test is much more powerful than distance, kernel, graph-based and moment tests, and against nonlinear dependence it is competitive with the strongest graph-based and max-type tests. Single-cell expression data illustrate whole-mouse calibration, and weekly industry returns illustrate time-shift calibration for serial dependence.""",
     r"""In simulations with $d = 20$ in which the populations share their means and covariances but differ in the shape of a few coordinates, the test is much more powerful than distance, median-bandwidth kernel, graph-based and moment tests; against nonlinear dependence it is competitive with the strongest graph-based and max-type tests, and its advantage declines when the departure is rotated or the dimension grows. Single-cell expression data illustrate whole-mouse calibration, and an exploratory analysis of weekly industry returns illustrates time-shift calibration for serial dependence and its limits."""),
])

# ---------------------------------------------------------------- introduction (R6, R2, R4 text, acronyms)
edit("sec1.tex", [
    (r"""kernel tests based on the maximum mean discrepancy \citep{gretton2012}""",
     r"""kernel tests based on the maximum mean discrepancy (MMD) \citep{gretton2012}"""),
    (r"""its area under the ROC curve (AUC) is computed""",
     r"""its area under the receiver operating characteristic curve (AUC) is computed"""),
    (r"""Training once and permuting held-out labels is known for two samples \citep{clemencon2009,kubler2022,chakravarti2023}; here one set of permutations covers the $k$ one-versus-rest comparisons and all learners together.""",
     r"""Training once and permuting held-out labels is known for two samples \citep{clemencon2009,kubler2022,chakravarti2023}, and for $k = 2$ the test reduces to that of \citet{chakravarti2023}; calibrating the maximum of a family of statistics by common permutations is standard in multiple testing \citep{westfall1993} and kernel aggregation \citep{schrab2023}. What is new is the family: the $k$ one-versus-rest comparisons and a library of learners, all fitted once on $\Dfit$, are covered by one set of permutations of the held-out labels."""),
    (r"""\item Comparing each population with its self-excluded reference loses no alternative,""",
     r"""\item The one-versus-rest scheme of multiclass classification, with the self-excluded mixture rather than the pooled sample as reference, loses no alternative,"""),
    (r"""\item In simulations, the test is much more powerful than distance, kernel, graph-based and moment tests when the populations differ in shape but share their means and covariances, and it is competitive with the strongest graph-based and max-type tests under nonlinear dependence. In single-cell expression data, with whole mice as units, it separates stages of tumour progression that distance and kernel tests do not separate.""",
     r"""\item In simulations with $d = 20$, the test is much more powerful than distance, median-bandwidth kernel, graph-based and moment tests when the populations differ in the shape of a few coordinates but share their means and covariances, and it is competitive with the strongest graph-based and max-type tests under nonlinear dependence; the advantage rests on learners that use the coordinates and declines after rotation and in higher dimension. In single-cell expression data, with whole mice as units, it rejects equality of the tumour stages where DISCO and MMD on the same test mice do not, although the combination of ten splits gives 0.054."""),
])

# ---------------------------------------------------------------- Section 2 (S1, S5)
edit("sec2.tex", [
    (r"""The four learners were selected in a pilot study with seeds not used elsewhere, from a larger set that also contained a feedforward network; Section~\ref{S-sec:S7} gives the study and the full specification.""",
     r"""The four learners were selected in a pilot study with seeds not used elsewhere, from a larger set that also contained a feedforward network (Section~\ref{S-sec:S7}); the pilot used the model families of Section~\ref{sec:sim}, so the simulations are not an independent assessment of this choice. In practice the library should be fixed before the test data are seen, with one learner per departure of interest: a sparse linear learner for location shifts, a quadratic learner for covariance changes, and a forest or nearest-neighbour learner for local nonlinear structure; a learner that duplicates a mechanism raises the critical value (Remark~\ref{rem:price}) without adding coverage."""),
])
s2 = io.open(L / "sec2.tex", encoding="utf-8").read()
a = r"""\textbf{Step 4 (Standardization).}"""
i = s2.index(a); j = s2.index("\n\n", i)
para = s2[i:j]
assert "Ties" not in para
io.open(L / "sec2.tex", "w", encoding="utf-8").write(
    s2[:j] + "\n\nTies, which occur for vote-based scores, reduce the null variance below $\\sigma_i^2$, so the standardization understates their statistics; this costs power, not validity, since Step 6 calibrates the maximum exactly." + s2[j:])

# ---------------------------------------------------------------- Section 4 (R2, R3, S3, acronyms)
edit("sec4.tex", [
    (r"""KMD \citep{huang2024}; MST \citep{friedman1979} and the statistics""",
     r"""the kernel measure of multi-sample dissimilarity (KMD) \citep{huang2024}; the minimum spanning tree (MST) test \citep{friedman1979} and the statistics"""),
    (r"""standardizing each coordinate by a univariate GARCH(1,1) filter""",
     r"""standardizing each coordinate by a univariate generalized autoregressive conditional heteroskedasticity (GARCH) filter of order (1,1)"""),
    (r"""The max-type covariance test stays at the level in (M3) and (M4), and Song--Chen has power at most 0.26 in (M5) with $s \le 4$.""",
     r"""The max-type covariance test stays at the level in (M3) and (M4), and Song--Chen has power at most 0.26 in (M5) with $s \le 4$. The advantage depends on the coordinates: after all populations are multiplied by one random orthogonal matrix, the power of MCA-kDET falls from 0.85 to 0.33 in (M3) with $\beta = 0.8$ and from 0.88 to 0.10 in (M5) with $s = 3$, because the forest and the quadratic learner use the coordinates, whereas the rotation-invariant Song--Chen statistic keeps 0.99 in (M3) (Supplementary Tables~\ref{S-tab:aucacc} and~\ref{S-tab:expanded-old})."""),
    (r"""In summary, MCA-kDET keeps its level and is most powerful against shape and dependence departures that moment, distance and kernel statistics do not see; a test built for one known alternative remains preferable for it.""",
     r"""In summary, MCA-kDET keeps its level, is the most powerful test against shape departures in a few coordinates with $d = 20$, and is close to the best competitor against dependence departures at strong signals; the specialized test of each model is more powerful at weak signals, and the advantage declines after rotation and for $d = 100$."""),
])

# ---------------------------------------------------------------- Section 5 (R4 text part, S2, acronyms)
edit("sec5.tex", [
    (r"""DISCO, MMD and MST use the same test cells and the same assignments.""",
     r"""DISCO, MMD and MST use the same test cells and the same assignments, so all tests are compared on the same eight test mice; since they need no fitting set, they could also use all 14 mice, which is not done here."""),
    (r"""Combined by twice the median (Remark~\ref{rem:split}), the ten splits, which reuse the same mice, give 0.054 for MCA-kDET.""",
     r"""Combined by twice the median (Remark~\ref{rem:split}), the ten splits, which reuse the same mice, give 0.054 for MCA-kDET, so once the dependence on the split is removed the evidence is not decisive at the 0.05 level."""),
    (r"""The regime of a week is set by the VIX \citep{cboe2026}""",
     r"""The regime of a week is set by the Cboe Volatility Index (VIX) \citep{cboe2026}"""),
    (r"""Read as exploratory evidence, the results show that,""",
     r"""Read as exploratory evidence, the results suggest that,"""),
])
print("edits applied")
