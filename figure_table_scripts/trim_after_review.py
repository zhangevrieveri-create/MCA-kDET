"""Keep 38 pages after the review edits: shorter wording only, no result removed."""
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


edit("sec1.tex", [
    (r"""\item In simulations with $d = 20$, the test is much more powerful than distance, median-bandwidth kernel, graph-based and moment tests when the populations differ in the shape of a few coordinates but share their means and covariances, and it is competitive with the strongest graph-based and max-type tests under nonlinear dependence; the advantage rests on learners that use the coordinates and declines after rotation and in higher dimension. In single-cell expression data, with whole mice as units, it rejects equality of the tumour stages where DISCO and MMD on the same test mice do not, although the combination of ten splits gives 0.054.""",
     r"""\item In simulations with $d = 20$, the test is much more powerful than distance, median-bandwidth kernel, graph-based and moment tests against shape departures with matched means and covariances, and competitive with the strongest graph-based and max-type tests under nonlinear dependence; the advantage declines after rotation and in higher dimension. In single-cell data, with whole mice as units, it rejects equality of tumour stages where DISCO and MMD do not, although ten combined splits give 0.054."""),
    (r"""calibrating the maximum of a family of statistics by common permutations is standard in multiple testing \citep{westfall1993} and kernel aggregation \citep{schrab2023}. What is new is the family: the $k$ one-versus-rest comparisons and a library of learners, all fitted once on $\Dfit$, are covered by one set of permutations of the held-out labels.""",
     r"""calibrating a maximum by common permutations is standard in multiple testing \citep{westfall1993} and kernel aggregation \citep{schrab2023}. What is new is the family: $k$ one-versus-rest comparisons and a library of learners, fitted once, calibrated by one set of permutations of the held-out labels."""),
])

edit("sec2.tex", [
    (r"""Ties, which occur for vote-based scores, reduce the null variance below $\sigma_i^2$, so the standardization understates their statistics; this costs power, not validity, since Step 6 calibrates the maximum exactly.""",
     r"""Ties, as for vote-based scores, reduce the null variance below $\sigma_i^2$; this costs power, not validity, since Step 6 calibrates the maximum exactly."""),
])

edit("sec5.tex", [
    (r"""The ten assignments that keep these three mice in the middle stage leave this comparison unchanged, so 0.018 is the smallest p-value the design allows when the middle-stage comparison attains the maximum.""",
     r"""The ten assignments that keep these three mice in the middle stage leave this comparison unchanged, so 0.018 is the smallest p-value the design allows here."""),
])

edit("sec6.tex", [
    (r"""Because Theorem~\ref{thm:valid} only uses that the statistic is a function of the test observations and their labels once $\Dfit$ is fixed, a standardized energy or kernel statistic computed on $\Dtest$ can enter the same maximum and be calibrated by the same permutations without any loss of validity; the cost is the larger critical value discussed in Remark~\ref{rem:price}.""",
     r"""Since Theorem~\ref{thm:valid} only needs the statistic to be a function of the test observations and labels given $\Dfit$, a standardized energy or kernel statistic on $\Dtest$ can enter the same maximum and permutations, at the cost of Remark~\ref{rem:price}."""),
])

edit("sec4.tex", [
    (r"""The advantage depends on the coordinates: after all populations are multiplied by one random orthogonal matrix, the power of MCA-kDET falls from 0.85 to 0.33 in (M3) with $\beta = 0.8$ and from 0.88 to 0.10 in (M5) with $s = 3$, because the forest and the quadratic learner use the coordinates, whereas the rotation-invariant Song--Chen statistic keeps 0.99 in (M3)""",
     r"""The advantage depends on the coordinates: after one random orthogonal rotation of all populations, MCA-kDET falls from 0.85 to 0.33 in (M3) with $\beta = 0.8$ and from 0.88 to 0.10 in (M5) with $s = 3$, since the forest and the quadratic learner use the coordinates, whereas the rotation-invariant Song--Chen statistic keeps 0.99 in (M3)"""),
])
print("trimmed")
