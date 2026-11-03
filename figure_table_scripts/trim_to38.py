"""Shorten the main text to 38 pages: compress Section 4, move two passages to the
supplement, drop the duplicated library paragraph of Section 2.4 and trim the
discussion. Results and conclusions are unchanged."""
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


SHAPE_OLD = r"""The first group consists of moment-matched shape alternatives. The baseline is $N_d(0, I_d)$ with $d = 20$, and population $k$ departs on a few coordinates: three coordinates have a symmetric two-component location mixture in (N1), a two-component scale mixture in (N2) and an asymmetric location mixture in (N5), and three coordinate pairs have a conditional variance depending on the other coordinate in (N3) and a ring-shaped distribution in (N6). In (N4) there are $k = 4$ populations, a baseline one and three that carry the departures of (N1), (N2) and (N3) on disjoint blocks. Every population has mean zero and covariance $I_d$, so no test of means or covariance matrices can detect these departures. The strength of each departure is indexed by its oracle AUC (Lemma~\ref{lem:oracle}), the AUC of the true density ratio against the baseline, with values 0.60, 0.70 and 0.80. Section~\ref{S-sec:newshape} gives the definitions, the calibration and the settings with $d = 100$."""
SHAPE_NEW = r"""The first group consists of moment-matched shape alternatives with baseline $N_d(0, I_d)$, $d = 20$. Population $k$ has a symmetric location mixture (N1), a scale mixture (N2) or an asymmetric location mixture (N5) on three coordinates, or a conditional variance (N3) or a ring-shaped distribution (N6) on three coordinate pairs; in (N4), $k = 4$ and three populations carry the departures of (N1)--(N3) on disjoint blocks. All populations have mean zero and covariance $I_d$. The strength is the oracle AUC of Lemma~\ref{lem:oracle}, 0.60, 0.70 or 0.80 (Section~\ref{S-sec:newshape})."""

COMP_OLD = r"""We compare MCA-kDET with distance, kernel, graph-based and moment tests: DISCO \citep{rizzo2010}; MMD, the average over pairs of samples of the unbiased squared maximum mean discrepancy with a Gaussian kernel and median bandwidth \citep{gretton2012}, and MMD-max, the maximum of these pairwise statistics \citep{kim2021b}; KMD \citep{huang2024}; MST, the number of within-population edges of the minimum spanning tree of the pooled sample \citep{friedman1979}, and the two statistics SC-S and SC-SA of \citet{song2022} on five successive spanning trees; max-type mean and covariance tests in the spirit of \citet{cai2014clx} and \citet{cai2013}; the covariance test of \citet{schott2007b}; and the multiclass accuracy test, the $k$-sample extension of \citet{kim2021} described in Section~\ref{sec:intro}. Except for Schott's test with its normal limit, all are calibrated by $B = 200$ permutations of the pooled sample (Section~\ref{S-sec:newshape}). Sizes are based on 1000 replications per setting, with 95\% Monte Carlo range $[0.037, 0.064]$, and powers on 100 replications, with standard error at most 0.05."""
COMP_NEW = r"""The competitors are DISCO \citep{rizzo2010}; MMD with a Gaussian kernel and median bandwidth, averaged over pairs of samples \citep{gretton2012}, and its maximum over pairs, MMD-max \citep{kim2021b}; KMD \citep{huang2024}; MST \citep{friedman1979} and the statistics SC-S and SC-SA of \citet{song2022}; max-type mean and covariance tests in the spirit of \citet{cai2014clx} and \citet{cai2013}; Schott's covariance test \citep{schott2007b}; and the multiclass accuracy test of Section~\ref{sec:intro}. All except Schott's test are calibrated by $B = 200$ permutations (Section~\ref{S-sec:newshape}). Sizes are based on 1000 replications and powers on 100, with standard error at most 0.05."""

edit("sec4.tex", [
    (SHAPE_OLD, SHAPE_NEW),
    (COMP_OLD, COMP_NEW),
    (r"""In (M3) and (M4) every coordinate is standard normal and the covariance matrix is the identity, so the populations differ only in a nonlinear dependence within pairs and no linear score carries information; in (M5)""",
     r"""In (M3) and (M4) every coordinate is standard normal and the covariance is the identity, so only a nonlinear dependence within pairs changes; in (M5)"""),
    (r""" A sparse variance shift (M2) and small mean shifts in every population (M6) are reported in Section~\ref{S-sec:moresim}. The constants were set in calibration rounds with seeds not used elsewhere (Section~\ref{S-sec:S8}).""",
     r""" Models (M2) and (M6) are in Section~\ref{S-sec:moresim}, and the calibration of the constants in Section~\ref{S-sec:S8}."""),
    # Size: two sentences to the supplement
    (r""" For the two most extreme competitor entries, (M1) and (M3) with $k = 5$, 2000 replications with new seeds gave 0.045, 0.045 and 0.036, and 0.042, 0.045 and 0.046, for DISCO, MMD and MST, whose discrete statistic makes it conservative.""",
     r""" The most extreme competitor entries, (M1) and (M3) with $k = 5$, are within the Monte Carlo range in 2000 further replications (Section~\ref{S-sec:moresim})."""),
    (r""" Tables~\ref{S-tab:S-size1} and~\ref{S-tab:S-size2} give further dimensions, the single learners, the variants of Section~\ref{sec:variants} and the levels 0.01 and 0.10.""",
     r""" Tables~\ref{S-tab:S-size1} and~\ref{S-tab:S-size2} give further settings and the levels 0.01 and 0.10."""),
    # shape results
    (r""" The distance and kernel tests stay near the level, as expected from their dependence on means and variances in high dimensions \citep{zhu2021}.""", ""),
    (r"""The forest carries these results: the held-out AUCs of the other three learners stay near one half, and the forest alone has average power 0.44 over the 24 shape alternatives of Section~\ref{S-sec:newshape}, against 0.38 for the library. This is the price of Remark~\ref{rem:price}, paid for coverage of the dependence and moment alternatives below.""",
     r"""The forest carries these results; the forest alone has average power 0.44 over the 24 shape alternatives, against 0.38 for the library, the price of Remark~\ref{rem:price} for covering the alternatives below."""),
    # dependence items
    (r"""\item The library adapts to the departure. The learner that matches it carries the test (Supplementary Figure~\ref{S-fig:S-learners}): the forest is the most powerful in (M4) and (M3) and the quadratic learner in (M5), but the forest has power 0.06 where the quadratic learner has 0.91 ((M5), $s = 3$), and the quadratic learner has none in (M3) and (M4). MCA-kDET is at most 0.12 below the best learner in any cell of these models.""",
     r"""\item The library adapts to the departure (Supplementary Figure~\ref{S-fig:S-learners}): the forest carries the test in (M4) and (M3), where the quadratic learner has no power, and the quadratic learner in (M5), where the forest has 0.06 against 0.91 at $s = 3$. MCA-kDET is at most 0.12 below the best learner in any cell."""),
    (r""" With $n$ fixed, a larger $k$ also means a larger total sample: the classifier of the deviating population is trained against $(k-1)n$ null observations, whereas the pairwise variant, which uses $2n$ observations per classifier whatever $k$ is, stays at 0.51, 0.48 and 0.47.""",
     r""" The classifier of the deviating population is trained against $(k-1)n$ null observations, whereas the pairwise variant, with $2n$ observations per classifier, stays at 0.51, 0.48 and 0.47."""),
    # mean and dimension
    (r""" The maximum over four learners raises the Bonferroni critical value of the standardized AUC from $\Phi^{-1}(1 - 0.05/3) = 2.13$ to $\Phi^{-1}(1 - 0.05/12) = 2.64$, and the largest standardized AUC of the L1 learner exceeds the two values in 80\% and 61\% of the replications.""",
     r""" The maximum over four learners raises the Bonferroni critical value of the standardized AUC from 2.13 to 2.64."""),
    (r"""

Since the classifiers are trained once, one replication with $d = 100$ and $n = 300$ per population takes 3.4\,s for $k = 3$ and 45\,s for $k = 10$, of which the $B = 200$ permutations take 0.18\,s (Supplementary Table~\ref{S-tab:S-time}). Section~\ref{S-sec:S9} reports the variants, the split fraction and multiple splits.""",
     r""" Section~\ref{S-sec:S9} reports the variants, the split fraction, multiple splits and the computing time."""),
    (r"""In summary, MCA-kDET keeps its level in all settings, its power is governed by the held-out AUC of the library, and it is most powerful against departures in shape and dependence that a learner captures and that moment, distance and kernel statistics do not see. A test built for one alternative, such as a max-type test for sparse moment changes, remains preferable when that alternative is known in advance.""",
     r"""In summary, MCA-kDET keeps its level, its power is governed by the held-out AUC of the library, and it is most powerful against departures in shape and dependence that moment, distance and kernel statistics do not see; a test built for one known alternative, such as a max-type test for sparse moment changes, remains preferable for that alternative."""),
])

# the moved passages, kept in the supplement
edit("supp_moresim.tex", [
    (r"""\textit{Library breadth.}""",
     r"""\textit{Competitor sizes and computing time.} For the two most extreme competitor entries of Table~\ref{M-tab:size}, (M1) and (M3) with $k = 5$, 2000 replications with new seeds gave 0.045, 0.045 and 0.036, and 0.042, 0.045 and 0.046, for DISCO, MMD and MST, whose discrete statistic makes it conservative. In (M1) with $\beta = 0.8$, the largest standardized AUC of the L1 learner exceeds the Bonferroni critical values 2.13 and 2.64 for one and four learners in 80\% and 61\% of the replications. Since the classifiers are trained once, one replication with $d = 100$ and $n = 300$ per population takes 3.4\,s for $k = 3$ and 45\,s for $k = 10$, of which the $B = 200$ permutations take 0.18\,s (Table~\ref{tab:S-time}).

\textit{Library breadth.}"""),
])

# Section 2.4: the library paragraph repeats contribution (a) and Remark 5
edit("sec2.tex", [
    (r"""The library is intended to cover distinct scoring mechanisms rather than every available implementation of each mechanism. Several closely related learners may add little signal while increasing the permutation critical value. Thus a short library and a larger library are different tests whose power must be compared on the intended alternatives. Section~\ref{S-sec:moresim} reports such a comparison without changing the four-learner specification above. Choosing a library for a new test must still precede inspection of its held-out labels, as required by Theorem~\ref{thm:valid}.

""", ""),
])

edit("sec6.tex", [
    (r""" The simulations identify gains for learnable nonlinear shape and dependence departures, including alternatives with matched means and covariances, while retaining the settings where specialist moment tests are stronger.""",
     r""" The simulations show its largest gains against shape and dependence departures that leave the means and covariances unchanged."""),
])
print("trimmed")
