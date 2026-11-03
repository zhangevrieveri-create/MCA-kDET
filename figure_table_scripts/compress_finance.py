"""Review item S2: compress Section 5.2 (finance) to keep 38 pages. All facts the reviewers
asked for are kept (exploratory status, calibration chosen in steps, failed TTS assumption,
Holm 0.050, lag-two results); details stay in Supplement S10. Also shortens the library
guidance added to Section 2.4."""
import io
from pathlib import Path

L = Path(__file__).resolve().parents[1] / "latex"
p = L / "sec5.tex"
s = io.open(p, encoding="utf-8").read()

i0 = s.index("Correlations among equity portfolios")
i1 = s.index("\\begin{table}", i0)
head = r"""Correlations among equity portfolios rise in volatile periods \citep{longin2001, ang2002, forbes2002}. We ask whether, once the volatility of each of 30 industries is removed, the joint distribution of weekly returns is the same in calm, normal and turbulent weeks. We use $T = 1395$ weekly log returns of the 30 value-weighted industry portfolios of \citet{french2026}, January 2000 to September 2026. The regime of a week is set by the Cboe Volatility Index (VIX) \citep{cboe2026} at the close of the previous week, with the 1990--1999 terciles 14.81 and 19.80 as thresholds (433, 433 and 529 weeks). A GARCH(1,1) filter \citep{bollerslev1986}, fitted without the regime labels, gives innovations $z_t$, which remain more dispersed in high-VIX weeks (average standard deviations 0.85, 0.96 and 1.13). Two hypotheses, Holm-adjusted, concern $z_t$ and its direction $u_t = z_t/\|z_t\|$, which removes a common scale: if $z_t = s_t A w_t$ with $w_t$ spherically symmetric and independent of $s_t > 0$, then $u_t$ depends on $A$ only.

Random halves of the 26-week blocks, separated by 4-week buffers, form $\Dfit$ (702 fitting and 585 test weeks). Because the dependence among industries drifts over 2000--2026, the week-label permutation of Step 6 lets flexible learners separate eras, and we calibrate by the truncated time-shift (TTS) test of \citet{yuan2024}: on the middle $T - 2r = 837$ weeks, $r = 279$ (388 fitting and 381 test weeks), week $t$ receives the regime of week $t + \delta$ for $|\delta| \le r$, the whole procedure is rerun with the split and seeds fixed, and $p = \min\{1, B/(r+1)\}$ with $B = \#\{\delta\colon T_\delta \ge T_0\}$. Its level guarantee requires the regime sequence to be stationary and independent of the whole return series, which is stronger than the hypothesis of interest and fails here because the VIX responds to returns. The calibration was also chosen in steps: after the week-label permutation failed a placebo check (Table~\ref{S-tab:S_fin_cal}), a rotation test was designated as primary, and the TTS test replaced it after the rotation and within-era results had been seen, though before any TTS statistic was computed (Section~\ref{S-sec:S10}). The p-values of this section are therefore exploratory.

"""
s = s[:i0] + head + s[i1:]

j0 = s.index(r"Table~\ref{tab:fin} reports the results.")
tail = r"""Table~\ref{tab:fin} reports the results. Under the TTS test MCA-kDET gives $p = 11/280 = 0.039$ for the innovations and $p = 7/280 = 0.025$ for the directions; both Holm-adjusted values equal $14/280 = 0.050$, the boundary of the 0.05 level, at which the rule $p \le \alpha$ rejects. On the directions DISCO and MMD give 0.946 and 1.000 and the linear learner 0.807, whereas the forest gives 0.018. The within-era and rotation checks, which have no proved level, give p-values of at most 0.038, and Schott's statistic gives 0.032 on the directions, 0.064 after the Holm adjustment. The largest statistics occur in the high-VIX regime (Remark~\ref{rem:local}, Table~\ref{S-tab:finT}). With the regime set by the VIX two weeks earlier, the TTS test gives $p = 0.207$ and 0.104 (Table~\ref{S-tab:S_fin_tts}); feedback from returns to the VIX can enter only the shifts with $\delta \ge 1$ (Section~\ref{S-sec:S10}).

Read as exploratory evidence, the results suggest that the joint distribution of the GARCH-standardized innovations differs across regimes of the lagged VIX and that the difference is not confined to a common scale; from low- to high-VIX weeks the average pairwise correlation of the innovations rises from 0.47 to 0.62 (Table~\ref{S-tab:S_fin_eco}). The design does not test whether the VIX adds information given the past returns, and because the regimes are not assigned at random it describes the difference rather than its cause.
"""
s = s[:j0] + tail
io.open(p, "w", encoding="utf-8").write(s)

q = L / "sec2.tex"
t = io.open(q, encoding="utf-8").read()
a = r"""In practice the library should be fixed before the test data are seen, with one learner per departure of interest: a sparse linear learner for location shifts, a quadratic learner for covariance changes, and a forest or nearest-neighbour learner for local nonlinear structure; a learner that duplicates a mechanism raises the critical value (Remark~\ref{rem:price}) without adding coverage."""
b = r"""In practice the library should be fixed in advance with one learner per departure of interest, since a learner that duplicates a mechanism raises the critical value (Remark~\ref{rem:price}) without adding coverage."""
assert t.count(a) == 1
io.open(q, "w", encoding="utf-8").write(t.replace(a, b))
print("compressed")
