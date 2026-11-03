"""Round 3: absorb the G2 sentence and keep 38 pages."""
import io
from pathlib import Path
p = Path(__file__).resolve().parents[1] / "latex" / "sec4.tex"
s = io.open(p, encoding="utf-8").read()
pairs = [
    (r"""Table~\ref{tab:shape} reports the power against the shape alternatives, departures that a classifier can learn but that leave the means and covariances unchanged. At oracle""",
     r"""Table~\ref{tab:shape} reports the power against the shape alternatives. At oracle"""),
    (r"""With the same fits and permutations, the accuracy of the forest is less powerful than its AUC, for instance 0.26 against 0.51 in (N1) and 0.50 against 0.79 in (N5) at oracle AUC 0.70 (Supplementary Table~\ref{S-tab:aucacc}), as the signal scales of Section~\ref{sec:why} suggest.""",
     r"""On the same fits and permutations, the accuracy of the forest is less powerful than its AUC, 0.26 against 0.51 in (N1) at oracle AUC 0.70 (Supplementary Table~\ref{S-tab:aucacc})."""),
    (r"""In summary, MCA-kDET keeps its level, its power is governed by the held-out AUC of the library, and it is most powerful against departures in shape and dependence that moment, distance and kernel statistics do not see; a test built for one known alternative, such as a max-type test for sparse moment changes, remains preferable for that alternative.""",
     r"""In summary, MCA-kDET keeps its level and is most powerful against shape and dependence departures that moment, distance and kernel statistics do not see; a test built for one known alternative remains preferable for it."""),
]
for a, b in pairs:
    assert s.count(a) == 1, a[:60]
    s = s.replace(a, b)
io.open(p, "w", encoding="utf-8").write(s)
print("ok")
