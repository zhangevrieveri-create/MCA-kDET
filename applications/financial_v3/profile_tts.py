"""Amendment 2, descriptive: profile of the TTS statistics T_delta against the shift delta.

    python profile_tts.py [--res=<results dir>]

The feedback of returns on the VIX can only enter alignments with delta >= 1 (the return of
week t paired with the regime of a later week), so a higher profile for delta > 0 than for
delta < 0 indicates feedback and a conservative TTS p-value. Descriptive only: no p-value.
Reads <res>/pvalues_jobs_tts.csv; writes <res>/tts_profile.csv and <res>/tts_profile.pdf.
"""
import os, sys
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.abspath(__file__))
RES = next((a.split("=", 1)[1] for a in sys.argv[1:] if a.startswith("--res=")), os.path.join(ROOT, "results"))
BANDS = [("near", 1, 52), ("far", 53, None)]                 # |delta| in 1..52 and 53..r

t = pd.read_csv(os.path.join(RES, "pvalues_jobs_tts.csv"))
t["variant"] = t.label_col.str.extract(r"^(\w+?)_d")[0]
t["delta"] = t.label_col.str.extract(r"_d(-?\d+)$")[0].astype(int)
t = t[t.method == "MCA"].dropna(subset=["stat_obs"])
r = int(t.delta.abs().max())

rows = []
for (var, rep), g in t.groupby(["variant", "representation"]):
    t0 = g.loc[g.delta == 0, "stat_obs"]
    for name, lo, hi in BANDS:
        hi = r if hi is None else hi
        pos = g.loc[g.delta.between(lo, hi), "stat_obs"]
        neg = g.loc[g.delta.between(-hi, -lo), "stat_obs"]
        rows.append(dict(variant=var, representation=rep, band=f"{lo}..{hi}",
                         T_obs=t0.iloc[0] if len(t0) == 1 else np.nan,
                         n_pos=len(pos), mean_pos=pos.mean(), median_pos=pos.median(),
                         n_neg=len(neg), mean_neg=neg.mean(), median_neg=neg.median(),
                         share_pos_ge_obs=(pos >= t0.iloc[0] - 1e-10).mean() if len(t0) == 1 else np.nan,
                         share_neg_ge_obs=(neg >= t0.iloc[0] - 1e-10).mean() if len(t0) == 1 else np.nan))
s = pd.DataFrame(rows)
s.to_csv(os.path.join(RES, "tts_profile.csv"), index=False)
print(s.to_string(index=False))

combos = sorted(t.groupby(["variant", "representation"]).groups)
fig, axes = plt.subplots(len(combos), 1, figsize=(7, 2.2 * len(combos)), sharex=True, squeeze=False)
for ax, (var, rep) in zip(axes[:, 0], combos):
    g = t[(t.variant == var) & (t.representation == rep)].sort_values("delta")
    ax.plot(g.delta, g.stat_obs, lw=0.6, color="0.45")
    ax.plot(g.delta, g.stat_obs.rolling(25, center=True, min_periods=5).mean(), lw=1.4, color="C0")
    t0 = g.loc[g.delta == 0, "stat_obs"]
    if len(t0) == 1:
        ax.axhline(t0.iloc[0], ls="--", lw=0.8, color="C3")
    ax.axvline(0, lw=0.5, color="k")
    ax.set_ylabel(r"$T_\delta$")
    ax.set_title(f"{var}, {rep}", fontsize=9, loc="left")
axes[-1, 0].set_xlabel(r"shift $\delta$ (weeks)")
fig.tight_layout()
fig.savefig(os.path.join(RES, "tts_profile.pdf"))
print("wrote", os.path.join(RES, "tts_profile.csv"), "and tts_profile.pdf")
