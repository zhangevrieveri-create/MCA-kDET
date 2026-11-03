"""Supplementary figure: MCA-kDET statistic T_delta against the shift delta (TTS, main regimes).

    python fig_tts_profile.py <output.pdf>

Reads results/pvalues_jobs_tts.csv (written by run_v3.R for jobs_tts.csv).
"""
import os, sys
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.abspath(__file__))
out = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "results", "fig_tts_profile.pdf")
t = pd.read_csv(os.path.join(ROOT, "results", "pvalues_jobs_tts.csv"))
t = t[(t.method == "MCA") & t.label_col.str.startswith("main_d")].copy()
t["delta"] = t.label_col.str.extract(r"_d(-?\d+)$")[0].astype(int)

fig, axes = plt.subplots(1, 2, figsize=(7.5, 2.6), sharey=True)
for ax, (rep, name) in zip(axes, [("innov", "Innovations"), ("dir", "Directions")]):
    g = t[t.representation == rep].sort_values("delta")
    t0 = g.loc[g.delta == 0, "stat_obs"].iloc[0]
    ax.plot(g.delta, g.stat_obs, lw=0.5, color="0.5")
    ax.plot(g.delta, g.stat_obs.rolling(25, center=True, min_periods=5).mean(), lw=1.2, color="black")
    ax.axhline(t0, ls="--", lw=0.8, color="black")
    ax.axvline(0, lw=0.4, color="black")
    ax.set_title(name, fontsize=9)
    ax.set_xlabel(r"shift $\delta$ (weeks)", fontsize=8)
    ax.tick_params(labelsize=7)
axes[0].set_ylabel(r"$T_\delta$", fontsize=8)
fig.tight_layout()
fig.savefig(out)
print("wrote", out)
