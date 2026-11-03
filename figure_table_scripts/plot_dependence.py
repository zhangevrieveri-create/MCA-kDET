"""Main-text power figure for the dependence models (M3)-(M5), with the added
competitors (Song--Chen, max-type covariance, KMD). Data: existing summaries
sim/results_v2_quick/E2_summary.csv and E9_summary.csv (E9 = new competitors on
the same E2 data). No simulation is run.

Run:  python -I scripts/plot_dependence.py   (from 正文重排稿_20261102)
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SIM = HERE.parents[2]                      # .../sim
sys.path.insert(0, str(SIM))
import pandas as pd
import matplotlib.pyplot as plt
import plot_publication as pp

pp.STYLE.update({
    "songchen": ("Song–Chen", "#D55E00", "s", (0, (5, 2))),
    "clx_cov": ("Max-type covariance", "#882255", "v", (0, (7, 2))),
    "kmd": ("KMD", "#56B4E9", "D", (0, (1.5, 1.5))),
})
METHODS = ["mca", "songchen", "clx_cov", "mst", "kmd", "disco"]
# ordered by how MCA-kDET compares with its strongest competitor: (M4), (M5), (M3)
PANELS = [("C2", "(a) M4: Asymmetric dependence", r"$\beta$"),
          ("D", "(b) M5: Sparse correlated block", "$s$"),
          ("C", "(c) M3: Pair dependence", r"$\beta$")]


def main():
    e2 = pd.read_csv(SIM / "results_v2_quick" / "E2_summary.csv")
    e9 = pd.read_csv(SIM / "results_v2_quick" / "E9_summary.csv")
    keys = ["model", "beta", "d", "k"]
    df = e2.merge(e9[keys + ["songchen", "clx_cov", "kmd", "reps"]], on=keys, suffixes=("", "_e9"))
    assert (df["reps"] == df["reps_e9"]).all()
    pp.configure_style()
    fig, axes = pp.make_canvas(1, 3, METHODS, "", height=3.15)
    rows = []
    for ax, (model, title, xlabel) in zip(axes[0], PANELS):
        sub = df[(df.model == model) & (df.beta > 0)].copy()
        sub["x"] = (sub.beta * 100).round().astype(int).clip(lower=1) if model == "D" else sub.beta
        if model == "D":
            sub["x"] = [int(-(-b * 100 // 1)) for b in sub.beta]   # s = ceil(100 beta)
        long = sub.melt(id_vars=["x"], value_vars=METHODS, var_name="method", value_name="power")
        pp.draw_curves(ax, long, "x", METHODS, title, xlabel, logx=(model == "D"))
        rows.append(long.assign(model=model))
    out = HERE.parent / "latex" / "figures" / "power_review" / "fig1_power_dependence.pdf"
    fig.savefig(out)
    fig.savefig(HERE / "fig1_power_dependence_preview.png", dpi=160)
    pd.concat(rows).to_csv(HERE / "fig1_power_dependence_data.csv", index=False)
    print("written", out)


if __name__ == "__main__":
    main()
