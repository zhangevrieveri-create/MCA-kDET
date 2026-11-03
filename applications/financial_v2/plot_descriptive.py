"""Fixed-order descriptive MMD map; no p-value or significance annotation."""
from pathlib import Path
import argparse
import json
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

ROOT = Path(__file__).resolve().parent
QA = Path.home()/".codex/skills/nature-figure/scripts"
sys.path.insert(0, str(QA))
from audit_panel_alignment import require_matplotlib_panel_alignment


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--representation",choices=["raw","innovations"],default="raw")
    args=parser.parse_args()
    innovations=args.representation=="innovations"
    out = ROOT/("innovations/descriptive" if innovations else "descriptive")
    panel = pd.read_csv(ROOT/"design_two_tests/main_panel.csv")
    matrix = pd.read_csv(out/"distance_matrix_d5.csv", index_col=0)
    assert matrix.index.tolist() == panel.order_book_id.tolist()
    D = matrix.to_numpy()
    contract = dict(backend="python", representation=args.representation,archetype="quantitative grid", panels=1,
                    question="How are empirical stock-return distribution distances arranged within and between the fixed three sectors?",
                    claim="Descriptive heterogeneity only; no claim of equality, significance or calibrated uncertainty.",
                    sorting="Historical sector selection and pre-test index-weight order; no outcome-based clustering",
                    observations="311 common complete five-trading-day vectors per stock; 15 selected stocks, 2018-2024 test calendar",
                    uncertainty="Not shown: calibration validation is incomplete; 105 stock pairs are dependent",
                    scale="One pooled fitting-only standardization and Gaussian bandwidth for every pair",
                    export=dict(width_mm=162.56, height_mm=144.78, text_floor_pt=7, formats=["pdf","svg","png"]),
                    data="distance_matrix_d5.csv, stock_pair_distances_d5.csv, distance_summary_d5.json")
    (out/"figure_contract.json").write_text(json.dumps(contract, indent=2)+"\n", encoding="utf-8")
    plt.rcParams.update({"font.family":"Arial", "font.size":8, "pdf.fonttype":42,
                         "svg.fonttype":"none", "axes.linewidth":.6,
                         "xtick.major.width":.5, "ytick.major.width":.5})
    fig = plt.figure(figsize=(6.4,5.7))
    ax = fig.add_axes([.125,.16,.67,.752])
    cax = fig.add_axes([.84,.16,.024,.752])
    vmax = np.ceil(D.max()*20)/20
    image = ax.imshow(D, cmap="Blues", vmin=0, vmax=vmax, interpolation="nearest", aspect="equal")
    labels = panel.order_book_id.str.split(".").str[0].tolist()
    ax.set_xticks(np.arange(15), labels, rotation=60, ha="right", rotation_mode="anchor", fontsize=8)
    ax.set_yticks(np.arange(15), labels, fontsize=8)
    ax.tick_params(length=0, pad=4)
    for line in (4.5,9.5):
        ax.axhline(line, color="white", lw=1.6)
        ax.axvline(line, color="white", lw=1.6)
    for start in (0,5,10):
        ax.add_patch(Rectangle((start-.5,start-.5),5,5,fill=False,edgecolor="#526477",lw=.7))
    for center, name in zip((2,7,12),("Banks","Pharmaceuticals","Machinery")):
        ax.text(center,-1.03,name,ha="center",va="bottom",fontsize=8,clip_on=False)
    cb = fig.colorbar(image,cax=cax)
    cb.set_label("Empirical MMD (unsquared)",fontsize=9,labelpad=8)
    cb.outline.set_linewidth(.5)
    cb.ax.tick_params(labelsize=8,length=2,pad=3)
    fig.text(.125,.022,"Descriptive estimated-innovation distances only." if innovations else "Descriptive distances only; stock pairs are dependent.",fontsize=8,color="#444444")
    fig.canvas.draw()
    require_matplotlib_panel_alignment(fig,json_out=out/"stock_distance_map.alignment.json",
                                      exclude_axes=[cax],strict=True)
    stem = out/"stock_distance_map"
    fig.savefig(stem.with_suffix(".pdf"))
    fig.savefig(stem.with_suffix(".svg"))
    fig.savefig(stem.with_suffix(".png"),dpi=300)
    plt.close(fig)
    caption = ("Empirical distribution distances for the fixed 15-stock panel. Rows and columns follow the three "
               "historically selected sectors, with five stocks per sector; boxes mark within-sector comparisons. "
               "Each observation is a nonoverlapping five-trading-day vector of adjusted log returns. "
               "The test period contributes 311 common complete calendar blocks per stock. All pairs use the same "
               "pooled 2015–2017 fitting-only coordinate scale and median-distance Gaussian bandwidth. Color denotes "
               "the norm of the empirical kernel mean difference, including diagonal terms. The map is descriptive: "
               "finite-sample distances can be positive under equal distributions, pairs are dependent, and no "
               "statistical significance or equivalence is asserted. The analysis targets the selected stocks on "
               "jointly observed dates, rather than all stocks in each sector or missing suspension periods.\n")
    if innovations:
        caption=caption.replace("Each observation is a nonoverlapping five-trading-day vector of adjusted log returns.",
            "Each observation is a nonoverlapping five-trading-day vector of estimated standardized GARCH residuals. Each stock's mean and GARCH parameters are fitted only on 2015–2017, then held fixed during forward filtering.")
        caption=caption.replace("pooled 2015–2017 fitting-only coordinate scale", "pooled residual fitting-only coordinate scale")
        caption+="The representation and fitted kernel scale differ from the original-return map; distance changes are not a directly comparable improvement measure. Estimated residuals are not assumed to be known or independent innovations.\n"
    (out/"stock_distance_map_caption.txt").write_text(caption,encoding="utf-8")
    print("Exported fixed-order descriptive map (PDF, SVG, PNG); no inference annotations.")


if __name__ == "__main__":
    main()
