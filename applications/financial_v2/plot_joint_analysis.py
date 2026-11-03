"""Describe the new sector-level order-statistic observation, without inference."""
from pathlib import Path
import json
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(Path.home()/".codex/skills/nature-figure/scripts"))
from audit_panel_alignment import require_matplotlib_panel_alignment


def main():
    out=ROOT/"joint_analysis/figures";out.mkdir(exist_ok=True)
    source=pd.read_csv(ROOT/"joint_analysis/sorted_sector_vectors.csv")
    data=source[(source.split=="test")&source.complete_all_stocks.astype(str).str.lower().eq("true")]
    assert len(data)==311*3
    palette=["#235A82","#A87541","#4F8A80"];names=["Banks","Pharmaceuticals","Machinery"]
    plt.rcParams.update({"font.family":"Arial","font.size":8,"pdf.fonttype":42,"svg.fonttype":"none",
      "axes.spines.top":False,"axes.spines.right":False,"axes.linewidth":.65,"legend.frameon":False})
    fig,axes=plt.subplots(1,2,figsize=(7.2,3.9))
    fig.subplots_adjust(left=.09,right=.98,bottom=.24,top=.80,wspace=.30)
    for g in range(1,4):
        a=data[data.sector==g][[f"x{i}" for i in range(1,6)]].to_numpy()
        assert a.shape==(311,5)
        boxes=axes[0].boxplot(a,positions=np.arange(1,6)+(g-2)*.24,widths=.20,patch_artist=True,
          manage_ticks=False,showfliers=True,whiskerprops={"lw":.6,"color":palette[g-1]},
          capprops={"lw":.6,"color":palette[g-1]},medianprops={"lw":.9,"color":"#283B47"},
          flierprops={"marker":".","markersize":1.5,"alpha":.28,"markeredgecolor":palette[g-1]})
        for box in boxes["boxes"]:box.set(facecolor=palette[g-1],edgecolor=palette[g-1],alpha=.45,linewidth=.7)
        spread=np.sort(a[:,-1]-a[:,0]);assert (spread>=0).all()
        axes[1].step(spread,np.arange(1,312)/311,where="post",lw=1.1,color=palette[g-1])
    axes[0].set_xticks(range(1,6),["Min","2nd","Median","4th","Max"])
    axes[0].set_yscale("symlog",linthresh=1);axes[0].set_ylabel("Weekly residual score",labelpad=5)
    axes[0].set_xlabel("Ordered stock position",labelpad=6)
    axes[0].set_title("a  Five-coordinate sector profiles",loc="left",fontsize=9,pad=10)
    axes[1].set_xscale("symlog",linthresh=1);axes[1].set_ylim(0,1.02)
    axes[1].set_ylabel("Empirical CDF",labelpad=5);axes[1].set_xlabel("Within-week range (max − min)",labelpad=6)
    axes[1].set_title("b  Cross-stock dispersion",loc="left",fontsize=9,pad=10)
    for ax in axes:ax.grid(axis="y",lw=.45,color="#E9EDF0")
    fig.legend([Patch(facecolor=c,edgecolor=c,alpha=.6)for c in palette],names,loc="upper center",bbox_to_anchor=(.54,.98),ncol=3,fontsize=8)
    fig.text(.09,.09,"311 common calendar blocks per sector; descriptive dependent observations, not test sample sizes.",fontsize=7.4,color="#4B5C69")
    fig.text(.09,.035,"The primary test uses 25 vectors per sector. Sorting removes stock identity; no equality claim.",fontsize=7.4,color="#4B5C69")
    fig.canvas.draw();require_matplotlib_panel_alignment(fig,json_out=out/"sector_joint_profiles.alignment.json",strict=True)
    fig.savefig(out/"sector_joint_profiles.pdf");fig.savefig(out/"sector_joint_profiles.svg");fig.savefig(out/"sector_joint_profiles.png",dpi=600)
    plt.close(fig)
    contract=dict(backend="python",archetype="quantitative grid",claim="The order-statistic representation exposes selected-sector cross-sectional profiles and dispersion without averaging stocks.",
      panels=["all five ranked coordinate distributions","distribution of cross-sectional range"],
      observations="311 common complete test blocks per sector; dependent descriptive observations. Primary test samples 25 vectors per sector.",
      statistics="Boxplot median/IQR, whiskers 1.5 IQR, all fliers retained; full-range empirical CDF, no confidence bands or inference",
      scope="Not stock-identity-preserving joint laws; differences may reflect individual marginals and cross-stock dependence. Neither is isolated by this figure.",
      export=dict(width_mm=182.88,minimum_glyph_pt=5,formats=["pdf","svg","png"]))
    (out/"figure_contract.json").write_text(json.dumps(contract,indent=2)+"\n",encoding="utf-8")
    caption=("Sector-level joint-profile representation on the fixed stock panel. Each stock's weekly residual score is the sum of five daily estimated standardized GARCH residuals divided by the square root of five; the five stock scores within a sector are sorted ascending. "
      "Panel a displays each ordered coordinate over all 311 common observed test blocks per sector: center lines are medians, boxes IQRs, whiskers extend to 1.5 IQR, and all outliers are retained. Panel b displays the full empirical distribution of the maximum-minus-minimum stock score in each block. "
      "Symmetric logarithmic axes retain the full tails. The figure uses all jointly observable blocks for description; the primary test uses the fixed sparse sample of 25 vectors per sector, not 311 independent test observations. "
      "The representation removes stock identity and daily order within the weekly sum. Differences can arise from stock marginals or dependence; the plot does not isolate either mechanism, prove equality or establish calibrated significance. "
      "The new exploratory representation was proposed after the earlier mixture analysis. Stocks, periods and the new primary schedule were fixed before viewing new p-values.\n")
    (out/"sector_joint_profiles_caption.txt").write_text(caption,encoding="utf-8")
    print("Exported sector joint-profile figure; descriptive only.")


if __name__=="__main__":main()
