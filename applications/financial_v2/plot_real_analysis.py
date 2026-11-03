"""Observed-data figures with fixed sampling and explicit inferential limits."""
from pathlib import Path
import sys
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parent
OUT=ROOT/"real_analysis/figures"
sys.path.insert(0,str(Path.home()/".codex/skills/nature-figure/scripts"))
from audit_panel_alignment import require_matplotlib_panel_alignment

plt.rcParams.update({"font.family":"Arial","font.size":8,"pdf.fonttype":42,"svg.fonttype":"none",
    "axes.spines.top":False,"axes.spines.right":False,"axes.linewidth":.65,"xtick.major.width":.6,
    "ytick.major.width":.6,"legend.frameon":False,"savefig.facecolor":"white"})
Q=["within_bank","within_pharma","within_machinery","between_sectors"]
QL=["Within banks","Within pharmaceuticals","Within machinery","Between sectors"]
METHODS=["MCA-kDET","DISCO","MMD","MST"]
COLORS=["#235A82","#9B713A","#4F8A80","#897197"]


def export(fig,name):
    fig.canvas.draw()
    require_matplotlib_panel_alignment(fig,json_out=OUT/f"{name}.alignment.json",strict=True)
    fig.savefig(OUT/f"{name}.pdf")
    fig.savefig(OUT/f"{name}.svg")
    fig.savefig(OUT/f"{name}.png",dpi=600)
    plt.close(fig)


def pvalues():
    x=pd.read_csv(ROOT/"real_analysis/results.csv")
    primary=x[(x.draw==1)&(x.stride==4)]
    assert np.all(x.p_holm_exploratory.to_numpy() > 0)
    assert min(primary.p_holm_exploratory.min(),x[(x.method=="MCA-kDET")&(x.stride==4)].p_holm_exploratory.min())>=.01
    fig,axes=plt.subplots(1,2,figsize=(7.2,3.75),sharey=True)
    fig.subplots_adjust(left=.235,right=.98,bottom=.23,top=.80,wspace=.15)
    for ax,title in zip(axes,("a  Fixed primary sample","b  MCA-kDET: all 10 draws")):
        ax.set_xscale("log");ax.set_xlim(.01,1.28);ax.set_ylim(3.5,-.5)
        ax.set_yticks(range(4),QL);ax.tick_params(axis="y",length=0,pad=7)
        ax.set_xticks([.01,.02,.05,.1,.5,1],[".01",".02",".05",".1",".5","1"])
        ax.axvline(.05,color="#B5A492",lw=.8,ls=(0,(3,3)),zorder=0)
        ax.grid(axis="x",which="major",lw=.4,color="#E7EBEE",zorder=0)
        ax.set_title(title,loc="left",fontsize=9,pad=10)
        ax.set_xlabel("Exploratory Holm-adjusted p",labelpad=7)
    for j,method in enumerate(METHODS):
        values=[float(primary[(primary.question==q)&(primary.method==method)].p_holm_exploratory.iloc[0]) for q in Q]
        axes[0].scatter(values,np.arange(4)+(j-1.5)*.13,s=23,marker=["o","s","^","D"][j],
                        color=COLORS[j],edgecolor="white",linewidth=.4,label=method,zorder=3)
    for i,q in enumerate(Q):
        g=x[(x.question==q)&(x.method=="MCA-kDET")&(x.stride==4)].sort_values("draw")
        axes[1].scatter(g.p_holm_exploratory,i+np.linspace(-.17,.17,10),s=19,color=COLORS[0],alpha=.7,edgecolor="white",lw=.3,zorder=3)
        first=g[g.draw==1]
        axes[1].scatter(first.p_holm_exploratory,[i-.17],s=40,facecolors="none",edgecolors="#203441",linewidth=.85,zorder=4)
    handles,labels=axes[0].get_legend_handles_labels()
    fig.legend(handles,labels,loc="upper center",bbox_to_anchor=(.60,.98),ncol=4,columnspacing=1.1,handletextpad=.4,fontsize=8)
    fig.text(.235,.085,"Dashed line: nominal .05 threshold. Open ring: prespecified first draw.",fontsize=7.5,color="#4B5C69")
    fig.text(.235,.033,"Exploratory only: residual calibration failed stress tests; draws share the same data.",fontsize=7.5,color="#4B5C69")
    export(fig,"real_method_comparison")
    caption=("Exploratory comparisons on the observed stock panel. Panel a retains the prespecified first draw and calendar stride four, "
      "using five-day estimated GARCH-residual vectors, five stocks per within-sector hypothesis, and three equal-weight stock-mixture sector populations. "
      "Each within-sector population has 15 test vectors and each sector mixture has 25. All four methods use the same vectors and 4,999 label permutations. "
      "Holm adjustment covers four hypotheses separately within each method; it does not cover all 16 method-hypothesis combinations jointly. "
      "Panel b shows all ten predeclared sampling/model-randomness draws for MCA-kDET at the same stride; these are reuses of one real dataset, not independent studies or power estimates. "
      "The open ring marks draw one. Points are not confidence intervals, p-values are not averaged, and no favorable draw is selected. "
      "The nominal .05 line is descriptive: stress tests found calibration failure, and comparator calibration in this financial setting is unverified. "
      "No threshold crossing establishes a confirmed discovery, and non-rejection does not establish equality. Source: results.csv and protocol.json.\n")
    (OUT/"real_method_comparison_caption.txt").write_text(caption,encoding="utf-8")


def distributions():
    panel=pd.read_csv(ROOT/"design_two_tests/main_panel.csv")
    calendar=pd.read_csv(ROOT/"design_two_tests/block_calendar.csv")
    valid=calendar[(calendar.dimension==5)&(calendar.split=="test")]
    blocks=valid.loc[valid.complete_all_main_stocks.astype(str).str.lower().eq("true"),"block"].tolist()
    rows=[];arrays=[]
    for label,path in (("Raw returns",ROOT/"design_two_tests/return_vectors.csv.gz"),
                       ("GARCH residuals",ROOT/"innovations/innovation_vectors.csv.gz")):
        data=pd.read_csv(path,float_precision="round_trip")
        data=data[(data.dimension==5)&(data.split=="test")]
        group=[]
        for stock in panel.order_book_id:
            a=data[data.stock==stock].set_index("block").loc[blocks,["r1","r2","r3","r4","r5"]].to_numpy().ravel()
            assert len(a)==1555 and np.isfinite(a).all()
            group.append(a)
            rows.extend(dict(representation=label,stock=stock,rank=i+1,value=value,cdf=(i+1)/len(a)) for i,value in enumerate(np.sort(a)))
        arrays.append(group)
    pd.DataFrame(rows).to_csv(OUT/"daily_coordinate_ecdf.csv.gz",index=False,compression="gzip")
    fig,axes=plt.subplots(2,3,figsize=(7.2,6.2),sharey=True)
    fig.subplots_adjust(left=.08,right=.98,bottom=.25,top=.935,hspace=.47,wspace=.20)
    palette=["#235A82","#B1743C","#4F8A80","#897197","#687786"]
    for row in range(2):
        lo=min(a.min() for a in arrays[row]);hi=max(a.max() for a in arrays[row])
        for col in range(3):
            ax=axes[row,col]
            for j in range(5):
                stock=col*5+j;a=np.sort(arrays[row][stock])
                ax.step(a,np.arange(1,len(a)+1)/len(a),where="post",lw=.9,color=palette[j],
                        label=panel.order_book_id.iloc[stock].split(".")[0])
            ax.set_xscale("symlog",linthresh=1);ax.set_xlim(lo,hi);ax.set_ylim(0,1.02)
            ax.set_xticks([-10,-2,-.5,0,.5,2,10],["−10","−2","−.5","0",".5","2","10"])
            ax.set_xlim(lo,hi)
            ax.set_yticks([0,.25,.5,.75,1]);ax.grid(axis="y",color="#E8ECEF",lw=.45)
            title=("Banks","Pharmaceuticals","Machinery")[col]
            ax.set_title(f"{'abcdef'[row*3+col]}  {title}",loc="left",fontsize=9,pad=8)
            ax.set_xlabel("Daily log return (%)" if row==0 else "Estimated standardized residual",labelpad=5,fontsize=8)
            if col==0:ax.set_ylabel("Empirical CDF",labelpad=7)
            if row==1:ax.legend(loc="upper center",bbox_to_anchor=(.5,-.27),ncol=2,fontsize=7.5,columnspacing=.8,handlelength=1.5)
    fig.text(.08,.015,"Daily coordinates from 311 common five-day blocks per stock; no confidence bands. Full tails retained on symlog axes.",fontsize=7.4,color="#4B5C69")
    export(fig,"stock_distribution_profiles")
    caption=("Within-sector empirical distribution profiles for the fixed stocks. Top panels show daily adjusted log-return coordinates in percent; "
      "bottom panels show estimated standardized GARCH residual coordinates. Five curves per panel correspond to historically selected stocks, with codes shown below. "
      "Each curve pools five coordinates from 311 jointly observed five-trading-day test blocks (1,555 daily values); these coordinates are dependent and are not treated as independent replications. "
      "The full tails are shown on symmetric logarithmic axes with a linear central region. No confidence bands or distributional-equality claims are attached. "
      "The inferential calculations compare full five-dimensional vectors, whereas these profiles describe pooled daily marginals only and cannot reveal every form of joint-distribution difference. "
      "GARCH means and parameters use only 2015–2017; test dates are 2018-01-30 through 2024 year end, conditional on common observability. "
      "Stocks and periods are unchanged after examining outcomes. Source: fixed return/innovation vectors and block calendar.\n")
    (OUT/"stock_distribution_profiles_caption.txt").write_text(caption,encoding="utf-8")


def main():
    OUT.mkdir(exist_ok=True)
    contract=dict(backend="python",archetype="quantitative grid",question="What do the observed stocks show, and how stable are exploratory method comparisons?",
        claim="Results are descriptive/exploratory; the prespecified MCA between-sector comparison does not cross .05 after Holm adjustment.",
        panels={"real_method_comparison":["prespecified four-method evidence","complete ten-draw MCA sensitivity"],
                "stock_distribution_profiles":["raw marginal heterogeneity in three sectors","GARCH-adjusted marginal profiles in the same three sectors"]},
        statistical_scope="No validated size or equivalence claim. Four hypotheses adjusted separately within method; draws are dependent uses of one dataset.",
        export={"width_mm":182.88,"minimum_text_pt":7.4,"formats":["pdf","svg","png"]})
    (OUT/"figure_contract.json").write_text(json.dumps(contract,indent=2)+"\n",encoding="utf-8")
    pvalues();distributions()
    print("Exported two observed-data figure sets with alignment manifests and captions.")


if __name__=="__main__":main()
