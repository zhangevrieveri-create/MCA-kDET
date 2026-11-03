"""Descriptive ECDFs: all five prespecified daily features, all fixed stocks."""
from pathlib import Path
import sys
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

ROOT=Path(__file__).resolve().parent;OUT=ROOT/'daily_features_v1/figures'
sys.path.insert(0,str(Path.home()/'.codex/skills/nature-figure/scripts'))
from audit_panel_alignment import require_matplotlib_panel_alignment
from daily_feature_trial import FEATURES


def main():
    OUT.mkdir(exist_ok=True)
    x=pd.read_csv(ROOT/'daily_features_v1/complete_daily_features.csv.gz')
    x=x[x.split.eq('test')]
    panel=pd.read_csv(ROOT/'design_two_tests/main_panel.csv')
    sectors=['银行','医药','机械'];names=['Bank stocks','Healthcare stocks','Machinery stocks']
    titles=['Overnight change','Intraday change','High/low range','Close position','Turnover']
    units=['Log change (%)','Log change (%)','Log ratio (%)','Position (0 to 1)','Vendor units']
    colors=['#0072B2','#D89000','#009E73','#B46DA0','#444444']
    plt.rcParams.update({'font.family':'Arial','font.size':8.5,'pdf.fonttype':42,'svg.fonttype':'none',
        'axes.linewidth':.6,'axes.spines.top':False,'axes.spines.right':False,'legend.frameon':False})
    fig,axes=plt.subplots(3,5,figsize=(10.8,7.5))
    fig.subplots_adjust(left=.075,right=.985,bottom=.16,top=.91,hspace=.94,wspace=.42)
    for row,(sector,name) in enumerate(zip(sectors,names)):
        stocks=panel[panel.first_industry_name.eq(sector)].order_book_id.tolist()
        handles=[];labels=[]
        for col,f in enumerate(FEATURES):
            ax=axes[row,col];g=x[x.sector.eq(sector)]
            for j,stock in enumerate(stocks):
                values=np.sort(g[g.stock.eq(stock)][f].to_numpy());assert np.isfinite(values).all()
                line,=ax.step(values,np.arange(1,len(values)+1)/len(values),where='post',color=colors[j],lw=1.05)
                if col==0:handles.append(line);labels.append(f'{stock.split(".")[0]} (n={len(values)})')
            ax.set_ylim(0,1.02);ax.set_yticks([0,.5,1]);ax.xaxis.set_major_locator(MaxNLocator(3))
            if col==3:ax.set_xlim(-.03,1.03)
            ax.set_xlabel(units[col],fontsize=8)
            ax.tick_params(labelsize=7.7,length=2.5)
            ax.set_title(f'{chr(97+5*row+col)}  {titles[col]}',loc='left',fontsize=8.5,pad=6)
            if col==0:ax.set_ylabel(name+'\nCumulative fraction',fontsize=8.2)
            else:ax.set_yticklabels([])
        fig.legend(handles,labels,ncol=5,loc='upper left',bbox_to_anchor=(.085,[.685,.395,.105][row]),
                   fontsize=8,columnspacing=1.8,handlelength=1.7,borderaxespad=0)
    fig.suptitle('Daily trading characteristics within fixed sectors | 2024 held-out observations',x=.075,y=.972,ha='left',fontsize=11)
    fig.text(.075,.035,'Source: RQData. Each line is one stock; all valid test days retained. Five features; no GARCH or five-day aggregation.',fontsize=8,color='#4B5C69')
    fig.text(.075,.014,'Descriptive empirical distributions only. Serial and same-day cross-stock dependence remain; no significance or equivalence claim.',fontsize=8,color='#4B5C69')
    fig.canvas.draw()
    require_matplotlib_panel_alignment(fig,json_out=str(OUT/'daily_feature_ecdf.alignment.json'),
        overlay_svg=str(OUT/'daily_feature_ecdf.alignment.svg'),tolerance_pt=1.5,gutter_tolerance_pt=1.5,strict=True)
    fig.savefig(OUT/'daily_feature_ecdf.pdf')
    fig.savefig(OUT/'daily_feature_ecdf.svg')
    fig.savefig(OUT/'daily_feature_ecdf.png',dpi=600)
    contract=dict(question='What differs in daily price shape and activity across the five fixed stocks within each sector?',
        core_claim='Descriptive marginal shapes differ across some stocks; all five daily features remain visible, without significance claims.',
        archetype='quantitative grid',backend='python',evidence='Three sectors by all five prespecified features; each stock identity and actual held-out n labelled.',
        observations=len(x),excluded_from_this_figure='Fitting/buffer days are outside the held-out descriptive target. All valid test observations plotted; no trimming.',
        uncertainty='No inferential bands: serial/cross-stock dependence uncalibrated; these are empirical CDFs, not confidence intervals.',
        scope='2024 daily feature trial; fixed historical cohort, not whole-sector population or event/country comparison.')
    (OUT/'figure_contract.json').write_text(json.dumps(contract,indent=2)+'\n',encoding='utf-8')
    print('Exported descriptive 15-panel daily feature ECDFs.')


if __name__=='__main__':main()
