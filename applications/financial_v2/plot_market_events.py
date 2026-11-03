"""Paired display of fixed-module distribution screens and raw downside profiles."""
from pathlib import Path
import sys
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter

ROOT=Path(__file__).resolve().parent;OUT=ROOT/'market_events/figures'
sys.path.insert(0,str(Path.home()/'.codex/skills/nature-figure/scripts'))
from audit_panel_alignment import require_matplotlib_panel_alignment

NAMES={'银行':'Banks','非银行金融':'Nonbank finance','房地产':'Real estate','建筑':'Construction','食品饮料':'Food and beverage',
 '金融体系':'Financial system','科技产业':'Technology','防御消费':'Defensive consumer','工业周期':'Industrial cyclicals','能源材料':'Energy and materials'}


def main():
    OUT.mkdir(exist_ok=True)
    data=pd.read_csv(ROOT/'market_events/results.csv');risk=pd.read_csv(ROOT/'market_events/risk_diagnostics.csv')
    modules=json.loads((ROOT/'market_events/modules.json').read_text(encoding='utf-8'))
    downside_vmax=max(.5,float(np.ceil(risk.joint_downside_frequency.max()*10)/10))
    plt.rcParams.update({'font.family':'Arial','font.size':8,'pdf.fonttype':42,'svg.fonttype':'none','axes.linewidth':.6})
    fig,axes=plt.subplots(3,2,figsize=(7.2,7.2))
    fig.subplots_adjust(left=.205,right=.975,bottom=.205,top=.925,hspace=.42,wspace=.30)
    cases=['China2015','China2024','US2020'];titles=['China: 2015','China: 9/24 policy','US: 2020 policy']
    for row,(case,title)in enumerate(zip(cases,titles)):
        names=[m['module']for m in modules if m['event']==case]
        g=data[data.event==case].pivot(index='module',columns='method',values='p_holm_exploratory').loc[names,['MCA-kDET','DISCO','MMD','MST']]
        P=g.to_numpy();assert np.all(P>0)
        a=axes[row,0];im=a.imshow(-np.log10(P),cmap='Blues',vmin=0,vmax=3,aspect='auto',interpolation='nearest')
        a.set_xticks(range(4),['MCA','DISCO','MMD','MST']);a.set_yticks(range(5),[NAMES[n]for n in names]);a.tick_params(length=0,pad=5)
        a.set_title(f"{'ace'[row]}  {title}: distribution screen",loc='left',fontsize=8.5,pad=9)
        for j in range(5):
            for k in range(4):
                if P[j,k]<.05:a.text(k,j,f'{P[j,k]:.3f}',ha='center',va='center',fontsize=7.5,color='white'if P[j,k]<.01 else '#173A53')
        R=risk[risk.event==case].pivot(index='module',columns='phase',values='joint_downside_frequency').loc[names,[1,2,3]].to_numpy()
        b=axes[row,1];im2=b.imshow(R,cmap='YlOrBr',vmin=0,vmax=downside_vmax,aspect='auto',interpolation='nearest')
        b.set_xticks(range(3),['Phase 1','Phase 2','Phase 3']);b.set_yticks(range(5),[]);b.tick_params(length=0,pad=5)
        b.set_title(f"{'bdf'[row]}  {title}: raw downside",loc='left',fontsize=8.5,pad=9)
        for j in range(5):
            for k in range(3):b.text(k,j,f'{100*R[j,k]:.0f}%',ha='center',va='center',fontsize=8,color='white'if R[j,k]>.33 else '#513B24')
        for ax in (a,b):
            ax.set_xticks(np.arange(-.5,len(ax.get_xticks())-.5,1),minor=True)
            ax.set_yticks(np.arange(-.5,5,1),minor=True);ax.grid(which='minor',color='white',lw=.8);ax.tick_params(which='minor',length=0)
    c1=fig.add_axes([.205,.11,.335,.018]);c2=fig.add_axes([.64,.11,.335,.018])
    cb=fig.colorbar(im,cax=c1,orientation='horizontal');cb.set_ticks([0,1,-np.log10(.05),2,3],labels=['1','.1','.05','.01','.001'])
    cb.set_label('Exploratory Holm p (15 modules per method)',fontsize=7.5,labelpad=3)
    cb2=fig.colorbar(im2,cax=c2,orientation='horizontal');cb2.ax.xaxis.set_major_formatter(PercentFormatter(1));cb2.set_label('Days with at least half below baseline left tails',fontsize=7.5,labelpad=3)
    for cb in (cb,cb2):cb.ax.tick_params(labelsize=7.5,length=2,pad=2);cb.outline.set_linewidth(.5)
    fig.text(.205,.041,'Left: filtered-vector tests. Right: all observed phase days, descriptive only; module dimensions differ.',fontsize=7.3,color='#4B5C69')
    fig.text(.205,.016,'No confirmed discoveries, causal policy effects, or resilience rankings. Exact windows are in the caption.',fontsize=7.3,color='#4B5C69')
    fig.canvas.draw();require_matplotlib_panel_alignment(fig,json_out=OUT/'event_module_screen.alignment.json',exclude_axes=[c1,c2],strict=True)
    fig.savefig(OUT/'event_module_screen.pdf');fig.savefig(OUT/'event_module_screen.svg');fig.savefig(OUT/'event_module_screen.png',dpi=600)
    plt.close(fig)
    contract=dict(backend='python',archetype='quantitative grid',claim='Distributional screening and observed downside changes answer different questions; neither alone establishes resilience.',
      panels=['five fixed modules by four tests','same-module downside by three stages'],
      statistics='Holm across 15 planned modules separately within method; downside frequency is descriptive, no uncertainty band',
      scale=f'Left minus-log10 p displayed with original p colorbar; right same 0-{100*downside_vmax:.0f}% range, no clipping of observed values',
      scope='Different module dimensions and country constructions; compare phases within a module, not resilience ranks across modules',
      export=dict(width_mm=182.88,minimum_text_pt=7.3,formats=['pdf','svg','png']))
    (OUT/'figure_contract.json').write_text(json.dumps(contract,indent=2)+'\n',encoding='utf-8')
    caption=('Fixed event-module distribution screens and descriptive downside profiles. Each row uses five economically specified modules; Chinese stock baskets have 5, 5, 3, 3 and 2 coordinates, while US modules contain 4, 4, 5, 5 and 5 value-weighted industry portfolios. '
      'Left panels show exploratory Holm-adjusted p-values across all 15 modules separately within each method, using the same test vectors and 4,999 label permutations. Numeric annotations are limited to p below nominal .05; this is not confirmed size control. '
      'Right panels use all jointly observed phase days and show the fraction on which at least half of coordinates fall below their own pre-event raw-return fifth percentiles. For the two-coordinate Chinese food module this includes one-coordinate breaches; it is not a pure co-crash or tail-dependence coefficient and is not used to rank modules of different dimensions. '
      'China 2015 phases: Jan 5-Mar 31, Apr 1-Jun 12, and Jun 15-Sep 30. China 2024: Mar 28-Jun 27, Jun 28-Sep 23, and Sep 24-Dec 23. US 2020: Sep 30-Dec 23 2019, Dec 24 2019-Mar 20 2020, and Mar 23-Jun 16 2020. '
      'Classifiers use the first 40% of each original phase calendar, followed by a five-day buffer and held-out tests; fit and test masks are applied without compressing calendar gaps. Daily GARCH parameters use only an earlier baseline and are shared across stages. '
      'These are historical exploratory event windows, not causal treatment assignments or point-in-time investment backtests. Existing calibration limitations remain; non-rejection cannot demonstrate safety, stability or equality. No confidence intervals are shown for the descriptive frequencies.\n')
    (OUT/'event_module_screen_caption.txt').write_text(caption,encoding='utf-8')
    print('Exported paired event-screen and downside-profile figure.')


if __name__=='__main__':main()
