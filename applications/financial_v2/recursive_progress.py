"""Monitor completed checkpoints without conflating inner and outer repetitions."""
from pathlib import Path
import argparse
import json
from datetime import datetime,timezone
import pandas as pd
from summarize_validation import interval

ROOT=Path(__file__).resolve().parent


def main():
    p=argparse.ArgumentParser();p.add_argument('--B',type=int,default=39);a=p.parse_args()
    folder=ROOT/'recursive_innovation_pilot_d5'
    request=pd.read_csv(folder/f'run_request_B{a.B:04d}.csv').iloc[0]
    cases=request.scenarios.split(',');planned=int(request.replications_per_case);first=int(request.first_replication)
    records=[];tables=[]
    for case in cases:
        started=completed=inner=0
        for r in range(first,first+planned):
            cell=folder/case/f'r{r:05d}'
            if (cell/'observed.rds').exists():started+=1
            files=list(cell.glob('inner_b*.csv')) if cell.exists() else []
            if files:
                x=pd.concat([pd.read_csv(f) for f in files],ignore_index=True)
                if x.b.duplicated().any():raise ValueError('Overlapping inner checkpoints')
                inner+=int((x.b<=a.B).sum())
            final=cell/f'result_B{a.B:04d}.csv'
            if final.exists():
                t=pd.read_csv(final)
                assert len(t)==4 and not t.question.duplicated().any()
                tables.append(t);completed+=1
        records.append(dict(scenario=case,outer_planned=planned,outer_started=started,outer_complete=completed,
                            inner_complete=inner,inner_planned=planned*a.B))
    result=dict(updated_utc=datetime.now(timezone.utc).isoformat(),B=a.B,status='development pilot; not confirmed calibration',
                per_scenario=records,real_pvalues_computed=False,pvalue_resolution=1/(a.B+1),
                holm_global_null_has_sufficient_resolution=1/(a.B+1)<=.05/4,
                restricted_common_stock_null_not_general_mixture_null=True)
    (folder/f'progress_B{a.B:04d}.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    if tables:
        x=pd.concat(tables,ignore_index=True);x.to_csv(folder/f'completed_outer_B{a.B:04d}.csv',index=False)
        summary=[]
        for (case,q),g in x.groupby(['scenario','question']):
            for alpha in (.01,.05,.10):
                k=int((g.p<=alpha+1e-12).sum());lo,hi=interval(k,len(g))
                summary.append(dict(scenario=case,question=q,alpha=alpha,n=len(g),rejections=k,rate=k/len(g),ci_lower=lo,ci_upper=hi,
                                    level_below_pvalue_resolution=alpha<1/(a.B+1)))
        pd.DataFrame(summary).to_csv(folder/f'outer_rates_B{a.B:04d}.csv',index=False)
    text=['# 全流程重估校准开发试跑进度','',f'内层B={a.B}，外层按场景分别计数。内层伪数据不是独立实证重复次数。',
          '只根据完整落盘块计数；尚未完成的外层重复不解释为未拒绝。部分完成的外层样本也不用于认定校准已通过。','']
    for r in records:
        text.append(f"{r['scenario']}：外层已完成 {r['outer_complete']}/{r['outer_planned']}，已启动 {r['outer_started']}；内层已完成 {r['inner_complete']}/{r['inner_planned']}。")
    text+=['','B=39时最小p=0.025，不能用于四项Holm在全局原假设下的5%家族水平检查，也不能检查1%单项水平。',
           '本候选施加全部股票创新同分布，只作机制诊断；不代表已覆盖较弱的行业混合原假设。没有真实股票检验p值。']
    (folder/f'开发试跑进度_B{a.B:04d}.md').write_text('\n'.join(text)+'\n',encoding='utf-8')
    print(json.dumps(result))


if __name__=='__main__':main()
