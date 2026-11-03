"""Assemble all fixed event modules and attach descriptive risk consequences."""
from pathlib import Path
import sys
import gzip
import hashlib
import json
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parent;OUT=ROOT/'market_events'
sys.path.insert(0,str(ROOT.parents[1]/'tests/out/financial_rdata_runtime'))
import rdata

LABELS={'China2015':'A股2015','China2024':'A股9·24','US2020':'美股2020'}
METHODS=['MCA-kDET','DISCO','MMD','MST']


def holm(p):
    p=np.asarray(p);order=np.argsort(p);result=np.empty(len(p))
    result[order]=np.minimum(1,np.maximum.accumulate(p[order]*np.arange(len(p),0,-1)))
    return result


def main():
    paths=sorted((OUT/'jobs').glob('*.rds'));assert len(paths)==15
    fp=rdata.read_rds(OUT/'fingerprint.rds')
    assert [hashlib.md5(Path(p).read_bytes()).hexdigest()for p in fp['paths']]==list(fp['hashes'])
    calendar=pd.read_csv(OUT/'phase_calendar.csv',parse_dates=['date'])
    rows=[];risks=[];filter_checks=[];curves=[];quality=[]
    for path in paths:
        with gzip.open(path,'rb')as f:
            while f.read(1024*1024):pass
        a=rdata.read_rds(path);r=a['rows'].copy();event=r.event.iloc[0];module=r.module.iloc[0]
        assert len(r)==4 and set(r.method)==set(METHODS)
        for method,v in a['audits']['1']['statistics'].items():
            p=(1+np.count_nonzero(v[1:]>=v[0]-1e-10))/len(v)
            assert abs(p-r.loc[r.method==method,'p_exploratory'].iloc[0])<1e-14
        rows.append(r)
        raw=np.ma.asarray(a['raw'],dtype=float).filled(np.nan)
        Z=np.ma.asarray(a['Z'],dtype=float).filled(np.nan)
        H=np.ma.asarray(a['H'],dtype=float).filled(np.nan)
        par=np.ma.asarray(a['parameters'],dtype=float).filled(np.nan);base=np.asarray(a['baseline'])
        assert not np.isinf(raw).any() and not np.isinf(Z).any()
        # Canonicalize R's signalling NA payload before NumPy arithmetic, retaining the same missingness mask.
        raw[np.isnan(raw)]=np.nan;Z[np.isnan(Z)]=np.nan
        dates=pd.DatetimeIndex(a['dates']);d=raw.shape[1]
        spec=json.loads((OUT/'modules.json').read_text(encoding='utf-8'))[int(a['q'][0])-1]
        for t,j in np.argwhere(np.abs(raw)>30):
            quality.append(dict(event=event,module=module,date=str(dates[t].date()),asset=spec['assets'][j],
                log_return_percent=float(raw[t,j]),in_filter_baseline=bool(base[t]),
                action='Preserved; extreme historical return requires source/corporate-action review, not outcome-based trimming'))
        assert par.shape==(d,11) and np.all(par[:,6]==0) and np.isfinite(par).all() and np.isfinite(H).all()
        mu=np.nanmean(raw[base],axis=0);v0=np.nanmean((raw[base]-mu)**2,axis=0)
        np.testing.assert_allclose(mu,par[:,0],atol=2e-12)
        np.testing.assert_allclose(v0,par[:,5],atol=2e-12)
        h=v0.copy();last2=v0.copy();expected=np.empty_like(H)
        for t in range(len(raw)):
            h=par[:,1]+par[:,2]*last2+par[:,3]*h
            expected[t]=h
            centered=raw[t]-mu;last2=np.where(np.isfinite(centered),centered**2,h)
        err=float(np.max(np.abs(expected-H)/np.maximum(1,H)))
        assert np.isfinite(err) and err<2e-10
        np.testing.assert_allclose((raw-mu)/np.sqrt(H),Z,atol=1e-11,equal_nan=True)
        filter_checks.append(dict(event=event,module=module,max_variance_relative_error=err,
            near_boundary=int(par[:,-1].sum()),failed_optimizer_starts=int(par[:,-2].sum())))
        threshold=np.nanquantile(raw[base],.05,axis=0)
        for phase in range(1,4):
            c=calendar[(calendar.event==event)&(calendar.phase==phase)]
            ids=dates.get_indexer(c.date);assert (ids>=0).all()
            complete=np.isfinite(raw[ids]).all(axis=1)
            X=raw[ids][complete];n=len(X)
            coprob=float(((X<threshold).sum(axis=1)>=int(np.ceil(d/2))).mean())
            for part,key in [('classifier_fit','fitting_dates'),('test','testing_dates')]:
                expect=c[(c.split==part)&complete].date.dt.strftime('%Y-%m-%d').tolist()
                assert expect==list(a[key][phase-1])
            # One-day equal-weight portfolio simple return; tested coordinates remain individual log returns.
            score=100*np.expm1(X/100).mean(axis=1)
            worst=max(1,int(np.ceil(.10*n)));loss=-float(np.sort(score)[:worst].mean())
            risks.append(dict(event=event,module=module,phase=phase,label=c.label.iloc[0],n_observed=n,dimension=d,
              joint_downside_frequency=coprob,portfolio_worst10pct_mean_loss=loss,tail_observations=worst,
              portfolio_mean=float(score.mean()),portfolio_sd=float(score.std(ddof=1))))
            # Coordinate-level empirical CDF remains descriptive; no sector averaging enters it.
            for j in range(d):
                curves.extend(dict(event=event,module=module,phase=phase,coordinate=j+1,value=v,cdf=(i+1)/n)
                  for i,v in enumerate(np.sort(X[:,j])))
    x=pd.concat(rows,ignore_index=True)
    assert len(x)==60 and not x.duplicated(['event','module','method']).any()
    x['p_holm_exploratory']=np.nan;x['p_holm_event_exploratory']=np.nan
    for method,g in x.groupby('method'):x.loc[g.index,'p_holm_exploratory']=holm(g.p_exploratory.to_numpy())
    for _,g in x.groupby(['event','method']):x.loc[g.index,'p_holm_event_exploratory']=holm(g.p_exploratory.to_numpy())
    x.to_csv(OUT/'results.csv',index=False)
    pd.DataFrame(risks).to_csv(OUT/'risk_diagnostics.csv',index=False)
    pd.DataFrame(filter_checks).to_csv(OUT/'filter_audit.csv',index=False)
    pd.DataFrame(quality,columns=['event','module','date','asset','log_return_percent','in_filter_baseline','action']).to_csv(OUT/'data_quality_flags.csv',index=False)
    pd.DataFrame(curves).to_csv(OUT/'raw_coordinate_ecdf.csv.gz',index=False,compression='gzip')
    m=x[x.method=='MCA-kDET'];significant=m[m.p_holm_exploratory<.05]
    text=['# 中外市场事件：风险模块分布筛查（探索性）','',
      '本轮回答同一模块跨阶段是否变化，不再把不同行业间存在差异作为主要发现。三个案例、每例五个固定模块，四方法共60个p值，均完整保留。',
      'A股2015按用户选择的上涨—高位见顶—回撤三阶段；2024以9月24日政策发布日为锚，分较早前期/临近政策前/政策后三个60交易日窗；美股2020以3月23日政策公告为锚，同样取两个前期窗和一个后期窗。',
      '2015的峰值是事后历史分期，两个政策事件也不是随机处理，均不据此前后对比作因果归因。2024不直接定义为永久牛熊转换。',
      '观测为固定资产坐标的日收益向量；主检验为每坐标一个共同的事件前拟合GARCH过滤器得到的残差向量。滤波参数不按阶段分别估计。',
      '每阶段原日历前40%用于分类器拟合，中间五个交易日缓冲，其余留出测试；完整坐标掩码在分割后应用，缺口不压缩。正式计算前检查全部组的拟合/测试完整样本均至少12条。',
      '四方法使用同一测试向量与4999次置换标签。主调整为每方法跨全部15个事件模块的Holm；事件内五项调整仅作并列敏感性，不替代主调整。',
      '中国使用2014年事前固定历史股票篮子；美国使用最新版本公开行业组合，组合成五个事前经济模块。两国模块不一一对应，不做中美名次或收益能力比较。',
      '这是依赖时间序列上的探索比较；已有校准失败记录未删除，新事件表示未获得精确水平保证。p较小不代表损害更大，未拒绝不等于风险稳定或抗跌。','',
      '## 主MCA结果','']
    for r in m.itertuples():
        text.append(f'- {LABELS[r.event]} / {r.module}：d={r.d}，各阶段测试数{r.test_n}，原始探索p={r.p_exploratory:.4f}，全15项Holm={r.p_holm_exploratory:.4f}。')
    text.extend(['','## 风险后果的描述证据','',
      '共同下跌频率：各坐标使用事件前拟合期原收益5%分位数为固定阈值，计算至少一半坐标同时低于各自阈值的日比例。',
      '注意：食品饮料只有两个坐标，至少一半包括一个坐标越界；该指标不等于纯共跌/尾部依赖系数。模块维数不同，只用于同一模块跨阶段比较，不按其绝对比例跨模块排名。',
      '另列等权组合最差10%观察日的平均损失，用于经济解释，明确不同于用于检验的完整向量；该项才对坐标求组合平均，不用它替代分布检验。',
      '这些指标使用各阶段全部共同可观测日、没有置信区间，是描述性后果，不是筛选后独立验证；尾部日数小，不能据此给出稳健投资建议或长期韧性排名。',
      '分布检验未拒绝的模块也完整显示风险指标。是否“经得起考验”，还需要独立时期的风险预测和实质差异评估。',''])
    for row in risks:
        if row['phase']==3:
            first=next(r for r in risks if r['event']==row['event']and r['module']==row['module']and r['phase']==1)
            middle=next(r for r in risks if r['event']==row['event']and r['module']==row['module']and r['phase']==2)
            text.append(f"- {LABELS[row['event']]} / {row['module']}：三阶段共同下跌比例{first['joint_downside_frequency']:.3f}/{middle['joint_downside_frequency']:.3f}/{row['joint_downside_frequency']:.3f}；最差10%日平均组合损失{first['portfolio_worst10pct_mean_loss']:.3f}%/{middle['portfolio_worst10pct_mean_loss']:.3f}%/{row['portfolio_worst10pct_mean_loss']:.3f}%，第三阶段仅{row['tail_observations']}个尾部观察日。")
    text.extend(['','## 当前边界与完整数据','',
      'results.csv保留四方法原始及两种调整数值；risk_diagnostics.csv保留所有阶段，包括与预期不一致的结果；feasibility.csv和phase_calendar.csv给出样本量及实际窗口。',
      '不存在按新p值更换股票、窗口、阈值或模块。未修改论文草稿和冻结模型/置换引擎，也未追加大型模拟。'])
    if quality:text.extend(['',f'数据质量记录：有{len(quality)}条绝对日对数收益超过30%的历史记录，原样保留，详见data_quality_flags.csv；需核实历史复牌/公司行动等来源，不自动删去或缩尾。'])
    (OUT/'事件风险筛查结果.md').write_text('\n'.join(text)+'\n',encoding='utf-8')
    result=dict(status='EXPLORATORY_EVENT_SCREEN_UNVALIDATED',event_modules=15,method_results=60,
      all_raw_pvalues_recomputed=True,source_fingerprints_match=True,all_phase_splits_verified=True,
      garch_recursion_independently_verified=True,max_variance_relative_error=max(r['max_variance_relative_error']for r in filter_checks),
      extreme_historical_return_flags=len(quality),
      MCA_global_Holm_nominal_crossings=len(significant),no_confirmatory_or_causal_claim=True)
    (OUT/'report_verification.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(m[['event','module','p_exploratory','p_holm_exploratory','test_n']].to_string(index=False));print(json.dumps(result))


if __name__=='__main__':main()
