"""Independently verify the daily trial and report predictive/descriptive evidence."""
from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
import sys
import hashlib
import json
import numpy as np
import pandas as pd
from scipy.stats import rankdata
from scipy.spatial.distance import cdist
from scipy.sparse.csgraph import minimum_spanning_tree

ROOT=Path(__file__).resolve().parent;OUT=ROOT/'daily_features_v1'
sys.path.insert(0,str(ROOT.parents[1]/'tests/out/financial_rdata_runtime'))
import rdata
from daily_feature_trial import FEATURES, sha


def holm(p):
    order=np.argsort(p);ans=np.empty(len(p));ans[order]=np.minimum(1,np.maximum.accumulate(np.asarray(p)[order]*np.arange(len(p),0,-1)));return ans


def prim(D):
    inside=np.zeros(len(D),dtype=bool);inside[0]=True;key=D[0].copy();key[0]=np.inf
    parent=np.zeros(len(D),dtype=int);edges=[]
    for _ in range(len(D)-1):
        key[inside]=np.inf;v=int(np.argmin(key));edges.append([parent[v],v]);inside[v]=True
        upd=(~inside)&(D[v]<key);key[upd]=D[v,upd];parent[upd]=v
    return np.asarray(edges)


def main():
    fp=rdata.read_rds(OUT/'analysis/fingerprint.rds')
    assert [hashlib.md5(Path(p).read_bytes()).hexdigest() for p in fp['paths']]==list(fp['hashes'])
    protocol=json.loads((OUT/'protocol.json').read_text(encoding='utf-8'))
    assert all(sha(ROOT.parents[1]/p)==h for p,h in protocol['frozen_sources'].items())
    features=pd.read_csv(OUT/'complete_daily_features.csv.gz');panel=pd.DataFrame(protocol['cohort'])
    methods=['MCA-kDET','DISCO','MMD','MST'];rows=[];aucs=[];checks=[];saved={}
    paths=sorted((OUT/'analysis/jobs').glob('*.rds'));assert len(paths)==6
    for path in paths:
        a=rdata.read_rds(path);r=a['rows'];sector=r.sector.iloc[0];rep=r.representation.iloc[0];rec=a['records']
        assert len(r)==4 and not rec.duplicated(['date','stock']).any()
        test=features[(features.sector==sector)&(features.split=='test')].set_index(['stock','date'])
        raw=test.loc[list(zip(rec.stock,rec.date)),list(a['features'])].to_numpy()
        fitting=features[(features.sector==sector)&(features.split=='fit')][list(a['features'])].to_numpy()
        st=a['coordinate_scale'];mu=np.asarray(st['mu']);sd=np.asarray(st['sd']);X=np.asarray(a['X'])
        np.testing.assert_allclose(fitting.mean(axis=0),mu,atol=1e-13)
        np.testing.assert_allclose(fitting.std(axis=0,ddof=1),sd,atol=1e-13)
        np.testing.assert_allclose((raw-mu)/sd,X,atol=1e-12)
        S=np.asarray(a['S']);L=np.asarray(a['Lall']);targets=np.asarray(a['targets']);N=len(X)
        np.testing.assert_array_equal(L[:,0],rec.population.to_numpy())
        n=np.array([(L[:,0]==g).sum() for g in range(1,6)])
        for g,count in enumerate(n,1):assert np.all((L==g).sum(axis=0)==count)
        ranks=rankdata(S,axis=0,method='average');T=np.empty_like(np.asarray(a['Tm']))
        D=cdist(X,X);D2=D**2;h2=float(np.asarray(a['h2'])[0]);K=np.exp(-D2/(2*h2))
        np.testing.assert_allclose(np.median(D2[np.triu_indices(N,1)]),h2,atol=1e-12)
        disco=np.zeros(L.shape[1]);within=[];products=[];I=[]
        for i in range(5):
            ind=(L==i+1).astype(float);I.append(ind);cols=np.flatnonzero(targets==i+1)
            auc=(ranks[:,cols].T@ind-n[i]*(n[i]+1)/2)/(n[i]*(N-n[i]))
            T[cols]=(auc-.5)/np.sqrt((N+1)/(12*n[i]*(N-n[i])))
            disco-=(ind*(D@ind)).sum(axis=0)/n[i]
            products.append(K@ind)
            within.append(((ind*products[-1]).sum(axis=0)-n[i])/(n[i]*(n[i]-1)))
        mmd=sum(within[i]+within[j]-2*(I[i]*products[j]).sum(axis=0)/(n[i]*n[j]) for i in range(5) for j in range(i+1,5))/10
        tree=minimum_spanning_tree(D).tocoo();assert tree.nnz==N-1
        mst=(L[tree.row]==L[tree.col]).sum(axis=0)
        scipy_mst_difference=float(np.max(np.abs(mst-np.asarray(a['statistics']['MST']))))
        if scipy_mst_difference:
            edges=prim(D);mst=(L[edges[:,0]]==L[edges[:,1]]).sum(axis=0)
        v={'MCA-kDET':T.max(axis=0),'DISCO':disco,'MMD':mmd,'MST':mst}
        errors={m:float(np.max(np.abs(v[m]-a['statistics'][m]))) for m in methods}
        assert max(errors.values())<1e-8,errors
        for m in methods:
            p=(1+np.count_nonzero(v[m][1:]>=v[m][0]-1e-10))/len(v[m])
            assert abs(p-r.loc[r.method==m,'iid_reference_p_unvalidated'].iloc[0])<1e-14
        b=a['auc'].copy();b['sector']=sector;b['representation']=rep
        for j,t in enumerate(targets):
            mi=n[t-1];u=(ranks[L[:,0]==t,j].sum()-mi*(mi+1)/2)/(mi*(N-mi))
            assert abs(u-b.auc.iloc[j])<1e-14
        rows.append(r);aucs.append(b);saved[(sector,rep)]=L
        checks.append(dict(sector=sector,representation=rep,all_1000_labelings_recomputed=True,
            max_statistic_errors=errors,scipy_MST_tie_difference=scipy_mst_difference,
            MST_engine='SciPy sparse MST' if not scipy_mst_difference else 'Literal Prim tie policy matching frozen R',
            warnings=list(a['warnings'])))
        print(path.name+' verified',flush=True)
    for sector in ['银行','医药','机械']:np.testing.assert_array_equal(saved[(sector,'full5')],saved[(sector,'price4')])
    result=pd.concat(rows,ignore_index=True);result['holm6_iid_reference_unvalidated']=np.nan
    for _,g in result.groupby('method'):result.loc[g.index,'holm6_iid_reference_unvalidated']=holm(g.iid_reference_p_unvalidated.to_numpy())
    result.to_csv(OUT/'reference_results.csv',index=False)
    A=pd.concat(aucs,ignore_index=True);A.to_csv(OUT/'stock_learner_auc.csv',index=False)
    macro=A.groupby(['sector','representation','learner']).auc.mean().reset_index(name='macro_auc')
    macro.to_csv(OUT/'macro_auc.csv',index=False)
    all_features=pd.read_csv(OUT/'daily_features.csv.gz');dependence=[];correlations=[]
    for sector,subset in all_features[all_features.split.eq('test')].groupby('sector'):
        for f in FEATURES:
            pv=subset.assign(value=subset[f].where(subset.valid)).pivot(index='date',columns='stock',values='value')
            for stock in pv.columns:
                dependence.append(dict(sector=sector,stock=stock,feature=f,lag1_correlation=float(pv[stock].autocorr())))
            cc=pv.corr();vals=cc.to_numpy()[np.triu_indices(len(cc),1)]
            correlations.append(dict(sector=sector,feature=f,median_absolute_same_day_correlation=float(np.median(abs(vals))),max_absolute_same_day_correlation=float(np.max(abs(vals)))))
    pd.DataFrame(dependence).to_csv(OUT/'temporal_dependence.csv',index=False)
    pd.DataFrame(correlations).to_csv(OUT/'same_day_dependence.csv',index=False)
    audit=dict(source_hashes_match=True,frozen_sources_unchanged=True,independent_reference_pvalues=24,
        independent_statistic_labelings=6000,feature_matrices_rebuilt_from_daily_values=True,
        all_fit_only_scales_verified=True,ablation_uses_identical_permutation_labels=True,
        no_claim_of_real_data_calibration=True,checks=checks)
    (OUT/'verification.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    lines=['# 日频五特征：行业内小试','',
      '一只股票为一个总体，一天为一条五维观测。固定银行、医药、机械各五股历史名单；2023拟合，2024前20交易日缓冲，其余逐日测试。每股拟合242日，测试220–222日。主表示不用GARCH，不构造五日向量，不跨股票平均，不稀疏抽样。',
      '', '五个坐标为隔夜对数涨跌、日内对数涨跌、对数高低价振幅、收盘区间位置、当日换手率。价格采用统一后复权，昨收从同一复权序列的上一市场交易日取得，不使用不复权prev_close混算；停牌及相邻停牌造成不可比的记录保留标记而不补值。高低相等的有成交日收盘位置定为0.5。换手率保留数据商today原单位，不擅自改标为小数比例。',
      '', '## 留出区分能力（描述性）','',
      '| 行业 | 学习器 | 五维宏平均AUC | 仅四个价格特征AUC | 差值 |','|---|---|---:|---:|---:|']
    for sector in ['银行','医药','机械']:
        for learner in ['l1','quad','knn','rf']:
            g=macro[(macro.sector==sector)&(macro.learner==learner)].set_index('representation').macro_auc
            lines.append(f'| {sector} | {learner} | {g["full5"]:.4f} | {g["price4"]:.4f} | {g["full5"]-g["price4"]:+.4f} |')
    lines+=['','AUC为每只股票对其余股票的留出AUC再对五股等权平均，四种学习器完整展示；不是概率同分布检验结论。全五维与价格四维分别按预先固定种子拟合，AUC差也可能包含学习器随机性，没有做提升显著性检验。',
      '', '## 置换参考数值的边界','',
      '六个配置（3行业×2表示）各四方法、B999，全部24个普通标签置换参考数值保留在reference_results.csv，另按每方法六配置作Holm的算术参考。三行业五维下四方法参考值均为0.001，达到B999分辨率；不写成p=0，更不据此称真实市场差异已显著。价格四维也有区分信息；部分MST参考值较大，同样保留。',
      '', '这些记录存在同日跨股票相关与时间相关，普通行标签置换的交换性没有验证；较小参考数值不能当正式p值。原独立样本有限样本定理不能直接套用。后续若需正式推断，须与作者核对依赖校准的假设和模拟表现；本轮没有改置换引擎或宣称分块重抽样已经有效。',
      '', '## 科学解释与限制','',
      '本轮说明新日频多特征在留出年份具有一定股票身份区分能力。换手率在银行与机械的部分学习器上增加区分，但价格特征本身也可能有信息；不把差异都归因于某一个变量。完整边际分布、均值/标准差、中位数以及所有学习器均保留。不按预想的“行业内应相同”改变特征或股票。',
      '', '当前试验只比较行业内部个股，尚未开展新特征的行业间、2015/2024跨阶段或海外分析。历史名单不代表2024全行业，康美药业等个股身份/状态与其他公司不同，均原样保留。不是因果分析、韧性排名、同分布等效证明或MCA优于全部对照的证据。',
      '', '所有24个参考p、六配置×1000标签四方法统计量、AUC及输入矩阵均独立复算；冻结源码与源指纹不变。首次受限环境PSOCK启动卡住，未出结果，保留日志后同种子恢复；完整任务每配置约15秒。论文未改，旧结果未删。']
    (OUT/'日频五特征试验说明.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    print(macro.to_string(index=False),flush=True)


if __name__=='__main__':main()
