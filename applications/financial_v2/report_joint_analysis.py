"""Combine retained within-stock tests with the new sector-joint hypothesis."""
from pathlib import Path
import json
import hashlib
import gzip
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parent
METHODS=["MCA-kDET","DISCO","MMD","MST"]
Q={"within_bank":"银行内部","within_pharma":"医药内部","within_machinery":"机械内部","between_sectors":"行业间·周横截面"}


def main():
    out=ROOT/"joint_analysis"
    original=pd.read_csv(ROOT/"real_analysis/results.csv")
    between=pd.read_csv(out/"between_results.csv")
    within=original[original.question!="between_sectors"].copy()
    assert len(within)==360 and len(between)==120
    within["target"]="individual_stock_five_day_residual_vector"
    between["target"]="sorted_sector_cross_section_weekly_residual_scores"
    x=pd.concat([within,between],ignore_index=True)
    assert len(x)==480 and not x.duplicated(["draw","stride","question","method"]).any()
    for _,g in x.groupby(["draw","stride","method"]):
        assert set(g.question)==set(Q)
        order=np.argsort(g.p_exploratory.to_numpy());p=np.empty(4)
        p[order]=np.minimum(1,np.maximum.accumulate(g.p_exploratory.to_numpy()[order]*np.arange(4,0,-1)))
        x.loc[g.index,"p_holm_exploratory"]=p
    x.to_csv(out/"results.csv",index=False)
    primary=x[(x.draw==1)&(x.stride==4)].copy()
    primary.to_csv(out/"primary_results.csv",index=False)
    counts=[]
    for (q,m,s),g in x.groupby(["question","method","stride"]):
        counts.append(dict(question=q,method=m,stride=int(s),draws=len(g),nominal_crossings_Holm=int((g.p_holm_exploratory<.05).sum()),
          min_p_Holm=float(g.p_holm_exploratory.min()),max_p_Holm=float(g.p_holm_exploratory.max())))
    pd.DataFrame(counts).to_csv(out/"sampling_stability.csv",index=False)
    for f in (out/"jobs").glob("*.rds"):
        with gzip.open(f,"rb")as r:
            while r.read(1024*1024):pass
    lines=["# 行业内个股与行业间联合结构：新探索分析","",
      "本分析在前一版行业混合分布结果之后提出，属于新的探索性表示，不冒称为事前注册。新表示的股票、日期、特征和主抽样在看新p值前固定。",
      "行业内：保留原来每股一个总体的五日标准化残差向量；银行、医药、机械各五股，三个原始p值不重新拟合、不挑抽样。",
      "行业间：同一行业每只股票在一周的五个日标准化残差求和、除以sqrt(5)，形成每股一个周分数，再对五股周分数升序排列，得到五维行业横截面。没有跨股票平均。该周分数不是重新拟合周GARCH后的创新。",
      "排序的坐标含义依次为最低、第二低、中位、第二高、最高，保留横截面的联合结构信息，主动去掉股票身份。比较的是三个行业的排序向量分布，不是个股等权混合分布，也不等于识别完整带股票身份的联合分布。",
      "沿用2015–2017拟合期、20交易日缓冲、2018-01-30至2024年底测试；行业间训练每组116条，每个测试日期只选一个行业的整条五维向量。主抽样每组25条，其他九次及间隔1/2全部保留。",
      "主抽样draw1、间隔4；四方法同观测、同4999次标签置换。每方法将三项既有行业内p值与新的行业间p值合成四项家族，重新作Holm。不得沿用旧行业混合的Holm。","",
      "## 固定主抽样的真实结果（全部仍为探索性）",""]
    for q in Q:
        lines.extend([f"### {Q[q]}",""])
        for method in METHODS:
            r=primary[(primary.question==q)&(primary.method==method)].iloc[0]
            lines.append(f"- {method}：原始探索p={r.p_exploratory:.4f}，四问题Holm探索p={r.p_holm_exploratory:.4f}。")
        lines.append("")
    lines.extend(["## 十次抽样的稳定性","","十次是同一真实数据的抽样/学习器随机性，不能解释成十个独立实证或功效。固定首个结果始终保留，不选最小p。"])
    for r in counts:
        if r["stride"]==4 and (r["question"]=="between_sectors" or r["method"]=="MCA-kDET"):
            lines.append(f"- {Q[r['question']]} / {r['method']}：Holm探索p范围[{r['min_p_Holm']:.4f},{r['max_p_Holm']:.4f}]；10次中{r['nominal_crossings_Holm']}次低于名义0.05。")
    lines.extend(["","## 写作边界","",
      "行业内与行业间现在回答不同问题：前者比较个股时间向量的边际分布，后者比较一周行业横截面状态的排序向量分布。两种结果可以共存，但不能通过改变表示来声称同一原假设得到更优p值。",
      "既有行业内原始p值保持完全相同；四项家族调整值可因新的行业间问题而变化。行业内未拒绝不证明相同；行业间即使数值低于0.05也不直接成为已保证水平的正式发现。",
      "原压力校准未覆盖新排序表示，不能宣称换表示已经解决过滤估计误差、非平稳、停牌选择或剩余时间依赖；全部真实p值仍标为探索性，原结果完整保留，论文草稿未改。"])
    lines.extend(["","## 其他固定间隔的完整对照（不是重新选择主分析）",""])
    for stride in (1,2,4):
        for method in METHODS:
            r=x[(x.draw==1)&(x.stride==stride)&(x.question=="between_sectors")&(x.method==method)].iloc[0]
            lines.append(f"- 间隔{stride} / {method}：每组{r.test_n_per_population}条，原始探索p={r.p_exploratory:.4f}，Holm探索p={r.p_holm_exploratory:.4f}。")
    (out/"行业联合结构结果说明.md").write_text("\n".join(lines)+"\n",encoding="utf-8")
    result=dict(status="EXPLORATORY_NEW_JOINT_TARGET",new_between_method_results=120,retained_within_method_results=360,
       original_within_pvalues_unchanged=True,Holm_recomputed_for_new_four_hypothesis_family=True,
       primary_draw=1,primary_stride=4,B=4999,source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in
          [ROOT/"real_analysis/results.csv",out/"between_results.csv",out/"protocol.json",out/"sorted_sector_vectors.csv",out/"fixed_schedules.csv"]})
    (out/"report_audit.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print(primary[["question","method","p_exploratory","p_holm_exploratory"]].to_string(index=False))


if __name__=="__main__":main()
