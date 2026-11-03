"""Audit and explain observed financial results without selecting favorable seeds."""
from pathlib import Path
import gzip
import json
import hashlib
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parent
Q={"within_bank":"银行内部","within_pharma":"医药内部","within_machinery":"机械内部","between_sectors":"行业之间"}
METHODS=["MCA-kDET","DISCO","MMD","MST"]


def main():
    out=ROOT/"real_analysis"
    x=pd.read_csv(out/"results.csv")
    assert len(x)==480 and not x.duplicated(["draw","question","stride","method"]).any()
    assert set(x.draw)==set(range(1,11)) and set(x.stride)=={1,2,4} and set(x.B)=={4999}
    assert x.p_exploratory.between(.0002,1).all() and np.allclose(x.p_exploratory*5000,np.rint(x.p_exploratory*5000))
    for _,g in x.groupby(["draw","stride","method"]):
        assert set(g.question)==set(Q)
        order=np.argsort(g.p_exploratory.to_numpy());expected=np.empty(4)
        expected[order]=np.minimum(1,np.maximum.accumulate(g.p_exploratory.to_numpy()[order]*np.arange(4,0,-1)))
        np.testing.assert_allclose(expected,g.p_holm_exploratory,atol=2e-14,rtol=0)
    files=list((out/"jobs").glob("*.rds"));assert len(files)==40
    for f in files:
        with gzip.open(f,"rb") as stream:
            while stream.read(1024*1024):pass
    audit=pd.read_csv(out/"independent_method_audit.csv");assert len(audit)==12
    python_audit=json.loads((out/"independent_python_verification.json").read_text(encoding="utf-8"))
    assert python_audit["all_480_pvalues_recomputed"] and python_audit["all_primary_pvalues_identical"]
    primary=x[(x.draw==1)&(x.stride==4)].copy();assert len(primary)==16
    sensitivity=[]
    for (question,method,stride),g in x.groupby(["question","method","stride"]):
        sensitivity.append(dict(question=question,method=method,stride=int(stride),draws=len(g),
            nominal_crossings_unadjusted=int((g.p_exploratory<.05).sum()),
            nominal_crossings_Holm=int((g.p_holm_exploratory<.05).sum()),
            min_p_Holm=float(g.p_holm_exploratory.min()),max_p_Holm=float(g.p_holm_exploratory.max())))
    pd.DataFrame(sensitivity).to_csv(out/"sampling_stability.csv",index=False)
    text=["# 真实股票实证结果（探索性）","",
          "这份结果来自真实米筐股票数据，不是模拟面板。15股名单及日期沿用历史信息规则，未按p值换股。",
          "训练2015–2017，缓冲20交易日，测试2018-01-30至2024年底；逐股训练期GARCH参数固定后向前过滤。检验对象是五日估计标准化残差向量，不是原始收益。",
          "固定主分析：第1次预定随机抽样、每隔4个五日块取样；每个被选块只取一股。行业内k=5，每股15条；行业间k=3，每行业25条，行业总体为所选五股等权混合分布。",
          "四方法使用完全相同的测试向量和4999次置换标签；DISCO/MMD/MST坐标用同一拟合段公共尺度，MMD带宽为置换不变的混合测试距离中位数。",
          "每种方法内对四问题作Holm；不声称四方法共16项的联合错误控制。","",
          "## 固定主分析的实际数字","",
          "以下p值是现有设计下的探索性计算值。压力模拟未通过，不能当成已保证0.05水平的正式显著性结论。"]
    for question in Q:
        text.extend(["",f"### {Q[question]}",""])
        for method in METHODS:
            r=primary[(primary.question==question)&(primary.method==method)].iloc[0]
            text.append(f"- {method}：原始探索p={r.p_exploratory:.4f}；四项Holm探索p={r.p_holm_exploratory:.4f}。")
    text += ["","## 十次抽样的稳定性","",
             "十次重复的是同一份真实数据的随机抽样与分类器随机性，不是十套独立实证、也不是功效实验。完整展示，不取最小p、不平均p值；第1次的主结果始终保留。"]
    for r in sensitivity:
        if r["stride"]==4:
            text.append(f"- {Q[r['question']]} / {r['method']}：Holm探索p范围[{r['min_p_Holm']:.4f},{r['max_p_Holm']:.4f}]；10次中{r['nominal_crossings_Holm']}次低于名义0.05。")
    text += ["","## 可以写与不能写的结论","",
             "固定主抽样下，MCA三项行业内问题和行业间问题均未跨过名义0.05阈值；MST的行业间数值跨过阈值，其他三个方法没有。不能概括为所有方法或MCA检出行业间差异，也不能以未拒绝证明行业内同分布。",
             "原收益描述性距离显示行业间平均距离较大，主要由银行与另外两行业的距离贡献；医药与机械更接近。该点估计不是显著性检验。",
             "目前未得到稳定、方法一致的行业间检出证据，不能确定原因是实际差异弱还是样本量和方法敏感性不足。由于核尺度在两种表示间不同，不以MMD距离下降幅度定量声称过滤改善。",
             "GARCH估计误差、剩余时间依赖、近边界参数、停牌选择及较小测试样本仍限制推断；9000次验证的失败结果完整保留。对比方法也未获得本金融场景下独立的校准保证。",
             "主结果来自事前固定的抽样，不因MCA不显著改选样本。间隔1/2仅作额外探索，完整数据在results.csv。","",
             "## 正文段落建议（不直接修改稿件）","",
             "We compared the five-day distributions of estimated GARCH-standardized residuals for a fixed panel of five stocks in each of three sectors. In the prespecified primary sample, the exploratory MCA-kDET p-value for the between-sector comparison was 0.1254 (Holm-adjusted: 0.5016); the three within-sector comparisons also did not cross the nominal 0.05 threshold. The corresponding between-sector Holm-adjusted values were 0.6176 for DISCO, 0.4656 for MMD, and 0.0232 for MST. These results provide exploratory method comparisons rather than confirmed discoveries: calibration stress tests revealed residual size distortion, and non-rejection cannot establish equality of distributions.","",
             "复核：全部480行p值/Holm独立复算；40个RDS的gzip校验通过。R复读曾出现项目已知的解压/绑定内部错误，日志保留；另用纯Python读取原RDS，对全部480个p值重算，并对固定主抽样12个问题/间隔单元、每格5000列标签的四方法统计量全部复核，最大误差5.7e-13，p值完全一致，避免依赖不稳定R复核进程。"]
    (out/"真实股票结果说明.md").write_text("\n".join(text)+"\n",encoding="utf-8")
    metadata=dict(real_data=True,status="EXPLORATORY_UNVALIDATED_CALIBRATION",rows=480,primary_rows=16,draws=10,
                  gzip_crc_passed_files=40,independent_method_audit_cells=12,Holm_independently_verified=True,
                  pure_Python_RDS_and_all_primary_permutation_statistics_verified=True,
                  primary_draw=1,primary_stride=4,B=4999,paper_modified=False,
                  source_sha256={str(f.relative_to(ROOT)):hashlib.sha256(f.read_bytes()).hexdigest() for f in
                      (out/"protocol.json",ROOT/"innovations/filters.rds",ROOT/"design_two_tests/main_panel.csv")})
    (out/"verification.json").write_text(json.dumps(metadata,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(metadata,ensure_ascii=False))


if __name__=="__main__":main()
