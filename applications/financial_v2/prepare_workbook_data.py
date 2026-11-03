"""Collect verified research outputs for one reader-facing workbook."""
from pathlib import Path
import json
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parent
Q={"within_bank":"银行内部","within_pharma":"医药内部","within_machinery":"机械内部","between_sectors":"行业之间"}
CASES={"null_iid":"独立正态","null_ar50":"AR(0.5)","null_ar90":"AR(0.9)","null_garch95":"GARCH(0.95)",
       "null_garch99":"GARCH(0.99)","null_t5_ar":"t5+AR(0.5)","null_heterogeneous":"不同GARCH参数、相同创新",
       "null_garch95_t5":"GARCH(0.95)、t5","null_garch999":"GARCH(0.999)",
       "alt_sector_mean":"行业均值差异","alt_sector_scale":"行业尺度差异","alt_sector_dependence":"行业依赖差异","alt_one_bank":"一只银行股偏离",
       "alt_sector_shape":"行业创新形状不同","null_equal_mixtures":"行业混合相同、个股成分不同"}


def calibration(folder,method,annotate_target=False):
    rate_path=folder/"calibration_rates.csv"
    if not rate_path.exists():
        return []
    rates=pd.read_csv(rate_path)
    rows=[]
    for r in rates[(rates.alpha==.05)&~rates.holm].itertuples():
        question=Q[r.question]+("（功效）" if annotate_target and r.target=="power" else "")
        rows.append([CASES.get(r.scenario,r.scenario),method,question,int(r.block_length),int(r.replications),int(r.rejections),r.rate,r.ci_lower,r.ci_upper])
    fwer=pd.read_csv(folder/"calibration_fwer.csv")
    for r in fwer[(fwer.alpha==.05)&fwer.holm].itertuples():
        rows.append([CASES.get(r.scenario,r.scenario),method,"Holm家族误报",int(r.block_length),int(r.replications),int(r.false_rejection_families),r.fwer,r.ci_lower,r.ci_upper])
    return rows


def main():
    panel=pd.read_csv(ROOT/"design_two_tests/main_panel.csv")
    params=pd.read_csv(ROOT/"innovations/garch_parameters.csv")
    acf=pd.read_csv(ROOT/"innovations/autocorrelation.csv")
    acf=acf[(acf.power==2)&(acf.lag==1)].pivot(index="stock",columns="representation",values="autocorrelation")
    merged=panel.merge(params,left_on="order_book_id",right_on="stock",validate="one_to_one").join(acf,on="stock")
    model=[]
    for r in merged.itertuples():
        model.append([r.order_book_id.split('.')[0],r.symbol,r.first_industry_name,r.mean,r.omega,r.alpha,r.beta,r.persistence,r.raw,r.filtered,int(r.near_boundary),int(r.fitting_observations)])
    pairs=[]
    summaries=[]
    for representation,folder in (("原始收益",ROOT/"descriptive"),("标准化残差",ROOT/"innovations/descriptive")):
        for d in (5,20):
            x=pd.read_csv(folder/f"stock_pair_distances_d{d}.csv")
            for r in x.itertuples():
                pairs.append([representation,d,r.stock_a.split('.')[0],r.stock_b.split('.')[0],r.sector_a,r.sector_b,
                              "行业内" if r.type=="within" else "行业间",r.mmd_unsquared,int(r.complete_test_calendar_blocks)])
            s=json.loads((folder/f"distance_summary_d{d}.json").read_text(encoding="utf-8"))
            summaries.append(dict(representation=representation,d=d,n=s["common_test_blocks"],W=s["W"],B=s["B"],Delta=s["Delta"]))
    rawcal=[]
    for method,folder in (("初轮移动块",ROOT/"validation_v1_d5"),("长块移动诊断",ROOT/"validation_longblock_pilot_d5"),
                          ("环形块诊断",ROOT/"validation_longblock_pilot_d5/circular")):
        rawcal+=calibration(folder,method)
    innovationcal=[]
    for phase,label in (("oracle","真创新基准"),("estimated","估计残差")):
        innovationcal+=calibration(ROOT/"innovation_validation_pilot_d5"/phase,label)
    samplingcal=calibration(ROOT/"sampling_fullfit_confirmation_d5/estimated","估计残差·确认",annotate_target=True)
    case_order={case:i for i,case in enumerate(("null_iid","null_garch95","null_garch99","null_heterogeneous","null_garch95_t5","null_garch999","null_equal_mixtures","alt_sector_shape","alt_one_bank"))}
    translated_order={CASES[case]:i for case,i in case_order.items()}
    samplingcal.sort(key=lambda r:(r[3]!=4,translated_order[r[0]],r[3],r[2]))
    for phase,label in (("oracle","真创新·试跑"),("estimated","估计残差·试跑")):
        samplingcal+=calibration(ROOT/"sampling_fullfit_pilot_d5"/phase,label,annotate_target=True)
    real=pd.read_csv(ROOT/"real_analysis/results.csv")
    methods=["MCA-kDET","DISCO","MMD","MST"]
    real_primary=[]
    for question in Q:
        for method in methods:
            r=real[(real.draw==1)&(real.stride==4)&(real.question==question)&(real.method==method)].iloc[0]
            real_primary.append([Q[question],method,int(r.populations),int(r.test_n_per_population),
                                 r.p_exploratory,r.p_holm_exploratory,r.statistic,int(r.B),"探索性；校准未通过"])
    real_stability=[]
    for (draw,stride,question),g in real.groupby(["draw","stride","question"]):
        by=g.set_index("method")
        real_stability.append([int(draw),int(stride),Q[question],int(g.test_n_per_population.iloc[0]),
          by.loc["MCA-kDET","p_exploratory"],by.loc["MCA-kDET","p_holm_exploratory"],
          by.loc["DISCO","p_holm_exploratory"],by.loc["MMD","p_holm_exploratory"],by.loc["MST","p_holm_exploratory"],
          "主结果" if draw==1 and stride==4 else "完整保留的敏感性"])
    joint=pd.read_csv(ROOT/"joint_analysis/results.csv")
    joint_primary=[]
    for question in Q:
        for method in methods:
            r=joint[(joint.draw==1)&(joint.stride==4)&(joint.question==question)&(joint.method==method)].iloc[0]
            joint_primary.append(["行业间·周横截面" if question=="between_sectors" else Q[question],method,
              int(r.populations),int(r.test_n_per_population),r.p_exploratory,r.p_holm_exploratory,r.statistic,int(r.B),"探索性；水平未保证"])
    joint_draws=[]
    for draw,g in joint[(joint.question=="between_sectors")&(joint.stride==4)].groupby("draw"):
        by=g.set_index("method")
        joint_draws.append([int(draw),int(g.test_n_per_population.iloc[0]),by.loc["MCA-kDET","p_exploratory"],
          by.loc["MCA-kDET","p_holm_exploratory"],by.loc["DISCO","p_holm_exploratory"],by.loc["MMD","p_holm_exploratory"],by.loc["MST","p_holm_exploratory"]])
    events=pd.read_csv(ROOT/"market_events/results.csv")
    risk=pd.read_csv(ROOT/"market_events/risk_diagnostics.csv")
    event_names={"China2015":"A股2015","China2024":"A股9·24","US2020":"美股2020"}
    event_screen=[]
    for (event,module),g in events.groupby(["event","module"],sort=False):
        by=g.set_index("method");m=by.loc["MCA-kDET"]
        frequencies=risk[(risk.event==event)&(risk.module==module)].sort_values("phase").joint_downside_frequency
        event_screen.append([event_names[event],module,int(m.d),m.test_n,m.p_exploratory,m.p_holm_exploratory,
          by.loc["DISCO","p_holm_exploratory"],by.loc["MMD","p_holm_exploratory"],by.loc["MST","p_holm_exploratory"],
          " / ".join(f"{100*v:.1f}%"for v in frequencies)])
    event_risk=[]
    for r in risk.itertuples():
        event_risk.append([event_names[r.event],r.module,r.label,int(r.n_observed),r.joint_downside_frequency,
          r.portfolio_worst10pct_mean_loss,int(r.tail_observations),r.portfolio_mean,r.portfolio_sd])
    precision=pd.read_csv(ROOT/"market_events/mc_precision/results.csv")
    event_precision=[[event_names[r.event],r.module,r.p_original,r.p_refined,r.holm_original,r.holm_refined,
                      r.mc_tail_lower,r.mc_tail_upper,int(r.B)] for r in precision.itertuples()]
    daily_ref=pd.read_csv(ROOT/"daily_features_v1/reference_results.csv")
    daily_macro=pd.read_csv(ROOT/"daily_features_v1/macro_auc.csv")
    daily_count=pd.read_csv(ROOT/"daily_features_v1/sample_counts.csv")
    daily_data=pd.read_csv(ROOT/"daily_features_v1/complete_daily_features.csv.gz")
    daily_trial=[];daily_auc=[]
    for sector in ["银行","医药","机械"]:
        for rep in ["full5","price4"]:
            g=daily_ref[(daily_ref.sector==sector)&(daily_ref.representation==rep)].set_index("method")
            rf=daily_macro[(daily_macro.sector==sector)&(daily_macro.representation==rep)&(daily_macro.learner=="rf")].macro_auc.iloc[0]
            daily_trial.append([sector,"五维含换手率" if rep=="full5" else "仅四个价格特征",g.fitting_n.iloc[0],g.test_n.iloc[0],
                *[g.loc[m,"iid_reference_p_unvalidated"] for m in methods],rf,"独立样本参考；非正式p值"])
        for learner in ["l1","quad","knn","rf"]:
            a=daily_macro[(daily_macro.sector==sector)&(daily_macro.learner==learner)].set_index("representation").macro_auc
            daily_auc.append([sector,learner,a["full5"],a["price4"],a["full5"]-a["price4"]])
    daily_counts=[]
    for r in pd.read_csv(ROOT/"design_two_tests/main_panel.csv").itertuples():
        c=daily_count[daily_count.stock==r.order_book_id].set_index("split")
        daily_counts.append([r.first_industry_name,r.order_book_id.split('.')[0],r.symbol,int(c.loc["fit","valid_days"]),int(c.loc["test","valid_days"])])
    fcols=["overnight_log_return","intraday_log_return","log_high_low_range","close_position","turnover_today"]
    daily_samples=[[r.date,r.stock.split('.')[0],r.name,r.sector,"拟合" if r.split=="fit" else "测试",*[getattr(r,f)for f in fcols]]for r in daily_data.itertuples()]
    spec=dict(summaries=summaries,
              daily_trial=daily_trial,daily_auc=daily_auc,daily_counts=daily_counts,daily_samples=daily_samples,
              daily_trial_headers=["行业","表示","各股拟合日数","各股测试日数","MCA参考值","DISCO参考值","MMD参考值","MST参考值","RF宏平均AUC","解释边界"],
              daily_auc_headers=["行业","学习器","五维宏平均AUC","价格四维AUC","差值"],
              daily_count_headers=["行业","股票代码","名称","拟合有效日","测试有效日"],
              daily_sample_headers=["交易日","股票代码","名称","行业","用途","隔夜对数涨跌(%)","日内对数涨跌(%)","对数振幅(%)","收盘位置[0,1]","换手率原单位"],
              event_screen=event_screen,event_risk=event_risk,
              event_precision=event_precision,
              event_precision_headers=["案例","模块","原p","加密p","原Holm15","新Holm15","MC下界","MC上界","置换次数B"],
              event_headers=["案例","模块","维数d","三阶段测试数","MCA探索p","MCA Holm15","DISCO Holm15","MMD Holm15","MST Holm15","三阶段半数越界比例"],
              event_risk_headers=["案例","模块","阶段","观察日数","半数越界比例","最差10%日均损失(%)","尾部日数","组合日均收益(%)","组合日波动(%)"],
              joint_primary=joint_primary,joint_draws=joint_draws,
              joint_draw_headers=["抽样编号","每行业测试数","MCA探索p","MCA Holm","DISCO Holm","MMD Holm","MST Holm"],
              real_primary_headers=["检验问题","方法","总体数k","每组测试数","探索p","Holm探索p","统计量","置换次数B","结果性质"],
              real_primary=real_primary,
              real_stability_headers=["抽样编号","日历间隔","检验问题","每组测试数","MCA探索p","MCA Holm","DISCO Holm","MMD Holm","MST Holm","展示规则"],
              real_stability=real_stability,
              overview_headers=["表示","维数d","描述面板完整段","行业内平均W","行业间平均B","差值B-W"],
              model_headers=["股票代码","名称","行业","训练均值(%)","omega","alpha","beta","持续性","原收益平方ACF1","残差平方ACF1","边界标记","训练有效日数"],
              models=model,
              pair_headers=["表示","维数d","股票a","股票b","行业a","行业b","类型","经验MMD距离","共同测试段"],pairs=pairs,
              calibration_headers=["模拟场景","流程","检验或家族","块长(向量)","重复次数","拒绝次数","拒绝率","区间下界","区间上界"],
              raw_calibration=rawcal,innovation_calibration=innovationcal,
              sampling_calibration=samplingcal,
              sampling_headers=["模拟场景","流程","检验或家族","日历间隔(5日块)","重复次数","拒绝次数","拒绝率","区间下界","区间上界"],
              innovation_calibration_available=bool(innovationcal),
              notes=["各行业五股，名单由2017年底历史信息确定。训练2015–2017，缓冲20交易日，测试2018–2024。",
                     "距离为描述性非平方经验MMD，股票对相关。两种表示使用各自训练尺度与带宽，不能直接比较改善幅度。",
                     "真实四问题、四方法已作探索计算；固定主抽样不替换。创新候选压力校准未通过，不能作正式显著性或等效结论。",
                     "区间为各场景拒绝率的逐格Clopper–Pearson 95%区间，不能称作同时区间。"],
              source="米筐RQData历史沪深300名单、复权收盘价、成交量、停牌和中信2019历史分类；拟合与校准由本项目计算。",
              source_url="https://www.ricequant.com/doc/rqdata/python/stock-mod")
    folder=ROOT/"workbook_tools";folder.mkdir(exist_ok=True)
    path=folder/"workbook_input.json"
    path.write_text(json.dumps(spec,ensure_ascii=False,allow_nan=False,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"stock_models":len(model),"distance_pairs":len(pairs),"raw_calibration_rows":len(rawcal),"innovation_calibration_rows":len(innovationcal),"sampling_calibration_rows":len(samplingcal)}))


if __name__=="__main__":
    main()
