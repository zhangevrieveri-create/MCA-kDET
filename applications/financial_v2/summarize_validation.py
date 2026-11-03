"""Summarize only atomically completed calibration chunks; never pool scenarios."""
from pathlib import Path
import argparse
import json
from datetime import datetime, timezone
import numpy as np
import pandas as pd
from scipy.stats import beta

QUESTIONS = ["within_bank", "within_pharma", "within_machinery", "between_sectors"]
NAMES = dict(zip(QUESTIONS, ["银行内部", "医药内部", "机械内部", "行业之间"]))


def interval(k, n):
    return (0.0 if k == 0 else beta.ppf(.025, k, n-k+1),
            1.0 if k == n else beta.ppf(.975, k+1, n-k))


def is_null(case, question):
    if case == "null_equal_mixtures":
        return question == "between_sectors"
    if case.startswith("null_"):
        return True
    if case == "alt_one_bank":
        return question in ("within_pharma", "within_machinery")
    return question != "between_sectors"


def summarize(root, main_length):
    files = sorted((root / "chunks").glob("*.csv"))
    if not files:
        raise SystemExit("No completed chunks yet")
    x = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
    keys = ["scenario", "rep", "d", "block_length", "question"]
    if x.duplicated(keys).any():
        raise ValueError("Overlapping chunks: duplicate replication keys")
    if not np.isfinite(x[["p", "p_holm", "statistic", "critical95"]]).all().all():
        raise ValueError("Nonfinite statistical output")
    if ((x.p < 0) | (x.p > 1) | (x.p_holm < x.p - 1e-12)).any():
        raise ValueError("Invalid p-values")
    for _, g in x.groupby(["scenario", "rep", "d", "block_length"]):
        if sorted(g.question) != sorted(QUESTIONS):
            raise ValueError("Incomplete four-question family")
        p = g.p.to_numpy(); order = np.argsort(p)
        holm = np.empty(4)
        holm[order] = np.minimum(1, np.maximum.accumulate(p[order] * np.arange(4, 0, -1)))
        if not np.allclose(holm, g.p_holm, atol=2e-14, rtol=0):
            raise ValueError("Holm mismatch against independent recomputation")
    rows = []
    for (case, d, length, question), g in x.groupby(["scenario", "d", "block_length", "question"]):
        for alpha in (.01, .05, .10):
            for adjusted, col in ((False, "p"), (True, "p_holm")):
                k = int((g[col] <= alpha + 1e-12).sum()); n = len(g)
                lo, hi = interval(k, n)
                rows.append(dict(scenario=case, d=d, block_length=length, question=question,
                                 target="size" if is_null(case, question) else "power",
                                 alpha=alpha, holm=adjusted, replications=n, rejections=k,
                                 rate=k/n, ci_lower=lo, ci_upper=hi))
    summary = pd.DataFrame(rows)
    families = []
    for (case, d, length), g in x.groupby(["scenario", "d", "block_length"]):
        nulls = [q for q in QUESTIONS if is_null(case, q)]
        for alpha in (.01, .05, .10):
            for adjusted, col in ((False, "p"), (True, "p_holm")):
                h = g[g.question.isin(nulls)].groupby("rep")[col].min()
                k = int((h <= alpha + 1e-12).sum()); n = len(h)
                lo, hi = interval(k, n)
                families.append(dict(scenario=case, d=d, block_length=length, alpha=alpha,
                                     holm=adjusted, replications=n, false_rejection_families=k,
                                     fwer=k/n, ci_lower=lo, ci_upper=hi))
    fwer = pd.DataFrame(families)
    summary.to_csv(root / "calibration_rates.csv", index=False)
    fwer.to_csv(root / "calibration_fwer.csv", index=False)
    counts = x.groupby("scenario").rep.nunique().to_dict()
    sampling = "sampling" in str(root)
    null_count = sum(all(is_null(case,q) for q in QUESTIONS) for case in counts)
    alt_count = len(counts)-null_count
    main = summary[(summary.block_length == main_length) & (summary.alpha == .05) & ~summary.holm]
    text = ["# 股票面板校准验证（按已完成块更新）", "",
            f"更新时间：{datetime.now(timezone.utc).isoformat()}", "",
            f"本目录包含 {null_count} 类四项全为原假设的场景、{alt_count} 类包含备择或部分原假设的场景；分别汇报，不合并重复次数。",
            "置信区间为逐格 Clopper–Pearson 95% 区间，不是同时区间。各场景次数见下文；100次只作开发试跑，1000次独立确认也不能替代真实数据假设审核。",
            "每次重复均重新生成相关面板并训练原四学习器。股票选择对应的共同缺失日历固定；不验证信息性缺失，也不证明真实市场平稳。",
            f"{'当前展示的日历抽样间隔' if sampling else '主块长'}：{main_length} 个收益向量。样本量逐场景：{counts}。", "",
            f"## {'固定主抽样间隔' if sampling else '主块长'}，α=0.05（未作Holm调整的单项拒绝率）", "",
            "|场景|问题|性质|次数|拒绝率|95%区间|", "|---|---|---|---:|---:|---|"]
    for row in main.itertuples():
        text.append(f"|{row.scenario}|{NAMES[row.question]}|{'误报率' if row.target == 'size' else '功效'}|{row.replications}|{row.rate:.3f}|[{row.ci_lower:.3f}, {row.ci_upper:.3f}]|")
    text += ["", "## 四项Holm调整后至少误拒一个真原假设的概率", "",
             "|场景|次数|FWER|95%区间|", "|---|---:|---:|---|"]
    for row in fwer[(fwer.block_length == main_length) & (fwer.alpha == .05) & fwer.holm].itertuples():
        text.append(f"|{row.scenario}|{row.replications}|{row.fwer:.3f}|[{row.ci_lower:.3f}, {row.ci_upper:.3f}]|")
    flags = main[(main.target == "size") & (main.ci_lower > .05)]
    text += ["", "## 执行审计", "",
             f"完整块 {len(files)}；结果行 {len(x)}；学习器警告最大值 {x.learner_warning_count.max()}；不收敛标志最大值 {x.fit_nonconvergence_flags.max()}。",
             "每个重复的四项Holm值已用Python独立复算，未发现不一致。原始统计量在R运行中逐次与冻结ovr_T核对。",
             f"主块长下真原假设单项区间下界超过0.05的格数：{len(flags)}（逐格筛查，须结合多重性与扩展试跑判断）。",
             "完整的 α=0.01/0.05/0.10、全部块长结果见 calibration_rates.csv；家族误报见 calibration_fwer.csv。",
             "未产生真实股票实证结论；未修改论文草稿、冻结模型或置换引擎。"]
    if sampling:
        text += ["", "本目录的block_length字段表示原日历抽样间隔，不是bootstrap块长。每个被选时间块只取一股；行业间每次独立均匀选股，不能强制各股出现次数相等。",
                 "null_equal_mixtures中行业内股票不同，但三行业各包含同样五种创新分布，行业混合分布确实相同；该场景只有行业间问题为真原假设。",
                 "使用冻结标签置换引擎；oracle对应独立创新块，estimated包含训练期GARCH估计误差。三种间隔并行展示，不能按真实p值选择。"]
    (root / "校准验证进度.md").write_text("\n".join(text) + "\n", encoding="utf-8")
    status = dict(updated_utc=datetime.now(timezone.utc).isoformat(), completed_chunks=len(files),
                  rows=len(x), replications_by_scenario=counts, main_block_length=main_length,
                  main_null_cells_lower_ci_above_nominal=len(flags),
                  max_learner_warnings=int(x.learner_warning_count.max()),
                  max_nonconvergence_flags=int(x.fit_nonconvergence_flags.max()))
    (root / "summary_status.json").write_text(json.dumps(status, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(status, ensure_ascii=False))
    print(main[["scenario", "question", "target", "replications", "rate", "ci_lower", "ci_upper"]].to_string(index=False))


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("directory", type=Path)
    p.add_argument("--main-length", type=int, default=8)
    a = p.parse_args()
    summarize(a.directory, a.main_length)
