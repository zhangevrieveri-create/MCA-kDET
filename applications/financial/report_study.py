"""Verify and report every prespecified financial analysis, including null findings."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from scipy.sparse.csgraph import minimum_spanning_tree
from scipy.stats import rankdata

ROOT = Path(__file__).resolve().parent
METHODS = ["MCA", "DISCO", "MMD", "MST"]
STYLE = {"MCA": ("MCA-kDET", "#171717", "o"), "DISCO": ("DISCO", "#0072B2", "^"),
         "MMD": ("MMD", "#B97600", "s"), "MST": ("MST", "#008B67", "D")}
LABEL = {"primary": "Main: 10 stocks", "dimension_20": "20 stocks", "dimension_30": "30 stocks",
         "block_13": "13-week blocks", "block_52": "52-week blocks", "buffer_8": "8-week buffer",
         "raw_blocks": "Unfiltered returns", "filtered_random": "Random week split", "vix_15_25": "VIX thresholds 15 / 25"}
LABEL.update({f"thin4_offset{i}": f"Every 4th test week, offset {i}" for i in range(4)})
LABEL.update({f"split_{i:02d}": f"Additional split {i}" for i in range(1, 6)})


def audit_filter():
    meta = pd.read_csv(ROOT / "prepared/weeks.csv")
    x = pd.read_csv(ROOT / "prepared/weekly_returns.csv", float_precision="round_trip").set_index("week_end")
    z = pd.read_csv(ROOT / "prepared/filtered_returns.csv", float_precision="round_trip").set_index("week_end")
    saved_h = pd.read_csv(ROOT / "prepared/conditional_variance.csv", float_precision="round_trip").set_index("week_end")
    parameters = pd.read_csv(ROOT / "prepared/garch_parameters.csv")
    baseline = meta.period.eq("baseline").to_numpy()
    assert meta.loc[baseline, "week_end"].max() < "2015-01-01"
    assert parameters.convergence.eq(0).all()
    for _, row in parameters.iterrows():
        values = x[row.stock].to_numpy()
        assert np.isclose(np.nanmean(values[baseline]), row["mean"], rtol=1e-12, atol=1e-12)
        centered = values - row["mean"]
        previous_h = previous_x2 = row.initial_variance
        h = []
        for residual in centered:
            variance = row.omega + row.alpha * previous_x2 + row.beta * previous_h
            h.append(variance); previous_h = variance
            previous_x2 = residual**2 if np.isfinite(residual) else variance
        np.testing.assert_allclose(h, saved_h[row.stock], rtol=1e-10, atol=1e-10)
        np.testing.assert_allclose(centered / np.sqrt(h), z[row.stock], rtol=1e-10, atol=1e-10, equal_nan=True)
    return parameters


def verify(spec):
    folder = ROOT / "results" / spec["name"]
    scores = pd.read_csv(folder / "test_scores_exact.csv", float_precision="round_trip")
    perm = pd.read_csv(folder / "permutation_statistics.csv")
    labels = pd.read_csv(folder / "permutation_labels.csv.gz").to_numpy().T
    assert labels.shape == (len(scores), spec["permutations"] + 1)
    assert (labels[:, 0] == scores.group_number).all()
    saved = pd.read_csv(folder / "pvalues.csv").set_index("method")
    components = []
    for learner in ["l1", "quad", "knn", "rf"]:
        current = []
        for target in (1, 2, 3):
            ranks = rankdata(scores[f"{learner}_target{target}"].to_numpy(), method="average")
            positive = labels == target
            n = positive[:, 0].sum(); m = len(ranks)
            assert (positive.sum(axis=0) == n).all()
            auc = (ranks @ positive - n * (n + 1) / 2) / (n * (m - n))
            current.append((auc - .5) / np.sqrt((m + 1) / (12 * n * (m - n))))
        components.extend(current)
        np.testing.assert_allclose(np.max(current, axis=0), perm[f"MCA_{learner}"], rtol=1e-10, atol=1e-10)
    np.testing.assert_allclose(np.max(components, axis=0), perm.MCA, rtol=1e-10, atol=1e-10)
    np.testing.assert_allclose(np.sum(components, axis=0), perm.Sum, rtol=1e-10, atol=1e-10)
    for method in saved.index:
        if method == "Comb":
            combined = np.minimum(rankdata(-perm.MCA, method="max"), rankdata(-perm.Sum, method="max")) / len(perm)
            p = np.mean(combined <= combined[0] + 1e-10)
        else:
            statistic = perm[method].to_numpy()
            p = np.mean(statistic >= statistic[0] - 1e-10)
        assert abs(p - saved.loc[method, "p_value"]) < 1e-12
    allocation = pd.read_csv(folder / "week_allocation.csv")
    assert not set(allocation.loc[allocation.role.eq("fit"), "week_end"]) & set(scores.week_end)
    if spec["split"] == "blocks":
        fit_ix = np.flatnonzero(allocation.role.eq("fit"))
        test_ix = np.flatnonzero(allocation.role.eq("test"))
        assert np.min(np.abs(fit_ix[:, None] - test_ix[None, :])) > spec["buffer"]
    if spec["thin"] > 1:
        test_ix = np.flatnonzero(allocation.role.eq("test"))
        assert np.all(test_ix % spec["thin"] == spec["offset"])
    matrix_name = "filtered_returns.csv" if spec["representation"] == "filtered" else "weekly_returns.csv"
    matrix = pd.read_csv(ROOT / "prepared" / matrix_name, float_precision="round_trip").set_index("week_end").iloc[:, :spec["dimension"]]
    fit_matrix = matrix.loc[allocation.loc[allocation.role.eq("fit"), "week_end"]].to_numpy()
    test_matrix = matrix.loc[scores.week_end].to_numpy()
    scaled = (test_matrix - fit_matrix.mean(axis=0)) / fit_matrix.std(axis=0, ddof=1)
    D2 = cdist(scaled, scaled, "sqeuclidean"); D = np.sqrt(D2)
    K = np.exp(-D2 / (2 * np.median(D2[np.triu_indices(len(D2), 1)])))
    for b in [0, 1, 2, 17, len(perm) - 1]:
        groups = [np.flatnonzero(labels[:, b] == i) for i in (1, 2, 3)]
        disco = -sum(D[np.ix_(g, g)].sum() / len(g) for g in groups)
        within = [(K[np.ix_(g, g)].sum() - len(g)) / (len(g) * (len(g) - 1)) for g in groups]
        mmd = np.mean([within[i] + within[j] - 2 * K[np.ix_(groups[i], groups[j])].mean() for i, j in [(0, 1), (0, 2), (1, 2)]])
        assert np.isclose(disco, perm.DISCO.iloc[b], rtol=1e-9, atol=1e-8)
        assert np.isclose(mmd, perm.MMD.iloc[b], rtol=1e-8, atol=1e-10)
    tree = minimum_spanning_tree(D).tocoo()
    assert len(tree.data) == len(D) - 1
    assert int(np.sum(scores.group_number.to_numpy()[tree.row] == scores.group_number.to_numpy()[tree.col])) == perm.MST.iloc[0]
    warnings = (folder / "warnings.txt").read_text(encoding="utf-8").strip()
    return dict(configuration=spec["name"], verified=True, fitting_weeks=len(fit_matrix), testing_weeks=len(test_matrix), warnings=warnings)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preview-dir", type=Path, required=True)
    args = parser.parse_args(); args.preview_dir.mkdir(parents=True, exist_ok=True)
    plan = json.loads((ROOT / "study_plan.json").read_text(encoding="utf-8"))
    audit = json.loads((ROOT / "prepared/data_audit.json").read_text(encoding="utf-8"))
    parameters = audit_filter()
    checks = pd.DataFrame([verify(s) for s in plan["specifications"]])
    results = ROOT / "results"
    checks.to_csv(results / "verification.csv", index=False, encoding="utf-8-sig")
    long = pd.concat([pd.read_csv(results / s["name"] / "pvalues.csv").assign(configuration=s["name"]) for s in plan["specifications"]], ignore_index=True)
    long.to_csv(results / "all_pvalues.csv", index=False, encoding="utf-8-sig")
    wide = long.pivot(index="configuration", columns="method", values="p_value").reindex([s["name"] for s in plan["specifications"]])
    wide.to_csv(results / "comparison.csv", encoding="utf-8-sig")
    instruments = pd.read_csv(ROOT / "raw/instruments.csv")
    selection = pd.read_csv(ROOT / "raw/panel_selection.csv").dropna(subset=["selected_rank"]).sort_values("selected_rank")
    stocks = selection.merge(instruments, on="order_book_id", validate="one_to_one").merge(parameters, left_on="order_book_id", right_on="stock", validate="one_to_one")
    stocks.to_csv(results / "stocks_and_filter.csv", index=False, encoding="utf-8-sig")
    allocation = pd.read_csv(results / "primary/week_allocation.csv")
    primary_counts = pd.read_csv(results / "primary/sample_counts.csv")
    grouped = allocation.groupby("group_number").agg(calendar_weeks=("role", "size"),
                   missing=("role", lambda r: r.eq("missing").sum()), fitting=("role", lambda r: r.eq("fit").sum()),
                   testing=("role", lambda r: r.eq("test").sum()), buffer=("role", lambda r: r.eq("buffer").sum())).reset_index()
    grouped.to_csv(results / "sample_attrition.csv", index=False, encoding="utf-8-sig")
    acf = pd.read_csv(ROOT / "prepared/autocorrelation.csv")
    primary_stocks = set(selection.order_book_id.iloc[:10])
    acf10 = acf.loc[acf.stock.isin(primary_stocks)]
    square_acf = acf10.loc[acf10.power.eq(2) & acf10.lag.eq(1)].pivot(index="stock", columns="representation", values="autocorrelation")
    primary = long.loc[long.configuration.eq("primary")].set_index("method")
    rolecounts = allocation.role.value_counts().to_dict()
    test_scores = pd.read_csv(results / "primary/test_scores_exact.csv")
    actual_start = allocation.loc[allocation.role.isin(["fit", "test"]), "week_end"].min()
    actual_end = allocation.loc[allocation.role.isin(["fit", "test"]), "week_end"].max()
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9,
        "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": .7,
        "pdf.fonttype": 42, "xtick.labelsize": 8, "ytick.labelsize": 8})
    figures = []; figdir = results / "figures"; figdir.mkdir(exist_ok=True)

    def save(fig, name, caption):
        destination = figdir / f"{name}.pdf"
        fig.savefig(destination, facecolor="white", metadata={"Title": caption})
        fig.savefig(args.preview_dir / f"{name}.png", facecolor="white", dpi=220)
        figures.append(dict(file=str(destination), width_inches=float(fig.get_figwidth()), height_inches=float(fig.get_figheight()),
                            caption=caption, sha256=hashlib.sha256(destination.read_bytes()).hexdigest()))
        plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0), layout="constrained", gridspec_kw={"width_ratios": [1.4, 1]})
    ax = axes[0]; colors = ["#0072B2", "#B97600", "#8A559D"]
    for group, color in enumerate(colors, 1):
        sub = allocation.loc[allocation.group_number.eq(group) & allocation.role.ne("missing")]
        ax.scatter(pd.to_datetime(sub.week_end), sub.lagged_vix, s=7, color=color,
                   marker=["o", "s", "^"][group - 1], label=["Low", "Middle", "High"][group - 1])
    missing = allocation.loc[allocation.role.eq("missing")]
    ax.scatter(pd.to_datetime(missing.week_end), missing.lagged_vix, s=9, color=".65", marker="x", linewidths=.5, label="Stock vector missing")
    for threshold in audit["vix_baseline_tercile_thresholds"]:
        ax.axhline(threshold, color=".5", ls="--", lw=.7)
    ax.set(ylabel="Previous-week VIX close", xlabel="A-share return week", title="a  External risk states")
    ax.legend(frameon=False, fontsize=7, loc="upper right", handletextpad=.2)
    ax = axes[1]
    for j, method in enumerate(METHODS):
        value = primary.loc[method, "p_value"]
        ax.scatter(value, j, c=STYLE[method][1], marker=STYLE[method][2], s=42)
        ax.annotate(f"{value:.4f}", (value, j), xytext=(6, 0), textcoords="offset points", va="center", fontsize=8)
    ax.axvline(.05, color=".4", ls="--", lw=.8)
    ax.set(xlim=(0, 1.12), yticks=range(4), yticklabels=[STYLE[m][0] for m in METHODS], xlabel="Approximate permutation p value", title="b  Filtered, blocked main analysis")
    ax.invert_yaxis(); ax.grid(axis="x", color=".92", lw=.6)
    save(fig, "financial_main", "Main analysis of ten baseline-selected A shares. Previous-Friday VIX defines risk states using pre-2015 weekly terciles; gray crosses mark weeks without a complete stock vector. Right: 4999 test-week label permutations, approximate after GARCH filtering; dashed line 0.05. Planned 2015–2025 calendar shown; usable A-share vectors end 2025-01-31 because the historical panel retains later suspended/delisted stocks.")

    fig, ax = plt.subplots(figsize=(7.0, 6.0), layout="constrained")
    for offset, method in zip(np.linspace(-.25, .25, 4), METHODS):
        ax.scatter(wide[method], np.arange(len(wide)) + offset, s=24, color=STYLE[method][1], marker=STYLE[method][2])
    ax.axvline(.05, color=".4", ls="--", lw=.8)
    ax.set(xscale="log", xlim=(.002, 1.25), yticks=range(len(wide)), yticklabels=[LABEL[i] for i in wide.index],
           xlabel="Approximate permutation p value", title="All prespecified financial analyses")
    ax.invert_yaxis(); ax.grid(axis="x", color=".92", lw=.6)
    fig.legend(handles=[Line2D([], [], color=STYLE[m][1], marker=STYLE[m][2], linestyle="none", label=STYLE[m][0]) for m in METHODS],
               loc="outside lower center", ncol=4, frameon=False, fontsize=8)
    save(fig, "financial_sensitivity", "All 18 prespecified configurations. Every-fourth-week offsets use real calendar spacing and the same main fitting weeks. Adding stocks also changes missing-week intersections. P values are approximate and unadjusted exploratory sensitivity results; repeated splits are dependent and not power estimates.")

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0), layout="constrained")
    ax = axes[0]; all_square = acf.loc[acf.power.eq(2) & acf.lag.eq(1)].pivot(index="stock", columns="representation", values="autocorrelation")
    extra = all_square.loc[~all_square.index.isin(primary_stocks)]
    ax.scatter(extra.raw, extra.filtered, facecolors="none", edgecolors=".6", s=25, label="Additional 20 stocks")
    ax.scatter(square_acf.raw, square_acf.filtered, color=".15", s=25, label="Main 10 stocks")
    ax.plot([-.15, .35], [-.15, .35], color=".6", lw=.7, ls="--")
    ax.axhline(0, color=".8", lw=.7)
    ax.set(xlim=(-.15, .35), ylim=(-.15, .35), xlabel="Unfiltered squared-return ACF(1)", ylabel="Filtered squared-return ACF(1)", title="a  Volatility-clustering diagnostic")
    ax.legend(frameon=False, fontsize=7, loc="upper left")
    ax = axes[1]
    for representation, color, marker in [("raw", "#B97600", "s"), ("filtered", "#008B67", "o")]:
        values = acf10.loc[acf10.power.eq(2) & acf10.representation.eq(representation)].groupby("lag").autocorrelation.mean()
        ax.plot(values.index, values.values, color=color, marker=marker, label=representation.title())
    ax.axhline(0, color=".6", lw=.7, ls="--")
    ax.set(xticks=range(1, 5), xlabel="Calendar-week lag", ylabel="Mean squared-return ACF", title="b  Main ten stocks")
    ax.legend(frameon=False, fontsize=8)
    save(fig, "financial_filtering", "Autocorrelations use available observation pairs at actual calendar lags; gaps are not compressed. Parameters were fitted only on pre-2015 weeks. Left: each dot is a stock. Right: main-panel average, not a confidence interval. Reduced squared-return autocorrelation does not establish residual independence.")

    filtered = pd.read_csv(ROOT / "prepared/filtered_returns.csv").set_index("week_end").iloc[:, :10]
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.8), layout="constrained")
    correlations = []
    for group, ax in enumerate(axes, 1):
        dates = test_scores.loc[test_scores.group_number.eq(group), "week_end"]
        corr = filtered.loc[dates].corr().to_numpy()
        correlations.append(dict(group=group, weeks=len(dates), mean_pairwise_correlation=float(corr[np.triu_indices(10, 1)].mean())))
        im = ax.imshow(corr, vmin=-1, vmax=1, cmap="RdBu_r", interpolation="nearest")
        ax.set(title=f"{['Low', 'Middle', 'High'][group - 1]} VIX: n = {len(dates)} weeks", xticks=[0, 2, 4, 6, 8], xticklabels=[1, 3, 5, 7, 9],
               yticks=[0, 2, 4, 6, 8], yticklabels=[1, 3, 5, 7, 9], xlabel="Stock rank")
    axes[0].set_ylabel("Stock rank")
    fig.colorbar(im, ax=axes, label="Filtered-return correlation", shrink=.8)
    save(fig, "financial_correlations", "Descriptive correlation matrices on the identical held-out weeks used for testing, by risk state. Stock numbers follow the pre-2015 selection ranking. These are not separate covariance tests and do not isolate dependence-only changes.")
    pd.DataFrame(correlations).to_csv(results / "descriptive_correlations.csv", index=False)
    (results / "figure_manifest.json").write_text(json.dumps(figures, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    def split_records(frame):
        return json.loads(frame.to_json(orient="split", force_ascii=False, double_precision=15))
    workbook = dict(main=split_records(primary.reset_index()), comparison=split_records(wide.reset_index()), stocks=split_records(stocks),
                    attrition=split_records(grouped), checks=split_records(checks), counts=split_records(primary_counts))
    (results / "workbook_data.json").write_text(json.dumps(workbook, ensure_ascii=False, allow_nan=False), encoding="utf-8")
    rows = ["# A股周收益实证结果（2026-10-27）", "", "已完成预定 18 项分析及独立复核。本例主分析没有检出显著差异；20/30维对比方法的一些结果显著，MCA未显著，须如实呈现。", "",
            "## 数据和事前方案", "", "米筐授权 API 已下载后复权日收盘价、成交量、全天停牌记录及 2014-12-31 沪深300权重。前60权重股中，仅按2010–2014价格覆盖≥98%、正成交量≥95%筛选，取前30只；主分析前10只，20/30维取嵌套名单。未按后续是否退市或检验结果重选股票；前10股行业集中于金融，不能代表全部A股。", "",
            "每周用中国市场当周最后交易日收盘价，以相邻Friday日历周的复权价格对数差×100构造收益。两个端点均须可观测且未全天停牌；整周休市、缺失和停牌不前填，不把复牌的多周收益当成单周收益。计划窗口2015-01-09至2025-12-26；保留后来停牌/退市的历史股票，因此实际可用向量区间为 " + actual_start + " 至 " + actual_end + "。主分析573个日历周中，共212训练、211测试、29缓冲、121缺失；缺失分布见 sample_attrition.csv。", "",
            f"市场状态用上一周美国最后交易日的VIX收盘，不使用当周A股收益。2010–2014周VIX三分位阈值为 {audit['vix_baseline_tercile_thresholds'][0]:.6f} 和 {audit['vix_baseline_tercile_thresholds'][1]:.6f}，固定用于样本期。美国前周周五收盘在中国周一开盘前已知，避免当周信息进入状态。", "",
            "逐股GARCH(1,1)均值和Gaussian QML参数只用2015年前的数据估计，随后仅根据过去收益递推条件方差。缺失周使用模型条件二阶矩推进方差。分类器在2015后拟合周上训练，逐坐标尺度也只由其拟合周估计。主拆分采用26周连续块、一半块拟合、另一半测试，测试周与拟合块间隔至少4个日历周；可用训练/测试每状态至少30周，按固定种子进行仅依标签/可用性的合格拆分。", "",
            "## 主结果", "", "| 方法 | 近似置换 p 值 | 蒙特卡洛标准误（近似） |", "|---|---:|---:|"]
    for method in METHODS:
        rows.append(f"| {STYLE[method][0]} | {primary.loc[method, 'p_value']:.4f} | {primary.loc[method, 'MC_SE']:.4f} |")
    rows += ["", "各方法使用相同训练/测试周及4999组共同置换标签，分母5000包含观测分配一次。这里4999是置换标签次数，不是4999次独立实证重复。过滤与时间隔离仅改善工作假设，周标签置换仍是近似校准，不能引用iid定理声称有限样本精确水平。主检验均未拒绝，不等于证明三个状态分布相同。", "",
             "| 状态 | 训练周 | 测试周 |", "|---|---:|---:|"]
    for _, row in primary_counts.iterrows():
        rows.append(f"| {['低VIX','中VIX','高VIX'][int(row.regime) - 1]} | {int(row.fitting)} | {int(row.testing)} |")
    rows += ["", "## 全部敏感性", "", "| 设定 | MCA-kDET | DISCO | MMD | MST |", "|---|---:|---:|---:|---:|"]
    for config in wide.index:
        rows.append("| " + LABEL[config] + " | " + " | ".join(f"{wide.loc[config,m]:.4f}" for m in METHODS) + " |")
    rows += ["", "18个版本中MCA均未达0.05。20维时DISCO/MMD/MST分别为0.0254/0.0244/0.0094，30维时DISCO/MST为0.0414/0.0036。13周分块、随机拆分、抽稀偏移1时MST也有显著结果；其余不显著结果全部保留。这些是事前列出的探索性检查，未作多重校正，不能选择其中最小p作为新的主结果；加股票同时改变完整周交集，不能把维数变化当成唯一原因。", "",
             "## 过滤诊断、局限和文章用途", "",
             f"前10股平方收益的一阶自相关平均值由 {square_acf.raw.mean():.4f} 降至 {square_acf.filtered.mean():.4f}；全部30股由 {acf.loc[acf.power.eq(2)&acf.lag.eq(1)&acf.representation.eq('raw'),'autocorrelation'].mean():.4f} 降至 {acf.loc[acf.power.eq(2)&acf.lag.eq(1)&acf.representation.eq('filtered'),'autocorrelation'].mean():.4f}。所有过滤器收敛。用真实日历滞后计算可用观测对相关，不压缩缺失周；相关减弱不能证明独立、同分布或正确设定。",
             "", "主要局限：主面板金融业集中；VIX衡量美国市场隐含风险，并非中国市场自身状态；历史股票后续停牌/退市导致完整向量缺失，且缺失比例随状态不同；单一预样本参数跨较长样本期可能失配；剩余序列依赖、缓冲长度和共同因子均影响近似置换解释。不能把此例说成MCA有更高实证功效或过滤后已消除全部依赖。", "",
             "正文建议：简要交代数据、外部指标、真实可用日期、近似校准及四方法主结果；金融例作为无显著证据的对照。补充材料保留18项结果、股票与过滤参数、样本删失表、过滤诊断和按状态相关矩阵。现有sec3‘外部指标/周收益/GARCH/分块缓冲’承诺已落实，但应改写为没有独立性保证的工作处理。论文草稿本轮未改。", "",
             "## 独立核验与复现", "", "Python重新根据参数、过去收益和缺失传播规则复算30条条件方差与残差，核对均值仅来自2015年前。18项中逐一核对分块缓冲、抽稀偏移和训练/测试隔离；重新计算所有标签的学习器AUC/最大型/和型及各p值、Comb联合排名；直接距离/核子矩阵抽核5个标签分配的DISCO/MMD，并由SciPy独立求MST核对观测边数。均通过，见verification.csv；若warnings列非空应按原文核对，不凭程序退出状态宣称无警告。原文件、原引擎未修改。", "",
             "复现顺序：fetch_data.py（文件完整时不登录、不重下）、prepare_data.py、prepare_filter.R、run_study.py、report_study.py。授权通过环境变量RQDATA_LICENSE提供，不把密钥写到脚本、报告或manifest。R进程要移除继承的C.UTF-8语言设置，并载入本机用户包目录。", "",
             "来源：[米筐行情API](https://www.ricequant.com/doc/rqdata/python/generic-api)、[历史指数权重](https://www.ricequant.com/doc/rqdata/python/indices-mod)、[停牌接口](https://www.ricequant.com/doc/rqdata/python/stock-mod)、[Cboe VIX历史数据](https://www.cboe.com/tradable-products/vix/vix-historical-data/)。下载参数、URL及SHA256在raw/download_manifest.json。"]
    (results / "股票实证结果.md").write_text("\n".join(rows) + "\n", encoding="utf-8")
    print(wide[METHODS].to_string())
    print(f"Verified all {len(checks)} financial configurations and all 30 GARCH recursions; wrote four figures and report.")


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    main()
