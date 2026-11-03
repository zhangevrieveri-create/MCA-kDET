"""Independently verify R results and produce scientific tables, figures and report."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import platform

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from scipy.sparse.csgraph import minimum_spanning_tree
from scipy.spatial.distance import cdist
from scipy.stats import rankdata

from run_study import ROOT, prepared_path

METHODS = ["MCA", "DISCO", "MMD", "MST"]
STYLE = {"MCA": ("MCA-kDET", "#171717", "o"), "DISCO": ("DISCO", "#0072B2", "^"),
         "MMD": ("MMD", "#B97600", "s"), "MST": ("MST", "#008B67", "D")}
STAGES = ["12 weeks", "18–20 weeks", "30 weeks"]
STAGE_STYLE = [("#0072B2", "o"), ("#B97600", "s"), ("#8A559D", "^")]
CONFIG_LABELS = {"primary": "Main", "genes_50": "50 genes", "genes_200": "200 genes",
                 "cells_20": "20 cells/mouse", "cells_40": "40 cells/mouse",
                 "middle_20w_only": "Middle: 20w only", "one_tumor": "One late tumor/mouse"}


def independent_check(spec: dict) -> dict:
    folder = ROOT / "results" / spec["name"]
    source = prepared_path(spec)
    samples = pd.read_csv(source / "samples.csv")
    scores = pd.read_csv(folder / "test_scores_exact.csv", float_precision="round_trip")
    X = pd.read_csv(source / "expression.csv.gz").set_index("sampleID")
    assert X.index.tolist() == samples.sampleID.tolist()
    fit = samples.loc[samples.split.eq("fit")]
    assert not set(fit.verified_animal_block) & set(scores.verified_animal_block)
    assignments = pd.read_csv(source / "animal_label_assignments.csv")
    labels = assignments.loc[:, scores.verified_animal_block].to_numpy().T
    actual = scores.group_number.to_numpy()
    observed = np.flatnonzero((labels == actual[:, None]).all(axis=0))
    assert len(observed) == 1
    order = np.r_[observed, np.delete(np.arange(labels.shape[1]), observed)]
    labels = labels[:, order]
    assert not assignments.duplicated().any()
    perms = pd.read_csv(folder / "permutation_statistics.csv")
    assert len(perms) == labels.shape[1]
    values = pd.read_csv(folder / "pvalues.csv").set_index("method")
    # Independent, explicit Mann–Whitney rank sums for every learner/group/assignment.
    all_t = []
    for learner in ("l1", "quad", "knn", "rf"):
        learner_t = []
        for i in (1, 2, 3):
            score = scores[f"{learner}_target{i}"].to_numpy()
            ranks = rankdata(score, method="average")
            positive = labels == i
            n = positive[:, 0].sum(); m = len(score)
            auc = (ranks @ positive - n * (n + 1) / 2) / (n * (m - n))
            learner_t.append((auc - .5) / np.sqrt((m + 1) / (12 * n * (m - n))))
        all_t.extend(learner_t)
        np.testing.assert_allclose(np.max(learner_t, axis=0), perms[f"MCA_{learner}"], rtol=1e-11, atol=1e-11)
    np.testing.assert_allclose(np.max(all_t, axis=0), perms.MCA, rtol=1e-11, atol=1e-11)
    np.testing.assert_allclose(np.sum(all_t, axis=0), perms.Sum, rtol=1e-11, atol=1e-11)
    for method in values.index:
        if method == "Comb":
            rank_max = rankdata(-perms.MCA, method="max") / len(perms)
            rank_sum = rankdata(-perms.Sum, method="max") / len(perms)
            combined = np.minimum(rank_max, rank_sum)
            p = np.mean(combined <= combined[0] + 1e-10)
        else:
            v = perms[method].to_numpy()
            p = np.mean(v >= v[0] - 1e-10)
        assert abs(p - values.loc[method, "p_value"]) < 1e-12, (spec["name"], method)
    # Check competitors by direct within/between submatrices for all 560 assignments.
    training = X.loc[fit.sampleID].iloc[:, :spec["dimension"]].to_numpy()
    z = (X.loc[scores.sampleID].iloc[:, :spec["dimension"]].to_numpy() - training.mean(axis=0)) / training.std(axis=0, ddof=1)
    D2 = cdist(z, z, metric="sqeuclidean"); D = np.sqrt(D2)
    bandwidth2 = np.median(D2[np.triu_indices(len(D2), 1)])
    K = np.exp(-D2 / (2 * bandwidth2))
    disco = []; mmd = []
    for lab in labels.T:
        groups = [np.flatnonzero(lab == i) for i in (1, 2, 3)]
        disco.append(-sum(D[np.ix_(g, g)].sum() / len(g) for g in groups))
        within = [(K[np.ix_(g, g)].sum() - len(g)) / (len(g) * (len(g) - 1)) for g in groups]
        mmd.append(np.mean([within[i] + within[j] - 2 * K[np.ix_(groups[i], groups[j])].mean()
                            for i, j in ((0, 1), (0, 2), (1, 2))]))
    np.testing.assert_allclose(disco, perms.DISCO, rtol=1e-10, atol=1e-8)
    np.testing.assert_allclose(mmd, perms.MMD, rtol=1e-9, atol=1e-10)
    tree = minimum_spanning_tree(D).tocoo()
    assert len(tree.data) == len(D) - 1
    # A different valid MST may be chosen under ties. Verify the primary observed
    # edge count when no distance ties occur; never replace an R statistic silently.
    mst_observed = int(np.sum(actual[tree.row] == actual[tree.col]))
    assert mst_observed == perms.MST.iloc[0], (spec["name"], "MST tied topology requires review")
    warnings = (folder / "warnings.txt").read_text(encoding="utf-8").strip()
    return dict(configuration=spec["name"], observations=len(samples), test_animals=scores.verified_animal_block.nunique(),
                complete_assignments=len(perms), independently_verified=True, learner_warnings=warnings)


def save_figure(fig, stem: str, preview: Path, manifest: list, caption: str):
    destination = ROOT / "results" / "figures"
    destination.mkdir(exist_ok=True)
    pdf = destination / f"{stem}.pdf"
    fig.savefig(pdf, facecolor="white", metadata={"Title": caption, "Author": "MCA-kDET single-cell study"})
    preview.mkdir(parents=True, exist_ok=True)
    png = preview / f"{stem}.png"
    fig.savefig(png, dpi=220, facecolor="white")
    manifest.append(dict(file=str(pdf), sha256=hashlib.sha256(pdf.read_bytes()).hexdigest(),
                         width_inches=float(fig.get_figwidth()), height_inches=float(fig.get_figheight()),
                         caption=caption, preview=str(png)))
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preview-dir", type=Path, required=True)
    args = parser.parse_args()
    plan = json.loads((ROOT / "study_plan.json").read_text(encoding="utf-8"))
    specs = plan["specifications"]
    audit = [independent_check(spec) for spec in specs]
    report_dir = ROOT / "results"
    pd.DataFrame(audit).to_csv(report_dir / "verification.csv", index=False, encoding="utf-8-sig")
    tables = []
    for spec in specs:
        table = pd.read_csv(report_dir / spec["name"] / "pvalues.csv")
        for col in ("name", "role", "cells", "cohort"):
            table[col] = spec[col]
        tables.append(table)
    all_results = pd.concat(tables, ignore_index=True).rename(columns={"name": "configuration", "cells": "cells_per_mouse"})
    all_results.to_csv(report_dir / "all_pvalues.csv", index=False, encoding="utf-8-sig")
    wide = all_results.pivot(index="configuration", columns="method", values="p_value").reindex([s["name"] for s in specs])
    wide.to_csv(report_dir / "comparison.csv", encoding="utf-8-sig")
    primary = pd.read_csv(report_dir / "primary" / "pvalues.csv").set_index("method")
    samples = pd.read_csv(ROOT / "prepared" / "analysis" / "samples.csv")
    X = pd.read_csv(ROOT / "prepared" / "analysis" / "expression.csv.gz").set_index("sampleID").iloc[:, :100]
    fit = samples.split.eq("fit").to_numpy()
    mu = X.to_numpy()[fit].mean(axis=0); sd = X.to_numpy()[fit].std(axis=0, ddof=1)
    z = (X.to_numpy() - mu) / sd
    _, singular, vt = np.linalg.svd(z[fit], full_matrices=False)
    pc = z @ vt[:2].T
    explained = singular[:2] ** 2 / (singular ** 2).sum()
    pca = samples[["sampleID", "verified_animal_block", "group_number", "split"]].copy()
    pca["PC1"] = pc[:, 0]; pca["PC2"] = pc[:, 1]
    pca.to_csv(report_dir / "pca_scores.csv", index=False, encoding="utf-8-sig")
    auc = pd.read_csv(report_dir / "primary" / "auc.csv")
    scores = pd.read_csv(report_dir / "primary" / "test_scores_exact.csv", float_precision="round_trip")
    animal_auc_rows = []
    for learner in ("l1", "quad", "knn", "rf"):
        for target in (1, 2, 3):
            for animal in scores.loc[scores.group_number.eq(target), "verified_animal_block"].unique():
                positive = scores.verified_animal_block.eq(animal)
                other = scores.group_number.ne(target)
                a = scores.loc[positive, f"{learner}_target{target}"].to_numpy()
                b = scores.loc[other, f"{learner}_target{target}"].to_numpy()
                value = np.mean((a[:, None] > b[None, :]) + .5 * (a[:, None] == b[None, :]))
                animal_auc_rows.append(dict(learner=learner, target=target, animal=animal, auc=value))
    animal_auc = pd.DataFrame(animal_auc_rows)
    animal_auc.to_csv(report_dir / "animal_auc.csv", index=False, encoding="utf-8-sig")
    # All available cohort animals, including the independent train/test allocation.
    metadata = pd.read_csv(ROOT / "prepared" / "cell_metadata.csv")
    cohort = metadata.loc[metadata.author_stage.isin(["06_KP_12w_ND", "07_KP_20w_ND", "08_KP_30w_ND"])].copy()
    cohort["group_number"] = cohort.author_stage.map({"06_KP_12w_ND": 1, "07_KP_20w_ND": 2, "08_KP_30w_ND": 3})
    animal_table = cohort.groupby("verified_animal_block", sort=False).agg(
        group_number=("group_number", "first"), available_cells=("sampleID", "size"),
        weeks=("geo_sack_timepoint (weeks)", "first"), sex=("geo_Sex", "first"),
        tumors=("geo_tumor", "nunique"), plates=("plateID", "nunique")).reset_index()
    allocation = samples.groupby("verified_animal_block").agg(split=("split", "first"), selected_cells=("sampleID", "size")).reset_index()
    animal_table = animal_table.merge(allocation, on="verified_animal_block", validate="one_to_one")
    animal_table.to_csv(report_dir / "animal_structure.csv", index=False, encoding="utf-8-sig")
    # Published cluster labels are purely descriptive; they were not used in feature selection/tests.
    composition = pd.crosstab(cohort.verified_animal_block, cohort.clusterK12, normalize="index").reset_index()
    composition.to_csv(report_dir / "author_cluster_composition.csv", index=False, encoding="utf-8-sig")

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.titlesize": 10,
                         "axes.labelsize": 9, "axes.spines.top": False, "axes.spines.right": False,
                         "pdf.fonttype": 42, "axes.linewidth": .7, "xtick.labelsize": 8, "ytick.labelsize": 8})
    manifest = []
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.05), layout="constrained", gridspec_kw={"width_ratios": [1.25, 1]})
    ax = axes[0]
    held = pca.loc[pca.split.eq("test")]
    for i, (color, marker) in enumerate(STAGE_STYLE, 1):
        sub = held.loc[held.group_number.eq(i)]
        ax.scatter(sub.PC1, sub.PC2, s=8, c=color, marker=marker, alpha=.20, linewidths=0, rasterized=True)
        centers = sub.groupby("verified_animal_block")[["PC1", "PC2"]].mean()
        ax.scatter(centers.PC1, centers.PC2, s=45, c=color, marker=marker, edgecolors="white", linewidths=.7, label=STAGES[i - 1], zorder=4)
    ax.set(xlabel=f"PC1 ({explained[0]:.1%} fitting variance)", ylabel=f"PC2 ({explained[1]:.1%} fitting variance)", title="a  Held-out expression profiles")
    ax.legend(frameon=False, fontsize=8, loc="upper right", handletextpad=.3, borderpad=.2)
    ax = axes[1]
    for j, method in enumerate(METHODS):
        name, color, marker = STYLE[method]
        p = primary.loc[method, "p_value"]
        ax.scatter(p, j, color=color, marker=marker, s=48, zorder=3)
        ax.annotate(f"{p:.4f}", (p, j), xytext=(6, 0), textcoords="offset points", va="center", fontsize=8)
    ax.axvline(.05, color="#666666", ls="--", lw=.8)
    ax.set(xscale="log", xlim=(.005, 1.8), yticks=range(4), yticklabels=[STYLE[m][0] for m in METHODS],
           xlabel="Whole-mouse permutation p value", title="b  Main method comparison")
    ax.invert_yaxis(); ax.grid(axis="x", color=".92", lw=.6)
    save_figure(fig, "single_cell_main", args.preview_dir, manifest,
                "Main analysis. Test cells (faint) and mouse centroids (large) projected using PCA fitted only to training cells; 100 training-selected genes. Panel b compares exact whole-mouse permutation p values; dashed line 0.05. No cell-level uncertainty intervals.")

    fig, axes = plt.subplots(1, 2, figsize=(7.8, 3.8), layout="constrained", gridspec_kw={"width_ratios": [1.4, 1]})
    ax = axes[0]
    selected_configs = list(CONFIG_LABELS)
    for offset, method in zip(np.linspace(-.25, .25, 4), METHODS):
        name, color, marker = STYLE[method]
        ax.scatter(wide.loc[selected_configs, method], np.arange(7) + offset, s=28, color=color, marker=marker, label=name)
    ax.axvline(.05, color=".4", ls="--", lw=.8)
    ax.set(xscale="log", xlim=(.006, 1.05), yticks=range(7), yticklabels=list(CONFIG_LABELS.values()),
           xlabel="Exact permutation p value", title="a  Prespecified sensitivity checks")
    ax.invert_yaxis(); ax.grid(axis="x", color=".92", lw=.6)
    legend_handles = [Line2D([], [], color=STYLE[m][1], marker=STYLE[m][2], linestyle="none", label=STYLE[m][0]) for m in METHODS]
    fig.legend(handles=legend_handles, frameon=False, fontsize=8, loc="outside lower center", ncol=4)
    ax = axes[1]
    split = wide.loc[[s["name"] for s in specs if s["role"] == "split_stability"]]
    for i, method in enumerate(METHODS):
        name, color, marker = STYLE[method]
        vals = split[method].to_numpy()
        jitter = np.random.default_rng(197 + i).uniform(-.15, .15, len(vals))
        ax.scatter(i + jitter, vals, color=color, marker=marker, s=23, alpha=.85)
        ax.plot([i - .22, i + .22], [np.median(vals)] * 2, color=color, lw=1.5)
        ax.scatter(i, primary.loc[method, "p_value"], facecolors="none", edgecolors=color, marker="D", s=72, linewidths=1.0)
    ax.axhline(.05, color=".4", ls="--", lw=.8)
    ax.set(yscale="log", ylim=(.005, 1.05), xticks=range(4), xticklabels=[STYLE[m][0] for m in METHODS],
           ylabel="Exact permutation p value", title="b  Ten additional animal splits")
    ax.tick_params(axis="x", labelrotation=25); ax.grid(axis="y", color=".92", lw=.6)
    save_figure(fig, "single_cell_sensitivity", args.preview_dir, manifest,
                "Sensitivity of all four methods to genes, cells/mouse, omitting 18-week mice and one tumor per late-stage mouse. All versions retained. Ten additional splits use seeds 20261028–20261106. P values from overlapping data are dependent; displayed split fractions are stability summaries, not power or independent replications.")

    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.7), layout="constrained", sharey=True)
    learners = ["l1", "quad", "knn", "rf"]
    for target, ax in enumerate(axes, 1):
        color, marker = STAGE_STYLE[target - 1]
        for j, learner in enumerate(learners):
            sub = animal_auc.loc[animal_auc.target.eq(target) & animal_auc.learner.eq(learner)]
            offsets = np.linspace(-.12, .12, len(sub))
            ax.scatter(j + offsets, sub.auc, s=27, color=color, marker=marker, edgecolors="white", linewidths=.4)
            aggregate = auc.loc[auc.target.eq(target) & auc.learner.eq(learner), "auc"].iloc[0]
            ax.plot([j - .25, j + .25], [aggregate] * 2, color=".15", lw=1)
        ax.axhline(.5, color=".5", ls="--", lw=.8)
        ax.set(xticks=range(4), xticklabels=["L1", "Quad", "kNN", "RF"], ylim=(0, 1.03), xlabel="Fitted learner",
               title=f"{chr(96 + target)}  {STAGES[target - 1]} vs rest")
    axes[0].set_ylabel("Held-out AUC")
    save_figure(fig, "single_cell_auc", args.preview_dir, manifest,
                "One-vs-rest held-out AUCs. Each dot compares one positive-stage test mouse with all other-stage test cells; horizontal bars are pooled AUCs. Scores were fitted on separate mice. Dots are dependent descriptive summaries, not confidence intervals or post-hoc tests.")

    fig, axes = plt.subplots(2, 2, figsize=(6.8, 4.5), layout="constrained")
    perm = pd.read_csv(report_dir / "primary" / "permutation_statistics.csv")
    for ax, method in zip(axes.ravel(), METHODS):
        v = perm[method].to_numpy(); standardized = (v - v.mean()) / v.std(ddof=0)
        ax.hist(standardized, bins=25, color="#D5DADD", edgecolor="white", linewidth=.5)
        ax.axvline(standardized[0], color=STYLE[method][1], lw=1.7)
        ax.set(xlabel="Permutation-standardized statistic", ylabel="Label assignments",
               title=f"{STYLE[method][0]}: p = {primary.loc[method, 'p_value']:.4f}")
    save_figure(fig, "single_cell_calibration", args.preview_dir, manifest,
                "All 560 animal label assignments. Each statistic is centered/scaled by its complete permutation distribution for display only; vertical line is observed. P values use original untransformed statistics, including ties and the observed assignment exactly once.")
    (report_dir / "figure_manifest.json").write_text(json.dumps(dict(
        python=platform.python_version(), matplotlib=matplotlib.__version__, numpy=np.__version__,
        figures=manifest), ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    # One compact structured source for the Excel workbook, preserving numeric types.
    def records(frame):
        return json.loads(frame.to_json(orient="split", force_ascii=False, double_precision=15))
    workbook_data = {
        "main": records(primary.reset_index()), "comparison": records(wide.reset_index()),
        "animals": records(animal_table), "auc": records(auc), "animal_auc": records(animal_auc),
        "genes": records(pd.read_csv(ROOT / "prepared" / "analysis" / "genes.csv")),
        "all_results": records(all_results), "checks": records(pd.DataFrame(audit)),
    }
    (report_dir / "workbook_data.json").write_text(json.dumps(workbook_data, ensure_ascii=False, allow_nan=False), encoding="utf-8")
    rows = ["# 单细胞实证结果（2026-10-27）", "", "已完成 1 项主分析、6 项敏感性分析和 10 次额外整鼠拆分。以下是探索性实证结果；正式稿应先明确动物块原假设，不能直接引用现有细胞 iid 定理保证水平。", "",
            "## 主分析", "", "同一 KP 基因型的 12 / 18–20 / 30 周，14 鼠、2583 个可用细胞；每鼠抽 30 个。训练 6 鼠、180 细胞，测试 8 鼠、240 细胞；基因只用训练细胞筛选并标准化，主分析 100 维。训练和测试动物不重叠。", "",
            "| 方法 | 整鼠置换 p 值 | ≥观测统计量的分配数 / 560 |", "|---|---:|---:|"]
    for method in METHODS:
        p = primary.loc[method, "p_value"]
        rows.append(f"| {STYLE[method][0]} | {p:.6f} | {round(p * 560)} / 560 |")
    rows += ["", "主分析在 0.05 水平下，MCA-kDET 和 MST 拒绝共同动物块分布原假设，DISCO、MMD 未拒绝。单个数据集的 p 值大小不能当作方法功效或普遍优劣的证据。所有方法用相同测试细胞和标签空间，分类器只训练一次。", "",
             "## 敏感性和拆分稳定性", "", "| 设定 | MCA-kDET | DISCO | MMD | MST |", "|---|---:|---:|---:|---:|"]
    for name in CONFIG_LABELS:
        rows.append("| " + CONFIG_LABELS[name] + " | " + " | ".join(f"{wide.loc[name, method]:.4f}" for method in METHODS) + " |")
    rows += ["", "去掉 18 周后，中间组仅 20 周，训练/测试动物数改为 1/2，总测试动物 7 鼠，完整标签空间为 210；这是样本与拆分一起变化的检查，不能单独归因于去掉一个时间点。单肿瘤检查在细胞抽样前给每只晚期小鼠随机选一个肿瘤，所有已选肿瘤均有至少 30 个可用细胞。", "",
             "另 10 次拆分种子固定为 20261028–20261106，每次重新在其训练细胞中筛基因、估尺度并拟合，测试细胞不参与训练。", "",
             "| 方法 | p 值中位数 | p 值范围 | p ≤ 0.05 次数 / 10 |", "|---|---:|---:|---:|"]
    for method in METHODS:
        vals = split[method]
        rows.append(f"| {STYLE[method][0]} | {vals.median():.4f} | [{vals.min():.4f}, {vals.max():.4f}] | {sum(vals <= .05 + 1e-12)} / 10 |")
    rows += ["", "MCA 的额外拆分有 3 次未达到 0.05，不能写成对所有拆分都稳定显著。各拆分反复使用同一批动物，不能计算独立重复的二项置信区间，也不能将其拒绝比例称为功效。补充材料保留全部 p 值。", "",
             "## 分类器读法与样本限制", "", "12 / 18–20 / 30 周一对其余的 RF AUC 分别为 " + ", ".join(f"{auc.loc[auc.learner.eq('rf') & auc.target.eq(i), 'auc'].iloc[0]:.4f}" for i in (1, 2, 3)) + "。中间组最容易与其余测试样本区分，不能据此断言晚期变化最大或中间阶段存在因果机制。AUC 只描述当前样本和训练的判别能力，不是生物学效应量的通用尺度。",
             "", "独立动物数仅 6/5/3；晚期训练只有一只鼠，个体及来源批次可影响泛化。18 周组全雌、20 周组全雄，性别与来源部分混杂；30 周为不同肿瘤采样，早期与晚期取材方式也可能不同。保留动物/肿瘤/plate/性别/采样时间表，不能把阶段关联写成纯进展时间效应或仅依赖结构的变化。作者已有 clusterK12 比例只作描述，未用于选择基因或检验；未将这些使用全数据得到的聚类当成独立验证。", "",
             "整鼠校准假设测试动物独立、各抽样细胞块在原假设下可交换；等量细胞和整鼠留出不会让同鼠细胞独立。这里的联合动物块原假设强于单个细胞边际同分布，尚需与作者共同确定正文的正式口径。", "",
             "## 文章中的放置建议", "", "正文：数据来源与动物数、整鼠拆分和训练筛选规则；四方法主结果表；single_cell_main.pdf。补充：完整动物结构、200 基因排名、单学习器和 Sum/Comb/Pair、全部敏感性和逐拆分 p 值，以及 sensitivity、auc、calibration 三幅图。统计量的结果表只给预定主检验；补充检验与调参检查均不作为额外确认性结论，不作最小 p 值筛选。", "",
             "论文待修改清单：第 3.1 节来源改为混合遗传背景；删除‘近交所以细胞独立’的解释；应用方法写清动物块置换的新原假设；中间阶段写成 18–20 周；第 5 节加入结果并明确样本/取材混杂。遵守项目约束，本轮未直接改论文草稿。", "",
             "## 核验与复现", "", "小规模 5 基因/l1 测试通过后再跑 17 项。Python 逐项独立回算所有学习器的全部标签分配 AUC/最大型/和型统计量、所有 p 值及 Comb 排名；DISCO/MMD 用直接距离/核子矩阵复算；MST 用 SciPy 独立求树核对观测边数。17 项全部一致，见 verification.csv。学习器没有发出收敛警告。输入/代码指纹、模型输出与完整标签空间保存在各自 analysis.rds，R 版本及包版本在 sessionInfo.txt。运行时间在并行工作负载中取得，不用来做正式速度优劣比较。", "",
             "入口：run_study.py --workers=4；report_study.py --preview-dir=<临时图像目录>。前者先冻结 study_plan.json，已完成且指纹一致的结果直接复用。R 调用要移除继承的 LC_ALL/LC_CTYPE/LANG=C.UTF-8，并设置已安装的用户包目录；本次普通 Windows 进程可正常加载包，沙箱运行出现解压异常，未重装包或改引擎。", "",
             "原始数据：[GSE154989](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE154989)。原论文：[Marjanovic et al. (2020)](https://doi.org/10.1016/j.ccell.2020.06.012)。来源文件 URL/SHA256 在 raw/download_manifest.json；本轮使用现有下载，未更动原文件。", "",
             "金融应用：用户随后授权米筐API取数，第二个实证已在相邻financial目录完成18项分析，详见 ../financial/results/股票实证结果.md。密钥没有写入项目文件。"]
    (report_dir / "单细胞实证结果.md").write_text("\n".join(rows) + "\n", encoding="utf-8")
    print(wide[METHODS].to_string())
    print(f"Independently verified {len(audit)} configurations; wrote four PDF figures and report.")


if __name__ == "__main__":
    main()
