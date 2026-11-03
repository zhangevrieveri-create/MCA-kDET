"""Summarizes the v3 financial application.

    python report_v3.py

Reads results/pvalues_jobs_main.csv, results/pvalues_jobs_calibration.csv (whichever
exist), results/main/*/ovr_T.csv and prepared/<ds>/; writes results/summary_*.csv and
results/报告_v3.md. Nothing is re-estimated here except descriptive statistics.
"""
import json, os, glob
import numpy as np, pandas as pd
from scipy.stats import beta as Beta

ROOT = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(ROOT, "results")
P = json.load(open(os.path.join(ROOT, "protocol.json"), encoding="utf-8"))
METH = ["MCA", "MCA_l1", "MCA_quad", "MCA_knn", "MCA_rf", "Sum", "Comb", "DISCO", "MMD", "MST", "Schott"]
REG = {1: "low VIX", 2: "medium VIX", 3: "high VIX"}
lines = ["# 第 5 章金融实证 v3：结果汇总（自动生成）", "",
         f"协议：`protocol.json`（{P['frozen_on']}）。p 值为依赖数据下的近似置换 p 值，校准证据见第 2 节。", ""]


def fmt(x):
    return "" if pd.isna(x) else f"{x:.4f}"


def cp(k, n, a=0.05):
    lo = Beta.ppf(a / 2, k, n - k + 1) if k > 0 else 0.0
    hi = Beta.ppf(1 - a / 2, k + 1, n - k) if k < n else 1.0
    return lo, hi


# ---------------------------------------------------------------- 1 main results
fm = os.path.join(RES, "pvalues_jobs_main.csv")
if os.path.exists(fm):
    m = pd.read_csv(fm)
    m["setting"] = m["dataset"] + " / " + m["representation"] + np.where(m["label_col"] == "regime", "", " / VIX 15-25")
    tab = m.pivot_table(index="setting", columns="method", values="p_value", aggfunc="first")
    tab = tab[[c for c in METH if c in tab.columns]]
    order = [s for s in [f"{d} / {r}" for d in ("ind30", "ind10", "ind49") for r in ("raw", "innov", "dir", "innov_pre", "dir_pre")]
             + ["ind30 / innov / VIX 15-25", "ind30 / dir / VIX 15-25"] if s in tab.index]
    tab = tab.loc[order]
    tab.to_csv(os.path.join(RES, "summary_main_pvalues.csv"))
    lines += ["## 1 主结果（B = 4999）", "", "| 设定 | " + " | ".join(tab.columns) + " |", "|---" * (len(tab.columns) + 1) + "|"]
    for s, r in tab.iterrows():
        lines.append(f"| {s} | " + " | ".join(fmt(r[c]) for c in tab.columns) + " |")
    lines.append("")
    prim = [f"ind30 / {h}" for h in P["primary"]["hypotheses"]]
    if all(s in tab.index for s in prim):
        pv = np.array([tab.loc[s, "MCA"] for s in prim]); o = np.argsort(pv)
        adj = np.empty(2); adj[o[0]] = min(1, 2 * pv[o[0]]); adj[o[1]] = min(1, max(adj[o[0]], pv[o[1]]))
        lines += ["主检验（MCA-kDET，Holm，两项）：" + "；".join(f"{s}: p = {pv[i]:.4f}，Holm = {adj[i]:.4f}" for i, s in enumerate(prim)), ""]
    fit_test = m.drop_duplicates("setting").set_index("setting")[["fit", "test", "d"]].loc[order]
    lines += ["各设定的拟合 / 测试周数（低 / 中 / 高）：", ""] + [f"- {s}: fit {r.fit}, test {r.test}, d = {r.d}" for s, r in fit_test.iterrows()] + [""]

# ---------------------------------------------------------------- 2 calibration
fc = os.path.join(RES, "pvalues_jobs_calibration.csv")
if os.path.exists(fc):
    c = pd.read_csv(fc)
    rows = []
    for (kind, rep, meth), g in c.groupby(["kind", "representation", "method"]):
        n = len(g); k = int((g["p_value"] <= 0.05).sum()); lo, hi = cp(k, n)
        rows.append(dict(kind=kind, representation=rep, method=meth, n=n, rejections=k, rate=k / n, ci_low=lo, ci_high=hi))
    cal = pd.DataFrame(rows); cal.to_csv(os.path.join(RES, "summary_calibration.csv"), index=False)
    lines += ["## 2 依赖下的校准检查（每格 200 次，B = 199，α = 0.05）", "",
              "| 检查 | 表示 | 方法 | 次数 | 拒绝率 | 95% 区间 |", "|---|---|---|---|---|---|"]
    for _, r in cal[cal.method.isin(["MCA", "MCA_l1", "MCA_quad", "MCA_knn", "MCA_rf", "DISCO", "MMD", "MST", "Schott"])].iterrows():
        lines.append(f"| {r.kind} | {r.representation} | {r.method} | {r.n} | {r.rate:.3f} | [{r.ci_low:.3f}, {r.ci_high:.3f}] |")
    lines.append("")

# ---------------------------------------------------------------- 3 localization
tt = []
for f in sorted(glob.glob(os.path.join(RES, "main", "*", "ovr_T.csv"))):
    d = pd.read_csv(f); d["setting"] = os.path.basename(os.path.dirname(f)); tt.append(d)
if tt:
    tt = pd.concat(tt); tt.to_csv(os.path.join(RES, "summary_ovr_T.csv"), index=False)
    lines += ["## 3 一对其余的标准化统计量 T（观测标签；行为状态，列为学习器）", ""]
    for s, g in tt.groupby("setting"):
        pv = g.pivot(index="regime", columns="learner", values="T")
        au = g.pivot(index="regime", columns="learner", values="auc")
        lines += [f"**{s}**", "", "| 状态 | " + " | ".join(f"{a} T (AUC)" for a in pv.columns) + " |", "|---" * (len(pv.columns) + 1) + "|"]
        for i in pv.index:
            lines.append(f"| {REG[i]} | " + " | ".join(f"{pv.loc[i, a]:.2f} ({au.loc[i, a]:.3f})" for a in pv.columns) + " |")
        lines.append("")

# ---------------------------------------------------------------- 4 economic descriptives
desc = []
for ds in P["datasets"]:
    pdir = os.path.join(ROOT, "prepared", ds)
    if not os.path.exists(os.path.join(pdir, "weeks.csv")):
        continue
    wk = pd.read_csv(os.path.join(pdir, "weeks.csv")); reg = wk.regime.values
    Z = pd.read_csv(os.path.join(pdir, "X_innov.csv")).values; X = pd.read_csv(os.path.join(pdir, "X_raw.csv")).values
    d = Z.shape[1]; q = np.quantile(Z, 0.05, axis=0); low = Z < q
    iu = np.triu_indices(d, 1)
    R1 = np.corrcoef(Z[reg == 1].T); w = np.linalg.solve(R1, np.ones(d)); w /= w.sum()
    pred = w @ R1 @ w
    for g in (1, 2, 3):
        Zg = Z[reg == g]; R = np.corrcoef(Zg.T); ev = np.sort(np.linalg.eigvalsh(R))[::-1]
        L = low[reg == g].astype(float); co = (L.T @ L / len(L))[iu]
        # scale-free version: each industry's own 5% quantile within the regime,
        # so that a larger innovation variance in a regime does not inflate the count
        Lw = (Zg_ := Z[reg == g]) < np.quantile(Zg_, 0.05, axis=0)
        Lw = Lw.astype(float); cow = (Lw.T @ Lw / len(Lw))[iu]
        Sg = np.cov(Zg.T)
        desc.append(dict(dataset=ds, regime=REG[g], weeks=len(Zg),
                         raw_vol=float(np.mean(X[reg == g].std(0))),
                         innov_var=float(np.mean(Zg.var(0))),
                         mean_corr=float(R[iu].mean()),
                         absorption=float(ev[:int(np.ceil(d / 5))].sum() / d),
                         innov_sd=float(np.mean(Zg.std(0))),
                         co_exceed=float(co.mean()), co_exceed_ratio=float(co.mean() / 0.0025),
                         co_exceed_within=float(cow.mean()), co_exceed_within_ratio=float(cow.mean() / 0.0025),
                         ew_var_ratio=float(np.ones(d) @ Sg @ np.ones(d) / d ** 2 / np.mean(np.diag(Sg))),
                         mv_risk_ratio=float((w @ np.corrcoef(Zg.T) @ w) / pred)))
if desc:
    desc = pd.DataFrame(desc); desc.to_csv(os.path.join(RES, "summary_economic.csv"), index=False)
    lines += ["## 4 经济后果（描述性，全部周，按状态）", "",
              "列：原始周收益平均波动；创新平均方差；创新平均两两相关；吸收比率（前 ⌈d/5⌉ 个特征值占比）；"
              "下尾共同越界概率（两行业同时低于各自 5% 分位，独立时为 0.0025）及其相对独立的倍数；"
              "等权组合创新方差相对单行业方差之比；低 VIX 期估计的最小方差组合（相关矩阵）在各状态的风险与低 VIX 期之比。", "",
              "共同越界有两种口径：全样本阈值（受各状态创新方差不同的影响）和状态内阈值（每个行业在该状态内的 5% 分位，不受尺度影响，解释依赖结构时用这一列）。", "",
              "| 数据 | 状态 | 周数 | 原始波动 | 创新方差 | 平均相关 | 吸收比率 | 共同越界（全样本阈值） | 倍数 | 共同越界（状态内阈值） | 倍数 | 等权方差比 | 最小方差组合风险比 |",
              "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for _, r in desc.iterrows():
        lines.append(f"| {r.dataset} | {r.regime} | {r.weeks} | {r.raw_vol:.2f} | {r.innov_var:.3f} | {r.mean_corr:.3f} | "
                     f"{r.absorption:.3f} | {r.co_exceed:.4f} | {r.co_exceed_ratio:.1f} | {r.co_exceed_within:.4f} | {r.co_exceed_within_ratio:.1f} | "
                     f"{r.ew_var_ratio:.3f} | {r.mv_risk_ratio:.3f} |")
    lines.append("")
lines += ["说明：第 4 节的量只用于解释检验结果的经济量级，不是另外的检验；最小方差组合用相关矩阵构造，只反映依赖结构的变化。"]
os.makedirs(RES, exist_ok=True)
open(os.path.join(RES, "报告_v3.md"), "w", encoding="utf-8").write("\n".join(lines) + "\n")
print("\n".join(lines))
