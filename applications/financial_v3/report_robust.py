"""Amendment 1: dependence-robust results.

    python report_robust.py

Reads results/pvalues_jobs_rotation.csv, results/pvalues_jobs_strata_calib.csv and
results/pvalues_jobs_strata_main.csv (whichever exist) and writes
results/报告_v3_依赖稳健.md and results/summary_rotation.csv, summary_strata_*.csv.
"""
import os, json
import numpy as np, pandas as pd
from scipy.stats import beta as Beta

ROOT = os.path.dirname(os.path.abspath(__file__)); RES = os.path.join(ROOT, "results")
A = json.load(open(os.path.join(ROOT, "protocol_amendment1.json"), encoding="utf-8"))
SHOW = ["MCA", "MCA_l1", "MCA_quad", "MCA_knn", "MCA_rf", "Sum", "DISCO", "MMD", "MST", "Schott"]
out = ["# 第 5 章金融实证 v3：依赖稳健的推断（修订 1）", "",
       f"修订原因与时间见 `protocol_amendment1.json`：{A['timing']}", ""]


def cp(k, n):
    return (Beta.ppf(.025, k, n - k + 1) if k else 0.0, Beta.ppf(.975, k + 1, n - k) if k < n else 1.0)


def holm(p):
    p = np.asarray(p, float); o = np.argsort(p); adj = np.empty_like(p); run = 0.0
    for r, i in enumerate(o):
        run = max(run, min(1.0, (len(p) - r) * p[i])); adj[i] = run
    return adj


# ---------------------------------------------------------------- rotation test
f = os.path.join(RES, "pvalues_jobs_rotation.csv")
if os.path.exists(f):
    r = pd.read_csv(f)
    rows = []
    for (rep, m), g in r.groupby(["representation", "method"]):
        if m not in SHOW:
            continue
        obs = g[g.label_col == "s0"]["stat_obs"]
        ref = g[g.label_col != "s0"]["stat_obs"].dropna().values
        if len(obs) != 1 or not np.isfinite(obs.iloc[0]):
            continue
        t0 = obs.iloc[0]
        rows.append(dict(representation=rep, method=m, T_obs=t0, shifts=len(ref),
                         p_rotation=(1 + np.sum(ref >= t0 - 1e-10)) / (1 + len(ref)),
                         ref_median=np.median(ref), ref_q95=np.quantile(ref, .95), ref_q99=np.quantile(ref, .99),
                         ref_max=ref.max()))
    rot = pd.DataFrame(rows)
    if len(rot):
        prim = rot[rot.method == "MCA"].set_index("representation")
        if {"innov", "dir"} <= set(prim.index):
            adj = holm(prim.loc[["innov", "dir"], "p_rotation"].values)
            rot["holm_MCA"] = np.nan
            for rep, a in zip(["innov", "dir"], adj):
                rot.loc[(rot.method == "MCA") & (rot.representation == rep), "holm_MCA"] = a
        rot.to_csv(os.path.join(RES, "summary_rotation.csv"), index=False)
        n = int(rot.shifts.max())
        out += ["## 1 主推断：旋转（循环平移）检验", "",
                f"参照分布：真实 VIX 状态序列循环平移 52 到 T−52 周的全部 {n} 种平移，每种平移都重跑整个流程（拆分固定，重新训练分类器）。"
                f"p 值 = (1 + 参照统计量 ≥ 观测统计量的个数) / (1 + {n})，最小为 {1/(1+n):.4f}。", "",
                "| 表示 | 方法 | 观测统计量 | 参照中位数 | 参照 95% 分位 | 参照 99% 分位 | 参照最大值 | 旋转 p 值 | Holm（MCA） |",
                "|---|---|---|---|---|---|---|---|---|"]
        rot["_o"] = rot.method.map({m: i for i, m in enumerate(SHOW)})
        for _, x in rot.sort_values(["representation", "_o"]).iterrows():
            h = "" if pd.isna(x.get("holm_MCA", np.nan)) else f"{x.holm_MCA:.4f}"
            out.append(f"| {x.representation} | {x.method} | {x.T_obs:.3f} | {x.ref_median:.3f} | {x.ref_q95:.3f} | "
                       f"{x.ref_q99:.3f} | {x.ref_max:.3f} | {x.p_rotation:.4f} | {h} |")
        out.append("")

# ---------------------------------------------------------------- within-era permutation: placebo
f = os.path.join(RES, "pvalues_jobs_strata_calib.csv")
cal = None
if os.path.exists(f):
    c = pd.read_csv(f); rows = []
    for (L, rep, m), g in c.groupby(["strata", "representation", "method"]):
        if m not in SHOW:
            continue
        n = len(g); k = int((g.p_value <= .05).sum()); lo, hi = cp(k, n)
        rows.append(dict(L=L, representation=rep, method=m, n=n, rate=k / n, ci_low=lo, ci_high=hi, calibrated=lo <= .05 <= hi))
    cal = pd.DataFrame(rows); cal.to_csv(os.path.join(RES, "summary_strata_placebo.csv"), index=False)
    out += ["## 2 次要推断：时期内置换的校准检查（原 200 个安慰剂状态序列，B = 199）", "",
            "| L（周） | 表示 | 方法 | 次数 | 拒绝率 | 95% 区间 | 区间含 0.05 |", "|---|---|---|---|---|---|---|"]
    for _, x in cal.iterrows():
        out.append(f"| {x.L} | {x.representation} | {x.method} | {x.n} | {x.rate:.3f} | [{x.ci_low:.3f}, {x.ci_high:.3f}] | {'是' if x.calibrated else '否'} |")
    out.append("")

# ---------------------------------------------------------------- within-era permutation: real labels
f = os.path.join(RES, "pvalues_jobs_strata_main.csv")
if os.path.exists(f):
    mres = pd.read_csv(f)
    mres["setting"] = mres.dataset + " / " + mres.representation + " / L = " + mres.strata.astype(str)
    tab = mres.pivot_table(index="setting", columns="method", values="p_value", aggfunc="first")
    tab = tab[[m for m in SHOW if m in tab.columns]]; tab.to_csv(os.path.join(RES, "summary_strata_main.csv"))
    out += ["## 3 次要推断：时期内置换下的真实 p 值（B = 4999）", "",
            "| 设定 | " + " | ".join(tab.columns) + " |", "|---" * (len(tab.columns) + 1) + "|"]
    for s, x in tab.iterrows():
        flag = ""
        if cal is not None:
            L = int(s.split("= ")[1]); rep = s.split(" / ")[1]
            ok = cal[(cal.L == L) & (cal.representation == rep) & (cal.method == "MCA")]
            if len(ok):
                flag = "" if bool(ok.calibrated.iloc[0]) else "（该 L 的 MCA 安慰剂检查未通过）"
        out.append(f"| {s}{flag} | " + " | ".join("" if pd.isna(x[m]) else f"{x[m]:.4f}" for m in tab.columns) + " |")
    out.append("")
open(os.path.join(RES, "报告_v3_依赖稳健.md"), "w", encoding="utf-8").write("\n".join(out) + "\n")
print("\n".join(out))
