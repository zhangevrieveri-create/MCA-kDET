"""Amendment 2: report of the truncated time-shift (TTS) test.

    python report_tts.py

Reads results/pvalues_jobs_tts.csv and, if present, results/pvalues_jobs_tts_strata.csv,
results/summary_rotation.csv and results/summary_strata_main.csv; writes
results/summary_tts.csv and results/报告_v3_TTS.md.
"""
import os, json
import numpy as np, pandas as pd

ROOT = os.path.dirname(os.path.abspath(__file__)); RES = os.path.join(ROOT, "results")
A2 = json.load(open(os.path.join(ROOT, "protocol_amendment2.json"), encoding="utf-8"))
SHOW = ["MCA", "MCA_l1", "MCA_quad", "MCA_knn", "MCA_rf", "Sum", "DISCO", "MMD", "MST", "Schott"]


def holm(p):
    p = np.asarray(p, float); o = np.argsort(p); adj = np.empty_like(p); run = 0.0
    for k, i in enumerate(o):
        run = max(run, min(1.0, (len(p) - k) * p[i])); adj[i] = run
    return adj


def fmt(x):
    return "" if pd.isna(x) else ("<0.001" if x < 0.001 else f"{x:.3f}")


t = pd.read_csv(os.path.join(RES, "pvalues_jobs_tts.csv"))
t["variant"] = t.label_col.str.extract(r"^(\w+?)_d")[0]
t["delta"] = t.label_col.str.extract(r"_d(-?\d+)$")[0].astype(int)
r = int(t.delta.max()); assert t.delta.min() == -r
rows = []
for (var, rep, m), g in t.groupby(["variant", "representation", "method"]):
    if m not in SHOW:
        continue
    g = g.dropna(subset=["stat_obs"])
    obs = g[g.delta == 0].stat_obs
    if len(obs) != 1:
        continue
    t0 = obs.iloc[0]; ref = g.stat_obs.values
    B = int(np.sum(ref >= t0 - 1e-10))                       # counts delta = 0 itself
    rows.append(dict(variant=var, representation=rep, method=m, shifts=len(ref), complete=len(ref) == 2 * r + 1,
                     T_obs=t0, B=B, p_tts=min(1.0, B / (r + 1)),
                     p_naive=B / len(ref), ref_median=np.median(ref), ref_max=ref.max()))
s = pd.DataFrame(rows)
s["holm_MCA"] = np.nan
for var in s.variant.unique():
    pm = s[(s.variant == var) & (s.method == "MCA")].set_index("representation")
    if {"innov", "dir"} <= set(pm.index):
        for rep, a in zip(["innov", "dir"], holm(pm.loc[["innov", "dir"], "p_tts"].values)):
            s.loc[(s.variant == var) & (s.method == "MCA") & (s.representation == rep), "holm_MCA"] = a
s.to_csv(os.path.join(RES, "summary_tts.csv"), index=False)

out = ["# 第 5 章金融实证 v3：截断平移（TTS）检验（修订 2）", "",
       f"规则与时间见 `protocol_amendment2.json`。{A2['timing']}", "",
       f"截断半径 r = {r}，窗口 {2 * r + 1} 个平移（含 δ = 0）。p = min(1, B/(r+1))，B 为统计量不小于观测值的平移个数（含 δ = 0），最小可达 {1 / (r + 1):.4f}。"
       "在“状态序列严格平稳且与创新过程独立”的原假设下 P(p ≤ α) ≤ α（Yuan and Shou 2024）。", ""]
if not s.complete.all():
    out += ["**注意：有平移的结果不全，下列 p 值不是最终结果。**", ""]
for var in ["main", "lag2"]:
    sv = s[s.variant == var]
    if sv.empty:
        continue
    out += [f"## {'主分析（状态由上一周收盘 VIX 决定）' if var == 'main' else '敏感性分析：状态由两周前收盘 VIX 决定'}", "",
            "| 方法 | 创新 p_TTS | 创新 B | 方向 p_TTS | 方向 B |", "|---|---|---|---|---|"]
    for m in SHOW:
        a = sv[(sv.method == m) & (sv.representation == "innov")]; b = sv[(sv.method == m) & (sv.representation == "dir")]
        if a.empty and b.empty:
            continue
        out.append(f"| {m} | {fmt(a.p_tts.iloc[0]) if len(a) else ''} | {int(a.B.iloc[0]) if len(a) else ''} | "
                   f"{fmt(b.p_tts.iloc[0]) if len(b) else ''} | {int(b.B.iloc[0]) if len(b) else ''} |")
    pm = sv[sv.method == "MCA"].set_index("representation")
    if "holm_MCA" in pm and pm.holm_MCA.notna().all():
        out += ["", f"MCA-kDET 经 Holm 校正：创新 {fmt(pm.loc['innov', 'holm_MCA'])}，方向 {fmt(pm.loc['dir', 'holm_MCA'])}。", ""]
f = os.path.join(RES, "pvalues_jobs_tts_strata.csv")
if os.path.exists(f):
    st = pd.read_csv(f)
    out += ["## 敏感性分析：两周前 VIX，时期内置换（L = 52，B = 4999）", "", "| 方法 | 创新 | 方向 |", "|---|---|---|"]
    for m in SHOW:
        a = st[(st.method == m) & (st.representation == "innov")]; b = st[(st.method == m) & (st.representation == "dir")]
        if len(a) or len(b):
            out.append(f"| {m} | {fmt(a.p_value.iloc[0]) if len(a) else ''} | {fmt(b.p_value.iloc[0]) if len(b) else ''} |")
    out.append("")
f = os.path.join(RES, "summary_rotation.csv")
if os.path.exists(f):
    rot = pd.read_csv(f)
    out += ["## 对照：修订 1 的旋转检验（不变）", "", "| 方法 | 创新 | 方向 |", "|---|---|---|"]
    for m in SHOW:
        a = rot[(rot.method == m) & (rot.representation == "innov")]; b = rot[(rot.method == m) & (rot.representation == "dir")]
        if len(a) or len(b):
            out.append(f"| {m} | {fmt(a.p_rotation.iloc[0]) if len(a) else ''} | {fmt(b.p_rotation.iloc[0]) if len(b) else ''} |")
    out.append("")
open(os.path.join(RES, "报告_v3_TTS.md"), "w", encoding="utf-8").write("\n".join(out))
print("\n".join(out))
