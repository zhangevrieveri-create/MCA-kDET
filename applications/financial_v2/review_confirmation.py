"""Review a completed independent confirmation without inventing a validity pass."""
from pathlib import Path
import hashlib
import json
from datetime import datetime, timezone
import numpy as np
import pandas as pd
from scipy.stats import binomtest
from summarize_validation import QUESTIONS, NAMES, is_null, interval

ROOT = Path(__file__).resolve().parent
CASES = ["null_iid", "null_garch95", "null_garch99", "null_heterogeneous",
         "null_garch95_t5", "null_garch999", "alt_sector_shape", "alt_one_bank", "null_equal_mixtures"]


def holm(p):
    p = np.asarray(p, dtype=float)
    order = np.argsort(p, kind="stable")
    result = np.empty(len(p))
    result[order] = np.minimum(1, np.maximum.accumulate(p[order] * np.arange(len(p), 0, -1)))
    return result


def review():
    folder = ROOT / "sampling_fullfit_confirmation_d5" / "estimated"
    request = json.loads((folder.parent / "execution_request.json").read_text(encoding="utf-8"))
    frozen = json.loads((ROOT / "sampling_primary_candidate.json").read_text(encoding="utf-8"))
    assert request["first_rep"] == 101 and request["last_rep"] == 1100
    assert request["replications_per_case"] == 1000 and request["B"] == 999
    assert frozen["primary_stride"] == request["primary_stride"] == 4
    for name, expected in frozen["source_sha256"].items():
        assert hashlib.sha256((ROOT / name).read_bytes()).hexdigest() == expected, name
    files = sorted((folder / "chunks").glob("*.csv"))
    x = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
    keys = ["scenario", "rep", "block_length", "question"]
    assert len(files) == 900 and len(x) == 108000 and not x.duplicated(keys).any()
    assert set(x.scenario) == set(CASES) and set(x.block_length) == {1, 2, 4}
    assert set(x.question) == set(QUESTIONS) and set(x.B) == {999} and set(x.d) == {5}
    seeds = []
    for case, g in x.groupby("scenario"):
        assert set(g.rep) == set(range(101, 1101))
        assert g.groupby("rep").size().eq(12).all()
        case_id = CASES.index(case)
        expected = (((case_id * 1000000 + (g.rep.to_numpy()-1)*100) * 65537 + 987654321) % (2**31))
        assert np.array_equal(g.seed.to_numpy(), expected)
        seeds.extend(g.drop_duplicates("rep").seed.tolist())
    assert len(set(seeds)) == 9000
    assert x.population_test_n.gt(0).all() and np.isfinite(x[["p", "p_holm", "statistic", "critical95"]]).all().all()
    assert x.p.between(.001, 1).all() and np.allclose(x.p*1000, np.rint(x.p*1000), atol=1e-10)
    for _, g in x.groupby(["scenario", "rep", "block_length"]):
        assert np.allclose(holm(g.p), g.p_holm, atol=2e-14, rtol=0)
    # The final independent AUC audit must cover every case and all 12 cells.
    audit = pd.read_csv(folder / "independent_sampling_audit.csv")
    assert set(audit.scenario) == set(CASES) and len(audit) == 108
    assert audit.max_statistic_error.le(3e-13).all() and audit.records_and_labels_valid.all()
    rows = []
    main = x[x.block_length == 4]
    for (case, question), g in main.groupby(["scenario", "question"]):
        k = int((g.p <= .05+1e-12).sum()); lo, hi = interval(k, len(g))
        rows.append(dict(scenario=case, question=question, target="size" if is_null(case, question) else "power",
                         n=len(g), rejections=k, rate=k/len(g), ci_lower=lo, ci_upper=hi))
    rates = pd.DataFrame(rows)
    probes = []
    for r in rates[rates.target == "size"].itertuples():
        probes.append(dict(scenario=r.scenario, question=r.question, n=r.n, rejections=r.rejections,
                           rate=r.rate, ci_lower=r.ci_lower, ci_upper=r.ci_upper,
                           binomial_p_upper=binomtest(r.rejections, r.n, .05, alternative="greater").pvalue))
    for case, g in main.groupby("scenario"):
        true = g[g.question.map(lambda q: is_null(case, q))]
        rejected = true.groupby("rep").p_holm.min().le(.05+1e-12)
        k = int(rejected.sum()); lo, hi = interval(k, len(rejected))
        probes.append(dict(scenario=case, question="Holm_true_null_FWER", n=len(rejected), rejections=k,
                           rate=k/len(rejected), ci_lower=lo, ci_upper=hi,
                           binomial_p_upper=binomtest(k, len(rejected), .05, alternative="greater").pvalue))
    probes = pd.DataFrame(probes)
    probes["binomial_p_holm"] = holm(probes.binomial_p_upper)
    # Multiplicity here concerns Monte Carlo inflation diagnostics, not the four real-data tests.
    flags = probes[probes.binomial_p_holm <= .05]
    rates.to_csv(folder / "primary_confirmation_rates.csv", index=False)
    probes.to_csv(folder / "primary_inflation_diagnostics.csv", index=False)
    status = "NUMERICAL_INFLATION_DETECTED" if len(flags) else "NO_MULTIPLICITY_ADJUSTED_INFLATION_DETECTED_REQUIRES_ASSUMPTION_REVIEW"
    result = dict(updated_utc=datetime.now(timezone.utc).isoformat(), status=status,
                  independent_simulated_panels=9000, replications_per_case=1000, primary_stride=4,
                  B=999, source_fingerprints_match=True, seed_mapping_and_nonoverlap_verified=True,
                  Holm_recomputed=True, saved_AUC_audit_cells=len(audit),
                  multiplicity_adjusted_inflation_flags=len(flags),
                  unadjusted_CI_inflation_flags=int(probes.ci_lower.gt(.05).sum()),
                  min_test_n=int(main.population_test_n.min()), max_test_n=int(main.population_test_n.max()),
                  learner_warnings_max=int(x.learner_warning_count.max()),
                  nonconvergence_flags_max=int(x.fit_nonconvergence_flags.max()),
                  automatic_real_inference_authorized_by_this_script=False,
                  caveat="Finite DGP checks neither establish exchangeability of estimated innovations nor verify informative missingness or market stationarity.")
    (folder / "confirmation_review.json").write_text(json.dumps(result, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
    text = ["# 股票创新检验：独立确认审核", "", "主候选在查看真实p值前固定为日历间隔4；九场景各1000次，B=999，新编号101–1100。",
            "各场景分别统计，不合并成9000次同一原假设实验。", "",
            f"实现审核通过：9000个不同面板种子、108000行结果、四项Holm独立复算、{len(audit)}个保存问题/间隔的逐对AUC核对。",
            f"误报诊断：{len(probes)}个真原假设单项或家族的单侧二项检验，经Holm调整后有{len(flags)}格提示超过0.05。",
            "此处多重调整用于模拟误报诊断，不等于真实四项检验的Holm值。逐格Clopper–Pearson区间不是同时区间。",
            "没有检出超水平不等于证明有效，更不能证明实际残差满足iid、非信息性停牌或稳定分布条件。", ""]
    for case, g in rates.groupby("scenario", sort=False):
        values = []
        for q in QUESTIONS:
            r = g[g.question == q].iloc[0]
            values.append(f"{NAMES[q]}{'误报' if r.target == 'size' else '功效'}{r.rate:.3f}")
        text.append(f"- {case}："+"；".join(values)+"。")
    text += ["", "全部场景和三个间隔均保留在calibration_rates.csv；主候选明细见primary_confirmation_rates.csv，误报诊断见primary_inflation_diagnostics.csv。",
             "若误报有实质问题，候选不进入正式推断。若数值筛查未发现问题，仍须审核实际过滤误差、剩余序列依赖、缺失机制及功效。",
             "本脚本不计算真实股票p值、不换股票或间隔、不修改论文草稿。"]
    (folder / "独立确认审核.md").write_text("\n".join(text)+"\n", encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False))


if __name__ == "__main__":
    review()
