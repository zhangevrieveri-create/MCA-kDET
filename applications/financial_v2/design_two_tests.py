"""Build a result-blind stock-panel design and audit feasible return blocks.

No classifiers, hypothesis statistics or p values are computed here.
"""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import rqdatac as rq

ROOT = Path(__file__).resolve().parent
RAW = ROOT / "raw"
OUT = ROOT / "design_two_tests"
SECTORS = ["银行", "医药", "机械"]
SELECTION_DATE = "2017-12-29"
TRAIN_START = "2015-01-01"
TRAIN_END = "2017-12-31"
TEST_START = "2018-01-01"
TEST_END = "2024-12-31"


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    RAW.mkdir(exist_ok=True); OUT.mkdir(exist_ok=True)
    connected = False

    def connect():
        nonlocal connected
        if connected:
            return
        secret = os.environ.get("RQDATA_LICENSE")
        if not secret:
            raise SystemExit("Set RQDATA_LICENSE in the process environment for uncached data")
        try:
            rq.init(username="license", password=secret)
        except Exception as error:
            raise SystemExit("RQData connection failed: " + type(error).__name__)
        connected = True

    weights_path = RAW / "hs300_weights_20171229.csv"
    if not weights_path.exists():
        connect()
        rq.index_weights("000300.XSHG", date=SELECTION_DATE).rename("weight").rename_axis("order_book_id").to_csv(weights_path)
    weights = pd.read_csv(weights_path)
    industry_path = RAW / "hs300_industry_20171229.csv"
    if not industry_path.exists():
        connect()
        rq.get_instrument_industry(weights.order_book_id.tolist(), source="citics_2019", level=1,
                                   date=SELECTION_DATE).reset_index().to_csv(industry_path, index=False, encoding="utf-8-sig")
    industry = pd.read_csv(industry_path)
    universe = weights.merge(industry, on="order_book_id", validate="one_to_one")
    candidates = pd.concat([universe.loc[universe.first_industry_name.eq(sector)].sort_values(
        ["weight", "order_book_id"], ascending=[False, True]).head(10) for sector in SECTORS], ignore_index=True)
    assert not candidates.order_book_id.duplicated().any()
    ids = candidates.order_book_id.tolist()
    meta_path = RAW / "candidate_instruments_2017panel.csv"
    if not meta_path.exists():
        connect()
        fields = ["order_book_id", "symbol", "listed_date", "de_listed_date"]
        records = [{field: getattr(inst, field, None) for field in fields} for inst in rq.instruments(ids)]
        pd.DataFrame(records).to_csv(meta_path, index=False, encoding="utf-8-sig")
    metadata = pd.read_csv(meta_path)
    baseline_path = RAW / "candidate_baseline_2015_2017.csv.gz"
    if not baseline_path.exists():
        connect()
        rq.get_price(ids + ["000300.XSHG"], start_date="2014-12-15", end_date="2017-12-31",
                     frequency="1d", fields=["close", "volume"], adjust_type="post", skip_suspended=False,
                     expect_df=True).to_csv(baseline_path, float_format="%.17g")
    baseline = pd.read_csv(baseline_path, parse_dates=["date"], float_precision="round_trip")
    base_close = baseline.pivot(index="date", columns="order_book_id", values="close")
    base_volume = baseline.pivot(index="date", columns="order_book_id", values="volume")
    base_dates = base_close.index[(base_close["000300.XSHG"].notna()) & (base_close.index >= TRAIN_START)]
    base_close = base_close.reindex(base_dates); base_volume = base_volume.reindex(base_dates)
    candidates = candidates.merge(metadata, on="order_book_id", validate="one_to_one")
    candidates["baseline_close_coverage"] = candidates.order_book_id.map(base_close.notna().mean())
    candidates["baseline_traded_fraction"] = candidates.order_book_id.map(base_volume.gt(0).mean())
    candidates["eligible"] = (candidates.baseline_close_coverage.ge(.98)
                              & candidates.baseline_traded_fraction.ge(.95)
                              & pd.to_datetime(candidates.listed_date).lt(TRAIN_START))
    candidates["within_sector_rank"] = candidates.groupby("first_industry_name", sort=False).cumcount() + 1
    candidates["eligible_rank"] = np.nan
    for sector in SECTORS:
        eligible = candidates.index[candidates.first_industry_name.eq(sector) & candidates.eligible]
        candidates.loc[eligible, "eligible_rank"] = np.arange(1, len(eligible) + 1)
    candidates["main_selected"] = candidates.eligible_rank.le(5)
    candidates["reserve_8stock_panel"] = candidates.eligible_rank.le(8)
    counts = candidates.loc[candidates.main_selected].groupby("first_industry_name").size()
    candidates.to_csv(OUT / "selection_audit.csv", index=False, encoding="utf-8-sig")
    if len(counts) != len(SECTORS) or not counts.eq(5).all():
        raise ValueError("A prespecified sector has fewer than five baseline-eligible candidates; do not silently replace it")
    selected = candidates.loc[candidates.main_selected].copy()
    selected.to_csv(OUT / "main_panel.csv", index=False, encoding="utf-8-sig")
    all_selected = candidates.loc[candidates.reserve_8stock_panel, "order_book_id"].tolist()
    price_path = RAW / "two_tests_daily_2015_2024.csv.gz"
    if not price_path.exists():
        connect()
        rq.get_price(all_selected + ["000300.XSHG"], start_date="2014-12-15", end_date=TEST_END,
                     frequency="1d", fields=["close", "volume"], adjust_type="post", skip_suspended=False,
                     expect_df=True).to_csv(price_path, float_format="%.17g")
    suspended_path = RAW / "two_tests_suspended_2015_2024.csv.gz"
    if not suspended_path.exists():
        connect()
        rq.is_suspended(all_selected, start_date="2014-12-15", end_date=TEST_END).to_csv(suspended_path)
    data = pd.read_csv(price_path, parse_dates=["date"], float_precision="round_trip")
    close = data.pivot(index="date", columns="order_book_id", values="close")
    volume = data.pivot(index="date", columns="order_book_id", values="volume")
    suspended = pd.read_csv(suspended_path, index_col=0, parse_dates=True)
    dates = close.index[close["000300.XSHG"].notna()]
    close = close.reindex(dates); volume = volume.reindex(dates)
    suspended = suspended.reindex(index=dates, columns=all_selected)
    endpoint_valid = (close[all_selected].gt(0) & volume[all_selected].gt(0) & suspended.eq(False))
    returns = 100 * np.log(close[all_selected] / close[all_selected].shift(1))
    returns = returns.where(endpoint_valid & endpoint_valid.shift(1, fill_value=False))
    returns.rename_axis("date").to_csv(OUT / "daily_returns.csv.gz", float_format="%.17g",
                                       compression={"method": "gzip", "mtime": 0})
    trading_after_train = dates[(dates >= TEST_START) & (dates <= TEST_END)]
    gap = trading_after_train[:20]
    test_dates = trading_after_train[20:]
    periods = {"fit": dates[(dates >= TRAIN_START) & (dates <= TRAIN_END)], "test": test_dates}
    main_ids = selected.order_book_id.tolist()
    rows = []; block_records = []; group_records = []
    for dimension in [5, 20]:
        for split, period_dates in periods.items():
            nblocks = len(period_dates) // dimension
            for block in range(nblocks):
                block_dates = period_dates[block * dimension:(block + 1) * dimension]
                block_returns = returns.loc[block_dates, main_ids]
                complete = block_returns.notna().all(axis=0)
                valid_joint = bool(complete.all())
                block_records.append(dict(dimension=dimension, split=split, block=block,
                                          start=str(block_dates.min().date()), end=str(block_dates.max().date()),
                                          complete_all_main_stocks=valid_joint))
                for _, stock in selected.iterrows():
                    values = block_returns[stock.order_book_id].to_numpy()
                    rows.append(dict(dimension=dimension, split=split, block=block, stock=stock.order_book_id,
                                     industry=stock.first_industry_name, complete=bool(complete[stock.order_book_id]),
                                     complete_all_main_stocks=valid_joint,
                                     **{f"r{j + 1}": value for j, value in enumerate(values)}))
            grid = [r for r in block_records if r["dimension"] == dimension and r["split"] == split]
            joint = np.array([r["complete_all_main_stocks"] for r in grid], dtype=bool)
            longest = 0; current = 0
            for valid in joint:
                current = current + 1 if valid else 0
                longest = max(longest, current)
            group_records.append(dict(dimension=dimension, split=split, calendar_blocks=nblocks,
                                      common_complete_blocks=int(joint.sum()), missing_joint_blocks=int((~joint).sum()),
                                      longest_complete_run=longest, trailing_days_discarded=len(period_dates) % dimension))
    pd.DataFrame(block_records).to_csv(OUT / "block_calendar.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(rows).to_csv(OUT / "return_vectors.csv.gz", index=False, float_format="%.17g",
                              compression={"method": "gzip", "mtime": 0})
    pd.DataFrame(group_records).to_csv(OUT / "feasibility.csv", index=False, encoding="utf-8-sig")
    design = dict(version="two_tests_design_v1", no_classifiers_fitted=True, no_tests_computed=True,
                  sectors=SECTORS, stocks_per_sector=5, selection_date=SELECTION_DATE,
                  industry_source="citics_2019, historical-date query; retrospective taxonomy",
                  selection_rule="Top10 by pre-test CSI300 weight, baseline 2015-2017 coverage>=.98, positive volume>=.95, listed before 2015; first5 eligible each sector. No outcome/survival filter.",
                  fit_period=[TRAIN_START, TRAIN_END], test_planned=[TEST_START, TEST_END],
                  gap_trading_days=len(gap), gap_dates=[str(gap.min().date()), str(gap.max().date())],
                  test_first_trading_day=str(test_dates.min().date()), primary_dimension=5, sensitivity_dimension=20,
                  analysis_plan={
                      "status": "design specification; dependence calibration unvalidated; no real test results",
                      "main_questions": ["three within-sector equality tests, five stock populations each",
                                         "one between-sector equality test, three equal-weight five-stock mixtures"],
                      "library": ["l1", "quad", "knn", "rf"],
                      "alpha": .05, "primary_pvalue_adjustment": "Holm across the four main MCA tests",
                      "bootstrap_repetitions": 4999,
                      "main_bootstrap_trading_day_span": 40,
                      "bootstrap_span_sensitivity": [20, 60],
                      "bootstrap_feature_blocks": {"d5": [4, 8, 12], "d20": [1, 2, 3]},
                      "bootstrap_mask_handling": "Resample full-calendar panel and common-validity marker together; compute AUC on selected complete rows, keep original calendar gaps; fixed observed AUC scales.",
                      "mca_centering": "(AUC_boot_component - AUC_observed_component)/fixed_observed_scale, then maximum across populations and learners",
                      "distance": "Gaussian-kernel MMD, unsquared, with one pooled fitting-only scale and one fitting-only median bandwidth for all stock/sector pairs",
                      "contrast": "Equally average 75 cross-sector stock-pair distances minus equally average 30 within-sector stock-pair distances; same five stocks/sector.",
                      "distance_uncertainty": "Synchronously bootstrap centered RKHS mean-difference errors; common 95% radius across 105 stock pairs and three sector-mixture pairs; derive simultaneous distance bounds and conservative contrast interval.",
                      "different_inference_families": "Four MCA tests have Holm familywise .05; the distance family has simultaneous .95 coverage. Do not claim .95 coverage for the union of both families.",
                      "null_validation_replications": 1000,
                      "null_validation_bootstrap_repetitions": 999,
                      "validation_levels": [.01, .05, .10],
                      "pilot_rule": "First smoke/100-rep numerical pilot, then 1000-rep validation. Freeze changes from null simulations, never from real test p values.",
                      "null_scenarios": ["same-marginal Gaussian panels with asymmetric cross-stock correlation",
                                         "same-marginal autoregressive panels",
                                         "same-marginal GARCH persistence .95 and .99",
                                         "same-marginal heavy-tailed panels"],
                      "stress_scenarios": ["actual fixed calendar missingness mask",
                                           "common volatility breaks (outside stationary baseline assumptions)"],
                      "positive_controls": ["within-sector equality plus between-sector mean/scale/dependence shifts",
                                            "one genuinely deviating stock within a sector"],
                      "equivalence_claim": "Not included without an author-defined substantive discrepancy tolerance",
                      "comparators": "Do not transfer this AUC calibration mechanically to degenerate DISCO/MMD null statistics or MST graphs; their inferential calibrations require separate validation.",
                  },
                  missingness_policy="Keep original block indices and joint validity mask. No compression of gaps, no synthetic filling, no changing selected stocks from post-2017 missingness.",
                  inference_status="Calendar missingness and dependence calibration require validation before main inference. Equal mixtures refer to the selected fixed panel, not all sector stocks.",
                  sha256={p.name: digest(p) for p in [weights_path, industry_path, meta_path, baseline_path, price_path, suspended_path,
                                                     Path(__file__)]})
    (OUT / "design.json").write_text(json.dumps(design, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(selected[["first_industry_name", "order_book_id", "symbol", "weight"]].to_string(index=False))
    print(pd.DataFrame(group_records).to_string(index=False))
    print("Data-only design audit complete; no learned scores, test statistics or p values.")


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    main()
