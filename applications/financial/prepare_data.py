"""Freeze the financial study and construct calendar-aligned, suspension-aware weeks."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent


def specifications():
    common = dict(dimension=10, representation="filtered", split="blocks", block=26,
                  buffer=4, thin=1, offset=0, regime="baseline_terciles", seed=20261027,
                  permutations=4999, role="sensitivity")
    specs = [dict(common, name="primary", role="primary"), dict(common, name="dimension_20", dimension=20),
             dict(common, name="dimension_30", dimension=30), dict(common, name="block_13", block=13),
             dict(common, name="block_52", block=52), dict(common, name="buffer_8", buffer=8),
             dict(common, name="raw_blocks", representation="raw"), dict(common, name="filtered_random", split="random"),
             dict(common, name="vix_15_25", regime="fixed_15_25")]
    specs.extend(dict(common, name=f"thin4_offset{i}", thin=4, offset=i) for i in range(4))
    specs.extend(dict(common, name=f"split_{i:02d}", seed=20261027 + i, role="split_stability") for i in range(1, 6))
    return specs


def main():
    raw = ROOT / "raw"; output = ROOT / "prepared"; output.mkdir(exist_ok=True)
    plan = dict(version="Ashare_VIX_weekly_v1", analysis_week_endpoints=["2015-01-09", "2025-12-26"],
                baseline=["2010-01-01", "2014-12-31"], primary_dimension=10, alpha=.05,
                panel_rule="2014-12-31 CSI300 weight, baseline-only data quality screen; nested 10/20/30 stocks.",
                returns="100 times log ratio of post-adjusted last-mainland-trading-day prices in consecutive Friday calendar weeks. Both endpoints must be observed and not suspended; never forward-fill missing weeks.",
                regimes="Previous Friday VIX close (last US trading day in that week); thresholds from 2010-2014 weekly VIX terciles. No current A-share return used in labels.",
                filtering="Per-stock GARCH(1,1) Gaussian QML parameters and mean estimated on pre-2015 weeks only, then one-step variance recursion from past observations. Missing weekly residuals propagated by their model conditional second moment.",
                calibration="4999 ordinary test-week label permutations, approximate only after filtering; frozen iid theorem not asserted for dependent weeks.",
                sensitivity_policy="All 18 configurations retained, including unfiltered/random splits and four prespecified calendar-thinning offsets. No smallest-p choice.",
                specifications=specifications())
    planfile = ROOT / "study_plan.json"
    if planfile.exists() and json.loads(planfile.read_text(encoding="utf-8")) != plan:
        raise ValueError("Frozen financial plan differs")
    planfile.write_text(json.dumps(plan, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    selection = pd.read_csv(raw / "panel_selection.csv").dropna(subset=["selected_rank"]).sort_values("selected_rank")
    stocks = selection.order_book_id.tolist()
    prices = pd.read_csv(raw / "daily_prices.csv.gz", parse_dates=["date"], float_precision="round_trip")
    assert not prices.duplicated(["order_book_id", "date"]).any()
    close = prices.pivot(index="date", columns="order_book_id", values="close")
    volume = prices.pivot(index="date", columns="order_book_id", values="volume")
    suspended = pd.read_csv(raw / "daily_suspended.csv.gz", index_col=0, parse_dates=True)
    calendar = close.index[close["000300.XSHG"].notna()]
    calendar = pd.Series(calendar, index=calendar)
    last_trading_day = calendar.resample("W-FRI").max()
    valid_weeks = last_trading_day.notna()
    weekly_close = pd.DataFrame(np.nan, index=last_trading_day.index, columns=stocks)
    weekly_valid = pd.DataFrame(False, index=last_trading_day.index, columns=stocks)
    endpoints = last_trading_day.loc[valid_weeks].to_numpy()
    weekly_close.loc[valid_weeks] = close.reindex(endpoints)[stocks].to_numpy()
    is_suspended = suspended.reindex(endpoints)[stocks]
    weekly_valid.loc[valid_weeks] = (close.reindex(endpoints)[stocks].gt(0).to_numpy()
                                    & volume.reindex(endpoints)[stocks].gt(0).to_numpy()
                                    & is_suspended.eq(False).to_numpy())
    weekly_return = 100 * np.log(weekly_close / weekly_close.shift(1))
    weekly_return = weekly_return.where(weekly_valid & weekly_valid.shift(1, fill_value=False))
    # Keep absent calendar weeks as NaN, so recurrence, buffers and thinning keep real spacing.
    vix = pd.read_csv(raw / "VIX_History.csv")
    vix.columns = vix.columns.str.strip()
    vix["DATE"] = pd.to_datetime(vix.DATE, format="%m/%d/%Y")
    assert not vix.DATE.duplicated().any()
    vix_weekly = vix.set_index("DATE").CLOSE.resample("W-FRI").last()
    baseline_vix = vix_weekly.loc["2010-01-01":"2014-12-31"]
    thresholds = baseline_vix.quantile([1 / 3, 2 / 3]).to_numpy()
    lagged_vix = vix_weekly.shift(1).reindex(weekly_return.index)
    assert lagged_vix.notna().all()
    meta = pd.DataFrame({"week_end": weekly_return.index, "last_cn_trading_day": last_trading_day.to_numpy(),
                         "lagged_vix": lagged_vix.to_numpy()})
    meta["regime"] = np.digitize(meta.lagged_vix, thresholds, right=True) + 1
    meta["regime_fixed"] = np.digitize(meta.lagged_vix, [15., 25.], right=True) + 1
    meta["period"] = np.where(meta.week_end.le("2014-12-31"), "baseline",
                              np.where(meta.week_end.between("2015-01-09", "2025-12-26"), "analysis", "warmup_or_partial"))
    meta["week_index"] = np.arange(len(meta))
    for dimension in (10, 20, 30):
        meta[f"valid_{dimension}"] = weekly_return.iloc[:, :dimension].notna().all(axis=1).to_numpy()
    meta.to_csv(output / "weeks.csv", index=False, encoding="utf-8-sig")
    weekly_return.rename_axis("week_end").to_csv(output / "weekly_returns.csv", float_format="%.17g")
    summary = []
    for dimension in (10, 20, 30):
        for period in ("baseline", "analysis"):
            mask = meta.period.eq(period)
            for regime in (1, 2, 3):
                chosen = mask & meta.regime.eq(regime)
                summary.append(dict(dimension=dimension, period=period, regime=regime,
                                    calendar_weeks=int(chosen.sum()), usable_weeks=int((chosen & meta[f"valid_{dimension}"]).sum())))
    summary = pd.DataFrame(summary)
    summary.to_csv(output / "data_counts.csv", index=False)
    audit = dict(stock_count=len(stocks), stocks_in_rank_order=stocks,
                 vix_baseline_tercile_thresholds=thresholds.tolist(),
                 vix_baseline_weeks=len(baseline_vix), no_inferential_statistics_computed=True,
                 source_manifest_sha256=hashlib.sha256((raw / "download_manifest.json").read_bytes()).hexdigest(),
                 missingness="Endpoint suspension/missing prices and full-China-holiday weeks remain missing. Reopening multiweek returns are not silently treated as one-week returns.")
    (output / "data_audit.json").write_text(json.dumps(audit, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"Frozen {len(plan['specifications'])} financial configurations; VIX thresholds={thresholds}")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
