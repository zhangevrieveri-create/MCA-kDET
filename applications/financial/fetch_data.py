"""Fetch licensed A-share data and public Cboe VIX; credentials only via environment."""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import requests
import rqdatac as rq

ROOT = Path(__file__).resolve().parent
RAW = ROOT / "raw"
BASELINE_DATE = "2014-12-31"
START = "2009-12-28"
END = "2025-12-31"
VIX_URL = "https://cdn-api.cboe.com/api/global/us_indices/daily_prices/VIX_History.csv"


def main():
    RAW.mkdir(parents=True, exist_ok=True)
    required = ["daily_prices.csv.gz", "daily_suspended.csv.gz", "historical_weights.csv",
                "panel_selection.csv", "instruments.csv", "VIX_History.csv"]
    if all((RAW / name).exists() for name in required):
        print("All raw input files cached; no API login or download required.")
        return
    secret = os.environ.get("RQDATA_LICENSE")
    if not secret:
        raise SystemExit("Set RQDATA_LICENSE in the process environment; never put the license in source files.")
    try:
        rq.init(username="license", password=secret)
    except Exception as error:
        # Avoid accidental credential echo from third-party exception text.
        raise SystemExit("RQData authentication/connection failed: " + type(error).__name__)
    print("RQData connected.", flush=True)
    if (RAW / "historical_weights.csv").exists():
        weights = pd.read_csv(RAW / "historical_weights.csv").set_index("order_book_id").weight
    else:
        weights = rq.index_weights("000300.XSHG", date=BASELINE_DATE).sort_values(ascending=False)
        weights.rename("weight").rename_axis("order_book_id").to_csv(RAW / "historical_weights.csv")
    candidates = weights.sort_values(ascending=False).head(60).index.tolist()
    if not (RAW / "baseline_prices.csv.gz").exists():
        baseline = rq.get_price(candidates, start_date=START, end_date=BASELINE_DATE,
                                frequency="1d", fields=["close", "volume"], adjust_type="post", skip_suspended=False)
        baseline.to_csv(RAW / "baseline_prices.csv.gz", float_format="%.17g")
    else:
        baseline = pd.read_csv(RAW / "baseline_prices.csv.gz", index_col=[0, 1], parse_dates=[1])
    base_close = baseline.close.unstack("order_book_id")
    base_volume = baseline.volume.unstack("order_book_id")
    base_close = base_close.loc[base_close.index >= "2010-01-01"]
    base_volume = base_volume.reindex(base_close.index)
    selection = pd.DataFrame({"weight_20141231": weights.reindex(candidates),
                              "baseline_close_coverage": base_close.notna().mean(),
                              "baseline_positive_volume_fraction": base_volume.gt(0).mean()}).loc[candidates]
    selection["baseline_eligible"] = (selection.baseline_close_coverage >= .98) & (selection.baseline_positive_volume_fraction >= .95)
    selected = selection.loc[selection.baseline_eligible].head(30).index.tolist()
    if len(selected) != 30:
        raise ValueError(f"Only {len(selected)} candidates passed the prespecified baseline data check; do not silently change rule")
    selection["selected_rank"] = pd.Series({stock: rank for rank, stock in enumerate(selected, 1)})
    selection.rename_axis("order_book_id").to_csv(RAW / "panel_selection.csv")
    instruments = rq.instruments(selected)
    fields = ["order_book_id", "symbol", "listed_date", "de_listed_date", "sector_code", "sector_code_name", "industry_name"]
    pd.DataFrame([{f: getattr(inst, f, None) for f in fields} for inst in instruments]).to_csv(RAW / "instruments.csv", index=False, encoding="utf-8-sig")
    print(f"Frozen 30-stock panel using {BASELINE_DATE} weights and pre-2015 coverage only.", flush=True)
    if not (RAW / "daily_prices.csv.gz").exists():
        prices = rq.get_price(selected + ["000300.XSHG"], start_date=START, end_date=END,
                              frequency="1d", fields=["close", "volume"], adjust_type="post", skip_suspended=False)
        prices.to_csv(RAW / "daily_prices.csv.gz", float_format="%.17g")
        print(f"Saved {len(prices)} daily price rows.", flush=True)
    if not (RAW / "daily_suspended.csv.gz").exists():
        rq.is_suspended(selected, start_date=START, end_date=END).to_csv(RAW / "daily_suspended.csv.gz")
    if not (RAW / "VIX_History.csv").exists():
        response = requests.get(VIX_URL, timeout=45)
        response.raise_for_status()
        if "DATE" not in response.text[:100] or "CLOSE" not in response.text[:100]:
            raise ValueError("Unexpected VIX download schema")
        (RAW / "VIX_History.csv").write_bytes(response.content)
    manifest = {
        "provider": "RQData user-authorized API and Cboe public VIX", "rqdatac_version": rq.__version__,
        "stock_interval": [START, END], "index_weight_date": BASELINE_DATE,
        "price_arguments": {"frequency": "1d", "fields": ["close", "volume"], "adjust_type": "post", "skip_suspended": False},
        "panel_rule": "First 30 by 2014-12-31 CSI300 weight among top60, with 2010-2014 close coverage>=.98 and positive-volume fraction>=.95. No post-2014 outcome or survival screen.",
        "vix_url": VIX_URL,
        "files": [{"name": p.name, "bytes": p.stat().st_size, "sha256": hashlib.sha256(p.read_bytes()).hexdigest()}
                  for p in sorted(RAW.iterdir()) if p.is_file() and p.name != "download_manifest.json"]}
    (RAW / "download_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print("All requested source files saved with SHA256 manifest.", flush=True)


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    main()
