"""Audit domestic implied-volatility sources and historical industry metadata.

This does not fit learners or compute p values. Credentials are process-local.
"""
import hashlib
import json
import os
from pathlib import Path
import sys

import pandas as pd
import rqdatac as rq

ROOT = Path(__file__).resolve().parent
RAW = ROOT / "raw"
SOURCES = {
    "VX0004.RI": "米筐50ETF波动率指数（数据商编制，国内期权隐含波动率）",
    "VX0001.RI": "米筐沪深300波动率指数（数据商编制，国内期权隐含波动率）",
    "000188.XSHG": "上证50ETF波动率指数（交易所中国波指历史序列，覆盖待查询）",
}


def main():
    RAW.mkdir(parents=True, exist_ok=True)
    secret = os.environ.get("RQDATA_LICENSE")
    if not secret:
        raise SystemExit("RQDATA_LICENSE required in process environment; do not save it in source")
    try:
        rq.init(username="license", password=secret)
    except Exception as error:
        raise SystemExit("RQData connection failed: " + type(error).__name__)
    records = []
    for code, description in SOURCES.items():
        file = RAW / (code.replace(".", "_") + ".csv")
        try:
            if file.exists():
                data = pd.read_csv(file)
            else:
                response = rq.get_price(code, start_date="2015-01-01", end_date="2025-12-31",
                                        frequency="1d", fields=["close"], expect_df=True)
                if response is None or len(response) == 0:
                    raise ValueError("No historical rows")
                data = response.reset_index()
                data.to_csv(file, index=False, float_format="%.17g", encoding="utf-8-sig")
            dates = pd.to_datetime(data["date"])
            record = dict(code=code, description=description, available=True, rows=len(data),
                          start=str(dates.min().date()), end=str(dates.max().date()),
                          missing_close=int(data.close.isna().sum()), sha256=hashlib.sha256(file.read_bytes()).hexdigest())
            print(f"{code}: {record['rows']} rows, {record['start']} to {record['end']}", flush=True)
        except Exception as error:
            # Some SDKs echo request context; never persist auth exception text.
            record = dict(code=code, description=description, available=False, error_type=type(error).__name__)
            print(f"{code}: unavailable ({type(error).__name__})", flush=True)
        records.append(record)
    weights = pd.read_csv(ROOT.parent / "financial/raw/historical_weights.csv")
    ids = weights.order_book_id.tolist()
    industry_file = RAW / "hs300_industry_20141231.csv"
    try:
        if industry_file.exists():
            industry = pd.read_csv(industry_file)
        else:
            industry = rq.get_instrument_industry(ids, source="citics_2019", level=1, date="2014-12-31").reset_index()
            industry.to_csv(industry_file, index=False, encoding="utf-8-sig")
        if "first_industry_name" in industry:
            counts = industry.groupby("first_industry_name").size().sort_values(ascending=False)
            counts.rename("stocks").to_csv(RAW / "historical_industry_counts.csv", encoding="utf-8-sig")
            print(counts.to_string(), flush=True)
        industry_record = dict(available=True, rows=len(industry), source="citics_2019", query_date="2014-12-31",
                               note="Historical-date membership queried in a later taxonomy; do not claim the 2019 taxonomy itself existed in 2014.")
    except Exception as error:
        industry_record = dict(available=False, error_type=type(error).__name__)
        print("Historical industry query unavailable: " + type(error).__name__, flush=True)
    manifest = dict(volatility_sources=records, industry=industry_record,
                    no_tests_computed=True, no_us_vix_fallback=True,
                    documentation="https://www.ricequant.com/doc/rqdata/python/ricequant-index")
    (RAW / "source_audit.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    main()
