"""Amendment 2 (protocol_amendment2.json): truncated time-shift (TTS) test.

    python make_jobs_v3c.py

Writes
  prepared/ind30/tts_labels.csv   one column per (variant, shift); 0 = week not used
  jobs_tts.csv                    TTS jobs, variants main and lag2, innov and dir
                                  (B = 19; only the observed statistics are used)
  jobs_tts_strata.csv             within-era permutation (L = 52, B = 4999) for lag2
No statistic is computed here.
"""
import os, json
import numpy as np, pandas as pd

ROOT = os.path.dirname(os.path.abspath(__file__))
A2 = json.load(open(os.path.join(ROOT, "protocol_amendment2.json"), encoding="utf-8"))
MIN_PER_GROUP = 5          # run_v3.R stops below this


def job(rows, kind, rep, label_col, B, seed, out, strata=0):
    rows.append(dict(job=len(rows) + 1, kind=kind, dataset="ind30", representation=rep,
                     x_file=f"prepared/ind30/X_{rep}.csv", label_file="prepared/ind30/tts_labels.csv",
                     label_col=label_col, role_file="prepared/ind30/weeks.csv", B=B,
                     seed=seed, out=out, competitors=1, strata=strata))


def tts_columns(y, valid_from, r):
    """y: regime sequence (length T, 0 where undefined); the usable weeks are
    valid_from, ..., T-1 (0-based). Window: valid_from + r, ..., T-1-r. Returns {delta: labels}."""
    T = len(y); lo, hi = valid_from + r, T - 1 - r
    assert hi - lo + 1 > 0
    cols = {}
    for d in range(-r, r + 1):
        lab = np.zeros(T, dtype=int)
        t = np.arange(lo, hi + 1)
        lab[t] = y[t + d]
        assert (lab[t] > 0).all()
        cols[d] = lab
    return cols


def main():
    wk = pd.read_csv(os.path.join(ROOT, "prepared", "ind30", "weeks.csv"))
    T = len(wk); reg = wk["regime"].values.astype(int); role = wk["role"].values
    r = int(np.floor(0.2 * T))
    assert r == 279 or T != 1395, "truncation radius differs from the protocol"
    lag2 = np.r_[0, reg[:-1]]                       # VIX at the close of week t-2
    variants = {"main": (reg, 0), "lag2": (lag2, 1)}
    labels = {"lag2_full": lag2}
    rows = []; check = []
    for vi, (name, (y, v0)) in enumerate(variants.items()):
        cols = tts_columns(y, v0, r)
        for d, lab in cols.items():
            labels[f"{name}_d{d}"] = lab
            cnt = [int(((lab == g) & (role == ro)).sum()) for g in (1, 2, 3) for ro in ("fit", "test")]
            check.append(dict(variant=name, delta=d, min_count=min(cnt), counts="/".join(map(str, cnt))))
        for ri, rep in enumerate(("innov", "dir")):
            seed = 9100000 + 10 * vi + ri             # same seed for every shift
            for d in range(-r, r + 1):
                job(rows, "tts", rep, f"{name}_d{d}", 19, seed, f"results/tts/{name}/ind30_{rep}_d{d}")
    ck = pd.DataFrame(check)
    bad = ck[ck.min_count < MIN_PER_GROUP]
    assert bad.empty, f"groups too small:\n{bad}"
    pd.DataFrame(labels).to_csv(os.path.join(ROOT, "prepared", "ind30", "tts_labels.csv"), index=False)
    pd.DataFrame(rows).to_csv(os.path.join(ROOT, "jobs_tts.csv"), index=False)
    ck.to_csv(os.path.join(ROOT, "prepared", "ind30", "tts_group_counts.csv"), index=False)
    srows = []
    for ri, rep in enumerate(("innov", "dir")):
        job(srows, "strata_lag2", rep, "lag2_full", 4999, 9200000 + ri, f"results/tts_strata/L52/ind30_{rep}_lag2", strata=52)
    pd.DataFrame(srows).to_csv(os.path.join(ROOT, "jobs_tts_strata.csv"), index=False)
    print(f"T = {T}, r = {r}, window = {T - 2 * r} weeks, shifts per series = {2 * r + 1}")
    print("TTS jobs:", len(rows), " within-era lag2 jobs:", len(srows))
    print("smallest fit/test count of a regime over all shifts:",
          ck.groupby("variant").min_count.min().to_dict())


if __name__ == "__main__":
    main()
