"""Chapter 5, financial application (v3): data preparation.

Reads Kenneth French daily industry portfolios (value weighted) and the Cboe VIX
history, builds Friday-week log returns, ex-ante VIX regimes, GARCH(1,1)
innovations, their directions, and the block split fixed in protocol.json.
Writes prepared/<dataset>/{weeks.csv, X_<transform>.csv, garch_<fit>.csv, audit.json}.

    python prepare_v3.py                 # all datasets whose raw file is present
    python prepare_v3.py --dataset=ind30

No test statistic or p-value is computed here.
"""
import io, json, os, sys, zipfile, hashlib
import numpy as np, pandas as pd
from scipy.optimize import minimize
from scipy.signal import lfilter

ROOT = os.path.dirname(os.path.abspath(__file__))
P = json.load(open(os.path.join(ROOT, "protocol.json"), encoding="utf-8"))


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def find_raw(name):
    for d in P["raw_dirs"]:
        p = os.path.join(ROOT, d, name)
        if os.path.exists(p):
            return p
    return None


def read_french(path):
    """Value-weighted daily returns (%) from a French *_Daily_CSV.zip or .csv."""
    if path.endswith(".zip"):
        with zipfile.ZipFile(path) as z:
            txt = z.read(z.namelist()[0]).decode("latin-1")
    else:
        txt = open(path, encoding="latin-1").read()
    lines = txt.splitlines()
    start = [i for i, l in enumerate(lines) if "Average Value Weighted Returns -- Daily" in l][0]
    rows = []
    for l in lines[start + 1:]:
        if not l.strip():
            break
        rows.append(l)
    df = pd.read_csv(io.StringIO("\n".join(rows)), index_col=0)
    df.index = pd.to_datetime(df.index.astype(str).str.strip(), format="%Y%m%d")
    df.columns = [c.strip() for c in df.columns]
    return df.replace([-99.99, -999], np.nan).astype(float)


# ---------------------------------------------------------------- GARCH(1,1), Gaussian QML
def garch_var(xc, omega, alpha, beta, v0):
    """h_t = omega + alpha x_{t-1}^2 + beta h_{t-1}, with x_0^2 = h_0 = v0 (as prepare_filter.R)."""
    u = np.empty_like(xc)
    u[0] = omega + (alpha + beta) * v0
    u[1:] = omega + alpha * xc[:-1] ** 2
    # h_t - beta h_{t-1} = u_t for t >= 2, h_1 = u_1
    h, _ = lfilter([1.0], [1.0, -beta], u[1:], zi=[beta * u[0]])
    return np.concatenate([[u[0]], h])


def decode(th):
    pers = 1 / (1 + np.exp(-th[1]))
    alpha = pers / (1 + np.exp(-th[2]))
    return np.exp(th[0]), alpha, pers - alpha


def fit_garch(x):
    mu = float(np.mean(x)); xc = x - mu; v0 = float(np.mean(xc ** 2))

    def nll(th):
        o, a, b = decode(th)
        if not np.all(np.isfinite([o, a, b])):
            return 1e100
        h = garch_var(xc, o, a, b, v0)
        if not np.all(np.isfinite(h)) or np.any(h <= 0):
            return 1e100
        return 0.5 * np.sum(np.log(h) + xc ** 2 / h)

    best = None
    for p in (0.85, 0.95, 0.985):
        th0 = np.array([np.log((1 - p) * v0), np.log(p / (1 - p)), np.log((0.1 / p) / (1 - 0.1 / p))])
        r = minimize(nll, th0, method="Nelder-Mead",
                     options={"maxiter": 4000, "xatol": 1e-8, "fatol": 1e-10})
        if np.isfinite(r.fun) and (best is None or r.fun < best.fun):
            best = r
    o, a, b = decode(best.x)
    return dict(mean=mu, omega=o, alpha=a, beta=b, persistence=a + b, v0=v0,
                nll=float(best.fun), converged=bool(best.success))


def apply_garch(x, par):
    xc = x - par["mean"]
    h = garch_var(xc, par["omega"], par["alpha"], par["beta"], par["v0"])
    return xc / np.sqrt(h), h


# ---------------------------------------------------------------- split
def block_split(regime, seed, block, buffer, minobs):
    T = len(regime); rng = np.random.default_rng(seed)
    bid = np.arange(T) // block; nb = bid.max() + 1
    for attempt in range(1000):
        fitb = rng.choice(nb, size=int(round(nb / 2)), replace=False)
        role = np.where(np.isin(bid, fitb), "fit", "test").astype(object)
        isfit = role == "fit"
        for t in np.where(~isfit)[0]:
            lo, hi = max(0, t - buffer), min(T, t + buffer + 1)
            if isfit[lo:hi].any():
                role[t] = "buffer"
        ok = all(((regime == g) & (role == r)).sum() >= minobs for g in (1, 2, 3) for r in ("fit", "test"))
        if ok:
            return role, attempt + 1
    raise RuntimeError("no admissible block split")


def weekly(df):
    lr = 100 * np.log1p(df / 100)
    wk = lr.resample("W-FRI").sum(min_count=1)
    nd = lr.iloc[:, 0].resample("W-FRI").count()
    last = df.index.max()
    keep = wk.index <= last          # drop the trailing incomplete week
    return wk[keep], nd[keep]


def main(only=None):
    vix_path = find_raw(P["vix_file"])
    v = pd.read_csv(vix_path, parse_dates=["DATE"]).set_index("DATE")["CLOSE"]
    vw = v.resample("W-FRI").last()
    pre = vw[P["vix_threshold_period"][0]:P["vix_threshold_period"][1]]
    thr = [float(pre.quantile(1 / 3)), float(pre.quantile(2 / 3))]
    lagged = vw.shift(1)
    for ds, spec in P["datasets"].items():
        if only and ds != only:
            continue
        path = find_raw(spec["file"])
        if path is None:
            print(f"[{ds}] raw file {spec['file']} not found in {P['raw_dirs']}; skipped")
            continue
        df = read_french(path)
        wk, ndays = weekly(df)
        a0, a1 = P["analysis_period"]; b0, b1 = P["presample_period"]
        A = wk[a0:a1]; Bp = wk[b0:b1]
        assert A.notna().all().all() and Bp.notna().all().all(), "missing industry returns"
        lv = lagged.reindex(A.index)
        assert lv.notna().all(), "missing lagged VIX"
        reg = np.where(lv <= thr[0], 1, np.where(lv <= thr[1], 2, 3))
        alt = P["vix_alt_thresholds"]
        reg_alt = np.where(lv <= alt[0], 1, np.where(lv <= alt[1], 2, 3))
        role, attempts = block_split(reg, P["split"]["seed"], P["split"]["block"], P["split"]["buffer"], P["split"]["min_per_group"])
        out = os.path.join(ROOT, "prepared", ds); os.makedirs(out, exist_ok=True)
        weeks = pd.DataFrame({"week_end": A.index.strftime("%Y-%m-%d"), "lagged_vix": lv.values,
                              "regime": reg, "regime_alt": reg_alt, "role": role,
                              "trading_days": ndays.reindex(A.index).values})
        weeks.to_csv(os.path.join(out, "weeks.csv"), index=False)
        X = A.values
        Z = np.empty_like(X); Zp = np.empty_like(X); pars = []; parsp = []
        for j, c in enumerate(A.columns):
            pf = fit_garch(X[:, j]); Z[:, j], _ = apply_garch(X[:, j], pf)
            pp = fit_garch(Bp[c].values)
            full = np.concatenate([Bp[c].values, X[:, j]]); zz, _ = apply_garch(full, pp)
            Zp[:, j] = zz[len(Bp):]
            pars.append(dict(industry=c, **pf)); parsp.append(dict(industry=c, **pp))
        pd.DataFrame(pars).to_csv(os.path.join(out, "garch_full.csv"), index=False)
        pd.DataFrame(parsp).to_csv(os.path.join(out, "garch_presample.csv"), index=False)
        reps = {"raw": X, "innov": Z, "dir": Z / np.linalg.norm(Z, axis=1, keepdims=True),
                "innov_pre": Zp, "dir_pre": Zp / np.linalg.norm(Zp, axis=1, keepdims=True)}
        for name, M in reps.items():
            pd.DataFrame(M, columns=A.columns).to_csv(os.path.join(out, f"X_{name}.csv"),
                                                      index=False, float_format="%.17g")
        cnt = {f"regime{g}_{r}": int(((reg == g) & (role == r)).sum()) for g in (1, 2, 3) for r in ("fit", "test", "buffer")}
        audit = dict(dataset=ds, raw_file=os.path.basename(path), raw_sha256=sha256(path),
                     vix_file=os.path.basename(vix_path), vix_sha256=sha256(vix_path),
                     first_week=weeks.week_end.iloc[0], last_week=weeks.week_end.iloc[-1],
                     weeks=len(weeks), d=X.shape[1], vix_thresholds=thr,
                     regime_counts={str(g): int((reg == g).sum()) for g in (1, 2, 3)},
                     split_attempts=attempts, counts=cnt,
                     garch_full_converged=int(sum(p["converged"] for p in pars)),
                     garch_presample_converged=int(sum(p["converged"] for p in parsp)),
                     persistence_full=[round(p["persistence"], 4) for p in pars])
        json.dump(audit, open(os.path.join(out, "audit.json"), "w"), indent=1)
        print(f"[{ds}] weeks={len(weeks)} d={X.shape[1]} thresholds={np.round(thr, 3)} "
              f"regimes={audit['regime_counts']} test={[cnt[f'regime{g}_test'] for g in (1,2,3)]} "
              f"fit={[cnt[f'regime{g}_fit'] for g in (1,2,3)]}")


if __name__ == "__main__":
    only = None
    for a in sys.argv[1:]:
        if a.startswith("--dataset="):
            only = a.split("=", 1)[1]
    main(only)
