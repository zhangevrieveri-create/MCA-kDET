"""Builds the job lists of the v3 financial application from protocol.json.

    python make_jobs_v3.py [--calib-dataset=ind30]

Writes jobs_main.csv, jobs_calibration.csv and, for the calibration data set,
prepared/<ds>/placebo_labels.csv and prepared/<ds>/boot/X_{innov,dir}_rNNN.csv.gz.
"""
import json, os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from prepare_v3 import fit_garch, apply_garch, garch_var, P, ROOT

from concurrent.futures import ProcessPoolExecutor


def boot_files(calib, r):
    return {rep: f"prepared/{calib}/boot/X_{rep}_r{r:03d}.csv.gz" for rep in ("innov", "dir")}


def gen_boot(calib, r, burn=200):
    Z = pd.read_csv(os.path.join(ROOT, "prepared", calib, "X_innov.csv")).values
    par = pd.read_csv(os.path.join(ROOT, "prepared", calib, "garch_full.csv"))
    T, d = Z.shape
    rng = np.random.default_rng([20261030, r])
    e = Z[rng.integers(0, T, size=T + burn)]
    X = np.empty((T, d))
    for j in range(d):
        p = par.iloc[j]; h = p["v0"]; x2 = p["v0"]; xs = np.empty(T + burn)
        om, al, be = p["omega"], p["alpha"], p["beta"]
        for t in range(T + burn):
            h = om + al * x2 + be * h
            xs[t] = np.sqrt(h) * e[t, j]; x2 = xs[t] ** 2
        X[:, j] = p["mean"] + xs[burn:]
    Zs = np.column_stack([apply_garch(X[:, j], fit_garch(X[:, j]))[0] for j in range(d)])
    Ds = Zs / np.linalg.norm(Zs, axis=1, keepdims=True)
    for rep, M in (("innov", Zs), ("dir", Ds)):
        pd.DataFrame(M).to_csv(os.path.join(ROOT, boot_files(calib, r)[rep]), index=False,
                               float_format="%.12g", compression="gzip")
    return r


def main():
    global jobs
    calib = "ind30"
    for a in sys.argv[1:]:
        if a.startswith("--calib-dataset="):
            calib = a.split("=", 1)[1]

    B = P["permutations"]; jobs = []
    seed0 = 7100000


    def add(ds, rep, label_col, kind, B, out, label_file="weeks.csv", x_file=None, competitors=1):
        jobs.append(dict(job=len(jobs) + 1, kind=kind, dataset=ds, representation=rep,
                         x_file=x_file or f"prepared/{ds}/X_{rep}.csv",
                         label_file=f"prepared/{ds}/{label_file}", label_col=label_col,
                         role_file=f"prepared/{ds}/weeks.csv", B=B, seed=seed0 + len(jobs) + 1,
                         out=out, competitors=competitors))


    avail = [ds for ds in P["datasets"] if os.path.exists(os.path.join(ROOT, "prepared", ds, "weeks.csv"))]
    print("prepared data sets:", avail)
    for ds in avail:
        for rep in ("raw", "innov", "dir"):
            add(ds, rep, "regime", "main", B, f"results/main/{ds}_{rep}")
    if "ind30" in avail:
        for rep in ("innov_pre", "dir_pre"):
            add("ind30", rep, "regime", "main", B, f"results/main/ind30_{rep}")
        for rep in ("innov", "dir"):
            add("ind30", rep, "regime_alt", "main", B, f"results/main/ind30_{rep}_vix15-25")
    main = pd.DataFrame(jobs); main.to_csv(os.path.join(ROOT, "jobs_main.csv"), index=False)
    print("main jobs:", len(main))

    jobs.clear()
    if calib in avail:
        wk = pd.read_csv(os.path.join(ROOT, "prepared", calib, "weeks.csv"))
        T = len(wk); reg = wk["regime"].values
        rng = np.random.default_rng(20261029)
        shifts = np.sort(rng.choice(np.arange(52, T - 51), size=200, replace=False))
        lab = pd.DataFrame({f"s{s}": np.roll(reg, s) for s in shifts})
        lab.to_csv(os.path.join(ROOT, "prepared", calib, "placebo_labels.csv"), index=False)
        for rep in ("innov", "dir"):
            for s in shifts:
                add(calib, rep, f"s{s}", "placebo", 199, f"results/placebo/{calib}_{rep}_s{s}",
                    label_file="placebo_labels.csv")
        # bootstrap null: resample innovation rows, regenerate GARCH returns, re-filter
        os.makedirs(os.path.join(ROOT, "prepared", calib, "boot"), exist_ok=True)
        todo = [r for r in range(1, 201) if not all(os.path.exists(os.path.join(ROOT, boot_files(calib, r)[k])) for k in ("innov", "dir"))]
        print(f"bootstrap data sets to generate: {len(todo)}"); sys.stdout.flush()
        if todo:
            with ProcessPoolExecutor(max_workers=max(1, (os.cpu_count() or 2) - 1)) as ex:
                for n, r in enumerate(ex.map(gen_boot, [calib] * len(todo), todo), 1):
                    if n % 20 == 0:
                        print(f"  {n}/{len(todo)}"); sys.stdout.flush()
        for r in range(1, 201):
            for rep in ("innov", "dir"):
                add(calib, rep, "regime", "bootstrap", 199, f"results/bootstrap/{calib}_{rep}_r{r:03d}", x_file=boot_files(calib, r)[rep])
    cal = pd.DataFrame(jobs); cal.to_csv(os.path.join(ROOT, "jobs_calibration.csv"), index=False)
    print("calibration jobs:", len(cal))


if __name__ == "__main__":
    main()
