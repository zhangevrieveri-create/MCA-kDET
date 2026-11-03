"""Amendment 1 (protocol_amendment1.json): dependence-robust inference.

    python make_jobs_v3b.py

Writes
  jobs_rotation.csv       ind30, innov and dir: the whole procedure with the real regime
                          sequence and with every circular shift s = 52, ..., T-52
                          (B = 19; only the observed statistics are used)
  jobs_strata_calib.csv   within-era permutation, L = 26, 52, 104 weeks, on the 200
                          placebo labelings of the original check (B = 199)
  jobs_strata_main.csv    within-era permutation on the real labels (B = 4999)
and prepared/ind30/rotation_labels.csv.
"""
import os, json
import numpy as np, pandas as pd

ROOT = os.path.dirname(os.path.abspath(__file__))
A = json.load(open(os.path.join(ROOT, "protocol_amendment1.json"), encoding="utf-8"))


def job(rows, seed0, kind, ds, rep, label_file, label_col, B, out, strata=0, competitors=1):
    rows.append(dict(job=len(rows) + 1, kind=kind, dataset=ds, representation=rep,
                     x_file=f"prepared/{ds}/X_{rep}.csv", label_file=f"prepared/{ds}/{label_file}",
                     label_col=label_col, role_file=f"prepared/{ds}/weeks.csv", B=B,
                     seed=seed0 + len(rows) + 1, out=out, competitors=competitors, strata=strata))


def main():
    wk = pd.read_csv(os.path.join(ROOT, "prepared", "ind30", "weeks.csv"))
    T = len(wk); reg = wk["regime"].values
    lo, hi = A["rotation"]["min_shift"], T - A["rotation"]["min_shift"]
    shifts = np.arange(lo, hi + 1)
    rot = pd.DataFrame({"s0": reg, **{f"s{s}": np.roll(reg, s) for s in shifts}})
    rot.to_csv(os.path.join(ROOT, "prepared", "ind30", "rotation_labels.csv"), index=False)
    rows = []
    for rep in A["rotation"]["representations"]:
        for s in np.concatenate([[0], shifts]):
            job(rows, 8100000, "rotation", "ind30", rep, "rotation_labels.csv", f"s{s}",
                A["rotation"]["B"], f"results/rotation/ind30_{rep}_s{s}")
    pd.DataFrame(rows).to_csv(os.path.join(ROOT, "jobs_rotation.csv"), index=False)
    print("rotation jobs:", len(rows), f"(shifts {lo}..{hi})")

    pl = pd.read_csv(os.path.join(ROOT, "prepared", "ind30", "placebo_labels.csv"))
    rows = []
    for L in A["strata"]["lengths"]:
        for rep in ("innov", "dir"):
            for c in pl.columns:
                job(rows, 8600000, "strata_placebo", "ind30", rep, "placebo_labels.csv", c, 199,
                    f"results/strata_placebo/L{L}/ind30_{rep}_{c}", strata=L)
    pd.DataFrame(rows).to_csv(os.path.join(ROOT, "jobs_strata_calib.csv"), index=False)
    print("within-era placebo jobs:", len(rows))

    rows = []
    for L in A["strata"]["lengths"]:
        for rep in ("innov", "dir"):
            job(rows, 8900000, "strata_main", "ind30", rep, "weeks.csv", "regime", 4999,
                f"results/strata_main/L{L}/ind30_{rep}", strata=L)
    for ds in ("ind10", "ind49"):
        if os.path.exists(os.path.join(ROOT, "prepared", ds, "weeks.csv")):
            for rep in ("innov", "dir"):
                job(rows, 8900000, "strata_main", ds, rep, "weeks.csv", "regime", 4999,
                    f"results/strata_main/L{A['strata']['primary']}/{ds}_{rep}", strata=A["strata"]["primary"])
    pd.DataFrame(rows).to_csv(os.path.join(ROOT, "jobs_strata_main.csv"), index=False)
    print("within-era main jobs:", len(rows))


if __name__ == "__main__":
    main()
