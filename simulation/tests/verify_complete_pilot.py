"""Independent CSV-only verification; never rewrites raw simulation results."""
import csv
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path

root = Path("results_v3_pilot_20261101_cpu16")
def read(path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))
grid = {r["cell"]: r for r in read(root / "formal_grid.csv")}
assert len(grid) == 27
rows = defaultdict(list)
seen = set()
for path in sorted((root / "libpilot3").glob("*.csv")):
    content = read(path)
    assert len(content) == 1, path
    r = content[0]
    cell, rep = r["cell"], int(r["rep"])
    assert cell in grid and (cell, rep) not in seen
    seen.add((cell, rep))
    g = grid[cell]
    assert r["err"] == "" and int(r["B"]) == 200
    assert int(r["seed"]) == int(float(g["seed0"])) + rep
    for key in ("lib", "libsets", "model", "n", "groups"):
        assert r[key] == g[key], (cell, rep, key)
    assert path.name == f"{cell}__r{rep:05d}-{rep:05d}.csv"
    for key, value in r.items():
        if key == "mca" or key.startswith(("mca_", "macc_")):
            if value in ("", "NA"):
                continue
            p = float(value)
            assert math.isfinite(p) and 0 <= p <= 1
            assert abs(201 * p - round(201 * p)) < 1e-9
    rows[cell].append(r)
assert len(seen) == 5100
for cell, g in grid.items():
    target = int(g["reps"])
    assert len(rows[cell]) == target
    assert {int(r["rep"]) for r in rows[cell]} == set(range(1, target + 1))
summary = {r["cell"]: r for r in read(root / "libpilot3_summary.csv")}
assert set(summary) == set(grid)
single = ["mca_" + x for x in next(iter(grid.values()))["lib"].split(",")]
sets = next(iter(grid.values()))["libsets"].split(";")
library_cols = ["mca_set_" + x for x in sets]
def rate(cell, col):
    vals = [float(r[col]) for r in rows[cell] if r.get(col) not in (None, "", "NA")]
    return sum(p <= .05 for p in vals) / len(vals) if vals else None
for cell, s in summary.items():
    assert int(s["errors"]) == 0 and int(s["reps"]) == len(rows[cell])
    for col in ["mca"] + single + library_cols:
        assert abs(float(s[col]) - rate(cell, col)) < 1e-12
rules = read(root / "report_tables/library_rule.csv")
assert len(rules) == 8
alternatives = [c for c, g in grid.items() if float(g["beta"]) > 0]
for j, r in enumerate(rules):
    powers = [rate(c, library_cols[j]) for c in alternatives]
    gaps = [max(rate(c, s) for s in single) - rate(c, library_cols[j]) for c in alternatives]
    assert r["library"] == f"L{j}"
    assert int(r["size"]) == len(sets[j].split("+"))
    assert abs(float(r["max_gap"]) - max(gaps)) < 1e-12
    assert abs(float(r["mean_power"]) - sum(powers) / len(powers)) < 1e-12
    assert int(r["cells_gap_over_005"]) == sum(g > .05 + 1e-12 for g in gaps)
    assert (r["eligible"] == "TRUE") == all(g <= .05 + 1e-12 for g in gaps)
flags = read(root / "report_tables/monitored_null_flags.csv")
assert len(flags) == 18
for r in flags:
    assert int(r["n"]) == 500
    assert abs(float(r["rate"]) - rate(r["cell"], r["method"])) < 1e-12
    assert int(r["rejections"]) == round(500 * float(r["rate"]))
    assert (r["outside_original_band"] == "TRUE") == (float(r["rate"]) < .031 or float(r["rate"]) > .069)
expected_counts = {"power_components.csv": 21, "power_libraries.csv": 21,
    "power_competitors.csv": 21, "library_rule.csv": 8, "monitored_null_flags.csv": 18,
    "learner_diagnostics.csv": 27, "timing.csv": 27}
for name, n in expected_counts.items():
    assert len(read(root / "report_tables" / name)) == n, name
for name in ("null_rates_intervals.csv", "cells_below_competitors.csv"):
    assert (root / "report_tables" / name).stat().st_size > 0
report = Path("试点报告_第三轮_功效与场景.md")
assert report.read_bytes() == (root / report.name).read_bytes()
text = report.read_text(encoding="utf-8")
assert "没有候选库满足" in text and not any(r["eligible"] == "TRUE" for r in rules)
assert "0.028" in text and "0.024" in text
appendix = Path("results_v2_quick/review3_tables.txt")
assert appendix.read_text(encoding="utf-8").strip() in text
result = {"status": "VERIFIED", "raw_repetitions": len(seen), "cells": len(grid),
    "null_repetitions": sum(len(rows[c]) for c, g in grid.items() if float(g["beta"]) == 0),
    "power_repetitions": sum(len(rows[c]) for c in alternatives), "recorded_errors": 0,
    "eligible_libraries": [], "null_flags": [r for r in flags if r["outside_original_band"] == "TRUE"],
    "report_sha256": hashlib.sha256(report.read_bytes()).hexdigest(),
    "scope": "CSV integrity and independent summary/rule recomputation; not a new mathematical validity proof."}
(root / "independent_verification.json").write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
print(json.dumps(result, ensure_ascii=False, indent=2))
