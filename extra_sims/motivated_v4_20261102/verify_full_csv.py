"""Independent CSV-only verification of R's summaries; no simulation is rerun."""
import argparse
import collections
import csv
import hashlib
import json
import math
from pathlib import Path
from scipy.stats import beta

parser = argparse.ArgumentParser()
parser.add_argument('--out', default='motivated_v4_20261102/full_B200')
args = parser.parse_args()
out = Path(args.out)

def read(path):
    with Path(path).open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))

def cp(y, n, alpha=.05):
    return (0.0 if y == 0 else float(beta.ppf(alpha / 2, y, n - y + 1)),
            1.0 if y == n else float(beta.ppf(1 - alpha / 2, y + 1, n - y)))

grid_rows = read(out / 'grid.csv')
grid = {r['cell']: r for r in grid_rows}
assert len(grid) == len(grid_rows) == 28
expected = {(c, r) for c, g in grid.items() for r in range(1, int(g['reps']) + 1)}
assert len(list((out / 'rows').glob('*.rds'))) == len(expected)
for f in read(out / 'source_md5.csv'):
    assert hashlib.md5(Path(f['file']).read_bytes()).hexdigest() == f['md5'], f['file']

by_method = collections.defaultdict(dict)
seen = set()
for r in read(out / 'report_tables/pvalues.csv'):
    key = (r['cell'], int(r['rep']))
    g = grid[key[0]]
    assert key in expected and int(r['seed']) == int(g['seed0']) + key[1]
    fullkey = key + (r['method'],)
    assert fullkey not in seen and r['method'] != 'kmd_eta'
    seen.add(fullkey)
    p = float(r['p'])
    assert math.isfinite(p) and 0 <= p <= 1
    by_method[(key[0], r['method'])][key[1]] = p
assert {(c, r) for c, r, m in seen} == expected
for (c, m), values in by_method.items():
    assert set(values) == set(range(1, int(grid[c]['reps']) + 1))

rates = {}
for r in read(out / 'report_tables/all_method_rates_and_intervals.csv'):
    key = (r['cell'], r['method'])
    assert key not in rates
    values = by_method[key]
    y = sum(p <= .05 for p in values.values())
    n = len(values)
    lo, hi = cp(y, n)
    assert int(r['rejections']) == y and int(r['reps']) == n
    assert int(r['B']) == int(grid[key[0]]['B'])
    assert abs(float(r['rate']) - y / n) < 1e-12
    assert abs(float(r['lower95']) - lo) < 1e-10
    assert abs(float(r['upper95']) - hi) < 1e-10
    rates[key] = y / n
assert set(rates) == set(by_method)

alts = [c for c, g in grid.items() if float(g['target']) > .5]
components = ['mca_' + x for x in grid_rows[0]['lib'].split(',')]
sets = grid_rows[0]['libsets'].split(';')
librows = read(out / 'report_tables/fixed_library_rule.csv')
for i, row in enumerate(librows):
    assert row['library'] == 'L' + str(i) and row['components'] == sets[i]
    powers = [rates[(c, 'mca_set_' + sets[i])] for c in alts]
    gaps = [max(rates[(c, m)] for m in components) - p for c, p in zip(alts, powers)]
    assert abs(float(row['mean_power']) - sum(powers) / len(powers)) < 1e-12
    assert abs(float(row['worst_single_component_deficit']) - max(gaps)) < 1e-12
    assert int(row['failing_cells']) == sum(x > .05 + 1e-12 for x in gaps)
    assert (row['eligible'].upper() == 'TRUE') == all(x <= .05 + 1e-12 for x in gaps)

# A paired conservative interval remains informative even if every repetition
# has the same rejection difference. Bonferroni simultaneous binomial bounds
# for discordant events p10,p01 imply coverage >=95% for p10-p01.
paired = []
for row in read(out / 'report_tables/paired_power_differences.csv'):
    c = row['cell']; lib = sets[0 if row['library'] == 'L0' else 6]
    a = by_method[(c, 'mca_set_' + lib)]
    b = by_method[(c, row['comparator'])]
    n = len(a)
    n10 = sum(a[r] <= .05 < b[r] for r in a)
    n01 = sum(b[r] <= .05 < a[r] for r in a)
    diff = (n10 - n01) / n
    assert int(row['library_only']) == n10 and int(row['comparator_only']) == n01
    assert abs(float(row['power_difference']) - diff) < 1e-12
    l10, u10 = cp(n10, n, .025); l01, u01 = cp(n01, n, .025)
    paired.append(dict(cell=c, library=row['library'], comparator=row['comparator'],
                       reps=n, difference=diff, lower95_conservative=l10-u01,
                       upper95_conservative=u10-l01, library_only=n10, comparator_only=n01))
with (out / 'report_tables/paired_difference_conservative_intervals.csv').open('w', encoding='utf-8-sig', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(paired[0]))
    w.writeheader(); w.writerows(paired)

for row in read(out / 'report_tables/integrity.csv'):
    assert int(row['expected']) == int(row['saved']) == int(grid[row['cell']]['reps'])
    assert int(row['record_errors']) == 0 and row['seeds'].upper() == 'TRUE'
    assert row['libraries'].upper() == row['unique_keys'].upper() == 'TRUE'
    assert int(row['B']) == int(grid[row['cell']]['B'])

result = dict(status='independently_verified', cells=len(grid), repetitions=len(expected),
              pvalue_rows=len(seen), method_cell_summaries=len(rates),
              seeds_counts_B_libraries_unique_keys=True, rates_intervals_library_rule=True,
              source_md5_unchanged=True, conservative_paired_intervals=len(paired),
              note='CSV summary verification, not an independent reimplementation of every test statistic.')
with (out / 'independent_verification.json').open('w', encoding='utf-8') as f:
    json.dump(result, f, ensure_ascii=False, indent=2)
print(json.dumps(result, ensure_ascii=False))
