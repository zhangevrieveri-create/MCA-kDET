"""Fixed-count precision extension of all 15 MCA screens; never selects a window or asset."""
from pathlib import Path
import os
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[key] = '1'
import sys
import json
import hashlib
import argparse
import numpy as np
import pandas as pd
from scipy.stats import rankdata, beta

ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'market_events'
DEST = OUT / 'mc_precision'
sys.path.insert(0, str(ROOT.parents[1] / 'tests/out/financial_rdata_runtime'))
import rdata


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def holm(p):
    order = np.argsort(p)
    answer = np.empty(len(p))
    answer[order] = np.minimum(1, np.maximum.accumulate(np.asarray(p)[order] * np.arange(len(p), 0, -1)))
    return answer


def interval(success, total):
    return (0. if success == 0 else beta.ppf(.025, success, total-success+1),
            1. if success == total else beta.ppf(.975, success+1, total-success))


def statistics(ranks, labels, targets, counts):
    n = len(labels)
    result = np.full(labels.shape[1], -np.inf)
    for pop, mi in enumerate(counts, 1):
        cols = np.flatnonzero(targets == pop)
        auc = (ranks[:, cols].T @ (labels == pop).astype(float) - mi*(mi+1)/2)/(mi*(n-mi))
        t = (auc-.5)/np.sqrt((n+1)/(12*mi*(n-mi)))
        result = np.maximum(result, t.max(axis=0))
    return result


def verify_saved():
    plan = json.loads((DEST/'plan.json').read_text(encoding='utf-8'))
    assert all(digest(ROOT/p) == sha for p, sha in plan['sources'].items())
    result = pd.read_csv(DEST/'results.csv')
    checks = []
    for job, path in enumerate(sorted((OUT/'jobs').glob('event_m*.rds')), 1):
        a = rdata.read_rds(path)['audits']['1']
        S = np.asarray(a['S']); targets = np.asarray(a['targets']); oldL = np.asarray(a['Lall'])
        new = np.load(DEST/f'event_m{job:02d}_added.npz')
        L = new['labels']; n = len(L)
        # Pairwise comparison weights, independently constructed without scipy ranks.
        weight = np.stack([((S[:, j, None] > S[None, :, j]).sum(axis=1)
                          + .5*(S[:, j, None] == S[None, :, j]).sum(axis=1)) for j in range(S.shape[1])], axis=1)
        calc = np.full(L.shape[1], -np.inf)
        for pop in range(1, 4):
            mi = int((oldL[:, 0] == pop).sum())
            assert np.all((L == pop).sum(axis=0) == mi)
            auc = (weight[:, targets == pop].T @ (L == pop).astype(float)-mi*mi/2)/(mi*(n-mi))
            calc = np.maximum(calc, ((auc-.5)/np.sqrt((n+1)/(12*mi*(n-mi)))).max(axis=0))
        np.testing.assert_allclose(calc, new['statistics'], atol=2e-12, rtol=0)
        rng = np.random.default_rng(np.random.SeedSequence([plan['seed_namespace'], job]))
        first = rng.permuted(np.broadcast_to(oldL[:, 0, None], (n, plan['chunk'])), axis=0)
        assert np.array_equal(first, L[:, :plan['chunk']])
        old = np.asarray(a['statistics']['MCA-kDET'])
        hits = int((old[1:] >= old[0]-1e-10).sum() + (calc >= old[0]-1e-10).sum())
        assert hits == result.exceedances.iloc[job-1]
        assert abs((hits+1)/50000-result.p_refined.iloc[job-1]) < 1e-14
        checks.append(float(np.max(np.abs(calc-new['statistics']))))
    adjusted = np.empty(15); running = 0.
    for i, index in enumerate(np.argsort(result.p_refined.to_numpy())):
        running = max(running, (15-i)*result.p_refined.iloc[index])
        adjusted[index] = min(1., running)
    np.testing.assert_allclose(adjusted, result.holm_refined, atol=1e-14)
    report = dict(added_labelings_independently_recomputed=675000,
        all_pvalues_and_holm_match=True, all_source_hashes_unchanged=True,
        all_label_counts_correct=True, seed_replay_matches=True,
        maximum_statistic_error=max(checks), engine='Pairwise comparison weights without rankdata; saved NPZ read back')
    (DEST/'saved_results_verification.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(report))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--smoke', action='store_true')
    parser.add_argument('--verify-only', action='store_true')
    args = parser.parse_args()
    if args.verify_only:
        verify_saved()
        return
    files = sorted((OUT/'jobs').glob('event_m*.rds'))
    assert len(files) == 15
    calendar = pd.read_csv(OUT/'phase_calendar.csv')
    original = pd.read_csv(OUT/'results.csv')
    plan = dict(B_original=4999, B_total=49999, B_added=45000,
        jobs=15, seed_namespace=2026102801, chunk=500, bit_generator='PCG64',
        seed_rule='SeedSequence([2026102801, job_index]); fresh stream for each job',
        scope='All 15 MCA modules; original 4999 permutations retained. Comparators unchanged.',
        adjustment='Holm over the same 15 MCA hypotheses; no window, module, score, filter or learner changes.',
        interpretation='Monte Carlo precision only, not increased statistical power or validated time-series calibration.',
        provenance='Post hoc precision extension declared before new permutations; all results retained.',
        sources={str(p.relative_to(ROOT)): digest(p) for p in files+[OUT/'protocol.json', OUT/'results.csv']})
    if not args.smoke:
        DEST.mkdir(exist_ok=True)
        plan_path = DEST/'plan.json'
        if plan_path.exists():
            assert json.loads(plan_path.read_text(encoding='utf-8')) == plan
        else:
            plan_path.write_text(json.dumps(plan, indent=2)+'\n', encoding='utf-8')
    rows = []; split_rows = []; audits = []
    for job, path in enumerate(files, 1):
        obj = rdata.read_rds(path)
        a = obj['audits']['1']
        S = np.asarray(a['S']); oldL = np.asarray(a['Lall'], dtype=np.uint8)
        targets = np.asarray(a['targets']); labels = oldL[:, 0]
        counts = np.bincount(labels)[1:]
        ranks = rankdata(S, axis=0, method='average')
        saved = np.asarray(a['statistics']['MCA-kDET'])
        old = statistics(ranks, oldL, targets, counts)
        np.testing.assert_allclose(old, saved, atol=2e-12, rtol=0)
        rng = np.random.default_rng(np.random.SeedSequence([plan['seed_namespace'], job]))
        added = 500 if args.smoke else plan['B_added']
        newL = np.empty((len(labels), added), dtype=np.uint8)
        newT = np.empty(added)
        for start in range(0, added, plan['chunk']):
            end = min(start+plan['chunk'], added)
            L = rng.permuted(np.broadcast_to(labels[:, None], (len(labels), end-start)), axis=0)
            for pop, mi in enumerate(counts, 1):
                assert np.all((L == pop).sum(axis=0) == mi)
            newL[:, start:end] = L
            newT[start:end] = statistics(ranks, L, targets, counts)
        # Independent, direct pairwise AUC check, including ties; no rank-sum formula.
        error = 0.
        for c in (0, added//2, added-1):
            scores = []
            for j, pop in enumerate(targets):
                x = S[newL[:, c] == pop, j]; y = S[newL[:, c] != pop, j]
                auc = ((x[:, None] > y).sum() + .5*(x[:, None] == y).sum())/(len(x)*len(y))
                scores.append((auc-.5)/np.sqrt((len(labels)+1)/(12*len(x)*len(y))))
            error = max(error, abs(max(scores)-newT[c]))
        assert error < 2e-12
        first = obj['rows'].iloc[0]
        original_row = original[(original.event == first.event)&(original.module == first.module)&(original.method == 'MCA-kDET')].iloc[0]
        old_hits = int(np.count_nonzero(saved[1:] >= saved[0]-1e-10))
        hits = old_hits+int(np.count_nonzero(newT >= saved[0]-1e-10))
        total = plan['B_original']+added
        low, high = interval(hits, total)
        rows.append(dict(event=first.event, module=first.module, B=total, p_original=original_row.p_exploratory,
            holm_original=original_row.p_holm_exploratory, p_refined=(hits+1)/(total+1),
            mc_tail_lower=low, mc_tail_upper=high, exceedances=hits, statistic=float(saved[0])))
        if not args.smoke:
            np.savez_compressed(DEST/f'event_m{job:02d}_added.npz', labels=newL, statistics=newT)
        dates = np.asarray(obj['dates'])
        Z = np.asarray(obj['Z'], dtype=float).copy(); Z[np.isnan(Z)] = np.nan
        for phase in range(1, 4):
            fit = np.asarray(obj['fitting_dates'][phase-1]); test = np.asarray(obj['testing_dates'][phase-1])
            zfit = Z[np.isin(dates, fit)]; ztest = Z[np.isin(dates, test)]
            c = calendar[(calendar.event == first.event)&(calendar.phase == phase)]
            split_rows.append(dict(event=first.event, module=first.module, phase=phase, label=c.label.iloc[0],
                dimension=Z.shape[1], fitting_n=len(fit), test_n=len(test),
                fit_start=fit[0], fit_end=fit[-1], test_start=test[0], test_end=test[-1],
                residual_mean_distance=float(np.linalg.norm(zfit.mean(axis=0)-ztest.mean(axis=0))),
                residual_covariance_distance=float(np.linalg.norm(np.atleast_2d(np.cov(zfit, rowvar=False))-np.atleast_2d(np.cov(ztest, rowvar=False)), ord='fro'))))
        audits.append(dict(job=job, original_statistic_error=float(np.max(np.abs(old-saved))),
                           direct_pairwise_added_error=error, all_label_counts_preserved=True))
        print(f'{job}/15 {first.event} {first.module}: {rows[-1]["p_refined"]:.6f}', flush=True)
        if args.smoke:
            print('SMOKE PASSED: all 5000 original statistics and 500 new labelings checked.')
            return
    result = pd.DataFrame(rows)
    result['holm_refined'] = holm(result.p_refined.to_numpy())
    result.to_csv(DEST/'results.csv', index=False)
    pd.DataFrame(split_rows).to_csv(DEST/'split_diagnostics.csv', index=False)
    assert all(digest(ROOT/p) == sha for p, sha in plan['sources'].items())
    verification = dict(all_original_sources_unchanged=True, modules=15, B=49999,
        added_permutations=675000, original_statistics_checked=75000, audits=audits,
        interval_scope='95% pointwise binomial interval for conditional permutation exceedance probability only; not a simultaneous interval, not a CI for an economic parameter.')
    (DEST/'verification.json').write_text(json.dumps(verification, indent=2)+'\n', encoding='utf-8')
    print(result[['event','module','p_original','p_refined','holm_refined']].to_string(index=False))


if __name__ == '__main__':
    main()
