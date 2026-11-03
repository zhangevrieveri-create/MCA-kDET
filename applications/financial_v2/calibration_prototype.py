"""Score-only prototype for synchronous, centered time-block calibration.

This file does not train models, analyze the real panel, establish a theorem,
or estimate rejection rates. Its smoke checks verify the proposed mechanics.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.stats import rankdata

ROOT = Path(__file__).resolve().parent


def auc_components(scores, groups, targets, valid):
    """scores: full-calendar time x stock x classifier, including missing rows."""
    pool = scores[valid].reshape(-1, scores.shape[2])
    labels = np.tile(groups, int(valid.sum()))
    if not np.isfinite(pool).all():
        raise ValueError("A complete-panel row has missing classifier scores")
    ranks = rankdata(pool, axis=0, method="average")
    result = np.empty(len(targets)); scale = np.empty(len(targets))
    for c, target in enumerate(targets):
        positive = labels == target
        n1 = int(positive.sum()); n0 = len(labels) - n1
        if min(n1, n0) < 2:
            raise ValueError("Too few positive/negative held-out scores")
        result[c] = (ranks[positive, c].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)
        # This is the frozen statistic's scaling convention, not a dependence SE.
        scale[c] = np.sqrt((len(labels) + 1) / (12 * n1 * n0))
    return result, scale


def moving_block_indices(n, length, rng):
    if not 1 <= length <= n:
        raise ValueError("Invalid block length")
    starts = rng.integers(0, n - length + 1, size=int(np.ceil(n / length)))
    return (starts[:, None] + np.arange(length)[None, :]).ravel()[:n]


def centered_block_calibration(scores, groups, targets, valid, *, length=8, B=499, seed=20261027):
    scores = np.asarray(scores, dtype=float)
    groups = np.asarray(groups, dtype=int); targets = np.asarray(targets, dtype=int)
    valid = np.asarray(valid, dtype=bool)
    if scores.shape[:2] != (len(valid), len(groups)):
        raise ValueError("Calendar/stock axes do not match metadata")
    observed_auc, fixed_scale = auc_components(scores, groups, targets, valid)
    observed_components = (observed_auc - .5) / fixed_scale
    rng = np.random.default_rng(seed)
    maxima = np.empty(B); sample_sizes = np.empty(B, dtype=int)
    for b in range(B):
        index = moving_block_indices(len(valid), length, rng)
        resampled_valid = valid[index]
        auc, _ = auc_components(scores[index], groups, targets, resampled_valid)
        # Center each AUC first; subtracting maxima would be a different procedure.
        maxima[b] = np.max((auc - observed_auc) / fixed_scale)
        sample_sizes[b] = int(resampled_valid.sum())
    statistic = float(observed_components.max())
    return dict(observed_auc=observed_auc, observed_components=observed_components,
                observed_max=statistic, bootstrap_max=maxima,
                critical_95=float(np.quantile(maxima, .95, method="higher")),
                approximate_p=(1 + int(np.sum(maxima >= statistic - 1e-10))) / (B + 1),
                resampled_complete_rows=sample_sizes)


def embedding_error_radius(pair_grams, counts, valid):
    """Centered RKHS mean-difference error, used for simultaneous distance bands.

    pair_grams[p,t,u] is <phi(X_a,t)-phi(X_b,t), phi(X_a,u)-phi(X_b,u)>.
    Missing rows must have zero Gram rows/columns. Same counts for all pairs.
    """
    valid = np.asarray(valid, dtype=bool)
    if np.any(pair_grams[:, ~valid, :]) or np.any(pair_grams[:, :, ~valid]):
        raise ValueError("Missing rows must not contribute to feature means")
    n_observed = valid.sum(); n_star = (counts * valid).sum(axis=1)
    if np.any(n_star == 0):
        raise ValueError("A bootstrap draw contains no complete observations")
    weights = counts * valid[None, :] / n_star[:, None] - valid[None, :] / n_observed
    squared = np.einsum("bt,ptu,bu->bp", weights, pair_grams, weights, optimize=True)
    if squared.min() < -1e-8:
        raise ValueError("Feature-difference Gram matrix is not positive semidefinite")
    errors = np.sqrt(np.maximum(squared, 0))
    return np.quantile(errors.max(axis=1), .95, method="higher"), errors


def smoke():
    rng = np.random.default_rng(20261027)
    groups = np.array([0, 0, 1, 1, 2, 2])
    targets = np.array([0, 1, 2, 0, 1, 2])
    valid = np.ones(48, dtype=bool); valid[[7, 8, 23, 39]] = False
    scores = rng.normal(size=(48, 6, 6))
    scores[~valid] = np.nan
    actual_auc, _ = auc_components(scores, groups, targets, valid)
    pool = scores[valid].reshape(-1, 6); labels = np.tile(groups, valid.sum())
    direct_auc = []
    for c, target in enumerate(targets):
        p = pool[labels == target, c]; n = pool[labels != target, c]
        direct_auc.append(np.mean((p[:, None] > n[None, :]) + .5 * (p[:, None] == n[None, :])))
    np.testing.assert_allclose(actual_auc, direct_auc, rtol=1e-14, atol=1e-14)
    constant = np.full((48, 6, 6), .3); constant[~valid] = np.nan
    tie = centered_block_calibration(constant, groups, targets, valid, length=4, B=99)
    assert tie["approximate_p"] == 1 and tie["observed_max"] == 0
    output = centered_block_calibration(scores, groups, targets, valid, length=4, B=99)
    assert output["bootstrap_max"].shape == (99,)
    assert 0 < output["resampled_complete_rows"].min() <= output["resampled_complete_rows"].max() <= len(valid)
    ids = np.empty((99, len(valid)), dtype=int)
    for b in range(99):
        index = moving_block_indices(len(valid), 4, rng)
        # Every length-four source block is consecutive on the full calendar.
        assert np.all(np.diff(index.reshape(-1, 4), axis=1) == 1)
        ids[b] = np.bincount(index, minlength=len(valid))
    a = rng.normal(size=(48, 3)); b = rng.normal(size=(48, 3))
    difference = a - b; difference[~valid] = 0
    grams = np.stack([difference @ difference.T, np.zeros((48, 48))])
    radius, errors = embedding_error_radius(grams, ids, valid)
    assert radius >= 0 and np.all(errors[:, 1] == 0)
    observed_mean = difference.sum(axis=0) / valid.sum()
    manual = []
    for count in ids:
        weight = count * valid / (count * valid).sum()
        manual.append(np.linalg.norm(weight @ difference - observed_mean))
    np.testing.assert_allclose(errors[:, 0], manual, rtol=1e-12, atol=1e-12)
    result = dict(scope="mechanical smoke checks only, synthetic scores/features",
                  auc_matches_direct_pair_comparison=True, constant_scores_p_equals_one=True,
                  bootstrap_preserves_full_calendar_and_synchronous_stock_axis=True,
                  centered_embedding_errors_match_explicit_feature_means=True,
                  zero_feature_difference_handled_without_naive_squared_MMD_bootstrap=True,
                  real_panel_analyzed=False, empirical_size_or_power_estimated=False,
                  theoretical_validity_proved=False)
    destination = ROOT / "design_two_tests" / "calibration_smoke.json"
    destination.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    smoke()
