# Prototype for TASK_第三轮_功效与场景_20261101.md, section 2.6 (1 Nov 2026).
# Usage: python proto_screen_sim.py <A|B|D> <seed> <beta> <reps>
# Not part of the R pipeline; numpy + scipy only.
# Prototype: screening-type "max inside the learner", with and without cross-fitting
# under TEST-LABEL permutations (D_fit labels fixed). Paired with full-sample
# max-type permutation tests (maxmean / CLX covariance) on the same data.
import numpy as np, sys, time
from scipy.stats import rankdata
rng0 = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 1)

def ar1(n, d, rho, rng):
    Z = rng.standard_normal((n, d)); X = np.empty_like(Z); X[:, 0] = Z[:, 0]
    c = np.sqrt(1 - rho**2)
    for j in range(1, d): X[:, j] = rho * X[:, j-1] + c * Z[:, j]
    return X

def gen(model, beta, n, d, k, rng):
    out = []
    j = np.arange(1, d+1)
    for g in range(k):
        dev = g == k-1
        if model == 'A':
            X = ar1(n, d, .5, rng)
            if dev: X = X + 0.6*beta*j**-3.
        elif model == 'B':
            X = ar1(n, d, .5, rng)
            if dev: X = X + rng.standard_normal((n, d)) * np.sqrt(4*beta*j**-3.)
        elif model == 'D':
            S = rng.choice([-1., 1.], size=(n, 1)); E = rng.standard_normal((n, d))
            if dev:
                s = int(np.ceil(round(100*beta, 9)))
                if s > 0: E[:, :s] = ar1(n, s, .9, rng)
            X = E + S
        out.append(X)
    return out

def feats(X, kind):
    if kind == 'mean': return X
    if kind == 'sq': return X**2
    if kind == 'prod':  # squares + all products
        d = X.shape[1]; iu = np.triu_indices(d, 1)
        return np.hstack([X**2, X[:, iu[0]] * X[:, iu[1]]])

def tstats(F, lab, k):
    # one-vs-rest standardized mean differences for all labellings
    # F: n x p ; lab: n x nb  -> t: nb x k x p
    n, p = F.shape; nb = lab.shape[1]
    mu = F.mean(0); sd = F.std(0) + 1e-12
    Fz = (F - mu) / sd
    t = np.empty((nb, k, p))
    for i in range(k):
        I = (lab == i).astype(float)          # n x nb
        ni = I.sum(0)                          # nb
        si = I.T @ Fz                          # nb x p  sum over class i
        tot = Fz.sum(0)
        mi = si / ni[:, None]; mo = (tot - si) / (n - ni)[:, None]
        t[:, i, :] = (mi - mo) / np.sqrt(1/ni + 1/(n - ni))[:, None]
    return t, mu, sd

def std_auc(score, is_i):
    # score: n x nb, is_i: n x nb bool -> standardized AUC per column
    n = score.shape[0]
    R = np.apply_along_axis(rankdata, 0, score)
    m1 = is_i.sum(0); m0 = n - m1
    auc = ((R * is_i).sum(0) - m1*(m1+1)/2) / (m1*m0)
    return (auc - .5) / np.sqrt((n+1)/(12*m1*m0))

def topm_score(t, m, F, mu, sd):
    # t: p vector ; returns score vector on F using top-m |t|
    sel = np.argsort(-np.abs(t))[:m]
    return ((F[:, sel] - mu[sel]) / sd[sel]) @ t[sel]

def one_rep(model, beta, n, d, k, B, kind, ms, rng):
    data = gen(model, beta, n, d, k, rng)
    fit, te = [], []
    for g in range(k):
        idx = rng.permutation(n); f = n//2
        fit.append(data[g][idx[:f]]); te.append(data[g][idx[f:]])
    Xf = np.vstack(fit); gf = np.repeat(np.arange(k), n//2)
    Xt = np.vstack(te);  Lt = np.repeat(np.arange(k), n - n//2)
    Ff = feats(Xf, kind); Ft = feats(Xt, kind)
    Lall = np.column_stack([Lt] + [rng.permutation(Lt) for _ in range(B)])
    nb = B + 1
    # direction 1: select on D_fit (fixed labels), AUC on D_test (permuted labels)
    tf, muf, sdf = tstats(Ff, gf[:, None], k)
    # direction 2: select on D_test (permuted labels), AUC on D_fit (fixed labels)
    tt, mut, sdt = tstats(Ft, Lall, k)
    res = {}
    for m in ms:
        T1 = np.empty((k, nb)); T2 = np.empty((k, nb))
        for i in range(k):
            s1 = topm_score(tf[0, i], m, Ft, muf, sdf)               # fixed score on test
            T1[i] = std_auc(np.repeat(s1[:, None], nb, 1), Lall == i)
            # direction 2: one score per labelling, evaluated on fit with true labels
            S2 = np.empty((Ff.shape[0], nb))
            for b in range(nb):
                S2[:, b] = topm_score(tt[b, i], m, Ff, mut, sdt)
            T2[i] = std_auc(S2, np.repeat((gf == i)[:, None], nb, 1))
        st1 = T1.max(0); stc = ((T1 + T2) / np.sqrt(2)).max(0)
        for name, st in (('split', st1), ('cross', stc)):
            res[f'{name}_m{m}'] = (1 + (st[1:] >= st[0] - 1e-10).sum()) / nb
    # full-sample permutation max-type competitor on the SAME data
    X = np.vstack(data); L = np.repeat(np.arange(k), n)
    Lfull = np.column_stack([L] + [rng.permutation(L) for _ in range(B)])
    F = feats(X, kind)
    tF, _, _ = tstats(F, Lfull, k)
    stF = np.abs(tF).max(axis=(1, 2))
    res['fullmax'] = (1 + (stF[1:] >= stF[0] - 1e-10).sum()) / nb
    return res

if __name__ == '__main__':
    model, beta, reps = sys.argv[1], float(sys.argv[3]), int(sys.argv[4])
    kind = {'A': 'mean', 'B': 'sq', 'D': 'prod'}[model]
    d = 100; n = 300; k = 3; B = 200; ms = (1, 3, 10)
    rng = np.random.default_rng(int(sys.argv[2]))
    t0 = time.time(); P = []
    for r in range(reps): P.append(one_rep(model, beta, n, d, k, B, kind, ms, rng))
    keys = P[0].keys()
    print(f"{model} beta={beta} reps={reps} ({time.time()-t0:.0f}s)")
    for kk in keys:
        p = np.array([q[kk] for q in P])
        print(f"  {kk:10s} alone {np.mean(p <= .05):.2f}   with library price (Bonferroni x6) {np.mean(p <= .05/6):.2f}")
