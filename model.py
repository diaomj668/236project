"""Data generation, projections, clustering and exact exchangeable-pair theory."""
from functools import lru_cache
from time import perf_counter
import numpy as np
from scipy.linalg import eigh, svd
from scipy.optimize import linear_sum_assignment
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score

METHODS = ('count', 'jaccard', 'weighted', 'bipartite_svd')


def activity_counts(n, mean, heterogeneity):
    if n % 5 or mean not in (3, 6, 12, 24):
        raise ValueError('n must be divisible by 5; mean must be 3, 6, 12, or 24')
    if heterogeneity == 'low':
        return np.full(n, mean, dtype=int)
    # For calibration mean=3, integer counts 1/11 preserve the exact mean.
    low = mean // 2 if heterogeneity == 'medium' else max(1, mean // 6)
    if heterogeneity not in ('medium', 'high'):
        raise ValueError(heterogeneity)
    high = 5 * mean - 4 * low
    return np.r_[np.full(n * 4 // 5, low), np.full(n // 5, high)].astype(int)


def generate(n, p, k, mean, heterogeneity, preference, r, rng):
    if p % k or k < 2 or p // k < 2 or r < 1:
        raise ValueError('Require K>=2, p divisible by K, category size>=2, r>=1')
    d = rng.permutation(activity_counts(n, mean, heterogeneity))
    if d.max() > p:
        raise ValueError('Purchase count exceeds available distinct products')
    z = np.repeat(np.arange(k), p // k)
    preferred = rng.integers(k, size=n)  # independent of activity
    # Lexicographic tie breaking by independent random uniforms.
    active = np.zeros(n, dtype=bool)
    active[np.lexsort((rng.random(n), d))[-n // 5:]] = True
    ratios = np.full(n, float(r))
    if preference == 'broad_active':
        ratios[active] = 1.
    elif preference != 'focused':
        raise ValueError(preference)
    b = np.zeros((n, p), dtype=float)
    for u in range(n):
        weights = np.where(z == preferred[u], ratios[u], 1.)
        b[u, rng.choice(p, size=d[u], replace=False, p=weights / weights.sum())] = 1
    assert np.array_equal(b.sum(axis=1), d)
    return b, z, d, preferred, active, ratios


def projections(b):
    d, buyers = b.sum(axis=1), b.sum(axis=0)
    shared = b.T @ b
    union = buyers[:, None] + buyers[None, :] - shared
    j = np.divide(shared, union, out=np.zeros_like(shared), where=union > 0)
    weights = np.divide(1., d - 1, out=np.zeros_like(d), where=d >= 2)
    w = b.T @ (weights[:, None] * b)
    for a in (shared, j, w):
        np.fill_diagonal(a, 0)
    return dict(count=shared, jaccard=j, weighted=w)


def cluster_embedding(embedding, valid, k, seed):
    labels = np.full(len(valid), -1, dtype=int)
    norms = np.linalg.norm(embedding, axis=1)
    good = valid & (norms > 1e-12)
    if good.sum() >= k:
        labels[good] = KMeans(n_clusters=k, n_init=20, random_state=seed).fit_predict(
            embedding[good] / norms[good, None])
    return labels


def spectral(w, k, seed):
    strength = w.sum(axis=1)
    valid = strength > 0
    idx = np.flatnonzero(valid)
    embed = np.zeros((len(w), k))
    if len(idx) >= k:
        a = w[np.ix_(idx, idx)] / np.sqrt(strength[idx, None] * strength[None, idx])
        # Largest algebraic eigenvalues, not largest absolute eigenvalues.
        _, v = eigh(a, subset_by_index=(len(idx) - k, len(idx) - 1))
        embed[idx] = v
    return cluster_embedding(embed, valid, k, seed)


def bipartite_svd(b, k, seed):
    du, dp = b.sum(axis=1), b.sum(axis=0)
    rows, cols = du > 0, dp > 0
    embed = np.zeros((b.shape[1], k))
    if min(rows.sum(), cols.sum()) >= k:
        m = b[np.ix_(rows, cols)] / np.sqrt(du[rows, None] * dp[None, cols])
        _, _, vt = svd(m, full_matrices=False, lapack_driver='gesdd')
        embed[cols] = vt[:k].T
    return cluster_embedding(embed, cols, k, seed)


def metrics(truth, labels, k):
    good = labels >= 0
    contingency = np.zeros((k, k), dtype=int)
    np.add.at(contingency, (truth[good], labels[good]), 1)
    i, j = linear_sum_assignment(-contingency)
    return dict(error=1 - contingency[i, j].sum() / len(truth),
                coverage=good.mean(),
                ari=adjusted_rand_score(truth[good], labels[good]) if good.sum() >= 2 else np.nan)


def edge_stats(w, z):
    i, j = np.triu_indices(len(z), 1)
    same = z[i] == z[j]
    strength = w.sum(axis=1)
    within, between = w[i[same], j[same]].mean(), w[i[~same], j[~same]].mean()
    return dict(within=within, between=between, contrast=within-between,
                strength_cv=strength.std() / strength.mean() if strength.mean() else np.nan)


@lru_cache(maxsize=None)
def pair_moments(p, k, d, r):
    """Exact sequential PPSWOR DP; average over a uniform preferred category.

    State c is the number selected from the preferred category. Remaining
    nonpreferred products are exchangeable, conditional on c.
    """
    s, t = p // k, p - p // k
    prob = np.zeros(s + 1)
    prob[0] = 1
    for step in range(d):
        nxt = np.zeros_like(prob)
        for c in np.flatnonzero(prob):
            other = step - c
            denom = r * (s-c) + t-other
            a = r * (s-c) / denom
            if c < s:
                nxt[c+1] += prob[c] * a
            nxt[c] += prob[c] * (1-a)
        prob = nxt
    c = np.arange(s+1)
    rest = d-c
    preferred_pairs = np.sum(prob*c*(c-1))
    other_pairs = np.sum(prob*rest*(rest-1)) * s*(s-1)/(t*(t-1))
    within_total = preferred_pairs + (k-1)*other_pairs
    qin = within_total / (k*s*(s-1))
    qout = (d*(d-1)-within_total) / (k*(k-1)*s*s)
    return float(qin), float(qout)


def theoretical_weights(p, k, degrees, ratios):
    q = np.array([pair_moments(p, k, int(d), float(r)) for d, r in zip(degrees, ratios)])
    scale = np.divide(1., degrees-1., out=np.zeros(len(degrees)), where=degrees>=2)
    count, weighted = q.sum(axis=0), (q*scale[:, None]).sum(axis=0)
    # Unconditional over random preferred categories; NOT E[Jaccard].
    buyer_mean = degrees.sum()/p
    japprox = count / (2*buyer_mean-count)
    return dict(count=count, weighted=weighted, jaccard=japprox)


def generate_separation(n_signal, p, k, signal_degree, r, max_noise_users,
                        max_noise_degree, signal_rng, noise_rng):
    """Fixed informative purchases plus nested uniform background purchases.

    Returned noise rows are ordered samples: prefixes support nested background
    degrees and user counts. No information-user purchases change across cells.
    """
    if n_signal % 5 or not 2 <= signal_degree <= p:
        raise ValueError('n_signal must be divisible by 5; require 2 <= signal_degree <= p')
    if max_noise_users < 0 or not 1 <= max_noise_degree <= p:
        raise ValueError('Invalid background dimensions')
    # The established generator validates categories and implements the same PPSWOR.
    b, z, *_ = generate(n_signal, p, k, signal_degree, 'low', 'focused', r, signal_rng)
    noise_order = np.array([noise_rng.permutation(p)[:max_noise_degree]
                            for _ in range(max_noise_users)], dtype=int)
    return b, z, noise_order


def append_background(signal_b, noise_order, users, degree):
    """Append a nested uniform background without modifying the signal matrix."""
    if not 0 <= users <= len(noise_order) or not 1 <= degree <= noise_order.shape[1]:
        raise ValueError('Background selection outside generated dimensions')
    background = np.zeros((users, signal_b.shape[1]), dtype=float)
    if users:
        background[np.arange(users)[:, None], noise_order[:users, :degree]] = 1.
    return np.vstack([signal_b, background])
