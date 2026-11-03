"""Post hoc SVD ablation: remove single-purchase users before the bipartite SVD.

Motivation (recorded before running): in focused, high-heterogeneity cells at
mean 6, 80% of users buy one product. They create no product pairs, so Count and
Weighted ignore them, but they still enter Du, Dp, and the diagonal of M'M.
This ablation reuses the recorded datasets and k-means seeds, recomputes the
original SVD (checked against saved errors), and reruns SVD on users with d>=2.
The mean-12 cell has no single-purchase users and serves as a negative control.
"""
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import argparse
import json
import sys
import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from model import generate, bipartite_svd, metrics  # noqa: E402
from run_experiments import summarize  # noqa: E402
from supplement import seed as supp_seed  # noqa: E402

CELLS = [  # (source, mean, heterogeneity, preference, r)
    ('main', 6, 'high', 'focused', 2.),
    ('main', 12, 'high', 'focused', 2.),
    ('sensitivity_r3', 6, 'high', 'focused', 3.),
]
GRID = [(m, h, p) for m in [6, 12] for h in ['low', 'medium', 'high'] for p in ['focused', 'broad_active']]
KEYS = ['source', 'mean_purchases', 'heterogeneity', 'preference', 'r']


def regenerate(row):
    """Rebuild B and the SVD k-means seed exactly as the original run did."""
    if row['source'] == 'main':
        children = np.random.SeedSequence(int(row['dataset_seed'])).spawn(5)
        rng = np.random.default_rng(children[0])
    else:
        rng = np.random.default_rng(int(row['dataset_seed']))
    b, z, d, *_ = generate(1000, 300, 3, int(row['mean_purchases']), row['heterogeneity'],
                           row['preference'], float(row['r']), rng)
    return b, z, d


def task(row):
    with threadpool_limits(limits=1):
        b, z, d = regenerate(row)
        seed = int(row['svd_kmeans_seed'])
        full = metrics(z, bipartite_svd(b, 3, seed), 3)
        drop = metrics(z, bipartite_svd(b[d >= 2], 3, seed), 3)
    assert abs(full['error'] - row['recorded_svd_error']) < 1e-12, 'SVD did not reproduce saved result'
    return {k: row[k] for k in KEYS + ['rep', 'dataset_seed']} | dict(
        single_purchase_users=int((d == 1).sum()),
        count_error=row['recorded_count_error'], weighted_error=row['recorded_weighted_error'],
        svd_error=full['error'], svd_drop_single_error=drop['error'],
        svd_drop_single_coverage=drop['coverage'])


def rows():
    main = pd.read_csv(ROOT/'outputs/main/results.csv')
    sens = pd.read_csv(ROOT/'outputs/supplement/sensitivity_r3/results.csv')
    out = []
    for source, mean, h, pref, r in CELLS:
        frame = main if source == 'main' else sens
        cell = frame[(frame.mean_purchases == mean) & (frame.heterogeneity == h) &
                     (frame.preference == pref) & (frame.r == r)]
        for rep, g in cell.groupby('rep'):
            g = g.set_index('method')
            out.append(dict(source=source, mean_purchases=mean, heterogeneity=h, preference=pref, r=r,
                            rep=int(rep), dataset_seed=int(g.dataset_seed.iloc[0]),
                            svd_kmeans_seed=int(g.loc['bipartite_svd', 'kmeans_seed']),
                            recorded_svd_error=float(g.loc['bipartite_svd', 'error']),
                            recorded_count_error=float(g.loc['count', 'error']),
                            recorded_weighted_error=float(g.loc['weighted', 'error'])))
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workers', type=int, default=1)
    parser.add_argument('--output', type=Path, default=ROOT/'outputs/supplement/svd_ablation')
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output/'design.json').write_text(json.dumps(dict(
        analysis_status='Post hoc exploratory ablation, added after the supplementary analyses',
        design_record_date='2026-11-03', cells=CELLS,
        intervention='Drop users with exactly one purchase, then run the unchanged bipartite_svd',
        reuse='Recorded datasets regenerated from saved seeds; original SVD k-means seed reused; '
              'original SVD errors must reproduce exactly',
        uncertainty='95% Student-t intervals over 100 datasets per cell; no multiplicity adjustment'),
        indent=2)+'\n')
    tasks = rows()
    with ProcessPoolExecutor(args.workers) as ex:
        df = pd.DataFrame(list(ex.map(task, tasks)))
    df.to_csv(args.output/'results.csv', index=False)
    df['difference'] = df.svd_error - df.svd_drop_single_error
    values = ['count_error', 'weighted_error', 'svd_error', 'svd_drop_single_error',
              'svd_drop_single_coverage', 'difference']
    summary = pd.concat([summarize(df, KEYS, v) for v in values], ignore_index=True)
    summary.to_csv(args.output/'summary.csv', index=False)
    lines = ['|Cell|Single-purchase users|Count|Weighted|SVD|SVD without single-purchase users|SVD difference (pp, 95% CI)|',
             '|---|---|---|---|---|---|---|']
    for (source, mean, h, pref, r), g in df.groupby(KEYS, sort=False):
        d = summary[(summary.metric == 'difference') & (summary.source == source) &
                    (summary.mean_purchases == mean) & (summary.r == r)].iloc[0]
        lines.append(f'|Mean {mean}, high, focused, r={r:g}|{int(g.single_purchase_users.mean())}|'
                     f'{100*g.count_error.mean():.2f}|{100*g.weighted_error.mean():.2f}|'
                     f'{100*g.svd_error.mean():.2f}|{100*g.svd_drop_single_error.mean():.2f}|'
                     f'{100*d["mean"]:.2f} ({100*d.ci_low:.2f} to {100*d.ci_high:.2f})|')
    (ROOT/'report-assets/svd_ablation.md').write_text('\n'.join(lines)+'\n')
    print('\n'.join(lines))


if __name__ == '__main__':
    main()
