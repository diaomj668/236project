"""User influence figure: theoretical versus simulated contributions by activity group.

Panel A shows the exact total unordered pair weight a user contributes under
Count and Weighted. Panel B shows the share of the expected within-minus-between
edge contrast supplied by the most active 20% of users in the focused r=2 main
conditions. Bars use the exact pair moments; points regenerate the recorded main
datasets from their seeds and measure each user's realized contribution.
"""
from pathlib import Path
import os
import sys
os.environ.setdefault('MPLCONFIGDIR', '/tmp/stat236-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import t

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from model import generate, pair_moments  # noqa: E402

P, K = 300, 3
S = P // K
N_WITHIN, N_BETWEEN = K*S*(S-1)/2, K*(K-1)/2*S*S
COLORS = {'count': '#28658e', 'weighted': '#328269'}
NAMES = {'count': 'Count', 'weighted': 'Activity weighted'}
CELLS = [(m, h) for m in [6, 12] for h in ['low', 'medium', 'high']]


def user_contrast(b, z, d):
    """Each user's realized contribution to mean within minus mean between weight."""
    c = np.stack([b[:, z == k].sum(axis=1) for k in range(K)], axis=1)
    within = (c*(c-1)/2).sum(axis=1)
    between = d*(d-1)/2-within
    return within/N_WITHIN-between/N_BETWEEN


def active_share(contrast, d, active):
    scale = np.divide(1., d-1., out=np.zeros(len(d)), where=d >= 2)
    out = {}
    for method, w in [('count', contrast), ('weighted', contrast*scale)]:
        out[method] = w[active].sum()/w.sum() if w.sum() else np.nan
    return out


def theory_share(degrees, active, r=2.):
    q = np.array([pair_moments(P, K, int(x), r) for x in degrees])
    return active_share(q[:, 0]-q[:, 1], degrees.astype(float), active)


def simulated_shares():
    x = pd.read_csv(ROOT/'outputs/main/results.csv')
    x = x[(x.method == 'count') & (x.preference == 'focused')]
    rows = []
    for row in x.itertuples():
        rng = np.random.default_rng(np.random.SeedSequence(int(row.dataset_seed)).spawn(5)[0])
        b, z, d, _, active, _ = generate(1000, P, K, row.mean_purchases, row.heterogeneity,
                                         'focused', row.r, rng)
        d = d.astype(float)
        share = active_share(user_contrast(b, z, d), d, active)
        theory = theory_share(d, active)
        rows += [dict(mean_purchases=row.mean_purchases, heterogeneity=row.heterogeneity, method=m,
                      simulated=share[m], theory=theory[m]) for m in share]
    return pd.DataFrame(rows)


def main():
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False,
                         'figure.dpi': 180})
    fig, axes = plt.subplots(1, 2, figsize=(11, 3.8), gridspec_kw=dict(width_ratios=[1, 1.5]))
    ax = axes[0]
    d = np.arange(1, 61)
    ax.plot(d, d*(d-1)/2, color=COLORS['count'], label='Count: d(d-1)/2')
    ax.plot(d, np.where(d >= 2, d/2, 0), color=COLORS['weighted'], label='Weighted: d/2 for d≥2')
    marks = np.array([1, 2, 3, 6, 12, 18, 26, 36, 52])
    ax.scatter(marks, marks*(marks-1)/2, color=COLORS['count'], s=14, zorder=3)
    ax.scatter(marks, np.where(marks >= 2, marks/2, 0), color=COLORS['weighted'], s=14, zorder=3)
    ax.set(yscale='symlog', xlabel='Distinct products purchased (d)', ylabel='Total pair weight per user',
           title='A. Exact influence of one user')
    ax.legend(fontsize=8, frameon=False)
    ax = axes[1]
    df = simulated_shares()
    width = .36
    for j, method in enumerate(['count', 'weighted']):
        g = df[df.method == method].groupby(['mean_purchases', 'heterogeneity'])
        summary = g.agg(theory=('theory', 'mean'), mean=('simulated', 'mean'),
                        sd=('simulated', 'std'), n=('simulated', 'size'))
        summary = summary.reindex(CELLS)
        x = np.arange(len(CELLS))+(j-.5)*width
        ax.bar(x, 100*summary.theory, width, color=COLORS[method], alpha=.35,
               label=f'{NAMES[method]}: exact expectation')
        half = t.ppf(.975, summary.n-1)*summary.sd/np.sqrt(summary.n)
        ax.errorbar(x, 100*summary['mean'], yerr=100*half, fmt='o', color=COLORS[method], ms=4,
                    capsize=2, label=f'{NAMES[method]}: simulated mean (95% CI)')
    ax.axhline(20, color='#888888', ls='--', lw=1, label='Active users\' share of users (20%)')
    ax.set(xticks=range(len(CELLS)), xticklabels=[f'Mean {m}\n{h}' for m, h in CELLS],
           ylabel='Share of block contrast (%)', ylim=(0, 105),
           title='B. Contrast supplied by the most active 20% (focused, r=2)')
    ax.legend(fontsize=7.5, frameon=False, loc='upper left')
    fig.tight_layout()
    fig.savefig(ROOT/'report-assets/user_influence.png')
    plt.close(fig)
    table = df.groupby(['mean_purchases', 'heterogeneity', 'method'])[['theory', 'simulated']].mean()
    print((100*table).round(2).to_string())


if __name__ == '__main__':
    main()
