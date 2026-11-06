"""Figures for the LaTeX report (paper/report.tex).

All figures share one style: Computer Modern text to match the LaTeX body,
one colour per method, Jaccard dashed because it nearly coincides with Count,
and a grey band for the random-label reference. Figures are written as vector
PDFs to paper/figures/. Every number is read from the saved experiment outputs;
nothing is re-simulated except the user-influence shares, which regenerate the
recorded focused main datasets from their seeds and are cached in paper/data/.
"""
from pathlib import Path
import os
import sys

os.environ.setdefault('MPLCONFIGDIR', '/tmp/stat236-matplotlib')
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from scipy.stats import t

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'scripts'))
OUT = ROOT / 'paper' / 'figures'
DATA = ROOT / 'paper' / 'data'
SUP = ROOT / 'outputs' / 'supplement'

plt.rcParams.update({
    'font.family': 'serif', 'font.serif': ['cmr10'], 'mathtext.fontset': 'cm',
    'axes.formatter.use_mathtext': True, 'axes.unicode_minus': False,
    'pdf.fonttype': 42, 'font.size': 9, 'axes.titlesize': 9.5, 'axes.labelsize': 9,
    'xtick.labelsize': 8.5, 'ytick.labelsize': 8.5, 'legend.fontsize': 8.5,
    'axes.spines.top': False, 'axes.spines.right': False, 'axes.linewidth': 0.7,
    'xtick.major.width': 0.7, 'ytick.major.width': 0.7, 'lines.linewidth': 1.4,
    'lines.markersize': 4.2, 'legend.frameon': False, 'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.02,
})

METHODS = ['count', 'jaccard', 'weighted', 'bipartite_svd']
LABEL = {'count': 'Count', 'jaccard': 'Jaccard', 'weighted': 'Weighted', 'bipartite_svd': 'Bipartite SVD'}
COLOR = {'count': '#1f4e79', 'jaccard': '#6b9bd1', 'weighted': '#b2182b', 'bipartite_svd': '#6f6f6f'}
STYLE = {'count': '-', 'jaccard': '--', 'weighted': '-', 'bipartite_svd': '-'}
MARK = {'count': 'o', 'jaccard': 'o', 'weighted': 's', 'bipartite_svd': '^'}
LEVELS = ['low', 'medium', 'high']
WIDTH = 6.5  # \textwidth in inches for 1in margins on letter paper


def errors(path):
    s = pd.read_csv(path)
    return s[s.metric == 'error']


def reference_band():
    frames = [pd.read_csv(SUP / f'{k}/summary.csv') for k in ['sensitivity_r3', 'separation', 'sensitivity_r4']]
    v = pd.concat(frames)
    v = v[v.metric == 'permutation_error']['mean'] * 100
    return float(v.min()), float(v.max())


BAND = reference_band()


def band(ax, label=False):
    ax.axhspan(*BAND, color='#bdbdbd', alpha=0.55, lw=0, zorder=0)
    if label:
        ax.text(0.02, BAND[1] + 1.2, 'random labels', transform=ax.get_yaxis_transform(),
                fontsize=7.5, color='#555555', va='bottom')


def method_line(ax, x, frame, method, **kw):
    y = frame['mean'].to_numpy() * 100
    err = np.vstack([y - frame.ci_low.to_numpy() * 100, frame.ci_high.to_numpy() * 100 - y])
    ax.errorbar(x, y, yerr=err, color=COLOR[method], ls=STYLE[method], marker=MARK[method],
                capsize=1.8, elinewidth=0.8, label=LABEL[method],
                mfc='white' if method == 'jaccard' else COLOR[method], **kw)


def legend_handles(methods=METHODS):
    return [Line2D([], [], color=COLOR[m], ls=STYLE[m], marker=MARK[m],
                   mfc='white' if m == 'jaccard' else COLOR[m], label=LABEL[m]) for m in methods]


def save(fig, name):
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / f'{name}.pdf')
    plt.close(fig)


# ---------------------------------------------------------------- Figure 1
def schematic():
    """Toy example: four light focused users and one heavy broad user."""
    cat_color = ['#4c72b0', '#dd8452', '#55a868']
    products = [f'{c}{i}' for c in 'ABC' for i in (1, 2, 3)]
    cat = {p: 'ABC'.index(p[0]) for p in products}
    users = {
        '1': ['A1', 'A2', 'A3'], '2': ['B1', 'B2', 'B3'], '3': ['C1', 'C2', 'C3'],
        '4': ['A1', 'A2'], '5': ['A1', 'A3', 'B2', 'B3', 'C1', 'C2'],
    }
    fig, axes = plt.subplots(1, 3, figsize=(WIDTH, 1.9), gridspec_kw=dict(width_ratios=[1.05, 1, 1]))
    ax = axes[0]
    uy = {u: 1 - i / 4 for i, u in enumerate(users)}
    py = {p: 1 - i / 8 for i, p in enumerate(products)}
    for u, items in users.items():
        heavy = u == '5'
        for p in items:
            ax.plot([0, 1], [uy[u], py[p]], color='#222222' if heavy else '#9a9a9a',
                    lw=1.1 if heavy else 0.7, zorder=1)
    for u in users:
        heavy = u == '5'
        ax.scatter(0, uy[u], s=95 if heavy else 55, color='#222222' if heavy else 'white',
                   edgecolor='#222222', lw=0.8, zorder=3)
        ax.text(-0.09, uy[u], f'$u_{u}$', ha='right', va='center', fontsize=8.5)
    for p in products:
        ax.scatter(1, py[p], s=50, color=cat_color[cat[p]], zorder=3, edgecolor='white', lw=0.5)
        ax.text(1.09, py[p], p, ha='left', va='center', fontsize=7.5)
    ax.text(0, 1.17, 'users', ha='center', fontsize=8.5)
    ax.text(1, 1.17, 'products', ha='center', fontsize=8.5)
    ax.set(xlim=(-0.35, 1.35), ylim=(-0.08, 1.25))
    ax.axis('off')
    ax.set_title('(a) Purchases', loc='left')

    ang = {}
    for k, c in enumerate('ABC'):
        for i in range(3):
            a = np.pi / 2 - 2 * np.pi * k / 3 + (i - 1) * 0.42
            ang[f'{c}{i+1}'] = a
    xy = {p: (np.cos(a), np.sin(a)) for p, a in ang.items()}

    def projection(ax, weighted, title):
        w = {}
        for u, items in users.items():
            d = len(items)
            unit = 1 / (d - 1) if weighted else 1.0
            for i in range(d):
                for j in range(i + 1, d):
                    key = tuple(sorted((items[i], items[j])))
                    w[key] = w.get(key, 0) + unit
        for (a, b), val in w.items():
            same = cat[a] == cat[b]
            ax.plot([xy[a][0], xy[b][0]], [xy[a][1], xy[b][1]], zorder=1,
                    color=cat_color[cat[a]] if same else '#8c8c8c', alpha=0.95 if same else 0.75,
                    lw=0.4 + 2.2 * val)
        for p in products:
            ax.scatter(*xy[p], s=50, color=cat_color[cat[p]], zorder=3, edgecolor='white', lw=0.5)
        within = sum(v for (a, b), v in w.items() if cat[a] == cat[b])
        total = sum(w.values())
        ax.text(0, -1.42, f'within-category weight: {100*within/total:.0f}%',
                ha='center', fontsize=7.8)
        ax.set(xlim=(-1.3, 1.3), ylim=(-1.55, 1.25), aspect='equal')
        ax.axis('off')
        ax.set_title(title, loc='left')

    projection(axes[1], False, '(b) Count projection')
    projection(axes[2], True, '(c) Weighted projection')
    fig.subplots_adjust(wspace=0.08)
    save(fig, 'fig1_schematic')


# ---------------------------------------------------------------- Figure 2
def user_influence():
    from build_user_influence import simulated_shares
    cache = DATA / 'user_influence_shares.csv'
    if not cache.exists():
        DATA.mkdir(parents=True, exist_ok=True)
        simulated_shares().to_csv(cache, index=False)
    df = pd.read_csv(cache)
    fig, axes = plt.subplots(1, 2, figsize=(WIDTH, 2.0), gridspec_kw=dict(width_ratios=[1, 1.45]))
    ax = axes[0]
    d = np.arange(1, 61)
    ax.plot(d, d * (d - 1) / 2, color=COLOR['count'], label=r'Count, $d(d-1)/2$')
    ax.plot(d, np.where(d >= 2, d / 2, 0), color=COLOR['weighted'], label=r'Weighted, $d/2$')
    marks = np.array([2, 3, 6, 12, 18, 26, 36, 52])
    ax.plot(marks, marks * (marks - 1) / 2, 'o', color=COLOR['count'], ms=3)
    ax.plot(marks, marks / 2, 's', color=COLOR['weighted'], ms=3)
    ax.set(yscale='log', xlabel='products bought by one user, $d$', ylabel='total pair weight',
           ylim=(0.8, 3000), xlim=(0, 60))
    ax.legend(loc='upper left')
    ax.set_title('(a) Influence of a single user', loc='left')

    ax = axes[1]
    cells = [(m, h) for m in (6, 12) for h in LEVELS]
    x = np.arange(len(cells))
    for j, m in enumerate(['count', 'weighted']):
        g = df[df.method == m].groupby(['mean_purchases', 'heterogeneity'])
        s = g.agg(theory=('theory', 'mean'), mean=('simulated', 'mean'), sd=('simulated', 'std'),
                  n=('simulated', 'size')).reindex(cells)
        off = (j - 0.5) * 0.36
        ax.bar(x + off, 100 * s.theory, 0.34, color=COLOR[m], alpha=0.85, label=LABEL[m], lw=0)
        half = t.ppf(0.975, s.n - 1) * s.sd / np.sqrt(s.n)
        ax.errorbar(x + off, 100 * s['mean'], yerr=100 * half, fmt='o', color='black', ms=2.2,
                    elinewidth=0.7, capsize=1.2)
    ax.axhline(20, color='#555555', ls=':', lw=0.9, label='their share of users')
    ax.set_xticks(x, [f'{h}\n$m={m}$' for m, h in cells])
    ax.set(ylabel='share of block contrast (%)', ylim=(0, 122), yticks=range(0, 101, 20))
    ax.legend(loc='upper left', ncol=3, handlelength=1.2, columnspacing=1.0)
    ax.set_title('(b) Contrast supplied by the top 20% of users', loc='left')
    fig.subplots_adjust(wspace=0.3)
    save(fig, 'fig2_influence')
    return df


# ---------------------------------------------------------------- Figure 3
def recovery_grid(summary, name, rows=(6, 12), title_r=None):
    fig, axes = plt.subplots(len(rows), 2, figsize=(WIDTH, 1.62 * len(rows)), sharex=True, sharey=True)
    axes = np.atleast_2d(axes)
    for i, m in enumerate(rows):
        for j, pref in enumerate(['focused', 'broad_active']):
            ax = axes[i, j]
            band(ax, label=(i == 0 and j == 0))
            d = summary[(summary.mean_purchases == m) & (summary.preference == pref)]
            for k, meth in enumerate(METHODS):
                f = d[d.method == meth].set_index('heterogeneity').reindex(LEVELS)
                method_line(ax, np.arange(3) + (k - 1.5) * 0.04, f, meth)
            pref_txt = 'focused active users' if pref == 'focused' else 'broad active users'
            ax.set_title(f'({"abcd"[2*i+j]}) $m={m}$, {pref_txt}', loc='left')
            ax.set(ylim=(-2, 72), xticks=range(3))
            ax.set_xticklabels(['low', 'medium', 'high'])
            if j == 0:
                ax.set_ylabel('misclustering (%)')
            if i == len(rows) - 1:
                ax.set_xlabel('activity heterogeneity')
    fig.legend(handles=legend_handles(), loc='upper center', ncol=4, bbox_to_anchor=(0.5, 1.03))
    fig.subplots_adjust(hspace=0.38, wspace=0.08)
    save(fig, name)


# ---------------------------------------------------------------- Figure 4
def paired():
    p = pd.read_csv(ROOT / 'outputs/main/paired_differences.csv')
    p = p[(p.method_a == 'count') & (p.method_b == 'weighted')]
    fig, axes = plt.subplots(1, 2, figsize=(WIDTH, 1.9), sharey=True)
    rows = [(m, h) for m in (6, 12) for h in LEVELS]
    for ax, pref, title in zip(axes, ['focused', 'broad_active'],
                               ['(a) Focused active users', '(b) Broad active users']):
        d = p[p.preference == pref].set_index(['mean_purchases', 'heterogeneity']).reindex(rows)
        y = np.arange(len(rows))[::-1]
        val = d['mean'].to_numpy() * 100
        err = np.vstack([val - d.ci_low * 100, d.ci_high * 100 - val])
        col = [COLOR['weighted'] if v > 0 and lo > 0 else COLOR['count'] if v < 0 and hi < 0 else '#6f6f6f'
               for v, lo, hi in zip(val, d.ci_low, d.ci_high)]
        ax.axvline(0, color='#555555', lw=0.8)
        for yy, v, e0, e1, c in zip(y, val, err[0], err[1], col):
            ax.errorbar(v, yy, xerr=[[e0], [e1]], fmt='o', color=c, capsize=1.8, ms=3.6, elinewidth=0.9)
        ax.set_yticks(y, [f'$m={m}$, {h}' for m, h in rows])
        ax.set_title(title, loc='left')
        ax.set_xlabel('Count error $-$ Weighted error (pp)')
    axes[0].set_xlim(-40, 6)
    axes[1].set_xlim(-6, 6)
    axes[0].text(-39, -0.9, 'Count better', fontsize=7.5, color=COLOR['count'])
    axes[1].text(5.6, -0.9, 'Weighted better', fontsize=7.5, color=COLOR['weighted'], ha='right')
    axes[0].set_ylim(-1.3, 5.6)
    fig.subplots_adjust(wspace=0.06)
    save(fig, 'fig4_paired')


# ---------------------------------------------------------------- Figure 5
def preference_strength():
    tables = {2: errors(ROOT / 'outputs/main/summary.csv'),
              3: errors(SUP / 'sensitivity_r3/summary.csv'),
              4: errors(SUP / 'sensitivity_r4/summary.csv')}
    cells = [(12, 'high', 'focused', '(a) focused active users, high heterogeneity'),
             (12, 'medium', 'broad_active', '(b) broad active users, medium heterogeneity')]
    fig, axes = plt.subplots(1, 2, figsize=(WIDTH, 1.7), sharey=True)
    for ax, (m, h, pref, title) in zip(axes, cells):
        band(ax, label=ax is axes[0])
        for k, meth in enumerate(METHODS):
            rows = []
            for r, s in tables.items():
                rows.append(s[(s.mean_purchases == m) & (s.heterogeneity == h) &
                              (s.preference == pref) & (s.method == meth)].iloc[0])
            method_line(ax, np.array([2, 3, 4]) + (k - 1.5) * 0.03, pd.DataFrame(rows), meth)
        ax.set(xticks=[2, 3, 4], xlabel='preference ratio $r$', ylim=(-2, 72))
        ax.set_title(title, loc='left')
    axes[0].set_ylabel('misclustering (%)')
    fig.legend(handles=legend_handles(), loc='upper center', ncol=4, bbox_to_anchor=(0.5, 1.13))
    fig.subplots_adjust(wspace=0.08)
    save(fig, 'fig5_strength')


# ---------------------------------------------------------------- Figure 6
def separation():
    s = errors(SUP / 'separation/summary.csv')
    fig, axes = plt.subplots(1, 2, figsize=(WIDTH, 1.75), sharey=True)
    titles = {'m12_r2': '(a) informative users: $m=12$, $r=2$', 'm6_r3': '(b) informative users: $m=6$, $r=3$'}
    for ax, setting in zip(axes, ['m12_r2', 'm6_r3']):
        band(ax, label=ax is axes[0])
        for meth in ['count', 'weighted']:
            base = s[(s.setting == setting) & (s.noise_users == 0) & (s.method == meth)]
            for a, ls in [(50, '-'), (200, '--')]:
                f = s[(s.setting == setting) & (s.noise_users == a) & (s.method == meth)].sort_values('degree_ratio')
                f = pd.concat([base, f])
                y = f['mean'].to_numpy() * 100
                err = np.vstack([y - f.ci_low * 100, f.ci_high * 100 - y])
                ax.errorbar(np.arange(5), y, yerr=err, color=COLOR[meth], ls=ls, marker=MARK[meth],
                            capsize=1.8, elinewidth=0.8, mfc=COLOR[meth] if a == 50 else 'white')
        ax.set_xticks(range(5), ['none', '1', '2', '4', '8'])
        ax.set(xlabel=r'background purchases per user, relative to $m$', ylim=(-2, 72))
        ax.set_title(titles[setting], loc='left')
    axes[0].set_ylabel('misclustering (%)')
    handles = [Line2D([], [], color=COLOR['count'], marker='o', label='Count'),
               Line2D([], [], color=COLOR['weighted'], marker='s', label='Weighted'),
               Line2D([], [], color='black', ls='-', label='50 background users'),
               Line2D([], [], color='black', ls='--', label='200 background users')]
    fig.legend(handles=handles, loc='upper center', ncol=4, bbox_to_anchor=(0.5, 1.13))
    fig.subplots_adjust(wspace=0.08)
    save(fig, 'fig6_separation')


# ---------------------------------------------------------------- Appendix
def appendix_strength_grids():
    fig, axes = plt.subplots(2, 4, figsize=(WIDTH, 3.9), sharex=True, sharey=True)
    for i, (r, path) in enumerate([(3, SUP / 'sensitivity_r3/summary.csv'), (4, SUP / 'sensitivity_r4/summary.csv')]):
        s = errors(path)
        for j, (m, pref) in enumerate([(6, 'focused'), (6, 'broad_active'), (12, 'focused'), (12, 'broad_active')]):
            ax = axes[i, j]
            band(ax)
            d = s[(s.mean_purchases == m) & (s.preference == pref)]
            for k, meth in enumerate(METHODS):
                f = d[d.method == meth].set_index('heterogeneity').reindex(LEVELS)
                method_line(ax, np.arange(3) + (k - 1.5) * 0.05, f, meth)
            ax.set_title(f'$r={r}$, $m={m}$, {"focused" if pref == "focused" else "broad"}', fontsize=8.5)
            ax.set(ylim=(-2, 72), xticks=range(3))
            ax.set_xticklabels(['low', 'med.', 'high'])
            if j == 0:
                ax.set_ylabel('misclustering (%)')
    fig.legend(handles=legend_handles(), loc='upper center', ncol=4, bbox_to_anchor=(0.5, 1.04))
    fig.subplots_adjust(hspace=0.3, wspace=0.07)
    save(fig, 'figA1_strength_grids')


def appendix_theory():
    x = pd.read_csv(ROOT / 'outputs/main/results.csv')
    fig, axes = plt.subplots(1, 3, figsize=(WIDTH, 2.0))
    for ax, meth in zip(axes, ['count', 'jaccard', 'weighted']):
        a = x[x.method == meth].groupby(['mean_purchases', 'heterogeneity', 'preference']).agg(
            obs=('contrast', 'mean'), sd=('contrast', 'std'), n=('contrast', 'size'),
            exp=('theory_contrast', 'mean')).reset_index()
        half = t.ppf(0.975, a.n - 1) * a.sd / np.sqrt(a.n)
        lim = max(a.exp.max(), a.obs.max()) * 1.08
        ax.plot([0, lim], [0, lim], color='#888888', lw=0.8, ls='--')
        ax.errorbar(a.exp, a.obs, yerr=half, fmt='o', color=COLOR[meth], ms=3, capsize=1.5, elinewidth=0.8)
        ax.set(xlim=(-lim * 0.03, lim), ylim=(-lim * 0.03, lim))
        ax.ticklabel_format(style='plain', useOffset=False)
        ax.set_title(f'{LABEL[meth]}', loc='left')
        ax.set_xlabel('approximate expected contrast' if meth == 'jaccard' else 'expected contrast')
    axes[0].set_ylabel('simulated contrast')
    fig.subplots_adjust(wspace=0.35)
    save(fig, 'figA2_theory')


if __name__ == '__main__':
    schematic()
    infl = user_influence()
    recovery_grid(errors(ROOT / 'outputs/main/summary.csv'), 'fig3_recovery')
    paired()
    preference_strength()
    separation()
    appendix_strength_grids()
    appendix_theory()
    print('figures written to', OUT)
    print((100 * infl.groupby(['mean_purchases', 'heterogeneity', 'method'])[['theory', 'simulated']].mean()).round(2))
