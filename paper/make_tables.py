"""Appendix tables for paper/report.tex, written from the saved summaries."""
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SUP = ROOT / 'outputs' / 'supplement'
OUT = ROOT / 'paper' / 'tables'
METHODS = ['count', 'jaccard', 'weighted', 'bipartite_svd']
HEAD = r'\textsc{Count} & \textsc{Jaccard} & \textsc{Weighted} & \textsc{SVD}'


def errors(path):
    s = pd.read_csv(path)
    return s[s.metric == 'error'].pivot_table(
        index=['mean_purchases', 'heterogeneity', 'preference'], columns='method', values='mean') * 100


def grid_table(path, caption, label):
    e = errors(path)
    lines = [r'\begin{table}[H]', r'\centering\small', rf'\caption{{{caption}}}', rf'\label{{{label}}}',
             r'\begin{tabular}{@{}cllcccc@{}}', r'\toprule',
             r'$m$ & Heterogeneity & Active users & ' + HEAD + r' \\', r'\midrule']
    for m in (6, 12):
        for h in ('low', 'medium', 'high'):
            for pref in ('focused', 'broad_active'):
                row = e.loc[(m, h, pref)]
                lines.append(f'{m} & {h} & {"focused" if pref == "focused" else "broad"} & '
                             + ' & '.join(f'{row[k]:.1f}' for k in METHODS) + r' \\')
        if m == 6:
            lines.append(r'\addlinespace')
    lines += [r'\bottomrule', r'\end{tabular}', r'\end{table}', '']
    return '\n'.join(lines)


def separation_table():
    s = pd.read_csv(SUP / 'separation/summary.csv')
    s = s[s.metric == 'error']
    e = s.pivot_table(index=['setting', 'noise_users', 'degree_ratio'], columns='method', values='mean') * 100
    names = {'m12_r2': r'$m=12$, $r=2$', 'm6_r3': r'$m=6$, $r=3$'}
    lines = [r'\begin{table}[H]', r'\centering\small',
             r'\caption{Fixed-information experiment: mean misclustering rate (\%) over 100 data sets. '
             r'$A$ is the number of added preference-free users and $D/m$ their purchases relative to the informative users.}',
             r'\label{tab:separation}', r'\begin{tabular}{@{}lcccccc@{}}', r'\toprule',
             r'Informative users & $A$ & $D/m$ & ' + HEAD + r' \\', r'\midrule']
    for setting in ('m12_r2', 'm6_r3'):
        first = True
        for a, ratio in [(0, 1.0)] + [(a, q) for a in (50, 200) for q in (1.0, 2.0, 4.0, 8.0)]:
            row = e.loc[(setting, a, ratio)]
            lab = names[setting] if first else ''
            first = False
            lines.append(f'{lab} & {a} & {"--" if a == 0 else int(ratio)} & '
                         + ' & '.join(f'{row[k]:.1f}' for k in METHODS) + r' \\')
        if setting == 'm12_r2':
            lines.append(r'\addlinespace')
    lines += [r'\bottomrule', r'\end{tabular}', r'\end{table}', '']
    return '\n'.join(lines)


if __name__ == '__main__':
    OUT.mkdir(parents=True, exist_ok=True)
    text = grid_table(ROOT / 'outputs/main/summary.csv',
                      r'Results at $r=2$: mean misclustering rate (\%) over 100 data sets per condition.',
                      'tab:main')
    text += grid_table(SUP / 'sensitivity_r3/summary.csv',
                       r'Results at $r=3$: mean misclustering rate (\%).', 'tab:r3')
    text += grid_table(SUP / 'sensitivity_r4/summary.csv',
                       r'Results at $r=4$: mean misclustering rate (\%).', 'tab:r4')
    text += separation_table()
    (OUT / 'error_tables.tex').write_text(text)
    print('wrote', OUT / 'error_tables.tex')
