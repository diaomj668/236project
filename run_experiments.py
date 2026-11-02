"""Reproducible paired simulation CLI. Run --help for available commands."""
import argparse
from concurrent.futures import ProcessPoolExecutor
from itertools import product, combinations
from pathlib import Path
import json
import platform
import time
import importlib.metadata
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, t
from sklearn.metrics import adjusted_rand_score
from threadpoolctl import threadpool_limits
from model import (METHODS, generate, projections, spectral, bipartite_svd,
                   metrics, edge_stats, theoretical_weights)

KEYS = ['mean_purchases', 'heterogeneity', 'preference', 'r']


def one_dataset(task):
    config, cell, rep, seed = task
    mean, heterogeneity, preference, r = cell
    children = np.random.SeedSequence(seed).spawn(5)
    rng = np.random.default_rng(children[0])
    seeds = [int(x.generate_state(1)[0]) for x in children[1:]]
    with threadpool_limits(limits=1):
        start = time.perf_counter()
        b, z, d, pref, active, ratios = generate(config['n'], config['p'], config['k'],
                                               mean, heterogeneity, preference, r, rng)
        generation_seconds = time.perf_counter() - start
        ws, labels, rows = {}, {}, []
        common = dict(zip(KEYS, cell)) | dict(rep=rep, dataset_seed=seed)
        theory = theoretical_weights(config['p'], config['k'], d, ratios)
        for method, km_seed in zip(METHODS, seeds):
            start = time.perf_counter()
            if method == 'bipartite_svd':
                pred = bipartite_svd(b, config['k'], km_seed)
                extras = dict(within=np.nan, between=np.nan, contrast=np.nan, strength_cv=np.nan,
                              theory_within=np.nan, theory_between=np.nan, theory_contrast=np.nan)
            else:
                # Measure construction separately for each method, without charging for others.
                buyers = b.sum(axis=0)
                if method == 'weighted':
                    a = np.divide(1., d-1., out=np.zeros(len(d)), where=d>=2)
                    w = b.T @ (a[:, None]*b)
                else:
                    w = b.T @ b
                    if method == 'jaccard':
                        union = buyers[:, None]+buyers[None, :]-w
                        w = np.divide(w, union, out=np.zeros_like(w), where=union>0)
                np.fill_diagonal(w, 0)
                pred = spectral(w, config['k'], km_seed)
                ws[method] = w
                extras = {}  # edge diagnostics computed outside the timed algorithm
            elapsed = time.perf_counter()-start
            if method != 'bipartite_svd':
                extras = edge_stats(ws[method], z)
                q = theory[method]
                extras |= dict(theory_within=q[0], theory_between=q[1], theory_contrast=q[0]-q[1])
            labels[method] = pred
            rows.append(common | dict(method=method, kmeans_seed=km_seed, runtime_seconds=elapsed)
                        | metrics(z, pred, config['k']) | extras)
        i, j = np.triu_indices(config['p'], 1)
        x, y = ws['count'][i, j], ws['jaccard'][i, j]
        nonzero = x > 0
        good = (labels['count']>=0) & (labels['jaccard']>=0)
        diagnostic = common | dict(
            realized_mean=d.mean(), min_degree=int(d.min()), max_degree=int(d.max()),
            isolated_products=int((b.sum(axis=0)==0).sum()),
            preferred_category_counts=json.dumps(np.bincount(pref, minlength=config['k']).tolist()),
            active_preferred_category_counts=json.dumps(np.bincount(pref[active], minlength=config['k']).tolist()),
            count_pair_mass=float(np.triu(ws['count'],1).sum()),
            expected_count_pair_mass=float(np.sum(d*(d-1)/2)),
            weighted_pair_mass=float(np.triu(ws['weighted'],1).sum()),
            expected_weighted_pair_mass=float(np.sum(d[d>=2]/2)),
            count_jaccard_spearman=float(spearmanr(x,y).statistic) if np.ptp(x)>0 else np.nan,
            count_jaccard_spearman_nonzero=float(spearmanr(x[nonzero],y[nonzero]).statistic) if nonzero.sum()>1 and np.ptp(x[nonzero])>0 else np.nan,
            count_jaccard_ari=adjusted_rand_score(labels['count'][good],labels['jaccard'][good]) if good.sum()>1 else np.nan,
            agreement_coverage=good.mean(), generation_seconds=generation_seconds,
            algorithm_seconds=sum(row['runtime_seconds'] for row in rows))
        assert np.isclose(diagnostic['count_pair_mass'], diagnostic['expected_count_pair_mass'])
        assert np.isclose(diagnostic['weighted_pair_mass'], diagnostic['expected_weighted_pair_mass'])
    return rows, diagnostic


def summarize(df, keys, value):
    a = df.groupby(keys, dropna=False)[value].agg(['count','mean','std']).reset_index()
    a['mcse'] = a['std']/np.sqrt(a['count'])
    margin = t.ppf(.975, a['count']-1)*a['mcse']
    a['ci_low'], a['ci_high'] = a['mean']-margin, a['mean']+margin
    a['metric'] = value
    return a


def analyze(out):
    df = pd.read_csv(out/'results.csv')
    diag = pd.read_csv(out/'diagnostics.csv')
    values = ['error','ari','coverage','contrast','strength_cv','runtime_seconds','within','between']
    summary = pd.concat([summarize(df, KEYS+['method'], x) for x in values], ignore_index=True)
    summary.to_csv(out/'summary.csv',index=False)
    pairs = []
    wide = df.pivot(index=KEYS+['rep','dataset_seed'],columns='method',values='error')
    for a,b in combinations(METHODS,2):
        delta = (wide[a]-wide[b]).rename('difference').reset_index()
        pairs.append(summarize(delta, KEYS,'difference').assign(method_a=a,method_b=b))
    pd.concat(pairs).to_csv(out/'paired_differences.csv',index=False)
    med = df.groupby(KEYS+['method']).error.median().reset_index()
    informative = med.loc[med.error.between(.05,.40),KEYS].drop_duplicates()
    # CPU-time estimate: serial algorithm + generation time, conservative worker-independent reference.
    seconds = (diag.algorithm_seconds+diag.generation_seconds).mean()*1200
    checkpoint = dict(informative_cells=len(informative), difficulty_pass=bool(len(informative)>=2),
                      estimated_1200_dataset_serial_hours=seconds/3600,
                      budget_pass=bool(seconds<=24*3600),
                      minimum_coverage=float(df.coverage.min()),
                      note='Pilot difficulty is not a method-selection criterion. Freeze settings explicitly before main.')
    (out/'checkpoint.json').write_text(json.dumps(checkpoint,indent=2))
    from figures import make_figures
    make_figures(df, summary, out)
    print(json.dumps(checkpoint,indent=2),flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage',choices=['smoke','pilot','calibrate','main'],default='pilot')
    parser.add_argument('--replicates',type=int)
    parser.add_argument('--n',type=int,default=1000)
    parser.add_argument('--p',type=int,default=300)
    parser.add_argument('--k',type=int,default=3)
    parser.add_argument('--means',nargs='+',type=int)
    parser.add_argument('--r',nargs='+',type=float)
    parser.add_argument('--seed',type=int)
    parser.add_argument('--workers',type=int,default=1)
    parser.add_argument('--output',type=Path)
    parser.add_argument('--resume',action='store_true')
    parser.add_argument('--analyze-only',action='store_true')
    parser.add_argument('--frozen-config',type=Path,help='Use a reviewed pilot manifest as final parameter settings')
    args = parser.parse_args()
    out = args.output or Path(__file__).parent/'outputs'/args.stage
    out.mkdir(parents=True,exist_ok=True)
    if args.analyze_only:
        analyze(out)
        return
    reps = args.replicates or dict(smoke=1,pilot=10,calibrate=10,main=100)[args.stage]
    if reps<1 or args.workers<1:
        parser.error('Replicates and workers must be positive')
    means = args.means or ([3,6,12,24] if args.stage=='calibrate' else [6,12])
    ratios = args.r or ([1.25,2.,4.] if args.stage=='calibrate' else [2.])
    config = dict(n=args.n,p=args.p,k=args.k)
    if args.frozen_config:
        frozen = json.loads(args.frozen_config.read_text())
        config = frozen['config']
        means,ratios = frozen['means'],frozen['ratios']
    seed = args.seed if args.seed is not None else dict(smoke=23600,pilot=23601,calibrate=23602,main=23603)[args.stage]
    cells = list(product(means,['low','medium','high'],['focused','broad_active'],ratios))
    manifest = dict(stage=args.stage,config=config,means=means,ratios=ratios,replicates=reps,seed=seed,
                    cells=cells,python=platform.python_version(),
                    versions={m:importlib.metadata.version(m) for m in ['numpy','scipy','pandas','scikit-learn','matplotlib','threadpoolctl']},
                    sampling='sequential PPS without replacement (numpy choice)',
                    theory='unconditional over independent uniform preferred categories; Jaccard is ratio-of-expectations approximation',
                    intervals='95% Student-t Monte Carlo confidence intervals for means; datasets are independent units')
    mf = out/'manifest.json'
    if mf.exists():
        old = json.loads(mf.read_text())
        if not args.resume:
            parser.error('Output already exists; use --resume or a new --output directory')
        if old != json.loads(json.dumps(manifest)):
            parser.error('Resume configuration differs from existing manifest')
    mf.write_text(json.dumps(manifest,indent=2))
    children = np.random.SeedSequence(seed).spawn(len(cells)*reps)
    tasks = [(config,cell,rep,int(children[c*reps+rep].generate_state(1)[0]))
             for c,cell in enumerate(cells) for rep in range(reps)]
    checkpoint_dir = out/'datasets'
    checkpoint_dir.mkdir(exist_ok=True)
    tasks = [task for task in tasks if not (checkpoint_dir/f'{task[3]}.json').exists()]
    # Append per completed dataset, so interrupted jobs can resume.
    executor = ProcessPoolExecutor(args.workers) if args.workers>1 else None
    iterator = executor.map(one_dataset,tasks) if executor else map(one_dataset,tasks)
    try:
        for done,(rows,diag) in enumerate(iterator,1):
            path = checkpoint_dir/f"{diag['dataset_seed']}.json"
            temp = path.with_suffix('.tmp')
            temp.write_text(json.dumps(dict(rows=rows, diagnostic=diag)))
            temp.replace(path)
            if done%10==0 or done==len(tasks):
                print(f'{done}/{len(tasks)} new datasets complete',flush=True)
    finally:
        if executor:
            executor.shutdown()
    saved = [json.loads(path.read_text()) for path in sorted(checkpoint_dir.glob('*.json'))]
    pd.DataFrame([row for item in saved for row in item['rows']]).to_csv(out/'results.csv',index=False)
    pd.DataFrame([item['diagnostic'] for item in saved]).to_csv(out/'diagnostics.csv',index=False)
    analyze(out)


if __name__=='__main__':
    main()
