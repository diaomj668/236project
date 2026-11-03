"""Post hoc supplementary experiments, frozen before new verification simulations.

These analyses were motivated by inspection of the completed r=2 study and
external exploratory descriptions. They are not prespecified main analyses.
"""
from concurrent.futures import ProcessPoolExecutor
from itertools import product, combinations
from pathlib import Path
from time import perf_counter
import hashlib
import importlib.metadata
import json
import platform
import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits
from sklearn.metrics import adjusted_rand_score
from model import (METHODS, generate, projections, spectral, bipartite_svd,
                   metrics, edge_stats, theoretical_weights, generate_separation,
                   append_background)
from run_experiments import summarize

ROOT = Path(__file__).resolve().parent
SENS_KEYS = ['mean_purchases', 'heterogeneity', 'preference', 'r']
SEP_KEYS = ['setting', 'signal_degree', 'r', 'noise_users', 'noise_degree', 'degree_ratio']
DESIGN = {
    'analysis_status': 'Post hoc exploratory supplement; original proposal and main design unchanged',
    'design_record_date': '2026-11-03',
    'selection_basis': 'Completed main results and user-supplied external exploratory descriptions; settings not selected from these fresh runs',
    'n_signal': 1000, 'p': 300, 'k': 3,
    'sensitivity': {'root_seed': 23604, 'r': 3., 'means': [6, 12],
                    'heterogeneity': ['low', 'medium', 'high'],
                    'preference': ['focused', 'broad_active']},
    'separation': {'root_seed': 23605,
                   'settings': [{'name': 'm12_r2', 'degree': 12, 'r': 2.},
                                {'name': 'm6_r3', 'degree': 6, 'r': 3.}],
                   'noise_users': [0, 50, 200], 'degree_ratios': [1, 2, 4, 8],
                   'noise_preference_ratio': 1.,
                   'nesting': 'Same informative matrix within setting/replicate; uniform background users and purchase prefixes nested across cells; A=0 represented once'},
    'kmeans': {'n_init': 20, 'seed_policy': 'Same recorded initialization seed for all methods within a dataset; also shared across separation cells within replicate'},
    'null_reference': {'root_seed': 23606, 'permutations_per_dataset': 20,
                       'definition': 'Permute balanced true product labels, keep fitted labels fixed, and recompute optimal matching; mean permutation error per dataset'},
    'uncertainty': '95% Student-t intervals over independent replicate units within each cell; all six paired method contrasts; no multiplicity adjustment',
    'main_results_sha256': hashlib.sha256((ROOT/'outputs/main/results.csv').read_bytes()).hexdigest(),
    'original_proposal_sha256': hashlib.sha256((ROOT/'236 proposal.pdf').read_bytes()).hexdigest(),
}


def seed(root, *parts):
    return int(np.random.SeedSequence([root, *parts]).generate_state(1)[0])


def evaluate(b, z, ratios, common, km_seed, null_seed):
    """Use unchanged four algorithms and the same initialization across methods."""
    k = len(np.unique(z))
    degrees = b.sum(axis=1)
    theory = theoretical_weights(b.shape[1], k, degrees, ratios)
    predictions, rows = {}, []
    for method in METHODS:
        start = perf_counter()
        if method == 'bipartite_svd':
            pred = bipartite_svd(b, k, km_seed)
            w = None
        else:
            if method == 'weighted':
                scale = np.divide(1., degrees-1, out=np.zeros(len(degrees)), where=degrees>=2)
                w = b.T @ (scale[:, None]*b)
            else:
                w = b.T @ b
                if method == 'jaccard':
                    buyers = b.sum(axis=0)
                    union = buyers[:, None]+buyers[None, :]-w
                    w = np.divide(w, union, out=np.zeros_like(w), where=union>0)
            np.fill_diagonal(w, 0)
            pred = spectral(w, k, km_seed)
        elapsed = perf_counter()-start
        extra = edge_stats(w, z) if w is not None else dict(within=np.nan,between=np.nan,contrast=np.nan,strength_cv=np.nan)
        if w is not None:
            q = theory[method]
            extra |= dict(theory_within=q[0], theory_between=q[1], theory_contrast=q[0]-q[1])
        predictions[method] = pred
        rows.append(common | dict(method=method,kmeans_seed=km_seed,runtime_seconds=elapsed)
                    | metrics(z,pred,k) | extra)
    # Same permuted truths across methods; label matching is repeated each time.
    rng = np.random.default_rng(null_seed)
    null_errors = {m: [] for m in METHODS}
    for _ in range(DESIGN['null_reference']['permutations_per_dataset']):
        permuted = rng.permutation(z)
        for m in METHODS:
            null_errors[m].append(metrics(permuted,predictions[m],k)['error'])
    for row in rows:
        row['permutation_error'] = float(np.mean(null_errors[row['method']]))
        row['null_seed'] = null_seed
    diag = common | dict(purchases=int(b.sum()),users=len(b),isolated_products=int((b.sum(0)==0).sum()),
                         count_pair_mass=float(np.sum(degrees*(degrees-1)/2)),
                         weighted_pair_mass=float(np.sum(degrees[degrees>=2]/2)))
    return rows, diag


def sensitivity_task(task):
    cell, rep = task
    mean,h,pref = cell
    cell_index = list(product([6,12],['low','medium','high'],['focused','broad_active'])).index(cell)
    ds = seed(23604,cell_index,rep,0)
    with threadpool_limits(limits=1):
        b,z,_,_,_,ratios = generate(1000,300,3,mean,h,pref,3.,np.random.default_rng(ds))
        common = dict(mean_purchases=mean,heterogeneity=h,preference=pref,r=3.,rep=rep,dataset_seed=ds,
                      dataset_id=f'r3_{cell_index}_{rep}')
        return evaluate(b,z,ratios,common,seed(23604,cell_index,rep,1),seed(23606,0,cell_index,rep))


def separation_task(task):
    setting_index, rep = task
    setting = DESIGN['separation']['settings'][setting_index]
    m,r = setting['degree'],setting['r']
    signal_seed,noise_seed = seed(23605,setting_index,rep,0),seed(23605,setting_index,rep,1)
    km_seed = seed(23605,setting_index,rep,2)
    rows,diagnostics = [],[]
    with threadpool_limits(limits=1):
        informative,z,order = generate_separation(1000,300,3,m,r,200,m*8,
                                                 np.random.default_rng(signal_seed),np.random.default_rng(noise_seed))
        signal_hash = hashlib.sha256(informative.astype(np.uint8).tobytes()).hexdigest()
        cells = [(0,m)]+[(a,m*ratio) for a in [50,200] for ratio in [1,2,4,8]]
        for cell_index,(a,d) in enumerate(cells):
            b=append_background(informative,order,a,d)
            common=dict(setting=setting['name'],signal_degree=m,r=r,noise_users=a,noise_degree=d,
                        degree_ratio=d/m,rep=rep,signal_seed=signal_seed,noise_seed=noise_seed,
                        signal_sha256=signal_hash,dataset_id=f'sep_{setting_index}_{rep}_{a}_{d}')
            result,diag=evaluate(b,z,np.r_[np.full(1000,r),np.ones(a)],common,km_seed,
                                 seed(23606,1,setting_index,rep,cell_index))
            diag |= dict(signal_count_mass=1000*m*(m-1)/2,background_count_mass=a*d*(d-1)/2,
                         signal_weighted_mass=1000*m/2,background_weighted_mass=a*d/2)
            rows.extend(result);diagnostics.append(diag)
    return rows, diagnostics


def paired_table(df, keys):
    wide=df.pivot(index=keys+['rep'],columns='method',values='error')
    pieces=[]
    for a,b in combinations(METHODS,2):
        delta=(wide[a]-wide[b]).rename('difference').reset_index()
        pieces.append(summarize(delta,keys,'difference').assign(method_a=a,method_b=b))
    return pd.concat(pieces,ignore_index=True)


def analyze(folder, kind):
    saved=[json.loads(p.read_text()) for p in sorted((folder/'datasets').glob('*.json'))]
    rows=[row for item in saved for row in item['rows']]
    diagnostics=[d for item in saved for d in item['diagnostics']]
    df=pd.DataFrame(rows);diag=pd.DataFrame(diagnostics)
    keys=SENS_KEYS if kind.startswith('sensitivity') else SEP_KEYS
    df.to_csv(folder/'results.csv',index=False);diag.to_csv(folder/'diagnostics.csv',index=False)
    summary=pd.concat([summarize(df,keys+['method'],value) for value in
                       ['error','ari','coverage','permutation_error','contrast','strength_cv','runtime_seconds']],ignore_index=True)
    summary.to_csv(folder/'summary.csv',index=False)
    paired_table(df,keys).to_csv(folder/'paired_differences.csv',index=False)
    return df


def run(workers=1, replicates=100, resume=False, analyze_only=False, output=None):
    if workers < 1 or replicates < 2:
        raise ValueError('Require workers>=1 and at least two replicates')
    out=Path(output) if output else ROOT/'outputs/supplement'
    out.mkdir(parents=True,exist_ok=True)
    manifest=DESIGN | dict(replicates=replicates,python=platform.python_version(),
                          versions={m:importlib.metadata.version(m) for m in
                                    ['numpy','scipy','pandas','scikit-learn','matplotlib','threadpoolctl']})
    path=out/'design.json'
    if path.exists():
        if json.loads(path.read_text()) != manifest:
            raise ValueError('Supplement design differs; use a new output directory')
        if not resume and not analyze_only:
            raise ValueError('Supplement already exists; use --resume')
    elif analyze_only:
        raise ValueError('No saved supplementary design')
    else:
        # Written before launching any fresh supplementary simulations.
        path.write_text(json.dumps(manifest,indent=2)+'\n')
    for kind in ['sensitivity_r3','separation']:
        folder=out/kind;checkpoints=folder/'datasets';checkpoints.mkdir(parents=True,exist_ok=True)
        if kind=='sensitivity_r3':
            tasklist=[(cell,rep) for cell in product([6,12],['low','medium','high'],['focused','broad_active']) for rep in range(replicates)]
            fn=sensitivity_task
            name=lambda task:f'{list(product([6,12],["low","medium","high"],["focused","broad_active"])).index(task[0])}_{task[1]}'
        else:
            tasklist=list(product(range(2),range(replicates)));fn=separation_task
            name=lambda task:f'{task[0]}_{task[1]}'
        tasks=[t for t in tasklist if not (checkpoints/f'{name(t)}.json').exists()]
        if not analyze_only:
            executor=ProcessPoolExecutor(workers) if workers>1 else None
            iterator=executor.map(fn,tasks) if executor else map(fn,tasks)
            try:
                for done,(task,result) in enumerate(zip(tasks,iterator),1):
                    rows,diags=result
                    if isinstance(diags,dict):diags=[diags]
                    target=checkpoints/f'{name(task)}.json';temp=target.with_suffix('.tmp')
                    temp.write_text(json.dumps(dict(rows=rows,diagnostics=diags)))
                    temp.replace(target)
                    if done%20==0 or done==len(tasks):
                        print(f'{kind}: {done}/{len(tasks)} new independent units complete',flush=True)
            finally:
                if executor:executor.shutdown()
        df=analyze(folder,kind)
        print(f'{kind}: {df.dataset_id.nunique()} datasets, {len(df)} method records',flush=True)
    from supplement_figures import make_figures
    make_figures(out)
