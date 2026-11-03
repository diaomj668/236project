"""Post hoc sensitivity at r=4: rerun the full original 12-cell grid.

Motivation (recorded before running): the proposal listed r=1.25, 2, and 4 as
pilot values, but the pilot used only r=2 because it already met the difficulty
criterion. The r=3 supplement was added after the main results. This run adds
the proposal's upper value, r=4, with fresh seeds and the unchanged methods,
so the preference-strength sensitivity covers r=2, 3, and 4. It is exploratory
and does not replace the frozen r=2 main study.
"""
from concurrent.futures import ProcessPoolExecutor
from itertools import product
from pathlib import Path
import argparse
import importlib.metadata
import json
import platform
import sys
import numpy as np
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from model import generate  # noqa: E402
from supplement import DESIGN as SUPPLEMENT, seed, evaluate, analyze  # noqa: E402
from supplement_figures import sensitivity_figure  # noqa: E402

R, ROOT_SEED, NULL_NAMESPACE = 4., 23607, 2
GRID = list(product([6, 12], ['low', 'medium', 'high'], ['focused', 'broad_active']))
DESIGN = {
    'analysis_status': 'Post hoc exploratory sensitivity; original proposal and main design unchanged',
    'design_record_date': '2026-11-03',
    'selection_basis': 'Upper pilot value r=4 listed in the proposal but not run; added after the main and r=3 results',
    'n': 1000, 'p': 300, 'k': 3, 'r': R, 'root_seed': ROOT_SEED,
    'means': [6, 12], 'heterogeneity': ['low', 'medium', 'high'], 'preference': ['focused', 'broad_active'],
    'kmeans': SUPPLEMENT['kmeans'],
    'null_reference': SUPPLEMENT['null_reference'] | {'namespace': [23606, NULL_NAMESPACE]},
    'uncertainty': SUPPLEMENT['uncertainty'],
    'main_results_sha256': SUPPLEMENT['main_results_sha256'],
    'original_proposal_sha256': SUPPLEMENT['original_proposal_sha256'],
}


def task(args):
    cell_index, rep = args
    mean, h, pref = GRID[cell_index]
    ds = seed(ROOT_SEED, cell_index, rep, 0)
    with threadpool_limits(limits=1):
        b, z, _, _, _, ratios = generate(1000, 300, 3, mean, h, pref, R, np.random.default_rng(ds))
        common = dict(mean_purchases=mean, heterogeneity=h, preference=pref, r=R, rep=rep,
                      dataset_seed=ds, dataset_id=f'r4_{cell_index}_{rep}')
        return evaluate(b, z, ratios, common, seed(ROOT_SEED, cell_index, rep, 1),
                        seed(23606, NULL_NAMESPACE, cell_index, rep))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workers', type=int, default=1)
    parser.add_argument('--replicates', type=int, default=100)
    parser.add_argument('--output', type=Path, default=ROOT/'outputs/supplement')
    parser.add_argument('--analyze-only', action='store_true')
    args = parser.parse_args()
    folder = args.output/'sensitivity_r4'
    checkpoints = folder/'datasets'
    checkpoints.mkdir(parents=True, exist_ok=True)
    manifest = DESIGN | dict(replicates=args.replicates, python=platform.python_version(),
                             versions={m: importlib.metadata.version(m) for m in
                                       ['numpy', 'scipy', 'pandas', 'scikit-learn', 'matplotlib', 'threadpoolctl']})
    path = folder/'design.json'
    if path.exists():
        old = json.loads(path.read_text())
        if {k: v for k, v in old.items() if k not in ('python', 'versions')} != \
           {k: v for k, v in manifest.items() if k not in ('python', 'versions')}:
            parser.error('Recorded r=4 design differs; use a new --output directory')
    elif args.analyze_only:
        parser.error('No saved r=4 design')
    else:
        # Written before launching any r=4 simulations.
        path.write_text(json.dumps(manifest, indent=2)+'\n')
    tasks = [t for t in product(range(len(GRID)), range(args.replicates))
             if not (checkpoints/f'{t[0]}_{t[1]}.json').exists()]
    if not args.analyze_only and tasks:
        with ProcessPoolExecutor(args.workers) as ex:
            for done, (t, (rows, diag)) in enumerate(zip(tasks, ex.map(task, tasks)), 1):
                target = checkpoints/f'{t[0]}_{t[1]}.json'
                temp = target.with_suffix('.tmp')
                temp.write_text(json.dumps(dict(rows=rows, diagnostics=[diag])))
                temp.replace(target)
                if done % 100 == 0 or done == len(tasks):
                    print(f'sensitivity_r4: {done}/{len(tasks)} new datasets complete', flush=True)
    df = analyze(folder, 'sensitivity_r4')
    print(f'sensitivity_r4: {df.dataset_id.nunique()} datasets, {len(df)} method records', flush=True)
    sensitivity_figure(args.output, 4)


if __name__ == '__main__':
    main()
