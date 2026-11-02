"""Reconstruct a recorded main-experiment purchase matrix and its ground truth."""
from pathlib import Path
import argparse
import sys
import json
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from model import generate

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mean',type=int,default=12)
    parser.add_argument('--heterogeneity',choices=['low','medium','high'],default='high')
    parser.add_argument('--preference',choices=['focused','broad_active'],default='focused')
    parser.add_argument('--rep',type=int,default=0)
    parser.add_argument('--dataset-seed',type=int)
    parser.add_argument('--output',type=Path,default=ROOT/'data/example')
    args=parser.parse_args()
    records=pd.read_csv(ROOT/'outputs/main/diagnostics.csv')
    if args.dataset_seed is not None:
        selected=records[records.dataset_seed==args.dataset_seed]
    else:
        selected=records[(records.mean_purchases==args.mean)&(records.heterogeneity==args.heterogeneity)&(records.preference==args.preference)&(records.rep==args.rep)]
    if len(selected)!=1:parser.error('Selection must identify exactly one recorded dataset.')
    row=selected.iloc[0]; manifest=json.loads((ROOT/'outputs/main/manifest.json').read_text())
    config=manifest['config'];seed=int(row.dataset_seed)
    # Matches one_dataset: generation uses the first of five spawned substreams.
    stream=np.random.SeedSequence(seed).spawn(5)[0]
    b,z,d,pref,active,ratios=generate(config['n'],config['p'],config['k'],int(row.mean_purchases),row.heterogeneity,row.preference,float(row.r),np.random.default_rng(stream))
    args.output.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(args.output/'purchases.npz',B=b.astype(np.uint8),product_category=z,user_degree=d,user_preference=pref,active=active,preference_ratio=ratios)
    pd.DataFrame(b.astype(np.uint8),index=pd.Index(range(config['n']),name='user_id'),columns=[f'product_{i}' for i in range(config['p'])]).to_csv(args.output/'purchase_matrix.csv')
    pd.DataFrame({'product_id':range(config['p']),'true_category':z}).to_csv(args.output/'products.csv',index=False)
    pd.DataFrame({'user_id':range(config['n']),'purchase_count':d,'preferred_category':pref,'active':active,'preference_ratio':ratios}).to_csv(args.output/'users.csv',index=False)
    metadata={k:(int(row[k]) if k in ['mean_purchases','rep','dataset_seed'] else float(row[k]) if k=='r' else row[k]) for k in ['mean_purchases','heterogeneity','preference','r','rep','dataset_seed']}
    metadata|={'config':config,'source':'Reconstructed from the recorded main-experiment seed; not a new simulation replicate.','indexing':'User IDs, product IDs, and category labels are zero based.'}
    (args.output/'metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    print(f'Exported B with shape {b.shape}, {int(b.sum())} purchases, dataset seed {seed}.')

if __name__=='__main__':main()
