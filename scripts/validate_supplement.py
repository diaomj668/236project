"""Audit supplementary pairing, fixed signal matrices, null references, and intervals."""
from pathlib import Path
import hashlib,json
import numpy as np
import pandas as pd
from scipy.stats import t
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'outputs/supplement'
design=json.loads((OUT/'design.json').read_text())
assert hashlib.sha256((ROOT/'outputs/main/results.csv').read_bytes()).hexdigest()==design['main_results_sha256']
assert hashlib.sha256((ROOT/'236 proposal.pdf').read_bytes()).hexdigest()==design['original_proposal_sha256']
replicates=design['replicates']
for kind,keys,cells in [('sensitivity_r3',['mean_purchases','heterogeneity','preference','r'],12),
                        ('separation',['setting','signal_degree','r','noise_users','noise_degree','degree_ratio'],18)]:
    folder=OUT/kind;x=pd.read_csv(folder/'results.csv');d=pd.read_csv(folder/'diagnostics.csv')
    assert len(x)==cells*replicates*4 and len(d)==cells*replicates
    assert x.groupby(keys+['method']).size().eq(replicates).all()
    assert x.groupby('dataset_id').size().eq(4).all()
    assert x.error.between(0,1).all() and x.permutation_error.between(0,1).all()
    assert x.coverage.eq(1).all()
    assert x.groupby('dataset_id').kmeans_seed.nunique().eq(1).all()
    s=pd.read_csv(folder/'summary.csv')
    for metric in ['error','permutation_error']:
        expected=x.groupby(keys+['method'])[metric].agg(['mean','std','count']).sort_index()
        actual=s[s.metric==metric].set_index(keys+['method']).sort_index()
        for col in ['mean','std','count']:np.testing.assert_allclose(actual[col],expected[col],atol=1e-14)
    wide=x.pivot(index=keys+['rep'],columns='method',values='error')
    pairs=pd.read_csv(folder/'paired_differences.csv')
    for (a,b),frame in pairs.groupby(['method_a','method_b']):
        expected=(wide[a]-wide[b]).groupby(keys).agg(['mean','std','count']).sort_index()
        actual=frame.set_index(keys).sort_index()
        margin=t.ppf(.975,replicates-1)*expected['std']/np.sqrt(replicates)
        np.testing.assert_allclose(actual['mean'],expected['mean'],atol=1e-14)
        np.testing.assert_allclose(actual.ci_low,expected['mean']-margin,atol=1e-14)
        np.testing.assert_allclose(actual.ci_high,expected['mean']+margin,atol=1e-14)
    if kind=='separation':
        assert d.groupby(['setting','rep']).signal_sha256.nunique().eq(1).all()
        assert d.groupby(['setting','rep']).size().eq(9).all()
        np.testing.assert_allclose(d.purchases,1000*d.signal_degree+d.noise_users*d.noise_degree)
        np.testing.assert_allclose(d.count_pair_mass,d.signal_count_mass+d.background_count_mass)
        np.testing.assert_allclose(d.weighted_pair_mass,d.signal_weighted_mass+d.background_weighted_mass)
        # D=m means constant weights for every contributing user.
        same=x[x.degree_ratio==1].pivot(index=['setting','noise_users','rep'],columns='method',values='error')
        np.testing.assert_allclose(same['count'],same['weighted'],atol=1e-12)
        assert x[x.noise_users==0].groupby(['setting','rep']).size().eq(4).all()
print('Supplementary counts, unchanged main/proposal hashes, paired intervals, nested-signal hashes, weighting equivalence, and null references passed.')

script=(ROOT/'website/dist/data.js').read_text()
records=json.loads(script.removeprefix('window.RESULTS = ').rstrip(';\n'))
assert len(records)==96 and {r['r'] for r in records}=={2.,3.}
for r in records:
    folder=ROOT/'outputs/main' if r['r']==2 else OUT/'sensitivity_r3'
    frame=pd.read_csv(folder/'summary.csv')
    row=frame[(frame.metric=='error')&(frame.mean_purchases==r['mean_purchases'])&(frame.heterogeneity==r['heterogeneity'])&(frame.preference==r['preference'])&(frame.method==r['method'])].iloc[0]
    np.testing.assert_allclose([r['mean'],r['ci_low'],r['ci_high']],[row['mean'],row.ci_low,row.ci_high],atol=1e-14)
print('All 96 interactive website records match saved main or supplementary summaries.')

# Post hoc SVD ablation: original SVD errors must equal the saved records.
abl=pd.read_csv(OUT/'svd_ablation/results.csv')
assert len(abl)==3*replicates and abl.svd_drop_single_coverage.eq(1).all()
for source,frame in [('main',pd.read_csv(ROOT/'outputs/main/results.csv')),('sensitivity_r3',pd.read_csv(OUT/'sensitivity_r3/results.csv'))]:
    saved=frame[frame.method=='bipartite_svd'].set_index(['dataset_seed']).error
    a=abl[abl.source==source].set_index('dataset_seed')
    np.testing.assert_allclose(a.svd_error,saved.loc[a.index],atol=1e-14)
print('SVD ablation reproduces all saved SVD errors before removing single-purchase users.')
