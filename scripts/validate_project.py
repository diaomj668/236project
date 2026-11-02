"""Validate saved experiment integrity and portable webpage assets."""
from pathlib import Path
from html.parser import HTMLParser
from zipfile import ZipFile
import json
import numpy as np
import pandas as pd
from scipy.stats import t
ROOT=Path(__file__).resolve().parents[1]
KEYS=['mean_purchases','heterogeneity','preference','r']
x=pd.read_csv(ROOT/'outputs/main/results.csv')
d=pd.read_csv(ROOT/'outputs/main/diagnostics.csv')
pilot=pd.read_csv(ROOT/'outputs/pilot/diagnostics.csv')
assert len(d)==1200 and len(x)==4800
assert len(d.dataset_seed.unique())==1200
assert not set(d.dataset_seed)&set(pilot.dataset_seed)
assert x.groupby(KEYS+['method']).size().eq(100).all()
assert x.groupby('dataset_seed').size().eq(4).all()
assert x.error.between(0,1).all() and x.coverage.between(0,1).all()
assert np.allclose(d.realized_mean,d.mean_purchases)
for a,b in [('count_pair_mass','expected_count_pair_mass'),('weighted_pair_mass','expected_weighted_pair_mass')]:
    np.testing.assert_allclose(d[a],d[b])
summary=pd.read_csv(ROOT/'outputs/main/summary.csv')
s=summary[summary.metric=='error'].set_index(KEYS+['method']).sort_index()
expected=x.groupby(KEYS+['method']).error.agg(['mean','std','count']).sort_index()
for column in ['mean','std','count']:np.testing.assert_allclose(s[column],expected[column])
np.testing.assert_allclose(s.mcse,expected['std']/10)
pairs=pd.read_csv(ROOT/'outputs/main/paired_differences.csv')
wide=x.pivot(index=KEYS+['rep','dataset_seed'],columns='method',values='error')
for (a,b),group in pairs.groupby(['method_a','method_b']):
    delta=(wide[a]-wide[b]).groupby(KEYS).agg(['mean','std','count']).sort_index()
    actual=group.set_index(KEYS).sort_index()
    np.testing.assert_allclose(actual['mean'],delta['mean'],atol=1e-14)
    margin=t.ppf(.975,99)*delta['std']/10
    np.testing.assert_allclose(actual.ci_low,delta['mean']-margin,atol=1e-14)
    np.testing.assert_allclose(actual.ci_high,delta['mean']+margin,atol=1e-14)
site=ROOT/'website/dist'
class CheckAssets(HTMLParser):
    def handle_starttag(self,tag,attrs):
        for key,value in attrs:
            if key in ['href','src'] and value and not value.startswith(('#','https:','http:','mailto:')):
                assert (site/value.split('#')[0]).is_file(),value
CheckAssets().feed((site/'index.html').read_text())
assert len(json.loads((site/'assets/results.json').read_text()))==48
for archive in ['reproducibility.zip','website.zip']:
    with ZipFile(ROOT/'deliverables'/archive) as z:
        assert z.testzip() is None
        assert not any('/.git/' in n or '/.venv/' in n or '/.openai/' in n for n in z.namelist())
example=np.load(ROOT/'data/example/purchases.npz')
assert example['B'].shape==(1000,300)
np.testing.assert_array_equal(example['B'].sum(axis=1),example['user_degree'])
print('Saved results, paired confidence intervals, exported data, website assets, and archives passed validation.')
