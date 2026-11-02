"""Build report figures and a compact browser dataset from saved main results."""
from pathlib import Path
import os, json
os.environ.setdefault('MPLCONFIGDIR','/tmp/stat236-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'report-assets'; OUT.mkdir(exist_ok=True)
x=pd.read_csv(ROOT/'outputs/main/results.csv')
s=pd.read_csv(ROOT/'outputs/main/summary.csv')
d=pd.read_csv(ROOT/'outputs/main/diagnostics.csv')
p=pd.read_csv(ROOT/'outputs/main/paired_differences.csv')
methods=['count','jaccard','weighted','bipartite_svd']
colors=['#28658e','#c47b24','#328269','#9e5362']
names=['Count','Jaccard','Activity weighted','Bipartite SVD']
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'figure.dpi':180})
fig,axes=plt.subplots(1,3,figsize=(11,3.5))
for ax,method in zip(axes,methods[:3]):
 a=x[x.method==method].groupby(['mean_purchases','heterogeneity','preference']).agg(observed=('contrast','mean'),sd=('contrast','std'),expected=('theory_contrast','mean')).reset_index()
 ax.errorbar(a.expected,a.observed,yerr=1.984*a.sd/10,fmt='o',color=colors[methods.index(method)],capsize=2,ms=4)
 lim=max(a.expected.max(),a.observed.max())*1.08
 ax.plot([0,lim],[0,lim],ls='--',color='#888888',lw=1)
 ax.set(title=names[methods.index(method)],xlabel='Expected contrast',ylabel='Mean simulated contrast')
 if method=='jaccard': ax.set_xlabel('Approximate contrast')
 ax.ticklabel_format(axis='both',style='sci',scilimits=(-2,2))
fig.tight_layout();fig.savefig(OUT/'contrast_validation.png');plt.close(fig)
fig,axes=plt.subplots(1,2,figsize=(10.5,3.5),sharey=True)
for ax,pref in zip(axes,['focused','broad_active']):
 a=p[(p.method_a=='count')&(p.method_b=='weighted')&(p.preference==pref)].copy()
 a['h_order']=a.heterogeneity.map({'low':0,'medium':1,'high':2})
 a=a.sort_values(['mean_purchases','h_order'])
 y=np.arange(len(a)); val=a['mean'].to_numpy()*100
 err=np.vstack([a['mean']-a.ci_low,a.ci_high-a['mean']])*100
 ax.errorbar(val,y,xerr=err,fmt='o',capsize=3,color='#28658e')
 ax.axvline(0,color='#888888',ls='--')
 ax.set(title='Focused preferences' if pref=='focused' else 'Broad active users',yticks=y,yticklabels=[f"Mean {r.mean_purchases}, {r.heterogeneity}" for r in a.itertuples()],xlabel='Count error minus weighted error (pp)')
 ax.invert_yaxis()
fig.suptitle('Paired differences with 95% Monte Carlo confidence intervals',fontsize=11)
fig.tight_layout();fig.savefig(OUT/'paired_weighting.png');plt.close(fig)
records=[]
for row in s[s.metric=='error'].to_dict('records'):
 records.append({k:(float(v) if isinstance(v,(np.floating,float)) else v) for k,v in row.items() if k in ['mean_purchases','heterogeneity','preference','method','mean','ci_low','ci_high']})
(ROOT/'website/dist/assets/results.json').write_text(json.dumps(records,indent=2))
# Display-ready appendix table using the same ordering as the experiment design.
lines=['|Mean|Activity|Preferences|Count|Jaccard|Weighted|SVD|','|---|---|---|---|---|---|---|']
for mean in [6,12]:
 for h in ['low','medium','high']:
  for pref in ['focused','broad_active']:
   a=x[(x.mean_purchases==mean)&(x.heterogeneity==h)&(x.preference==pref)].groupby('method').error.mean()
   lines.append(f"|{mean}|{h.title()}|{'Focused' if pref=='focused' else 'Broad active'}|"+'|'.join(f'{100*a[m]:.2f}' for m in methods)+'|')
(OUT/'all_errors.md').write_text('\n'.join(lines))
print('Report assets built from recorded main results.')
