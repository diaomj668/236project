"""Build r=3, r=4, and fixed-information tables and website downloads from saved results."""
from pathlib import Path
import json
import shutil
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'outputs/supplement'
ASSETS=ROOT/'website/dist/assets/supplement';ASSETS.mkdir(parents=True,exist_ok=True)
for path in (OUT/'figures').glob('*.png'):shutil.copy2(path,ASSETS/path.name)
shutil.copy2(OUT/'design.json',ASSETS/'design.json')
shutil.copy2(OUT/'sensitivity_r4/design.json',ASSETS/'sensitivity_r4_design.json')
for kind,prefix in [('sensitivity_r3','sensitivity'),('sensitivity_r4','sensitivity_r4'),('separation','separation')]:
    for filename in ['results.csv','summary.csv','diagnostics.csv','paired_differences.csv']:
        shutil.copy2(OUT/kind/filename,ASSETS/f'{prefix}_{filename.replace("paired_differences","paired")}')
# One combined browser dataset covering preference strengths r=2, 3, and 4.
frames=[pd.read_csv(ROOT/'outputs/main/summary.csv'),pd.read_csv(OUT/'sensitivity_r3/summary.csv'),
        pd.read_csv(OUT/'sensitivity_r4/summary.csv')]
cols=['mean_purchases','heterogeneity','preference','r','method','mean','ci_low','ci_high']
records=pd.concat([f[f.metric=='error'][cols] for f in frames],ignore_index=True).to_dict('records')
(ROOT/'website/dist/data.js').write_text('window.RESULTS = '+json.dumps(records,indent=2)+';\n')
sep=pd.read_csv(OUT/'separation/summary.csv')
sens=pd.read_csv(OUT/'sensitivity_r3/summary.csv')
sens4=pd.read_csv(OUT/'sensitivity_r4/summary.csv')
methods=['count','jaccard','weighted','bipartite_svd']
def error(frame,method,**filters):
    f=frame[(frame.method==method)&(frame.metric=='error')]
    for key,value in filters.items():f=f[f[key]==value]
    return float(f.iloc[0]['mean'])*100
finding=f'''<p><strong>Fixed information, added background.</strong> With m=12 and r=2, 50 background users buying 48 products each raise Count error to {error(sep,'count',setting='m12_r2',noise_users=50,noise_degree=48):.2f}%, while Weighted gives {error(sep,'weighted',setting='m12_r2',noise_users=50,noise_degree=48):.2f}%. Without background, both give {error(sep,'count',setting='m12_r2',noise_users=0):.2f}%.</p>'''
p=ROOT/'website/dist/index.html';html=p.read_text()
start=html.index('<div id="supplement-findings">')+len('<div id="supplement-findings">')
end=html.index('</div><figure',start)
html=html[:start]+finding+html[end:];p.write_text(html)
for frame,name in [(sens,'sensitivity_errors.md'),(sens4,'sensitivity_r4_errors.md')]:
 lines=['|Mean|Activity|Preferences|Count|Jaccard|Weighted|SVD|','|---|---|---|---|---|---|---|']
 for mean in [6,12]:
  for h in ['low','medium','high']:
   for pref in ['focused','broad_active']:
    values=[error(frame,m,mean_purchases=mean,heterogeneity=h,preference=pref) for m in methods]
    lines.append(f'|{mean}|{h.title()}|{"Focused" if pref=="focused" else "Broad active"}|'+ '|'.join(f'{v:.2f}' for v in values)+'|')
 (ROOT/'report-assets'/name).write_text('\n'.join(lines)+'\n')
lines=['|Setting|A|D|Count|Jaccard|Weighted|SVD|','|---|---|---|---|---|---|---|']
for setting,m in [('m12_r2',12),('m6_r3',6)]:
 for a,d in [(0,m)]+[(a,m*q) for a in [50,200] for q in [1,2,4,8]]:
  values=[error(sep,method,setting=setting,noise_users=a,noise_degree=d) for method in methods]
  lines.append(f'|{setting}|{a}|{d}|'+ '|'.join(f'{v:.2f}' for v in values)+'|')
(ROOT/'report-assets/separation_errors.md').write_text('\n'.join(lines)+'\n')
print('Supplementary website data, downloads, findings, and report tables built.')
