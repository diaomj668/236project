"""Figures for post hoc supplementary analyses; no assumed detection threshold."""
from pathlib import Path
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/stat236-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from model import METHODS
from figures import random_label_reference, add_reference
NAMES=dict(zip(METHODS,['Count','Jaccard','Activity weighted','Bipartite SVD']))
COLORS=dict(zip(METHODS,['#28658e','#c47b24','#328269','#9e5362']))

def curve(ax,frame,x,label,color):
    y=frame['mean'].to_numpy()*100
    err=np.vstack([frame['mean']-frame.ci_low,frame.ci_high-frame['mean']])*100
    ax.errorbar(x,y,yerr=err,marker='o',capsize=2,lw=1.6,label=label,color=color,ms=4)

def sensitivity_figure(out,r):
    out=Path(out);figs=out/'figures';figs.mkdir(exist_ok=True)
    plt.rcParams.update({'figure.dpi':180,'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    reference=random_label_reference()
    s=pd.read_csv(out/f'sensitivity_r{r}/summary.csv')
    fig,axes=plt.subplots(2,2,figsize=(10.5,6.6),sharey=True)
    for ax,(mean,pref) in zip(axes.flat,[(m,p) for m in [6,12] for p in ['focused','broad_active']]):
        data=s[(s.metric=='error')&(s.mean_purchases==mean)&(s.preference==pref)]
        for method in METHODS:
            f=data[data.method==method].set_index('heterogeneity').reindex(['low','medium','high'])
            curve(ax,f,np.arange(3),NAMES[method],COLORS[method])
        ax.set(title=f'Mean {mean}; {"focused" if pref=="focused" else "broad active"}',
               xticks=np.arange(3),xticklabels=['Low','Medium','High'],ylabel='Misclustering rate (%)',xlabel='Activity heterogeneity',ylim=(-1,70))
        add_reference(ax,reference,scale=100)
    axes[0,0].legend(fontsize=8)
    fig.suptitle(f'Design 1 (matched total) at r={r}: mean misclustering rate with 95% intervals; dashed line: random-label reference',fontsize=12)
    fig.tight_layout();fig.savefig(figs/f'sensitivity_r{r}.png');plt.close(fig)

def make_figures(out):
    out=Path(out);figs=out/'figures'
    sensitivity_figure(out,3)
    reference=random_label_reference()
    s=pd.read_csv(out/'separation/summary.csv');s=s[s.metric=='error']
    fig,axes=plt.subplots(2,2,figsize=(10.5,6.8),sharey=True)
    for row,setting in enumerate(['m12_r2','m6_r3']):
        for col,a in enumerate([50,200]):
            ax=axes[row,col]
            for method in METHODS:
                base=s[(s.setting==setting)&(s.noise_users==0)&(s.method==method)]
                f=s[(s.setting==setting)&(s.noise_users==a)&(s.method==method)].sort_values('degree_ratio')
                frame=pd.concat([base,f])
                curve(ax,frame,np.arange(5),NAMES[method],COLORS[method])
            ax.set(title=f'{setting.replace("_", ", ")}; {a} background users',
                   xticks=np.arange(5),xticklabels=['None','1','2','4','8'],ylabel='Misclustering rate (%)',
                   xlabel='Background / informative purchase count',ylim=(-1,70))
            add_reference(ax,reference,scale=100)
    axes[0,0].legend(fontsize=8)
    fig.suptitle('Design 2 (fixed information): added preference-free users; dashed line: random-label reference',fontsize=12)
    fig.tight_layout();fig.savefig(figs/'separation_recovery.png');plt.close(fig)
    p=pd.read_csv(out/'separation/paired_differences.csv')
    p=p[(p.method_a=='count')&(p.method_b=='weighted')&(p.noise_users>0)].copy()
    vmax=max(abs(p['mean']).max()*100,1)
    fig,axes=plt.subplots(1,2,figsize=(9,3.1),layout='constrained')
    for ax,setting in zip(axes,['m12_r2','m6_r3']):
        data=p[p.setting==setting].pivot(index='noise_users',columns='degree_ratio',values='mean').reindex(index=[50,200],columns=[1,2,4,8])
        values=-data.to_numpy()*100  # Weighted minus Count: negative favors Weighted.
        im=ax.imshow(values,cmap='RdBu_r',vmin=-vmax,vmax=vmax,aspect='auto')
        for i in range(2):
            for j in range(4):
                ax.text(j,i,('0.0' if abs(values[i,j])<.05 else f'{values[i,j]:+.1f}'),ha='center',va='center',color='white' if abs(values[i,j])>.6*vmax else 'black',fontsize=11)
        ax.set(title=setting.replace('_',', '),xticks=range(4),xticklabels=['1','2','4','8'],yticks=range(2),yticklabels=['50','200'],xlabel='Background / informative purchase count',ylabel='Background users')
    fig.colorbar(im,ax=axes,label='Weighted minus Count error (pp)',shrink=.8)
    fig.savefig(figs/'separation_differences.png');plt.close(fig)
    # Tables for the report and auditable browser summary.
    errors=s.pivot(index=['setting','noise_users','degree_ratio'],columns='method',values='mean').reset_index()
    errors.to_csv(out/'separation/mean_errors.csv',index=False)
    s_all=pd.read_csv(out/'separation/summary.csv')
    s_all[s_all.metric=='permutation_error'].to_csv(out/'separation/permutation_reference.csv',index=False)
    print('Supplementary figures generated from saved results.',flush=True)
