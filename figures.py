"""Report figures; all uncertainty bars are dataset-level mean confidence intervals."""
import os
os.environ.setdefault('MPLCONFIGDIR', '/tmp/stat236-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from model import METHODS


def random_label_reference():
    """Range of mean label-permutation errors from the supplementary runs, or None.

    Same n, p, K and balanced truth as the main study; used only as a visual
    reference for 'no recovery', not as a theoretical boundary.
    """
    from pathlib import Path
    import pandas as pd
    root = Path(__file__).resolve().parent/'outputs/supplement'
    paths = [root/'sensitivity_r3/summary.csv', root/'separation/summary.csv']
    if not all(p.exists() for p in paths):
        return None
    s = pd.concat([pd.read_csv(p) for p in paths])
    v = s[s.metric == 'permutation_error']['mean']
    return (float(v.min()), float(v.max())) if len(v) else None


def add_reference(ax, reference, scale=1.0, label=True):
    if reference is None:
        return
    lo, hi = reference[0]*scale, reference[1]*scale
    ax.axhspan(lo, hi, color='#999999', alpha=.25, lw=0, zorder=0)
    ax.axhline((lo+hi)/2, color='#666666', ls='--', lw=1, zorder=0,
               label='Random-label reference' if label else None)


def make_figures(df, summary, out):
    folder = out/'figures'
    folder.mkdir(exist_ok=True)
    plt.rcParams.update({'figure.dpi':140,'axes.spines.top':False,'axes.spines.right':False})
    recovery_figures(df, summary, folder)
    figure_rest(df, out, folder)


NAMES = dict(zip(METHODS, ['Count', 'Jaccard', 'Weighted', 'SVD']))
COLORS = dict(zip(METHODS, ['#28658e', '#c47b24', '#328269', '#9e5362']))
# Count and Jaccard nearly coincide, so methods are dodged horizontally and Count is drawn on top.
OFFSETS = dict(zip(METHODS, [-.09, -.03, .03, .09]))
DRAW_ORDER = ['bipartite_svd', 'weighted', 'jaccard', 'count']


def label_line_ends(ax, ends, gap=4.2):
    """Label lines at their right end, nudging labels apart; merge them when lines end together."""
    ends = sorted(ends)
    if ends[-1][0]-ends[0][0] < gap:
        ax.text(2.2, np.mean([y for y, _, _ in ends]), 'All methods', color='#444444',
                va='center', fontsize=8.5, fontweight='bold')
        return
    placed = []
    for y, text, color in ends:
        if placed and y-placed[-1] < gap:
            y = placed[-1]+gap
        placed.append(y)
        ax.text(2.2, y, text, color=color, va='center', fontsize=8.5, fontweight='bold')


def recovery_figures(df, summary, folder):
    reference = random_label_reference()
    means = sorted(df.mean_purchases.unique(), reverse=True)
    for r in sorted(df.r.unique()):
        fig, axes = plt.subplots(len(means), 2, figsize=(10.5, 3.3*len(means)), sharey=True, squeeze=False)
        for ax, (mean, pref) in zip(axes.flat, [(m, p) for m in means for p in ['focused', 'broad_active']]):
            data = summary[(summary.metric=='error')&(summary['mean_purchases']==mean)&(summary.preference==pref)&(summary.r==r)]
            if reference is not None:
                ax.axhspan(100*reference[0], 100*reference[1], color='#bdbdbd', lw=0, zorder=0)
                ax.text(-.25, 100*reference[1]+1.2, 'random labels', ha='left', va='bottom', fontsize=8, color='#777777')
            ax.grid(axis='y', color='#eeeeee', zorder=0)
            ax.spines[['top', 'right']].set_visible(False)
            ends = []
            for method in DRAW_ORDER:
                s = data[data.method==method].set_index('heterogeneity').reindex(['low','medium','high'])
                y = 100*s['mean'].to_numpy()
                error = np.vstack([y-100*s.ci_low.to_numpy(), 100*s.ci_high.to_numpy()-y])
                x = np.arange(3)+OFFSETS[method]
                count = method == 'count'
                ax.plot(x, y, color=COLORS[method], lw=2.4 if count else 1.6, zorder=4 if count else 3)
                # Intervals need >=2 replicates; smoke runs show points only.
                ax.errorbar(x, y, yerr=error if np.isfinite(error).all() else None, fmt='o', ms=5, lw=1,
                            capsize=0, color=COLORS[method], mfc='white' if method == 'jaccard' else COLORS[method], zorder=5)
                ends.append((y[-1], NAMES[method], COLORS[method]))
            label_line_ends(ax, ends)
            ax.set(title=f'Mean {mean} purchases · {"focused" if pref == "focused" else "broad active"}',
                   xticks=range(3), xticklabels=['Low', 'Medium', 'High'], xlim=(-.3, 2.75), ylim=(0, 70))
        for ax in axes[:, 0]:
            ax.set_ylabel('Misclustering rate (%)')
        for ax in axes[-1]:
            ax.set_xlabel('Activity heterogeneity')
        fig.suptitle(f'Recovery at r={r:g}: mean misclustering rate with 95% Monte Carlo intervals', fontsize=11.5)
        fig.tight_layout()
        fig.savefig(folder/f'recovery_r{r:g}.png', dpi=180); plt.close(fig)


def figure_rest(df, out, folder):
    s=df[df.method!='bipartite_svd']
    fig,ax=plt.subplots(figsize=(6,5))
    for method in METHODS[:3]:
        x=s[s.method==method]
        ax.scatter(x.theory_contrast,x.contrast,s=12,alpha=.4,label=method)
    lo=min(ax.get_xlim()[0],ax.get_ylim()[0]); hi=max(ax.get_xlim()[1],ax.get_ylim()[1])
    ax.plot([lo,hi],[lo,hi],color='black',ls='--',lw=1)
    ax.set(xlabel='Theoretical contrast (Jaccard: approximation)',ylabel='Realized within-minus-between weight')
    ax.legend(); fig.tight_layout(); fig.savefig(folder/'theory_vs_simulation.png'); plt.close(fig)
    d=np.arange(1,101)
    fig,ax=plt.subplots(figsize=(6,4))
    ax.plot(d,d*(d-1)/2,label='Count: d(d-1)/2')
    ax.plot(d,np.where(d>=2,d/2,0),label='Weighted: d/2 for d>=2')
    ax.set(xlabel='Distinct products purchased',ylabel='Total unordered pair weight',yscale='symlog')
    ax.legend(); fig.tight_layout(); fig.savefig(folder/'user_influence.png'); plt.close(fig)
    example_network(out)
    import pandas as pd
    paired=pd.read_csv(out/'paired_differences.csv')
    paired=paired[paired.method_a=='count'].copy()
    fig,ax=plt.subplots(figsize=(10,max(4,len(paired)*.17)))
    y=np.arange(len(paired))
    err=np.vstack([paired['mean']-paired.ci_low,paired.ci_high-paired['mean']])
    ax.errorbar(paired['mean'],y,xerr=err if np.isfinite(err).all() else None,fmt='o',capsize=2)
    ax.axvline(0,color='gray',ls='--')
    ax.set(yticks=y,yticklabels=[f"m={row['mean_purchases']} {row['heterogeneity']} {row['preference']} r={row['r']} vs {row['method_b']}" for row in paired.to_dict('records')],xlabel='Count error minus comparator error (95% MC CI)')
    fig.tight_layout(); fig.savefig(folder/'paired_differences.png'); plt.close(fig)


def example_network(out):
    """Illustrative data only; rendering subsamples edges, inference never does."""
    from model import generate, projections
    b,z,_,_,_,_=generate(40,30,3,6,'medium','broad_active',2,np.random.default_rng(23699))
    ws=projections(b)
    fig,axes=plt.subplots(1,2,figsize=(11,5))
    colors=np.array(['#3677b5','#ed9643','#3b9e79'])
    ax=axes[0]
    for u,i in zip(*np.nonzero(b[:12])):
        ax.plot([0,1],[u/11,i/29],color='gray',alpha=.15,lw=.5)
    ax.scatter(np.zeros(12),np.arange(12)/11,c='gray',s=20)
    ax.scatter(np.ones(30),np.arange(30)/29,c=colors[z],s=25)
    ax.set(title='Bipartite sample (12 of 40 users)',xticks=[0,1],xticklabels=['Users','Products'],yticks=[])
    ax=axes[1]
    angles=np.arange(30)*2*np.pi/30
    xy=np.column_stack([np.cos(angles),np.sin(angles)])
    w=ws['weighted']
    for i,j in zip(*np.where(np.triu(w,1)>0)):
        ax.plot(xy[[i,j],0],xy[[i,j],1],color='gray',alpha=.13,lw=.4+min(w[i,j],2))
    ax.scatter(xy[:,0],xy[:,1],c=colors[z],s=40,zorder=3)
    ax.set(title='Weighted product projection',aspect='equal',xticks=[],yticks=[])
    fig.tight_layout(); fig.savefig(out/'figures'/'network_example.png'); plt.close(fig)
