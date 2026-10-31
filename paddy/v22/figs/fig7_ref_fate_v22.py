# Fig 7 (v22): weathering and carbon fate over ten years at the reference site, 10 and 20 t/ha/yr (treatment minus
# year-matched control). Data: ../timeline_ref10.json, ../timeline_ref20.json (timeline22.py, --alk=pda).
import sys as _s; _s.path.insert(0,"/mnt/project-files/riceland-scepter/v22/figs"); import bigfont
import json,numpy as np,matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False})
C=dict(pda='#1baf7a',cc='#9b59d0',ret='#c9c7c0',cat='#52514e')
fig,ax=plt.subplots(3,1,figsize=(7,10.5))
for a,dose in zip(ax[:2],(10,20)):
    R=json.load(open(f'../timeline_ref{dose}.json')); y=np.array([r['y'] for r in R])
    cat=np.cumsum([r['cat'] for r in R]); pda=np.cumsum([r['pda'] for r in R]); cc=np.minimum(np.cumsum([r['cc'] for r in R]),np.maximum(cat-pda,0))
    pda=np.maximum(pda,0); cc=np.maximum(cc,0); ret=cat-pda-cc
    a.stackplot(y,pda,cc,ret,colors=[C['pda'],C['cc'],C['ret']],edgecolor='white',lw=.6)
    a.plot(y,cat,c=C['cat'],lw=1.6)
    a.annotate(f"{cat[-1]:.0f}: cation release\n(CDR potential, not removal)",(y[-1],cat[-1]),xytext=(-6,-2),textcoords='offset points',ha='right',va='top',fontsize=8)
    for v,b,l in ((pda[-1],0,'PDA'),(cc[-1],pda[-1],'calcite'),(ret[-1],pda[-1]+cc[-1],'retained')):
        a.annotate(f'{l} {v:.1f}',(y[-1],b+v/2),xytext=(6,0),textcoords='offset points',va='center',fontsize=8)
    a.set_xlim(1,11.6); a.set_xticks(range(1,11)); a.set_xlabel('Year'); a.grid(c='#ecebe7',lw=.6); a.set_axisbelow(True)
    a.set_title(f"{'ab'[dose//20]}  {dose} t ha$^{{-1}}$ yr$^{{-1}}$: cumulative (t CO$_2$ ha$^{{-1}}$)",loc='left',fontsize=9)
    ax[2].plot(y,[r['pda'] for r in R],marker='o',ms=3,c=C['pda'],ls='-' if dose==20 else '--',label=f'PDA export, {dose} t')
    ax[2].plot(y,[r['cc'] for r in R],marker='o',ms=3,c=C['cc'],ls='-' if dose==20 else '--',label=f'calcite change, {dose} t')
[x.set_ylabel('t CO$_2$ ha$^{-1}$') for x in ax[:2]]
ax[2].axhline(0,c='#888',lw=.8); ax[2].set_xticks(range(1,11)); ax[2].set_xlabel('Year'); ax[2].grid(c='#ecebe7',lw=.6)
ax[2].set_title('c  Annual inorganic C terms (t CO$_2$ ha$^{-1}$ yr$^{-1}$)',loc='left',fontsize=9); ax[2].legend(frameon=False,fontsize=8,ncol=2,loc='upper left')
from matplotlib.patches import Patch
fig.legend([Patch(color=C[k]) for k in ('pda','cc','ret')],['Potentially durable alkalinity exported at the base (PDA = ALK − 2Fe(II) − 2NH$_4^+$)',
  'Soil carbonate (calcite; counted only up to cations not exported)','Retained cations (exchange sites, porewater): potential only'],loc='lower center',ncol=1,frameon=False,fontsize=8)
fig.tight_layout(rect=(0,0.07,1,1),h_pad=1.2); fig.savefig('fig7_ref_fate_v22.png')
