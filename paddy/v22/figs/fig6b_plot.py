# Fig 6b (v22): structural test of the empirical Si feedback (w 0 / 0.2 / 0.4), reference site, 10 yr. Data: fig6b_si_v22.json
import sys as _s; _s.path.insert(0,"/mnt/project-files/riceland-scepter/v22/figs"); import bigfont
import json,numpy as np,matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False})
R=json.load(open('fig6b_si_v22.json')); W=['0','0.2','0.4']; CW={'0':'#888','0.2':'#2a78d6','0.4':'#0b3f8a'}
LW={'0':'w = 0 (Si feedback off)','0.2':'w = 0.2 (v22 reference)','0.4':'w = 0.4 (stronger)'}
fig,ax=plt.subplots(3,1,figsize=(7,11)); y=np.arange(1,11)
for w in W:
    ax[0].plot(y,R[w]['ctrl_ch4'],marker='o',ms=3,c=CW[w],label=LW[w])
    for d,ls in ((10,'--'),(20,'-')): ax[1].plot(y,R[w][str(d)]['pct'],marker='o',ms=3,ls=ls,c=CW[w],label=f'{LW[w]}, {d} t' )
ax[0].set_title('b1  Control CH$_4$ (kg ha$^{-1}$ season$^{-1}$)',loc='left',fontsize=9)
ax[1].set_title('b2  Season CH$_4$ change vs matched control (%)',loc='left',fontsize=9); ax[1].axhline(0,c='#888',lw=.8)
ax[0].legend(frameon=False,fontsize=7.5)
ax[1].legend([Line2D([],[],c=CW[w]) for w in W]+[Line2D([],[],c='k',ls='--'),Line2D([],[],c='k')],[LW[w] for w in W]+['10 t ha$^{-1}$ yr$^{-1}$','20 t ha$^{-1}$ yr$^{-1}$'],frameon=False,fontsize=7,ncol=2,loc='lower left'); ax[1].set_ylim(-112,3); fig.font_boost=1.05
C=dict(ch4='#2a78d6',n2o='#d6a10f',pda='#1baf7a',cc='#9b59d0'); a=ax[2]; xs=[]; lab=[]
for i,(d,w) in enumerate([(d,w) for d in (10,20) for w in W]):
    x=i+(0.6 if d==20 else 0); s=R[w][str(d)]; pos=neg=0
    for k,sg in (('ch4',1),('n2o',1),('pda',-1),('cc',-1)):
        v=sg*s[k]; l=pos if v>=0 else neg; a.bar(x,v,bottom=l,width=0.7,color=C[k],edgecolor='white',lw=.6)
        if v>=0: pos+=v
        else: neg+=v
    a.plot(x,s['net'],'D',c='#222',ms=6,zorder=5); a.annotate(f"{s['net']:.0f}",(x,s['net']),xytext=(0,-12),textcoords='offset points',ha='center',fontsize=7.5)
    xs.append(x); lab.append(f'w {w}\n{d} t')
a.set_xticks(xs,lab,fontsize=7.5); a.axhline(0,c='#888',lw=.8)
a.set_title('b3  Net GHG change over 10 yr (t CO$_2$e ha$^{-1}$)',loc='left',fontsize=9)
from matplotlib.patches import Patch; from matplotlib.lines import Line2D
a.legend([Patch(color=C[k]) for k in C]+[Line2D([],[],marker='D',c='#222',ls='')],['ΔCH$_4$','ΔN$_2$O','PDA export','calcite','net'],frameon=False,fontsize=7.5,loc='lower left')
for b in ax[:2]: b.set_xlabel('Year'); b.set_xticks(range(1,11)); b.grid(c='#ecebe7',lw=.6)
a.grid(axis='y',c='#ecebe7',lw=.6); a.set_ylim(bottom=min(a.get_ylim()[0],-66)); fig.tight_layout(h_pad=1.5); fig.savefig('fig6b_si_v22.png')
