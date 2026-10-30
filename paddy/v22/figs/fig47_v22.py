# Figs 4 and 7 (v22): first-year (season + fallow) GHG budget vs dose at four archetypes, and flooded vs AWD.
# Usage (from /home/user/pd): python3 fig47_v22.py <outdir> <figdir> --alk=pda
import sys,os,json; O=sys.argv[1]; FD=sys.argv[2]; sys.argv=['x','--days=500']+sys.argv[3:]
exec(open('paddy/ghg_budget.py').read().split('ctrl = season')[0])
import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False})
K=lambda n: n*44/28*0.273; AK='potentially durable ALK (t CO2/ha)'
ARCH=[('rf189','warm rainfed',1),('irr061','irrigated',0),('irr043','acid high OM',0),('rf184','hot calcareous',1)]  # 1 = AWD in its own run
def run(t,d,sw=False):
    if sw: return f'r22sw_{t}_{d}_y1'
    return f'r22_{t}_{d}_y1' if d in (0,20) else f'r22y1_{t}_{d}_y1'
def bud(c,r):
    a={x:season(f'{O}/flood_{c}{x}') for x in 'sf'}; b={x:season(f'{O}/flood_{r}{x}') for x in 'sf'}
    t=dict(ch4=sum(b[x]['CH4 (kg/ha)']-a[x]['CH4 (kg/ha)'] for x in 'sf')*27e-3,
           n2o=K(sum(b[x]['N2O (kg N/ha)']-a[x]['N2O (kg N/ha)'] for x in 'sf')),
           alk=sum(b[x][AK]-a[x][AK] for x in 'sf'), cat=sum(b[x]['cation release (t CO2/ha)'] for x in 'sf'),
           cc=sum(b[x]['calcite (t CO2/ha)']-a[x]['calcite (t CO2/ha)'] for x in 'sf'),
           c0=a['s']['CH4 (kg/ha)'], dpct=100*(b['s']['CH4 (kg/ha)']/a['s']['CH4 (kg/ha)']-1))
    t['cc']=min(t['cc'],max(t['cat']-t['alk'],0)); t['net']=t['ch4']+t['n2o']-t['alk']-t['cc']; t['pot']=t['ch4']+t['n2o']-t['cat']
    return t
ok=lambda n: os.path.isdir(f'{O}/flood_{n}f')
C=dict(ch4='#2a78d6',n2o='#d6a10f',alk='#1baf7a',cc='#9b59d0')
def stack(ax,x,t,w):
    pos=neg=0
    for k,s in (('ch4',1),('n2o',1),('alk',-1),('cc',-1)):
        v=s*t[k]; l=pos if v>=0 else neg; ax.bar(x,v,bottom=l,width=w,color=C[k],edgecolor='white',lw=.8)
        if v>=0: pos+=v
        else: neg+=v
    ax.plot(x,t['net'],'D',c='#222',ms=6,zorder=5); ax.plot(x,t['pot'],'o',mfc='white',mec='#222',ms=6,zorder=5)
from matplotlib.patches import Patch; from matplotlib.lines import Line2D
LEG=([Patch(color=C[k]) for k in C]+[Line2D([],[],marker='D',c='#222',ls=''),Line2D([],[],marker='o',mfc='white',mec='#222',ls='')],
     ['ΔCH$_4$ (season + fallow)','ΔN$_2$O (season + fallow)','Potentially durable alkalinity exported','Calcite stored','Net, realised','Net, potential (all cations)'])
res={}
# Fig 4
fig,axs=plt.subplots(1,4,figsize=(15,4.3),sharey=True,dpi=170)
for ax,(t,lab,awd) in zip(axs,ARCH):
    c=run(t,0); res[t]={}
    for i,d in enumerate((5,10,20,30)):
        r=run(t,d)
        if ok(c) and ok(r): res[t][d]=b=bud(c,r); stack(ax,i,b,0.62); ax.annotate('%+.0f%%'%b['dpct'],(i,0),xytext=(0,4),textcoords='offset points',ha='center',fontsize=7.5,color='#2a78d6')
    c0=res[t].get(20,next(iter(res[t].values()),{})).get('c0',np.nan)
    ax.set_title(f'{t} {lab} ({"AWD" if awd else "flooded"})\ncontrol {c0:.0f} kg CH$_4$/season',loc='left',fontsize=9)
    ax.set_xticks(range(4),['5','10','20','30']); ax.set_xlabel('Basalt (t ha$^{-1}$, applied once)'); ax.axhline(0,c='#888',lw=.8); ax.grid(axis='y',c='#e5e4df',lw=.6); ax.set_axisbelow(True)
axs[0].set_ylabel('Change in first year, season + fallow (t CO$_2$e ha$^{-1}$)')
fig.legend(*LEG,loc='lower center',ncol=6,frameon=False,fontsize=8); fig.tight_layout(rect=(0,0.07,1,1)); fig.savefig(f'{FD}/fig4_dose_v22.png')
# Fig 7: own water management vs switched
W={}
for t,lab,awd in ARCH:
    for sw in (False,True):
        wm='AWD' if awd!=sw else 'Flooded'; c,r=run(t,0,sw),run(t,20,sw)
        if ok(c) and ok(r): W[(t,wm)]=bud(c,r)
fig,ax=plt.subplots(1,3,figsize=(15,4.3),dpi=170); x=np.arange(4); col={'Flooded':'#2a78d6','AWD':'#e8622a'}
for j,wm in enumerate(('Flooded','AWD')):
    v=[W.get((t,wm),{}) for t,_,_ in ARCH]
    ax[0].bar(x+(j-.5)*.38,[b.get('c0',np.nan) for b in v],.36,color=col[wm],label=wm)
    ax[1].bar(x+(j-.5)*.38,[b.get('dpct',np.nan) for b in v],.36,color=col[wm],label=wm)
    for i,b in enumerate(v):
        if b: stack(ax[2],i+(j-.5)*.38,b,.36); ax[2].annotate(wm[0],(i+(j-.5)*.38,0),xytext=(0,3),textcoords='offset points',ha='center',fontsize=7.5)
rat=[W[(t,'AWD')]['c0']/W[(t,'Flooded')]['c0'] for t,_,_ in ARCH if (t,'AWD') in W and (t,'Flooded') in W]
for a,ti in zip(ax,('a  Control CH$_4$ (kg ha$^{-1}$ season$^{-1}$)','b  CH$_4$ change at 20 t ha$^{-1}$, first season (%)','c  Net change at 20 t ha$^{-1}$, first year (t CO$_2$e ha$^{-1}$; F flooded, A AWD)')):
    a.set_title(ti,loc='left',fontsize=9); a.set_xticks(x,[t for t,_,_ in ARCH]); a.grid(axis='y',c='#e5e4df',lw=.6); a.set_axisbelow(True); a.axhline(0,c='#888',lw=.8)
ax[0].legend(frameon=False); ax[1].legend(frameon=False,loc='lower left')
if rat: ax[0].annotate('AWD/flooded control CH$_4$: %.2f–%.2f'%(min(rat),max(rat)),(0.02,0.97),xycoords='axes fraction',va='top',fontsize=8)
fig.legend(*LEG,loc='lower center',ncol=6,frameon=False,fontsize=8); fig.tight_layout(rect=(0,0.07,1,1)); fig.savefig(f'{FD}/fig7_awd_v22.png')
json.dump(dict(dose={t:v for t,v in res.items()},awd={f'{k[0]} {k[1]}':v for k,v in W.items()},awd_ratio=rat),open(f'{FD}/fig47_v22.json','w'),indent=1)
for t in res: print(t,{d:{k:round(v,2) for k,v in b.items()} for d,b in res[t].items()})
for k,b in W.items(): print(k,{kk:round(v,2) for kk,v in b.items()})
