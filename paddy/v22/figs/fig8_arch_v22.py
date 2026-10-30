# v22 (frozen) 10-yr year-round (season + fallow) GHG and CDR budget per site, vs year-matched control
import sys,os,json; O=sys.argv[1]; out=sys.argv[2]; sys.argv=['x','--days=500']+sys.argv[3:]
exec(open('paddy/ghg_budget.py').read().split('ctrl = season')[0])
import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
SITES=[]
for code,lab in [('r22_ref','Reference'),('r22_rf189','rf189 warm rainfed (23%)'),('r22_irr061','irr061 irrigated (22%)'),('r22_irr043','irr043 acid high OM (19%)'),
                 ('r22_rf184','rf184 hot calcareous (14%)'),('r22_irr092','irr092 cool irrigated (13%)'),('r22_rf151','rf151 rainfed acid high OM (10%)')]:
    for d in (10,20): SITES.append((f'{lab}, {d} t','r10_ref_0' if code=='r22_ref' else f'{code}_0',f'{code}_{d}'))  # ref control: v21 (no rock, identical physics)
K=lambda n: n*44/28*0.273
AK='potentially durable ALK (t CO2/ha)' if '--alk=pda' in sys.argv else 'ALK export (t CO2/ha)'
rows=[]
SITES=[x for x in SITES if x[0].endswith('20 t') and not x[0].startswith('Reference') and os.path.isdir(f'{O}/flood_{x[2]}_y10f') and os.path.isdir(f'{O}/flood_{x[1]}_y10f')]
for lab,c,r in SITES:
    t=dict(ch4s=0,ch4f=0,n2o=0,alk=0,cc=0,cat=0,ny=0)
    for y in range(1,11):
        p={k:{x:season(f'{O}/flood_{n}_y{y}{x}') for x in 'sf'} for k,n in (('c',c),('b',r))}
        a,b=p['c'],p['b']
        t['ch4s']+=(b['s']['CH4 (kg/ha)']-a['s']['CH4 (kg/ha)'])*27e-3
        if y==1: t['ch4y1']=sum(b[x]['CH4 (kg/ha)']-a[x]['CH4 (kg/ha)'] for x in 'sf')*27e-3   # yr-1 CH4 difference (season+fallow)
        t['ch4f']+=(b['f']['CH4 (kg/ha)']-a['f']['CH4 (kg/ha)'])*27e-3
        t['n2o']+=K(sum(b[x]['N2O (kg N/ha)']-a[x]['N2O (kg N/ha)'] for x in 'sf'))
        dalk=sum(b[x][AK]-a[x][AK] for x in 'sf')
        cat=sum(b[x]['cation release (t CO2/ha)'] for x in 'sf')
        dcc=sum(b[x]['calcite (t CO2/ha)']-a[x]['calcite (t CO2/ha)'] for x in 'sf')
        t['alk']+=dalk; t['cat']+=cat; t['cc']+=dcc; t['ny']+=1
    # calcite counted as CDR only up to the rock cations not exported as alkalinity (cumulative cap)
    t['cc']=min(t['cc'],max(t['cat']-t['alk'],0)); t['ret']=t['cat']-t['alk']-t['cc']
    t['net']=t['ch4s']+t['ch4f']+t['n2o']-t['alk']-t['cc']
    # diagnostic sensitivity: CH4 difference held at its year-1 value (N2O and inorganic C unchanged); not a simulated Fe-ageing scenario
    t['netfix']=t['ch4y1']*t['ny']+t['n2o']-t['alk']-t['cc']; t['netpot']=t['ch4s']+t['ch4f']+t['n2o']-t['cat']
    rows.append((lab,t)); print(lab,{k:round(v,1) for k,v in t.items()})
json.dump(rows,open(out.replace('.png','.json'),'w'),indent=1)
fig,ax=plt.subplots(1,2,figsize=(13,2.2+0.45*len(SITES)),dpi=170,sharey=True); y=np.arange(len(rows))[::-1]
C=dict(ch4s='#2a78d6',ch4f='#8fb9ea',n2o='#d6a10f',alk='#1baf7a',cc='#9b59d0',ret='#c9c7c0')
for i,(lab,t) in zip(y,rows):
    pos=neg=0
    for k,s in (('ch4s',1),('ch4f',1),('n2o',1),('alk',-1),('cc',-1)):
        v=s*t[k]; l=pos if v>=0 else neg
        ax[0].barh(i,v,left=l,height=0.62,color=C[k],edgecolor='white',linewidth=1)
        if v>=0: pos+=v
        else: neg+=v
    ax[0].plot(t['net'],i,'D',color='#222',ms=6,zorder=5); ax[0].plot(t['netfix'],i,'|',color='#c0392b',ms=14,mew=2.2,zorder=6); ax[0].plot(t['netpot'],i,'o',mfc='white',mec='#222',ms=6,zorder=5)
    ax[0].annotate('%+.0f'%t['net'],(t['net'],i),xytext=(0,9),textcoords='offset points',ha='center',fontsize=7.5)
    l=0
    for k in ('alk','cc','ret'):
        ax[1].barh(i,t[k],left=l,height=0.62,color=C[k],edgecolor='white',linewidth=1); l+=t[k]
    ax[1].annotate('%.0f'%t['cat'],(t['cat'],i),xytext=(4,0),textcoords='offset points',va='center',fontsize=7.5)
ax[0].set_yticks(y,[l for l,_ in rows],fontsize=8); ax[0].axvline(0,color='#888',lw=.8)
ax[0].set_title('(a) Net GHG change, 10 yr, season + fallow (t CO$_2$e ha$^{-1}$)',fontsize=9,loc='left')
ax[1].set_title('(b) Where the released rock cations went, 10 yr (t CO$_2$ ha$^{-1}$)',fontsize=9,loc='left')
for a in ax:
    a.grid(axis='x',color='#e5e4df',lw=.6); a.set_axisbelow(True); a.tick_params(labelsize=8)
    for s in ['top','right']: a.spines[s].set_visible(False)
from matplotlib.patches import Patch; from matplotlib.lines import Line2D
ax[0].legend([Patch(color=C[k]) for k in ('ch4s','ch4f','n2o','alk','cc')]+[Line2D([],[],marker='D',color='#222',ls=''),Line2D([],[],marker='o',mfc='white',mec='#222',ls=''),Line2D([],[],marker='|',color='#c0392b',mew=2.2,ms=12,ls='')],
  ['ΔCH$_4$ season','ΔCH$_4$ fallow','ΔN$_2$O (all year)',('Potentially durable alkalinity exported' if 'potentially' in AK else 'Alkalinity exported (CDR)'),'Calcite stored (soil carbonate)','Net, full CH$_4$ trajectory','Net if all released cations were counted (potential)','Net with CH$_4$ difference fixed at year 1 (sensitivity, not a bound)'],fontsize=7.5,frameon=False,loc='upper center',bbox_to_anchor=(0.5,-0.06),ncol=3)
ax[1].legend([Patch(color=C[k]) for k in ('alk','cc','ret')],[('Exported as potentially durable alkalinity' if 'potentially' in AK else 'Exported as alkalinity'),'Stored as calcite','Still in the soil (exchange sites, porewater)'],fontsize=7.5,frameon=False,loc='upper center',bbox_to_anchor=(0.5,-0.06),ncol=1)
pass
fig.tight_layout(rect=(0,0.04,1,1)); fig.savefig(out)
