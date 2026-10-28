# Fig 3 draft: verification of the paddy model code.
# a  PHREEQC example 11 cation-exchange benchmark (outlet breakthrough), base SCEPTER code
# b  pH of model porewaters re-speciated in PHREEQC (phreeqc.dat, charge balance) vs the model's own pH
# c  mineral saturation, model log Omega vs PHREEQC SI (activity coefficients removed, as the model runs ideal)
# d  element mass-balance closure for season and fallow runs (irr061, flooded and AWD, 0 and 20 t/ha)
# Data: ex11.json (ex11_data.py), phreeqc_respeciation.json (phq.py), mass_balance.json (massbal2.py)
import json,os,numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
D=os.path.dirname(os.path.abspath(__file__))
ex=json.load(open(f'{D}/ex11.json')); ph=json.load(open(f'{D}/phreeqc_respeciation.json')); mb=json.load(open(f'{D}/mass_balance.json'))
INK='#52514e'; GRID='#ecebe7'; REF='#9a9992'
fig,ax=plt.subplots(2,2,figsize=(10.5,8.4),dpi=170)

# a: ex11
a=ax[0,0]; C={'na':'#2a78d6','k':'#eb6834','ca':'#1baf7a','cl':'#8a5cd0'}; N={'na':'Na','k':'K','ca':'Ca','cl':'Cl'}
for n,c in C.items():
    a.plot(ex['pv_s'],ex['s'][n],color=c,lw=2,label=N[n])
    p=np.array(ex['pv_p']); v=np.array(ex['p'][n]); k=slice(0,None,3)
    a.scatter(p[k],v[k],s=16,facecolor='none',edgecolor=c,linewidth=1,zorder=3)
a.set_xlim(0,2.6); a.set_xlabel('Pore volumes',fontsize=8); a.set_ylabel('Outlet concentration (mM)',fontsize=8)
pvg=np.linspace(0.3,2.5,100); rm=max(np.sqrt(np.mean((np.interp(pvg,ex['pv_s'],ex['s'][n])-np.interp(pvg,ex['pv_p'],ex['p'][n]))**2)) for n in C)
a.set_title('a  Cation-exchange benchmark (PHREEQC ex. 11)',fontsize=9,loc='left')
h=[Line2D([],[],color=C[n],lw=2) for n in C]+[Line2D([],[],color=INK,lw=2),Line2D([],[],marker='o',ls='',mfc='none',mec=INK)]
a.legend(h,[N[n] for n in C]+['SCEPTER','PHREEQC'],fontsize=7,frameon=False,ncol=2,loc='center left',bbox_to_anchor=(0.0,0.6))
a.annotate(f'RMSE ≤ {rm:.3f} mM',(0.03,0.95),xycoords='axes fraction',fontsize=7.5,color=INK,va='top')

# b, c: re-speciation
G={'Flooded':['flood_ctrl','flood_b20','flood_b20_late'],'AWD, drained':['awd_drained'],'AWD, reflooded':['awd_reflooded'],'Fallow':['fallow_b20','fallow_ctrl_end']}
GC=dict(zip(G,['#2a78d6','#eb6834','#1baf7a','#8a5cd0']))
grp=lambda tag: next(g for g,t in G.items() if tag in t)
b=ax[0,1]; b.plot([5,8.2],[5,8.2],color=REF,lw=.8,ls='--',zorder=1)
for g in G:
    r=[x for x in ph if grp(x['tag'])==g]
    b.scatter([x['ph_s'] for x in r],[x['ph_p'] for x in r],s=22,color=GC[g],edgecolor='#fcfcfb',linewidth=.8,label=g,zorder=3)
d=max(abs(x['ph_p']-x['ph_s']) for x in ph)
b.set_xlim(5.2,8); b.set_ylim(5.2,8); b.set_xlabel('Model pH',fontsize=8); b.set_ylabel('PHREEQC pH (charge balance)',fontsize=8)
b.set_title('b  Porewater pH re-speciated in PHREEQC',fontsize=9,loc='left'); b.legend(fontsize=7,frameon=False,loc='upper left')
b.annotate(f'n = {len(ph)}, max |ΔpH| = {d:.2f}',(0.97,0.04),xycoords='axes fraction',ha='right',fontsize=7.5,color=INK)

c=ax[1,0]; c.plot([-4,2.5],[-4,2.5],color=REF,lw=.8,ls='--',zorder=1)
M=[('Calcite','cc_s','cc_pi',0,'o'),('Siderite (−1)','sid_s','sid_pi',1,'s'),('Fe(OH)$_3$(a)','fh_s','fh_p',0,'^'),('SiO$_2$(a)','amsi_s','amsi_p',0,'D')]
MC=['#2a78d6','#eb6834','#1baf7a','#8a5cd0']
for (l,ks,kp,sh,mk),col in zip(M,MC):
    xs=np.array([x[ks] for x in ph]); ys=np.array([x[kp] for x in ph])-sh; ok=xs>-6
    c.scatter(xs[ok],ys[ok],s=18,marker=mk,color=col,edgecolor='#fcfcfb',linewidth=.6,label=f'{l}: mean Δ {np.mean(ys[ok]-xs[ok]):+.2f}',zorder=3)
c.set_xlim(-4,2.5); c.set_ylim(-4,2.5); c.set_xlabel('Model log Ω',fontsize=8); c.set_ylabel('PHREEQC SI (ideal solution)',fontsize=8)
c.set_title('c  Mineral saturation',fontsize=9,loc='left'); c.legend(fontsize=7,frameon=False,loc='upper left')
c.annotate('Siderite shifted by the model\'s +1 log-unit\napparent solubility. With activity\ncorrections, calcite SI is 0.15–0.26 lower.',
           (0.97,0.04),xycoords='axes fraction',ha='right',fontsize=6.8,color=INK)

# d: closure
e=ax[1,1]; EL=['C','Fe','N','Ca','Mg']; x=np.arange(len(EL)); w=0.36
rc=[max(abs(v['budget'][el]['internal'])/v['budget'][el]['gross'] for v in mb.values()) for el in EL]
sc=[max(abs(v['storage'][el]['dinv']-v['storage'][el]['dtflx'])/v['budget'][el]['gross'] for v in mb.values()) for el in EL]
e.bar(x-w/2-0.01,rc,w,color='#2a78d6',label='Reaction stoichiometry (internal terms)')
e.bar(x+w/2+0.01,sc,w,color='#eb6834',label='Storage vs profile inventory')
e.set_yscale('log'); e.set_ylim(1e-17,1e-1); e.set_xticks(x,EL,fontsize=8)
e.set_ylabel('Worst relative imbalance (of gross flux)',fontsize=8)
e.set_title('d  Element mass balance (8 season/fallow runs)',fontsize=9,loc='left'); e.legend(fontsize=7,frameon=False,loc='upper right')
e.annotate('C: untracked oxalate +\nacetate, 4% of\naerobic OM-C',(x[0]+0.02,2e-4),fontsize=6.8,color=INK,va='top')
for ax_ in ax.flat:
    ax_.tick_params(labelsize=7.5); ax_.grid(color=GRID,lw=.6); ax_.set_axisbelow(True)
    for s in ['top','right']: ax_.spines[s].set_visible(False)
fig.tight_layout(); fig.savefig(f'{D}/fig3_verification_draft.png')
print('ex11 rmse',rm,'max dpH',d,'closure',dict(zip(EL,rc)),'storage',dict(zip(EL,sc)))
