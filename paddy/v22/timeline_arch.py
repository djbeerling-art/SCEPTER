# Cumulative climate terms over 10 yr per archetype (20 t/ha/yr vs year-matched control): CH4 avoided (full, fixed at yr 1),
# realised inorganic CDR (PDA + calcite), cation release (potential). Metric 'gwp100' (CH4 x27) or 'co2we' (GWP*, Planavsky
# et al. 2026: [4E(t) - 3.75E(t-20)] x 32; illustrative, not CDR). Usage: timeline_arch.py <outdir> <png> gwp100|co2we
import sys,json; O=sys.argv[1]; out=sys.argv[2]; M=sys.argv[3]; sys.argv=['x','--days=500','--alk=pda']
exec(open('paddy/ghg_budget.py').read().split('ctrl = season')[0])
import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
AK='potentially durable ALK (t CO2/ha)'
S=[('Reference','r10_ref_0','r22_ref_20'),('rf189 warm rainfed (23%)','r22_rf189_0','r22_rf189_20'),('irr061 irrigated (22%)','r22_irr061_0','r22_irr061_20'),
   ('irr043 acid high OM (19%)','r22_irr043_0','r22_irr043_20'),('rf184 hot calcareous (14%)','r22_rf184_0','r22_rf184_20'),
   ('irr092 cool irrigated (13%)','r22_irr092_0','r22_irr092_20'),('rf151 rainfed acid high OM (10%)','r22_rf151_0','r22_rf151_20')]
def ch4co2(E):   # E: t CH4/ha/yr difference (avoided < 0) -> avoided t CO2e/ha/yr (positive = benefit)
    if M=='gwp100': return [-27*e for e in E]
    return [-(4*E[t]-3.75*(E[t-20] if t>=20 else 0))*32 for t in range(len(E))]
fig,axs=plt.subplots(2,4,figsize=(15,7),dpi=160,sharex=True); res={}
for ax,(lab,c,r) in zip(axs.flat,S):
    E=[];inorg=[];cat=[]
    for y in range(1,11):
        a={x:season(f'{O}/flood_{c}_y{y}{x}') for x in 'sf'}; b={x:season(f'{O}/flood_{r}_y{y}{x}') for x in 'sf'}
        d=lambda k: sum(b[x][k]-a[x][k] for x in 'sf')
        E.append(d('CH4 (kg/ha)')*1e-3); inorg.append(d(AK)+d('calcite (t CO2/ha)')); cat.append(sum(b[x]['cation release (t CO2/ha)'] for x in 'sf'))
    full=ch4co2(E); fix=ch4co2([E[0]]*10); y=np.arange(1,11)
    res[lab]=dict(ch4_full=full,ch4_fix=fix,inorganic=inorg,cations=cat)
    ax.plot(y,np.cumsum(full),'-o',ms=3,color='#2a78d6',label='CH$_4$ avoided (full)')
    ax.plot(y,np.cumsum(fix),'--',color='#2a78d6',label='CH$_4$ avoided (fixed at yr 1)')
    ax.plot(y,np.cumsum(inorg),'-o',ms=3,color='#1baf7a',label='realised inorganic CDR (PDA + calcite)')
    ax.plot(y,np.cumsum(cat),':',color='#888',label='cation release (potential CDR)')
    ax.set_title(lab,fontsize=9,loc='left'); ax.grid(color='#e5e4df',lw=.6); ax.tick_params(labelsize=8)
    for s in ['top','right']: ax.spines[s].set_visible(False)
axs.flat[-1].axis('off'); h,l=axs.flat[0].get_legend_handles_labels(); axs.flat[-1].legend(h,l,fontsize=8,frameon=False,loc='center')
unit='t CO$_2$e ha$^{-1}$ (GWP100)' if M=='gwp100' else 't CO$_2$-we ha$^{-1}$ (GWP*, illustrative)'
fig.suptitle(f'Cumulative terms, 20 t ha$^{{-1}}$ yr$^{{-1}}$ basalt vs year-matched control, {unit}',fontsize=10)
for a in axs[1]: a.set_xlabel('year',fontsize=8)
fig.tight_layout(); fig.savefig(out); json.dump(res,open(out.replace('.png','.json'),'w'),indent=1)
