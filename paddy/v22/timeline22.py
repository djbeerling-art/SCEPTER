# Year-by-year climate terms at one site (treatment minus year-matched control), t CO2e/ha/yr: CH4 avoided (GWP100 27,
# full and fixed at yr 1), N2O, exported PDA, calcite change, and cation release (potential); cumulative inorganic vs CH4.
import sys,json; O=sys.argv[1]; C=sys.argv[2]; T=sys.argv[3]; out=sys.argv[4]; sys.argv=['x','--days=500','--alk=pda']
exec(open('paddy/ghg_budget.py').read().split('ctrl = season')[0])
import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
AK='potentially durable ALK (t CO2/ha)'; K=lambda n: n*44/28*0.273
R=[]
for y in range(1,11):
    a={x:season(f'{O}/flood_{C}_y{y}{x}') for x in 'sf'}; b={x:season(f'{O}/flood_{T}_y{y}{x}') for x in 'sf'}
    d=lambda k: sum(b[x][k]-a[x][k] for x in 'sf')
    R.append(dict(y=y,ch4=-d('CH4 (kg/ha)')*27e-3,n2o=-K(d('N2O (kg N/ha)')),pda=d(AK),cc=d('calcite (t CO2/ha)'),
                  cat=sum(b[x]['cation release (t CO2/ha)'] for x in 'sf')))
for r in R: r['ch4fix']=R[0]['ch4']
json.dump(R,open(out.replace('.png','.json'),'w'),indent=1)
print(' yr  CH4avoid  fixed  N2O   PDA  calcite  cations')
for r in R: print(f"{r['y']:3d} {r['ch4']:7.2f} {r['ch4fix']:6.2f} {r['n2o']:5.2f} {r['pda']:5.2f} {r['cc']:7.2f} {r['cat']:7.2f}")
y=np.array([r['y'] for r in R]); f=lambda k: np.cumsum([r[k] for r in R])
fig,ax=plt.subplots(1,2,figsize=(11,4.2),dpi=170)
for k,c,l in (('ch4','#2a78d6','CH$_4$ avoided (full)'),('ch4fix','#8fb9ea','CH$_4$ avoided (fixed at yr 1)'),('pda','#1baf7a','PDA exported'),('cc','#9b59d0','calcite stored'),('cat','#888','cation release (potential CDR)')):
    ax[0].plot(y,[r[k] for r in R],'-o',ms=3,color=c,label=l); ax[1].plot(y,f(k),'-o',ms=3,color=c,label=l)
ax[1].plot(y,f('pda')+f('cc'),'--',color='#0b6b48',label='realised inorganic CDR (PDA + calcite)')
ax[0].set_title('(a) Annual, t CO$_2$e ha$^{-1}$ yr$^{-1}$',loc='left',fontsize=9); ax[1].set_title('(b) Cumulative, t CO$_2$e ha$^{-1}$',loc='left',fontsize=9)
for a in ax:
    a.set_xlabel('year'); a.grid(color='#e5e4df',lw=.6); a.tick_params(labelsize=8)
    for s in ['top','right']: a.spines[s].set_visible(False)
ax[1].legend(fontsize=7,frameon=False); fig.tight_layout(); fig.savefig(out)
