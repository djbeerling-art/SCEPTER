# Fig 6a (v22): Morris screen, 11 factors, r = 8 (96 points), 1 year after a 7-yr settle, basalt run minus year-matched control.
# Baseline CH4 (control) ranked separately from the basalt-induced change. Data: morris22/ (paddy/screen/morris22.py analyse)
import os,pandas as pd,numpy as np,matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import Patch; from matplotlib.lines import Line2D
D=os.path.dirname(os.path.abspath(__file__)); M=f'{D}/morris22'
plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False})
L={'soilfe3':'Soil ferrihydrite (0.4–3 wt%)','om':'OM input (250–900 g C m⁻² yr⁻¹)','irralk':'Irrigation alkalinity, soil pH (~5.0–6.7)','cec':'CEC (8–30 cmol kg⁻¹)',
   'temp':'Temperature (20–32 °C)','basfe3':'Rock Fe(III) oxide (0–5 wt%)','rock':'Rock type (5 basalts)','dose':'Basalt dose (5–30 t ha⁻¹)',
   'awd':'Water: flooded vs AWD','siw':'Si feedback strength w (0.1–0.4)','q':'Percolation (0.25–1.1 m yr⁻¹)'}
TR={'basfe3','rock','dose'}
P=[('base_ch4','a1  Baseline: control CH$_4$ (kg ha$^{-1}$ season$^{-1}$)'),('ch4_pct','a2  Basalt effect: season CH$_4$ change (%)'),('net','a3  Basalt effect: 1-yr net GHG (t CO$_2$e ha$^{-1}$)')]
fig,ax=plt.subplots(1,3,figsize=(15,4.6),dpi=170)
for a,(m,t) in zip(ax,P):
    s=pd.read_csv(f'{M}/morris_{m}.csv')
    if m=='base_ch4': s=s[~s.names.isin(TR)]
    s=s.sort_values('mu_star'); y=np.arange(len(s))
    a.barh(y,s.mu_star,xerr=s.mu_star_conf,color=np.where(s.mu<0,'#2a78d6','#eb6834'),height=.65,error_kw=dict(lw=.7,ecolor='#52514e'))
    a.plot(s.sigma,y,'o',mfc='none',mec='#222',ms=5)
    a.set_yticks(y,[L[n] for n in s.names],fontsize=7.5); a.set_title(t,loc='left',fontsize=9); a.grid(axis='x',c='#ecebe7',lw=.6); a.set_axisbelow(True)
    a.set_xlabel('μ* (bars, 95% CI) and σ (circles)',fontsize=8)
ax[0].annotate('Dose, rock type and rock Fe(III)\ndo not enter the control runs',(0.98,0.04),xycoords='axes fraction',ha='right',fontsize=7,color='#52514e')
r=pd.read_csv(f'{M}/results.csv')
ax[2].annotate(f'CH$_4$ fell in {(r.ch4_pct<0).sum()}/{len(r)} runs (median {r.ch4_pct.median():.0f}%);\nnet negative in {(r.net<0).sum()}/{len(r)} (median {r.net.median():.2f})',(0.98,0.04),xycoords='axes fraction',ha='right',fontsize=7,color='#52514e')
fig.legend([Patch(color='#eb6834'),Patch(color='#2a78d6'),Line2D([],[],marker='o',ls='',mfc='none',mec='#222')],['raises the metric (μ > 0)','lowers the metric (μ < 0)','σ: interactions / non-linearity'],
           frameon=False,fontsize=7.5,ncol=3,loc='lower center',bbox_to_anchor=(0.5,0))
fig.suptitle('Model version: scepter_paddy_v22 (frozen; K$_{Si}$ 2 mM)',x=0.99,ha='right',fontsize=7.5,color='#52514e')
fig.tight_layout(rect=(0,0.05,1,0.97)); fig.savefig(f'{D}/fig6a_morris_v22.png')
