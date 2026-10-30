# Fig 3 (v22): model against the five evaluation trials (numbers from Section 3, Table 3).
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt, numpy as np
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
B,S='#2a78d6','#e8622a'
fig,(a,b)=plt.subplots(1,2,figsize=(12,5.2),gridspec_kw={'width_ratios':[1,1.05]})
# (label, obs, model, colour, calibration-informed, production?)
P=[('Wang 2014 2 t/ha',-27,-35,S,1,0),('Wang 2014 4 t/ha',-43,-43,S,0,0),('Wang 2014 8 t/ha',-49,-52,S,0,0),
   ('Gong 36 t/ha',-18,-23,B,1,0),('Park 10 wt%',-25,-20,B,0,0),('Park 20 wt%',-37,-34,B,0,0),
   ('Wang 2014 production 2 t/ha',-16,-27,S,0,1),('Wang 2014 production 4 t/ha',-33,-33,S,0,1),('Wang 2014 production 8 t/ha',-41,-38,S,0,1)]
a.plot([-60,5],[-60,5],'--',c='0.55',lw=1); a.axhline(0,c='0.8',lw=1); a.axvline(0,c='0.8',lw=1)
off={'Wang 2014 2 t/ha':(5,-12),'Wang 2014 4 t/ha':(6,3),'Wang 2014 8 t/ha':(6,-10),'Gong 36 t/ha':(6,-4),'Park 10 wt%':(6,4),'Park 20 wt%':(-72,6),
     'Wang 2014 production 2 t/ha':(6,-4),'Wang 2014 production 4 t/ha':(6,4),'Wang 2014 production 8 t/ha':(-58,-3)}
for l,o,m,c,cal,pr in P:
    a.scatter(o,m,s=55,facecolor='white' if pr else c,edgecolor=c,lw=1.6,zorder=3)
    if cal: a.scatter(o,m,s=170,facecolor='none',edgecolor='k',lw=1,zorder=2)
    a.annotate(l.replace('Wang 2014 production','prod.'),(o,m),xytext=off[l],textcoords='offset points',fontsize=8.5,color='0.3')
a.errorbar(-18,-10,fmt='none'); a.annotate('',(-18,-10),(-18,-23),arrowprops=dict(arrowstyle='-',ls=':',color=B))
a.scatter(-18,-10,s=40,marker='x',c=B,zorder=3); a.annotate('Gong, no Si feedback',(-18,-10),xytext=(6,0),textcoords='offset points',fontsize=8.5,color='0.3')
from matplotlib.lines import Line2D as L
a.legend(handles=[L([],[],marker='o',ls='',mfc=B,mec=B,label='Basalt'),L([],[],marker='o',ls='',mfc=S,mec=S,label='Steel slag'),
    L([],[],marker='o',ls='',mfc='white',mec=S,label='CH$_4$ production (in-situ C$_2$H$_2$)'),
    L([],[],marker='o',ls='',ms=11,mfc='none',mec='k',label='calibration-informed')],loc='upper left',frameon=False,fontsize=9)
a.set(xlim=(-60,5),ylim=(-60,5),xlabel='Observed CH$_4$ change (%)',ylabel='Model CH$_4$ change (%)')
a.set_title('a  CH$_4$ response to amendments',loc='left'); a.grid(c='0.93')
# (b) ratios, log scale
R=[('Wang 2014 control CH$_4$\n(mg m$^{-2}$ h$^{-1}$)',2.95,3.11,None,0,'2.95 vs 3.11'),
   ('Gong control CH$_4$\n(kg ha$^{-1}$)',85,52,None,1,'85 vs 52'),
   ('Uchibayashi Na release\n100 t/ha',4.8,4.6,None,1,'4.8 vs 4.6'),
   ('Uchibayashi Mg release\n100 t/ha',12.1,6.5,None,1,'12.1 vs 6.5'),
   ('Wang 2024 inorganic C gain\n0–20 cm (t CO$_2$ ha$^{-1}$)',0.31,1.20,0.48,0,'0.31–0.48 vs 1.20'),
   ('Wang 2014 ferric Fe gain\n8 t/ha (g kg$^{-1}$)',0.5,5.6,None,0,'0.5 vs 5.6')]
y=np.arange(len(R))[::-1]
for yi,(l,m,o,m2,cal,t) in zip(y,R):
    r=m/o; c='0.55' if cal else B
    if m2: b.barh(yi,m2/o-r,left=r,height=0.5,color=c,alpha=0.35); 
    b.plot([1,r],[yi,yi],c=c,lw=2); b.scatter(r,yi,s=55,color=c,zorder=3)
    b.annotate(t,(max(1,r),yi),xytext=(8,-3),textcoords='offset points',fontsize=8.5,color='0.3')
b.axvline(1,ls='--',c='0.55',lw=1); b.set_xscale('log'); b.set_xlim(0.05,5)
b.set_xticks([0.1,0.25,0.5,1,2,4]); b.set_xticklabels(['0.1','0.25','0.5','1','2','4'])
b.set_yticks(y); b.set_yticklabels([r[0] for r in R],fontsize=9); b.set_xlabel('Model / observed (log scale)')
b.legend(handles=[L([],[],marker='o',ls='',c=B,label='prediction'),L([],[],marker='o',ls='',c='0.55',label='calibration-informed')],
    loc='upper left',frameon=False,fontsize=9)
b.set_title('b  Baselines, cation release and carbonate/Fe checks',loc='left'); b.grid(c='0.93',axis='x')
fig.tight_layout(); fig.savefig('fig3_eval_v22.png',dpi=200)
