# Fig 5 (v22): paddy biogeochemical dynamics in the reference control column (no rock), one year (season + fallow).
# Usage: python3 fig5_dynamics_v22.py <outdir> <run> <10-yr control run> <png>   e.g. ... flood_r22y1_ref_0 r10_ref_0 fig5.png
import sys as _s; _s.path.insert(0,"/mnt/project-files/riceland-scepter/v22/figs"); import bigfont
import sys,glob,numpy as np,matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
O,R,C10,out=sys.argv[1:5]
plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False})
def L(f): h=open(f).readline().split(); return h,np.loadtxt(f,skiprows=1)
def series(run,kind,t0):
    fs=[x for x in sorted(glob.glob(f'{O}/flood_{run}/prof/{kind}-[0-9][0-9][0-9].txt')) if not x.endswith('-000.txt')]; T=[];Z=None;D={}
    for f in fs:
        h,d=L(f); T.append(d[0,-1]+t0); Z=d[:,0]
        for i,k in enumerate(h[1:-1],1): D.setdefault(k,[]).append(d[:,i])
    return np.array(T)*365,Z,{k:np.array(v).T for k,v in D.items()}
def year(kind):
    t1,z,a=series(f'{R}_y1s',kind,0); t2,_,b=series(f'{R}_y1f',kind,0.31)
    return np.concatenate([t1,t2]),z,{k:np.hstack([a[k],b[k]]) for k in a}
def flux(run,sp,t0):
    h,d=L(f'{O}/flood_{run}/flx/flx_gas-{sp}.txt'); c={k:d[:,i] for i,k in enumerate(h)}
    e=c['dif']+c.get('ch4eb',0)+c.get('ch4pl',0)+c.get('n2opl',0)
    return (c['time']+t0)*365,e
fig=plt.figure(figsize=(7.4,10)); gs=fig.add_gridspec(3,2,hspace=0.5,wspace=0.55)
tg,z,G=year('prof_gas(aq)'); ta,_,A=year('prof_aq')
def heat(ax,t,v,title,cmap,lab,log=False,vmin=None,vmax=None):
    from matplotlib.colors import LogNorm
    m=ax.pcolormesh(t,z,v,shading='nearest',cmap=cmap,norm=LogNorm(vmin=vmin,vmax=vmax) if log else None,vmin=None if log else vmin,vmax=None if log else vmax)
    ax.set_ylim(0.5,0); ax.axhspan(0.20,0.24,color='w',alpha=0.0); [ax.axhline(y,c='w',lw=.8,ls='--') for y in (0.2,0.24)]
    ax.axvline(113,c='k',lw=.8); ax.set_xlabel('Day of year'); ax.set_ylabel('Depth (m)'); ax.set_title(title,loc='left',fontsize=9)
    cb=fig.colorbar(m,ax=ax,pad=0.02); cb.set_label(lab,fontsize=8)
heat(fig.add_subplot(gs[0,0]),tg,np.maximum(G['po2']*1e6,1e-3),'a  Dissolved O$_2$','Blues','µmol L$^{-1}$',True,1e-2,300)
heat(fig.add_subplot(gs[0,1]),ta,A['fe2']*1e3,'b  Porewater Fe$^{2+}$','Oranges','mmol L$^{-1}$')
heat(fig.add_subplot(gs[1,0]),ta,A['ph'],'c  Porewater pH','viridis','pH')
heat(fig.add_subplot(gs[1,1]),tg,np.maximum(G['pch4']*1e6,1e-2),'d  Dissolved CH$_4$','Reds','µmol L$^{-1}$',True,1,2000)
ax=fig.add_subplot(gs[2,0]); t1,e1=flux(f'{R}_y1s','pch4',0); t2,e2=flux(f'{R}_y1f','pch4',0.31)
ax.plot(np.r_[t1,t2],np.r_[e1,e2]*16.04e3/8760,c='#d6542a',lw=1.4,label='CH$_4$')
ax.set_ylabel('CH$_4$ emission (mg m$^{-2}$ h$^{-1}$)',color='#d6542a'); ax.set_xlabel('Day of year'); ax.axvline(113,c='k',lw=.8)
b=ax.twinx(); t1,n1=flux(f'{R}_y1s','pn2o',0); t2,n2=flux(f'{R}_y1f','pn2o',0.31)
b.plot(np.r_[t1,t2],np.r_[n1,n2]*28.01*10/365*1e3,c='#6a3fb0',lw=1.2); b.set_ylabel('N$_2$O emission (g N ha$^{-1}$ d$^{-1}$)',color='#6a3fb0'); b.spines['right'].set_visible(True)
ax.set_title('e  Emissions',loc='left',fontsize=9)
ax.annotate('season',(56,0.97),xycoords=('data','axes fraction'),ha='center',fontsize=8,va='top'); ax.annotate('fallow',(240,0.97),xycoords=('data','axes fraction'),ha='center',fontsize=8,va='top')
ax=fig.add_subplot(gs[2,1])
def fe3(run,last=True):
    fs=[x for x in sorted(glob.glob(f'{O}/flood_{run}/prof/prof_sld-[0-9][0-9][0-9].txt')) if not x.endswith('-000.txt')]; h,d=L(fs[-1] if last else fs[0])
    return d[:,0],d[:,h.index('amfe3')],d[:,h.index('gt')]
for run,last,lab,c,ls in ((f'{R}_y1s',False,'start of year','#1b1b1b','-'),(f'{R}_y1s',True,'end of season','#d6542a','-'),
                          (f'{R}_y1f',True,'end of fallow','#2a78d6','-'),(f'{C10}_y10f',True,'after 10 more years','#1b1b1b','--')):
    zz,fh,gt=fe3(run,last); ax.plot(fh,zz,c=c,ls=ls,lw=1.5,label=lab)
ax.axhspan(0.2,0.24,color='#8f7b5f',alpha=0.25,lw=0); ax.annotate('plough pan',(0.02,0.22),xycoords=('axes fraction','data'),ha='left',va='center',fontsize=8)
ax.set_ylim(0.5,0); ax.set_xlabel('Ferrihydrite (mol m$^{-3}$)'); ax.set_ylabel('Depth (m)'); ax.legend(frameon=False,fontsize=7.5,loc='lower right')
ax.set_title('f  Fe(III) across the plough pan',loc='left',fontsize=9)
fig.savefig(out,bbox_inches='tight')
