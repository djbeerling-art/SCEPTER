"""Goethite x silica test: absolute and relative As(III) mobilisation by Si, and the diagnostic split of sorbed As
between ferrihydrite and goethite sites (equilibrium shares of the single sorbed pool, as in as3sp/as5sp).
Usage: python3 paddy/arsenic/gt_si_partition.py <dir with flood_gt<g> and flood_gtSi<g> runs> [g ...]
Root zone 0-0.2 m; inventories in mmol m-2; Kd = sorbed/dissolved (dimensionless, bulk)."""
import sys
import glob,os,numpy as np,json
_trap=getattr(np,"trapezoid",None) or getattr(np,"trapz")
C={'as3':(5.1882,-0.0224,-0.0171,3.8028,0.0104,0.8685),'as5':(7.2189,-1.0735,-0.2022,9.3859,0.1815,1.9115)}
FGT,FS=0.02,1.0
def ld(d,p):
    out=[]
    for f in sorted(glob.glob(os.path.join(d,'prof',p))):
        h=open(f).readline().split(); x=np.loadtxt(f,skiprows=1)
        if out and len(h)!=x.shape[1]: h=out[-1][0]
        out.append((h,x))
    if out and len(out)>1 and len(out[0][0])!=out[0][1].shape[1]: out[0]=(out[1][0],out[0][1])
    return out
def run(d):
    aq,sl,bs=ld(d,'prof_aq-[0-9]*.txt'),ld(d,'prof_sld-[0-9]*.txt'),ld(d,'bsd-[0-9]*.txt')
    rows=[]
    for (ha,a),(hs,s),(hb,b) in zip(aq,sl,bs):
        z=a[:,0]; m=z<=0.2; dz=0.02
        g=lambda h,x,k:x[:,h.index(k)]
        ph=np.clip(g(ha,a,'ph'),5,8); si=np.clip(g(ha,a,'si'),2e-5,2e-3)
        th=g(hb,b,'poro')*g(hb,b,'sat')
        fh=FS*0.2*g(hs,s,'amfe3'); gt=FGT*g(hs,s,'gt')
        r={'t':a[0,ha.index('time')]}
        for sp,ad in (('as3','as3ad'),('as5','as5ad')):
            A,B,Cc,aS,bS,n=C[sp]; x=ph-7
            K=10**(A+B*x+Cc*x*x); D=1+10**(aS+bS*x)*si**n
            sfh=fh/D; sh=sfh/(sfh+gt)
            sorb=g(hs,s,ad); dis=g(ha,a,sp)*1e3*th   # mol/m3 bulk
            r[sp+'_dis']=(dis*m).sum()*dz; r[sp+'_fh']=(sorb*sh*m).sum()*dz; r[sp+'_gt']=(sorb*(1-sh)*m).sum()*dz
            r[sp+'_conc']=g(ha,a,sp)[m].mean()
            r[sp+'_Kfh']=(K*sfh/(1e3*th))[m].mean(); r[sp+'_Kgt']=(K*gt/(1e3*th))[m].mean()  # dimensionless sorbed/dissolved
        r['si']=g(ha,a,'si')[m].mean(); r['ph']=g(ha,a,'ph')[m].mean()
        rows.append(r)
    return rows
res={}
os.chdir(sys.argv[1])
for gtl in (sys.argv[2:] or ['0.0001','0.02','0.06']):
    c,b=run(f'flood_gt{gtl}'),run(f'flood_gtSi{gtl}')
    T=lambda rr,k: np.array([r[k] for r in rr])
    t=T(c,'t')
    def tm(rr,k):
        v=T(rr,k); return _trap(v,t)/(t[-1]-t[0])
    def at(rr,k,i=-1): return T(rr,k)[i]
    o={}
    for lab,rr in (('C',c),('Si',b)):
        o[lab]=dict(si_uM=tm(rr,'si')*1e6,pH=tm(rr,'ph'),as3_uM_mean=tm(rr,'as3_conc')*1e6,as3_uM_t0=at(rr,'as3_conc',0)*1e6,
          as5_uM_t0=at(rr,'as5_conc',0)*1e6,as3_uM_end=at(rr,'as3_conc')*1e6,
          Kd3_fh_end=at(rr,'as3_Kfh'),Kd3_gt_end=at(rr,'as3_Kgt'),Kd3_fh_t0=at(rr,'as3_Kfh',0),Kd3_gt_t0=at(rr,'as3_Kgt',0),
          **{f'{k}_end_mmolm2':at(rr,k)*1e3 for k in ('as3_dis','as3_fh','as3_gt','as5_dis','as5_fh','as5_gt')},
          **{f'{k}_t0_mmolm2':at(rr,k,0)*1e3 for k in ('as3_dis','as3_fh','as3_gt','as5_dis','as5_fh','as5_gt')})
    o['dAs3_mean_uM']=o['Si']['as3_uM_mean']-o['C']['as3_uM_mean']; o['R_As_mean']=o['Si']['as3_uM_mean']/o['C']['as3_uM_mean']
    o['dAs3_end_uM']=o['Si']['as3_uM_end']-o['C']['as3_uM_end']; o['R_As_end']=o['Si']['as3_uM_end']/o['C']['as3_uM_end']
    # capacity-only prediction: same total As(III), redistributed by Kd: R = (1+Kd_C)/(1+Kd_Si)
    kc=1+o['C']['Kd3_fh_end']+o['C']['Kd3_gt_end']; ks=1+o['Si']['Kd3_fh_end']+o['Si']['Kd3_gt_end']
    o['R_As_capacity_pred']=kc/ks
    res[gtl]=o
print(json.dumps(res,indent=1,default=float))
