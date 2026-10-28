# Element budgets and closure checks for SCEPTER paddy runs.
# (1) Reaction closure: in each element budget (sum over species of stoichiometry x column-integrated cumulative flux),
#     every internal reaction column should cancel; open terms are transport, inputs, storage and the plant/N2 terms.
# (2) Storage check: change of the element inventory in the profiles (solid + aqueous + exchangeable; gases excluded)
#     against the integrated storage term (tflx) of the same species over the same interval.
import sys,json,glob,numpy as np
FR_AN=0.95
E={'C':{'sld-g2':1,'sld-cc':1,'sld-sid':1,'sld-urea':1,'gas-pco2':1,'gas-pch4':1},
   'Fe':{'aq-fe2':1,'aq-fe3':1,'sld-amfe3':1,'sld-gt':1,'sld-sid':1,'sld-gr':3,'sld-fes':1,'sld-fa':2,'sld-ilmn':1,'sld-fer':1},
   'N':{'sld-urea':2,'sld-g2':'n2c','gas-pnh3':1,'aq-no3':1,'gas-pn2o':2},
   'Ca':{'aq-ca':1,'sld-cc':1,'sld-an':FR_AN,'sld-dp':1},
   'Mg':{'aq-mg':1,'sld-fo':2,'sld-dp':1,'sld-en':1}}
OPEN=['tflx','adv','dif','rain','rtco2','nh4up','siup','ch4pl','ch4eb','n2opl','g2n2p','g2n0','rtoxy']
def L(f):
    h=open(f).readline().split(); return h,np.loadtxt(f,skiprows=1)
def cum(run,sp):
    h,d=L(f'{run}/flx/int_flx_{sp}.txt'); return h,d[:,0],d[:,1:]*d[:,:1]   # cumulative mol/m2 at each output time
def budget(run,n2c=0.03):
    out={}
    for el,sp in E.items():
        tot={}; t=None
        for s,nu in sp.items():
            nu=n2c if nu=='n2c' else nu
            try: h,t,c=cum(run,s)
            except OSError: continue
            for i,k in enumerate(h[1:]): tot[k]=tot.get(k,0)+nu*c[-1,i]
        gross=sum(abs(v) for k,v in tot.items() if k!='tflx' and k!='res')
        internal={k:v for k,v in tot.items() if k not in OPEN}
        out[el]=dict(terms={k:v for k,v in tot.items() if abs(v)>1e-9},gross=gross,
                     internal=sum(internal.values()),worst=max(internal,key=lambda k:abs(internal[k])),t=float(t[-1]))
    return out
def storage(run,n2c=0.03):
    P=f'{run}/prof'; n=max(int(f[-7:-4]) for f in glob.glob(f'{P}/prof_sld-[0-9][0-9][0-9].txt'))
    hs0,s0=L(f'{P}/prof_sld-001.txt'); hs1,s1=L(f'{P}/prof_sld-{n:03d}.txt')
    hx0,x0=L(f'{P}/prof_ex(tot)-001.txt'); hx1,x1=L(f'{P}/prof_ex(tot)-{n:03d}.txt')
    t0,t1=s0[0,-1],s1[0,-1]; dz=np.diff(np.concatenate([[0],(s0[:-1,0]+s0[1:,0])/2,[2*s0[-1,0]-(s0[-1,0]+s0[-2,0])/2]]))
    out={}
    for el,sp in E.items():
        dinv=dtf=0.
        for s,nu in sp.items():
            nu=n2c if nu=='n2c' else nu
            kind,name=s.split('-')
            if kind=='gas': continue
            if kind=='sld': h0,a0,h1,a1=hs0,s0,hs1,s1
            else: h0,a0,h1,a1=hx0,x0,hx1,x1
            if name not in h0: continue
            dinv+=nu*((a1[:,h1.index(name)]-a0[:,h0.index(name)])*dz).sum()
            h,t,c=cum(run,s); tf=np.interp([t0,t1],t,c[:,h.index('tflx')-1]); dtf+=nu*(tf[1]-tf[0])
        out[el]=dict(dinv=dinv,dtflx=dtf)
    return out
if __name__=='__main__':
    res={}
    for run in sys.argv[2:]:
        tag=run.rstrip('/').split('flood_')[-1]; b=budget(run); s=storage(run); res[tag]=dict(budget=b,storage=s)
        for el in E:
            x=b[el]; y=s[el]
            print(f"{tag:16s} {el:2s} gross {x['gross']:8.3f} internal {x['internal']:+.2e} ({abs(x['internal'])/x['gross']:.1e}, worst {x['worst']}) "
                  f"storage dInv {y['dinv']:+8.4f} vs tflx {y['dtflx']:+8.4f} err {abs(y['dinv']-y['dtflx'])/max(x['gross'],1e-12):.1e}")
    json.dump(res,open(sys.argv[1],'w'),indent=0)
