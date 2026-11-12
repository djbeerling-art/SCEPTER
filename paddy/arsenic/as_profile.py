import sys,glob,os,numpy as np
d=sys.argv[1]; zsel=float(sys.argv[2]) if len(sys.argv)>2 else 0.1
def load(pat):
    fs=sorted(glob.glob(os.path.join(d,'prof',pat)))
    out=[]
    for f in fs:
        h=open(f).readline().split(); a=np.loadtxt(f,skiprows=1)
        out.append((h,a))
    return out
aq=load('prof_aq-[0-9]*.txt'); sl=load('prof_sld-[0-9]*.txt')
print('time   pH    Si_uM   As3_uM   As5_uM   amfe3   as3ad   as5ad  (mol/m3)  totAs_mg/kg-equiv   z=%.3f'%zsel)
for (ha,a),(hs,s) in zip(aq,sl):
    i=np.argmin(abs(a[:,0]-zsel))
    g=lambda h,arr,k: arr[i,h.index(k)] if k in h else np.nan
    print('%6.3f %5.2f %7.1f %8.3f %8.4f %8.1f %8.4f %8.4f'%(g(ha,a,'time'),g(ha,a,'ph'),g(ha,a,'si')*1e6,g(ha,a,'as3')*1e6,g(ha,a,'as5')*1e6,g(hs,s,'amfe3'),g(hs,s,'as3ad'),g(hs,s,'as5ad')))
