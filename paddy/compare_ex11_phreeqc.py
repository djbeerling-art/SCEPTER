# Compare the SCEPTER cation-exchange benchmark (test_phreeqc_ex11.py, init then rstrt)
# with PHREEQC example 11 run through phreeqpython, at the column outlet.
# Usage: SCEPTER_OUT=<outdir used for the ex11 runs> python3 paddy/compare_ex11_phreeqc.py
import glob,os
import numpy as np
import phreeqpython  # pip install phreeqpython
O=os.environ.get('SCEPTER_OUT','./scepter_output/')+'phreeqc_ex/ex11_test_rstrt/prof/'
rows=[]
for f in sorted(glob.glob(O+'prof_aq-[0-9][0-9][0-9].txt')):
    d=np.loadtxt(f,skiprows=1); rows.append(d[-1])
s=np.array(rows)  # z na k ca mg no3 cl ph time
pv_s=s[:,8]/0.25
# PHREEQC ex11 (standard, phreeqc.dat)
pp=phreeqpython.PhreeqPython(database='phreeqc.dat')
ip=pp.ip
inp='''
SOLUTION 0  CaCl2
        units            mmol/kgw
        temp             25.0
        pH               7.0     charge
        pe               12.5    O2(g)   -0.68
        Ca               0.6
        Cl               1.2
SOLUTION 1-40  Initial solution for column
        units            mmol/kgw
        temp             25.0
        pH               7.0     charge
        pe               12.5    O2(g)   -0.68
        Na               1.0
        K                0.2
        N(5)             1.2
END
EXCHANGE 1-40
        -equilibrate 1
        X                0.0011
END
TRANSPORT
        -cells           40
        -lengths         0.002
        -shifts          100
        -time_step       720.0
        -flow_direction  forward
        -boundary_conditions   flux  flux
        -diffusion_coefficient 0.0e-9
        -dispersivities  0.002
        -correct_disp    true
        -punch_cells     40
        -punch_frequency 1
        -print_cells     40
        -print_frequency 20
SELECTED_OUTPUT
        -reset false
        -step
        -totals Na Cl K Ca
END
'''
ip.run_string(inp)
out=ip.get_selected_output_array()
hdr=out[0]; a=np.array([r for r in out[1:] if r[0]>0],dtype=float)
pv_p=(a[:,0]+0.5)/40.
names=['na','k','ca','cl']; si={'na':1,'k':2,'ca':3,'cl':6}; pi={'na':1,'cl':2,'k':3,'ca':4}
print('PV   ' + '  '.join(f'{n}_SCEP  {n}_PHRQ' for n in names))
for pv in [0.5,0.8,1.0,1.2,1.5,1.8,2.0,2.2,2.5]:
    vals=[]
    for n in names:
        vs=np.interp(pv,pv_s,s[:,si[n]])*1e3; vp=np.interp(pv,pv_p,a[:,pi[n]])*1e3
        vals.append(f'{vs:6.3f} {vp:6.3f}')
    print(f'{pv:4.1f}  '+'   '.join(vals))
# RMSE over 0.3-2.5 PV
pvg=np.linspace(0.3,min(2.5,pv_s.max(),pv_p.max()),100)
for n in names:
    d=np.interp(pvg,pv_s,s[:,si[n]])*1e3-np.interp(pvg,pv_p,a[:,pi[n]])*1e3
    print(n,'RMSE mM',round(float(np.sqrt((d**2).mean())),4),'max',round(float(np.abs(d).max()),4))
