# Re-speciate SCEPTER porewaters (charge_balance-*.txt: every species, mol/L, ideal solution) in PHREEQC (phreeqc.dat)
# and compare pH (from charge balance) and mineral saturation (SCEPTER sat_sld-*.txt = Omega).
import sys,re,json,numpy as np,phreeqpython
BASE={'mg':'Mg','ca':'Ca','na':'Na','k':'K','al':'Al','fe2':'Fe(2)','fe3':'Fe(3)','si':'Si','so4':'S(6)','hs':'S(-2)',
      'no3':'N(5)','nh4':'N(-3)','nh3':'N(-3)','cl':'Cl','co3':'C(4)','hco3':'C(4)','co2':'C(4)'}
def parse(sp):
    # component list of a SCEPTER species name, e.g. 'mg(oh)1', 'fe2hco3', '(nh4)1so4', 'h1so4', 'hhs'
    out={}
    s=sp
    for m in re.finditer(r'\((\w+?)\)(\d)',s):
        out[m.group(1)]=out.get(m.group(1),0)+int(m.group(2))
    s=re.sub(r'\(\w+?\)\d','',s)
    s=re.sub(r'^h1','',s)
    if s=='hhs': s='hs'
    for b in ['fe2','fe3','hco3','co3','co2','nh4','nh3','so4','no3','mg','ca','na','al','si','cl','hs','k']:
        if s.startswith(b): out[b]=out.get(b,0)+1; s=s[len(b):]
        if s.startswith(b): out[b]=out.get(b,0)+1; s=s[len(b):]
    if s.startswith('hco3'): out['hco3']=out.get('hco3',0)+1
    elif s.startswith('co3'): out['co3']=out.get('co3',0)+1
    return out
def totals(h,row):
    t={}
    for k,v in zip(h,row):
        if k in ('z','h','oh','IS','tot_charge','o2','n2o','ch4'): continue
        for c,n in parse(k).items():
            if c in BASE: t[BASE[c]]=t.get(BASE[c],0)+n*v
    return t
def L(f):
    h=open(f).readline().split(); return h,np.loadtxt(f,skiprows=1)
TEMP=float(__import__('os').environ.get('TEMP_C','25'))
pp=phreeqpython.PhreeqPython(database='phreeqc.dat')
def run(prof,idx,tag):
    hc,c=L(f'{prof}/charge_balance-{idx}.txt'); hs,s=L(f'{prof}/sat_sld-{idx}.txt'); hg,g=L(f'{prof}/prof_gas-{idx}.txt')
    rows=[]
    for i in range(len(c)):
        cb=dict(zip(hc,c[i])); t=totals(hc,c[i])
        lines=['SOLUTION 1','  units mol/kgw',f'  temp {TEMP}','  pH 7 charge']+[f'  {k} {v:.6e}' for k,v in t.items() if v>1e-15]
        inp='\n'.join(lines)+'\nSELECTED_OUTPUT\n -reset false\n -pH\n -ionic_strength\n -si Calcite Siderite Fe(OH)3(a) SiO2(a)\n -activities Ca+2 CO3-2 Fe+2\n -molalities Ca+2 CO3-2 Fe+2\nEND\n'
        pp.ip.run_string(inp); o=pp.ip.get_selected_output_array(); r={k.split('(mol')[0]:v for k,v in zip(o[0],o[1])}
        # PHREEQC SI with activity coefficients set to 1 (SCEPTER's ideal-solution convention)
        lg=lambda sp: np.log10(r[f'la_{sp}'] if f'la_{sp}' in r else 1)
        g_ca=np.log10(10**r['la_Ca+2']/r['m_Ca+2']); g_co3=np.log10(10**r['la_CO3-2']/r['m_CO3-2']); g_fe=np.log10(10**r['la_Fe+2']/r['m_Fe+2'])
        sd=dict(zip(hs,s[i]))
        rows.append(dict(tag=tag,z=cb['z'],ph_s=-np.log10(cb['h']),ph_p=r['pH'],IS_s=cb['IS'],IS_p=r['mu'],
            cc_s=np.log10(max(sd['cc'],1e-30)),cc_p=r['si_Calcite'],cc_pi=r['si_Calcite']-g_ca-g_co3,
            sid_s=np.log10(max(sd['sid'],1e-30)),sid_p=r['si_Siderite'],sid_pi=r['si_Siderite']-g_fe-g_co3,
            fh_s=np.log10(max(sd['amfe3'],1e-30)),fh_p=r['si_Fe(OH)3(a)'],amsi_s=np.log10(max(sd['amsi'],1e-30)),amsi_p=r['si_SiO2(a)']))
    return rows
if __name__=='__main__':
    out=[]
    for spec in sys.argv[2:]:
        prof,idx,tag=spec.split(':'); out+=run(prof,idx,tag)
    json.dump(out,open(sys.argv[1],'w'),indent=0)
    for r in out[::3]: print({k:(round(v,3) if isinstance(v,float) else v) for k,v in r.items()})
