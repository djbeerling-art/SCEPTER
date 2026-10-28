# Whole-rock oxides (wt%) -> SCEPTER rock recipe, by a simplified CIPW norm.
# Order: apatite (P2O5, removes CaO only), ilmenite (TiO2 + FeO), Fe2O3 (all as reducible Fe(III), amfe3, as used
# for JB-1b against Park et al. 2026; --fe3mgn puts it in magnetite with FeO instead), orthoclase (K2O), albite (Na2O),
# anorthite (remaining CaO with Al2O3), diopside (remaining CaO with MgO), hypersthene (en/fer, remaining MgO/FeO),
# then olivine (fo/fa) replaces hypersthene if silica runs short, quartz takes any excess. Excess Al2O3 is inert.
# Usage: python3 basalt_norm.py <name> SiO2=.. TiO2=.. Al2O3=.. Fe2O3=.. FeO=.. MnO=.. MgO=.. CaO=.. Na2O=.. K2O=.. P2O5=..
#        -> writes data/dust_<name>.in (wt fractions normalised to 1) and prints the recipe
# FeO only (FeOT as FeO) -> Fe2O3 = 0.15*FeOT*1.111 by convention (Fe3+/Fetot 0.15, typical of fresh basalt) unless given.
import sys,os

MW = {'SiO2':60.08,'TiO2':79.87,'Al2O3':101.96,'Fe2O3':159.69,'FeO':71.84,'MnO':70.94,'MgO':40.30,'CaO':56.08,
      'Na2O':61.98,'K2O':94.20,'P2O5':141.94}
FW = {'kfs':278.33,'ab':262.22,'an':278.21,'dp':216.55,'en':100.39,'fer':131.93,'fo':140.69,'fa':203.77,'qtz':60.08,
      'ilmn':151.71,'mgn':231.53,'amfe3':159.69,'inrt':101.96}   # amfe3 counted as Fe2O3 mass; inrt = spare Al2O3

def norm(ox, fe3mgn=False):
    ox = dict(ox)
    if 'Fe2O3' not in ox:
        feot = ox.pop('FeOT', ox.get('FeO', 0.))
        ox['Fe2O3'] = 0.15*feot*1.1113; ox['FeO'] = 0.85*feot
    m = {k: ox.get(k, 0.)/MW[k] for k in MW}
    m['FeO'] += m['MnO']; m['MnO'] = 0.
    out = dict.fromkeys(FW, 0.)
    m['CaO'] -= min(m['CaO'], 10/3*m['P2O5'])                     # apatite
    t = min(m['TiO2'], m['FeO']); out['ilmn'] = t; m['FeO'] -= t   # ilmenite
    if fe3mgn:
        t = min(m['Fe2O3'], m['FeO']); out['mgn'] = t; m['FeO'] -= t; m['Fe2O3'] -= t
    out['amfe3'] = m['Fe2O3']
    si = m['SiO2']; al = m['Al2O3']
    t = min(m['K2O'], al); out['kfs'] = 2*t; al -= t; si -= 6*t       # KAlSi3O8 (2 per K2O)
    t = min(m['Na2O'], al); out['ab'] = 2*t; al -= t; si -= 6*t       # NaAlSi3O8
    t = min(m['CaO'], al); out['an'] = t; al -= t; si -= 2*t; ca = m['CaO'] - t
    out['inrt'] = al
    mg, fe = m['MgO'], m['FeO']; xmg = mg/(mg+fe) if mg+fe > 0 else 1.
    t = min(ca, mg); out['dp'] = t; mg -= t; si -= 2*t                 # CaMgSi2O6 (Fe-free diopside)
    hy = mg + fe                                                       # (Mg,Fe)SiO3
    if si >= hy:
        out['en'] = mg; out['fer'] = fe; si -= hy; out['qtz'] = si
    else:                                                              # olivine (Mg,Fe)2SiO4 for the silica shortfall
        ol = max(0., min(hy/2, hy - si)); hy2 = hy - 2*ol
        x = mg/hy if hy > 0 else xmg
        out['fo'], out['fa'] = ol*x, ol*(1-x); out['en'], out['fer'] = hy2*x, hy2*(1-x)
    w = {k: v*FW[k] for k, v in out.items() if v > 1e-9}
    s = sum(w.values())
    return {k: v/s for k, v in w.items()}

def write(name, r, note='CIPW-type norm of whole-rock oxides (paddy/global/basalt_norm.py)'):
    root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    with open(os.path.join(root, 'data', f'dust_{name}.in'), 'w') as f:
        f.write(f'** {name}: {note}\n')
        for k, v in sorted(r.items(), key=lambda kv: -kv[1]): f.write(f'{k:7s} {v:.4f}\n')

if __name__ == '__main__':
    name = sys.argv[1]
    ox = {a.split('=')[0]: float(a.split('=')[1]) for a in sys.argv[2:] if '=' in a}
    r = norm(ox, fe3mgn='--fe3mgn' in sys.argv)
    write(name, r)
    print(' '.join(f'{k} {v:.3f}' for k, v in sorted(r.items(), key=lambda kv: -kv[1])))
