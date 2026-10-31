# Compare a 1-yr irr061 check run on Stanage with the cloud v23 reference (paddy/v23/irr061_y1_v23.json).
# Usage (repo root): python3 stanage_handoff/compare_ref.py <run dir = $SCEPTER_OUT/paddy> [tag=g22_irr061]
# Pass = every CH4, N2O, CO2-degassing and cation-release value within 1% (calcite/PDA compared in absolute t CO2/ha:
# calcite in near-saturated runs differs by up to ~0.09 t between identical cloud runs).
import sys, json
R = sys.argv[1]; tag = sys.argv[2] if len(sys.argv) > 2 else 'g22_irr061'
sys.argv = ['x', '--days=500']
exec(open('paddy/ghg_budget.py').read().split('ctrl = season')[0])
ref = json.load(open('paddy/v23/irr061_y1_v23.json')); bad = 0
for d in (0, 20):
    for x in 'sf':
        got = season(f'{R}/flood_{tag}_{d}_y1{x}'); r = ref[f'f{d}_{x}']
        for k, v in r.items():
            g = got[k]
            if 'calcite' in k or 'ALK' in k: ok = abs(g-v) <= 0.1; dev = f'{g-v:+.3f} t'
            else: ok = abs(g-v) <= 0.01*max(abs(v), 1e-9) or abs(g-v) < 1e-3; dev = f'{100*(g/v-1) if v else 0:+.2f}%'
            bad += not ok
            print(f'{d:2d} t {"season" if x=="s" else "fallow"} {k:38s} ref {v:10.4f} stanage {g:10.4f} {dev:>10s} {"OK" if ok else "FAIL"}')
print('PASS' if bad == 0 else f'{bad} values outside tolerance')
