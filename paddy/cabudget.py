# Whole-column Ca budget from native flx files (time-integrated, mol m-2): aq Ca (incl. exchangeable), wollastonite, calcite.
# Closure: sum of d(inventory) = sum of rain inputs + transport (adv+dif) ; residual reported. Usage: cabudget.py <ctrl> <trt> ...
import sys, numpy as np, pandas as pd
def integ(d, f):
    import os
    if not os.path.exists(f'{d}/flx/int_flx_{f}.txt'): return {k: 0.0 for k in ('tflx','rain','res','adv','dif')}
    x = pd.read_csv(f'{d}/flx/int_flx_{f}.txt', sep=r'\s+').iloc[-1]   # model's own time-averaged fluxes
    return {c: float(x[c]*x['time']) for c in x.index if c != 'time'}
def budget(d):
    aq, ex, w, c = integ(d, 'aq-ca'), integ(d, 'aqex-ca'), integ(d, 'sld-wls'), integ(d, 'sld-cc')
    B = dict(d_aq=aq['tflx'], d_ex=ex['tflx'], d_wls=w['tflx'], d_cc=c['tflx'],
             rain_aq=aq['rain'], rain_wls=w['rain'], rain_cc=c['rain'], transport=aq['adv']+aq['dif'],
             wls_to_aq=aq.get('wls',0.0), cc_to_aq=aq['cc'], res_aq=aq['res'], res_wls=w['res'], res_cc=c['res'])
    other = {k: v for k, v in aq.items() if k not in ('tflx','adv','dif','rain','wls','cc','res') and abs(v) > 1e-6}
    B['other_aq'] = sum(other.values()); B['other_terms'] = other
    return B
a = sys.argv[1:]
for cd, td in zip(a[::2], a[1::2]):
    bc, bt = budget(cd), budget(td)
    print('==', td.split('flood_')[-1], 'minus', cd.split('flood_')[-1])
    for k in bc:
        if k == 'other_terms': print(f'  other aq terms trt: {bt[k]}'); continue
        print(f'  {k:10s} ctrl {bc[k]:+9.3f}  trt {bt[k]:+9.3f}  diff {bt[k]-bc[k]:+9.3f}')
    for lab, B in (('ctrl', bc), ('trt', bt)):
        dinv = B['d_aq'] + B['d_wls'] + B['d_cc']; src = B['rain_aq'] + B['rain_wls'] + B['rain_cc'] + B['transport'] + B['other_aq']
        print(f'  {lab}: d(inventory) {dinv:+.3f} vs net input {-src:+.3f}  (closure error {dinv+src:+.2e})')
