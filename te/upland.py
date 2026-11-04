"""Upland (aerated cropland) SCEPTER driver for the trace-element (TE) study.

Spin-up follows the upstream GMD24 inert-soil approach (inert matrix + OM + CEC,
Ca rain and N fertiliser setting soil pH), then basalt is applied as dust.
"""
import os, sys, shutil, argparse
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
import make_inputs as mi

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
EXE = os.path.join(ROOT, 'scepter_gmd24')
OUT = os.environ.get('TE_OUT', '/home/user/te_runs/')


def write_common(run, *, ztot=0.5, nz=30, ttot, temp, fdust=0., taudust=0.,
                 omrain, poro, moist, q, w, zdust=0.2, restart='false',
                 rstrt='self', mix=1, cec, alpha=3.4, ca_rain, n_rain,
                 sld=None, dust=None, psd=True, nstep=10, psdrain=None,
                 uptake=False, parent=None):
    """parent: optional [(mineral, wt fraction)] in the soil parent material besides
    the inert matrix (e.g. kaolinite and gibbsite for an acid-buffered tropical soil)."""
    pr = [('inrt', 1.0 - sum(f for _, f in parent))] + list(parent) if parent else [('inrt', 1.0)]
    if parent:
        sld = list(dict.fromkeys([m for m, _ in parent] + (sld or [])))
    mi.get_input_frame(outdir=OUT, runname=run, ztot=ztot, nz=nz, ttot=ttot,
                       temp=temp, fdust=fdust, fdust2=0, taudust=taudust,
                       omrain=omrain, zom=ztot, poro=poro, moistsrf=moist,
                       zwater=100000, zdust=zdust, w=w, q=q, p=10e-6,
                       nstep=nstep, rstrt=rstrt, runid=run)
    mi.get_input_switches(outdir=OUT, runname=run, w_scheme=1, mix_scheme=mix,
                          poro_iter='false', sldmin_lim='true', display=1,
                          report=0, restart=restart, rough=os.environ.get('TE_ROUGH', 'true'),
                          act_ON='false', dt_fix='false', cec_on='true',
                          dz_fix='true', close_aq='false', poro_evol=os.environ.get('TE_POROEVOL', 'true'),
                          sa_evol_1='true', sa_evol_2='false',
                          psd_bulk='true' if psd else 'false',
                          psd_full='true' if psd else 'false', season='false')
    sld_list = ['inrt', 'g2', 'amnt'] + (sld or [])
    aq_list = ['ca', 'k', 'mg', 'na', 'no3']
    if sld:
        aq_list += ['si', 'al', 'fe2', 'fe3']
    mi.get_input_tracers(outdir=OUT, runname=run, sld_list=sld_list,
                         aq_list=aq_list, gas_list=['pco2', 'po2'] if (sld or uptake) else ['pco2'],
                         exrxn_list=(['fe2o2'] if sld else [])
                         + (['g2k', 'g2ca', 'g2mg'] if uptake else []))
    rain = [('ca', ca_rain)] + [(sp, 1e-20) for sp in aq_list
                                     if sp not in ('ca', 'fe3')]
    if sld:
        rain += [('fe3', 1e-21)]
    mi.get_input_tracer_bounds(outdir=OUT, runname=run, pr_list=pr,
                               rain_list=rain)
    mi.get_input_sld_properties(outdir=OUT, runname=run, filename='dust.in',
                                sld_varlist=dust or [])
    mi.get_input_sld_properties(outdir=OUT, runname=run, filename='cec.in',
                                sld_varlist=[(sp, cec, 5.9, 4.8, 10.47, 10.786, 16.47, alpha)
                                             for sp in ('inrt', 'g2')])
    # N fertiliser as NH4NO3 normalised to OC rain (as GMD24)
    nr = n_rain / 14. * 80. / 2. / omrain
    mi.get_input_sld_properties(outdir=OUT, runname=run, filename='OM_rain.in',
                                sld_varlist=[('g2', 1), ('amnt', nr)])
    mi.get_input_sld_properties(outdir=OUT, runname=run, filename='2ndslds.in',
                                srcfile=os.path.join(ROOT, 'data/2ndslds_def.in'))
    if psdrain:   # (median radius m, sd of log10 radius, weight) tuples
        with open(os.path.join(OUT, run, 'psdrain.in'), 'w') as f:
            f.write('** mean radius [m], standard deviation in log10 [-], weight [-]\n')
            for r, sd, w_ in psdrain:
                f.write(f'{r:.4e}\t{sd:.4f}\t{w_}\n')
    else:
        mi.get_input_sld_properties(outdir=OUT, runname=run, filename='psdrain.in',
                                    srcfile=os.path.join(ROOT, 'data/psdrain_100um.in'))
    shutil.copy(EXE, os.path.join(OUT, run, 'scepter'))


def run(run):
    d = os.path.join(OUT, run)
    return os.system(f'cd {d} && ./scepter > log.txt 2>&1')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('name')
    ap.add_argument('--ttot', type=float, default=10000)
    ap.add_argument('--temp', type=float, default=10)
    ap.add_argument('--q', type=float, default=0.25)
    ap.add_argument('--w', type=float, default=1e-3)
    ap.add_argument('--moist', type=float, default=0.5)
    ap.add_argument('--cec', type=float, default=12)
    ap.add_argument('--om', type=float, default=300)
    ap.add_argument('--ca', type=float, default=5e-6)
    ap.add_argument('--n', type=float, default=15)  # g N/m2/yr
    ap.add_argument('--uptake', action='store_true')  # crop K, Ca, Mg uptake (g2k/g2ca/g2mg)
    ap.add_argument('--parent', default='')  # e.g. 'ka:0.2,gb:0.02'
    a = ap.parse_args()
    parent = [(x.split(':')[0], float(x.split(':')[1])) for x in a.parent.split(',') if x]
    write_common(a.name, ttot=a.ttot, temp=a.temp, omrain=a.om, poro=0.5,
                 moist=a.moist, q=a.q, w=a.w, cec=a.cec, ca_rain=a.ca,
                 n_rain=a.n, uptake=a.uptake, parent=parent or None)
    print('rc', run(a.name))


# Feedstock grain size as SCEPTER dust PSD (lognormal in radius).
# Pioneer Valley basalt (Zacharias et al. 2026 SI Table S2c sieve data, 75 um to
# 4.75 mm): lognormal fit D50 = 106 um, sigma_ln = 1.16 (sd log10 = 0.503),
# P80 = 281 um (283 um by interpolation); max misfit 0.6% passing.
# Blue Ridge metabasalt (Corn Belt; Beerling et al. 2024): P80 = 267 um. Spread
# not reported: sigma_ln 1.16 assumed (as Pioneer Valley) -> D50 = 101 um.
# Tawau basaltic andesite: <5 mm jaw/cone-crushed screenings (Lewis et al. 2021,
# Table 1); cumulative mass PSD read from the vector data of Lewis 2021 SI Fig S3
# (D50 ~0.9 mm, P80 ~1.7 mm), fitted with two lognormal modes (max misfit 2% passing).
PSD = {'pioneer': [(53e-6, 0.503, 1)],
       'blueridge': [(50e-6, 0.503, 1)],
       'tawau': [(170e-6, 0.527, 0.518), (662e-6, 0.191, 0.482)]}

ROCKS = {  # SCEPTER mineral recipes (wt fraction)
    'hillhouse': [('la', 0.346), ('agt', 0.327), ('fo', 0.147), ('anl', 0.071)],
    # Pioneer Valley metabasalt (Vermont; Zacharias et al. 2026 SI Table S2b,
    # thin-section point count): clinopyroxene 35.1 (as augite), plagioclase 33.7
    # (as labradorite), sericite 10.6 (as muscovite), chlorite 9.2, actinolite 6.7
    # (as tremolite), Fe-Ti oxides 3.8 (as magnetite), quartz 0.9 %.
    'pioneer': [('agt', 0.351), ('la', 0.337), ('mscv', 0.106), ('clch', 0.092),
                ('tm', 0.067), ('mgn', 0.038), ('qtz', 0.009)],
    # Blue Ridge metabasalt (Corn Belt), Dataset S1 XRD (Southampton): chlorite 34,
    # albite 23.3, epidote (piemontite) 17.8, actinolite 11.9 (as tremolite),
    # quartz 9, calcite 2.6, muscovite 1.4 wt%.
    'blueridge': [('clch', 0.34), ('ab', 0.233), ('ep', 0.178), ('tm', 0.119),
                  ('qtz', 0.09), ('cc', 0.026), ('mscv', 0.014)],
    # Tawau basaltic andesite (Sabah; Lewis et al. 2021 SI Table S2, SEM-EDS + XRD):
    # plagioclase An64 44.6 (as labradorite), quartz 19.1, sanidine 13.3 (as K-feldspar),
    # diopside 12.4, smectite 5.8, magnetite 2.0, epistilbite 1.0, calcite 0.8, biotite
    # 0.6, apatite 0.4 wt%. Smectite, epistilbite, biotite, apatite left out
    # (unreactive or absent from SCEPTER); the rest renormalised to 1.
    'tawau': [('la', 0.484), ('qtz', 0.207), ('kfs', 0.144), ('dp', 0.134),
              ('mgn', 0.022), ('cc', 0.009)],
}
SECONDARY = ['ka', 'gb', 'amsi', 'gt', 'cc']
# TE_SECONDARY="+amfe3,cabd,mgbd,amal" adds precipitating phases (ferrihydrite-like
# amorphous Fe(III), Ca/Mg beidellite smectite, amorphous Al hydroxide)
if os.environ.get('TE_SECONDARY', '').startswith('+'):
    SECONDARY = SECONDARY + os.environ['TE_SECONDARY'][1:].split(',')


def amend(spin, name, *, rock='hillhouse', dose=40., years=40., temp, q, moist,
          cec, ca_rain, omrain=300, w=0., zdust=0.2, nstep=10, pulse=None, psdrain=None,
          n_rain=0., uptake=False, parent=None):
    """Restart from spin-up and apply rock continuously at dose t/ha/yr, or,
    with pulse = tau (yr, < 1), once: dose t/ha spread over the first tau years.
    n_rain (g N/m2/yr, as NH4NO3) should match the spin-up so the fertiliser
    acidity carries on through the amendment run."""
    d = os.path.join(OUT, name)
    if os.path.exists(d):
        shutil.rmtree(d)
    shutil.copytree(os.path.join(OUT, spin), d)
    for sub in ('flx',):
        shutil.rmtree(os.path.join(d, sub), ignore_errors=True)
    minerals = [m for m, _ in ROCKS[rock]] if dose > 0 else []
    write_common(name, ttot=years, temp=temp, fdust=dose * 100.,
                 taudust=pulse if pulse else years,
                 omrain=omrain, poro=0.5, moist=moist, q=q, w=w, zdust=zdust,
                 restart='true', rstrt=os.path.join(OUT, spin), mix=3, cec=cec, ca_rain=ca_rain,
                 n_rain=n_rain, sld=minerals + SECONDARY, dust=ROCKS[rock] if dose > 0 else [],
                 nstep=nstep, psdrain=psdrain, uptake=uptake, parent=parent,
                 psd=os.environ.get('TE_PSD', 'true') == 'true')
    return run(name)
