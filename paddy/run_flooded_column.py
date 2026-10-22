# Flooded soil column test for the paddy Fe(III) reduction reaction (g2fe3).
# Usage (from the repository root, after `make --file=makefile_ex11`):
#   python3 paddy/run_flooded_column.py <tag> <extrxns comma list | none> [ttot_yr] [--fe2buf]
# e.g. python3 paddy/run_flooded_column.py fe3 fe2o2,g2fe3 0.5
# --fe2buf switches on cation exchange (10 cmol/kg on inrt, Fe2+ included) and lets siderite precipitate.
# --ch4 implies --fe2buf and adds SO42-/HS- and CH4; list g2so4,g2ch4,ch4o2 in the reactions to use them.
#   e.g. python3 paddy/run_flooded_column.py ch4 fe2o2,g2fe3,g2so4,g2ch4,ch4o2 0.5 --ch4
#   add fe2gr for the green rust pH-stat, ch4pl,ch4eb for plant transport and ebullition of CH4,
#   and fe2fs,fesox for FeS precipitation and oxidation
# --awd drives alternate wetting and drying through the seasonal wetness forcing (Wet_temp.in): flooded for
#   the first month, then repeated ~1 week drains (surface saturation 0.5, water table at 0.5 m) and ~2 week
#   reflooding until mid-year, then drained fallow. --season uses the same forcing machinery but stays flooded
#   (the control for --awd).
# Output goes to $SCEPTER_OUT (default ./scepter_output/) under paddy/flood_<tag>.
import os,sys
sys.path.insert(0,os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import spinup

outdir = os.environ.get('SCEPTER_OUT','./scepter_output/')
flags  = [a for a in sys.argv[1:] if a.startswith('--')]
args   = [a for a in sys.argv[1:] if not a.startswith('--')]
tag    = args[0]
rxns   = args[1].split(',') if args[1] != 'none' else []
ttot   = float(args[2]) if len(args) > 2 else 0.5
ch4    = '--ch4' in flags
fe2buf = '--fe2buf' in flags or ch4
awd    = '--awd' in flags
# optional --fe3=<wt fraction of amfe3, default 0.03> and --om=<OM input g C m-2 yr-1, default 300>
opt    = dict(a[2:].split('=') for a in flags if '=' in a)
fe3wt  = float(opt.get('fe3',0.03))
omin   = float(opt.get('om',300))
season = awd or '--season' in flags

def awd_schedule(t_flood=0.08,t_drain=0.02,t_reflood=0.04,t_end=0.5,s_dry=0.5,ramp=0.003):
    # surface saturation through the year (linear interpolation between points, repeated every year)
    t,v = [0.0],[1.0]
    tc = t_flood
    while tc + t_drain < t_end:
        t += [tc, tc+ramp, tc+t_drain, tc+t_drain+ramp]
        v += [1.0, s_dry, s_dry, 1.0]
        tc += t_drain + t_reflood
    t += [t_end, t_end+ramp, 0.99]
    v += [1.0, s_dry, s_dry]
    return [t, v if awd else [1.0]*len(v)]

spinup.run_a_scepter_run(
    f'paddy/flood_{tag}',outdir,
    # ---- frame.in ----
    ztot=0.3,nz=30,ttot=ttot,temp=25,omrain=omin,zom=0.15,poro=0.5,
    moistsrf=1.0,           # water saturation 1 throughout: flooded
    zwater=0.5 if awd else 1000,w=0,q=0.5,p=1e-5,nstep=10,rstrt='self',
    # ---- switches.in ----
    act_ON='false',cec_on='true' if fe2buf else 'false',display=1,report=1,season='true' if season else 'false',
    psd_bulk='false',psd_full='false',
    # ---- tracers ----
    sld_list=['inrt','g2','amfe3'] + (['sid','gr'] if fe2buf else []) + (['fes'] if ch4 else []),
    aq_list=(['fe2','fe3','ca','mg','na','k'] if fe2buf else ['fe2','fe3','ca'])
        + (['so4','hs'] if ch4 else []),
    gas_list=['pco2','po2'] + (['pch4'] if ch4 else []),
    exrxn_list=rxns,
    # ---- boundary values ----
    pr_list=[('inrt',1-fe3wt),('amfe3',fe3wt)],
    # rain inputs are free (primary) species with -Drain_input_as_primary:
    # free Fe3+ must be near Fe(OH)3 equilibrium, otherwise hydrolysis makes total Fe(III) huge
    rain_list=[('ca',2.05e-4 if ch4 else 5e-6),('fe3',1e-21),('fe2',1e-10)]
        + ([('mg',2e-6),('na',1e-5),('k',2e-6)] if fe2buf else [])
        + ([('so4',2e-4),('hs',1e-12)] if ch4 else []),  # irrigation water ~0.2 mM SO42-, balanced by Ca above
    atm_list=[('pco2',3.16e-4),('po2',0.21),('pnh3',1e-50),('pn2o',1e-50),('pch4',1.9e-6)],
    # cec [cmol/kg], log KH\\X for Na, K, Ca, Mg, Al, beta; Fe2+ follows Mg (+0.158) in the code
    sld_varlist_cec=[('inrt',10,5.9,4.8,10.47,10.786,16.47,3.4)] if fe2buf else [],sld_varlist_omrain=[('g2',1.0)],
    srcfile_2ndslds='./data/2ndslds_paddy.in' if fe2buf else './data/2ndslds_def.in',
    use_local_storage=False,exename_src='scepter_ex11',show_runtime_res=False,
    # ---- seasonal forcing (only read when season is on) ----
    moist_temp=awd_schedule(),T_temp=[[0.0,0.99],[25,25]],q_temp=[[0.0,0.99],[0.5,0.5]],
    dust_temp=[[0.0,0.99],[0,0]],
    )
