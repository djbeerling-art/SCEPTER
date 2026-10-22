# Flooded soil column test for the paddy Fe(III) reduction reaction (g2fe3).
# Usage (from the repository root, after `make --file=makefile_ex11`):
#   python3 paddy/run_flooded_column.py <tag> <extrxns comma list | none> [ttot_yr] [--fe2buf]
# e.g. python3 paddy/run_flooded_column.py fe3 fe2o2,g2fe3 0.5
# --fe2buf switches on cation exchange (10 cmol/kg on inrt, Fe2+ included) and lets siderite precipitate.
# --ch4 implies --fe2buf and adds SO42-/HS- and CH4; list g2so4,g2ch4,ch4o2 in the reactions to use them.
#   e.g. python3 paddy/run_flooded_column.py ch4 fe2o2,g2fe3,g2so4,g2ch4,ch4o2 0.5 --ch4
#   add fe2gr for the green rust pH-stat and ch4pl,ch4eb for plant transport and ebullition of CH4
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

spinup.run_a_scepter_run(
    f'paddy/flood_{tag}',outdir,
    # ---- frame.in ----
    ztot=0.3,nz=30,ttot=ttot,temp=25,omrain=300,zom=0.15,poro=0.5,
    moistsrf=1.0,           # water saturation 1 throughout: flooded
    zwater=1000,w=0,q=0.5,p=1e-5,nstep=10,rstrt='self',
    # ---- switches.in ----
    act_ON='false',cec_on='true' if fe2buf else 'false',display=1,report=1,season='false',
    psd_bulk='false',psd_full='false',
    # ---- tracers ----
    sld_list=['inrt','g2','amfe3'] + (['sid','gr'] if fe2buf else []),
    aq_list=(['fe2','fe3','ca','mg','na','k'] if fe2buf else ['fe2','fe3','ca'])
        + (['so4','hs'] if ch4 else []),
    gas_list=['pco2','po2'] + (['pch4'] if ch4 else []),
    exrxn_list=rxns,
    # ---- boundary values ----
    pr_list=[('inrt',0.97),('amfe3',0.03)],
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
    )
