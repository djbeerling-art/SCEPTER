# Flooded soil column test for the paddy Fe(III) reduction reaction (g2fe3).
# Usage (from the repository root, after `make --file=makefile_ex11`):
#   python3 paddy/run_flooded_column.py <tag> <extrxns comma list | none> [ttot_yr]
# e.g. python3 paddy/run_flooded_column.py fe3 fe2o2,g2fe3 0.5
# Output goes to $SCEPTER_OUT (default ./scepter_output/) under paddy/flood_<tag>.
import os,sys
sys.path.insert(0,os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import spinup

outdir = os.environ.get('SCEPTER_OUT','./scepter_output/')
tag    = sys.argv[1]
rxns   = sys.argv[2].split(',') if sys.argv[2] != 'none' else []
ttot   = float(sys.argv[3]) if len(sys.argv) > 3 else 0.5

spinup.run_a_scepter_run(
    f'paddy/flood_{tag}',outdir,
    # ---- frame.in ----
    ztot=0.3,nz=30,ttot=ttot,temp=25,omrain=300,zom=0.15,poro=0.5,
    moistsrf=1.0,           # water saturation 1 throughout: flooded
    zwater=1000,w=0,q=0.5,p=1e-5,nstep=10,rstrt='self',
    # ---- switches.in ----
    act_ON='false',cec_on='false',display=1,report=1,season='false',
    psd_bulk='false',psd_full='false',
    # ---- tracers ----
    sld_list=['inrt','g2','amfe3'],
    aq_list=['fe2','fe3','ca'],
    gas_list=['pco2','po2'],
    exrxn_list=rxns,
    # ---- boundary values ----
    pr_list=[('inrt',0.97),('amfe3',0.03)],
    # rain inputs are free (primary) species with -Drain_input_as_primary:
    # free Fe3+ must be near Fe(OH)3 equilibrium, otherwise hydrolysis makes total Fe(III) huge
    rain_list=[('ca',5e-6),('fe3',1e-21),('fe2',1e-10)],
    atm_list=[('pco2',3.16e-4),('po2',0.21),('pnh3',1e-50),('pn2o',1e-50)],
    sld_varlist_cec=[],sld_varlist_omrain=[('g2',1.0)],
    srcfile_2ndslds='./data/2ndslds_def.in',
    use_local_storage=False,exename_src='scepter_ex11',show_runtime_res=False,
    )
