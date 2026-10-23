# Flooded soil column test for the paddy Fe(III) reduction reaction (g2fe3).
# Usage (from the repository root, after `make --file=makefile_paddy`):
#   python3 paddy/run_flooded_column.py <tag> <extrxns comma list | none> [ttot_yr] [--fe2buf]
# e.g. python3 paddy/run_flooded_column.py fe3 fe2o2,g2fe3 0.5
# --fe2buf switches on cation exchange (10 cmol/kg on inrt, Fe2+ included) and lets siderite precipitate.
# --ch4 implies --fe2buf and adds SO42-/HS- and CH4; list g2so4,g2ch4,ch4o2 in the reactions to use them.
#   e.g. python3 paddy/run_flooded_column.py ch4 fe2o2,g2fe3,g2so4,g2ch4,ch4o2 0.5 --ch4
#   add fe2gr for the green rust pH-stat, ch4pl,ch4eb for plant transport and ebullition of CH4,
#   and fe2fs,fesox for FeS precipitation and oxidation
# --awd drives alternate wetting and drying through the seasonal wetness forcing (Wet_temp.in): flooded for
#   the first month, then repeated ~1 week drains (see --zdry, --sdry below) and ~2 week
#   reflooding until mid-year, then drained fallow. --season uses the same forcing machinery but stays flooded
#   (the control for --awd).
# --n adds nitrogen: NO3- (aq), NH3/NH4+ (pnh3) and N2O (pn2o) gases, atmospheric N2O 330 ppb; list
#   amo2o (nitrification), g2n21 (NO3- -> N2O) and g2n2p (paddy N2O -> N2) in the reactions to use them.
# --k=<file> copies a rate-constant override file (header line, then lines of `name krxn1 krxn2`, -999 keeps
#   the default) into the run directory as krxn_ext.in, so calibration runs need no recompiling.
# --kfw=<m/yr> sets the floodwater gas transfer velocity between the flooded surface layer and the atmosphere
#   (written to floodwater.in; default 36.5 m/yr = 0.1 m/d, limited by the diffusive boundary layer at the soil
#   surface); --kfw=0 restores SCEPTER's closed top boundary for dissolved gases.
# --fallow runs the column drained (surface saturation 0.5, water table at 0.5 m), e.g. to spin up the soil OM
#   and Fe(III) stocks before a rice season; --rstrt=<tag> starts from the end state of run flood_<tag>.
# --nfert=<kg N/ha> (with --n) adds urea, spread evenly over the first 0.31 yr (113 d, IPCC default season)
#   through the OM input (same depth distribution); urea hydrolyses to NH3/NH4+ in a few days.
# --zdry=<m> and --sdry=<sat> set the AWD drain depth (water table, default 0.15 m: 'safe AWD', ponding
#   allowed to fall to 15 cm below the surface) and the surface saturation reached in a drain (default 0.6).
# --basalt=<t/ha> applies crushed basalt (mineral mix of data/dust_basalt.in without hornblende, renormalised;
#   particle radius 10 um) over the first 0.01 yr, mixed into the top 0.15 m, and tracks Si, Al and the secondary
#   minerals calcite, kaolinite, gibbsite and amorphous silica. --basalt=0 tracks the same species without
#   applying basalt (for spin-ups and controls, so restarts carry the same species).
#   --basfe3=<wt fraction> adds Fe(III) oxide (as amfe3) to the basalt, for its Fe2O3 and oxidised grain surfaces;
#   --basmgn=<wt fraction> adds magnetite. --rock=<name> takes the mineralogy from data/dust_<name>.in instead
#   (e.g. Hillhouse: labradorite, augite, olivine, analcime; Oregon, Tichum: other basalts; Tawau: Malaysian andesite).
#   Magnetite is not allowed to precipitate (data/2ndslds_paddy.in): it would take over from green rust and lock up
#   reducible soil Fe(III).
# --refsoil (with --basalt=0 or a dose, so Si is tracked): reference paddy soil with CEC 15 cmol/kg, 1 wt% phytolith
#   (amorphous) silica and river-type irrigation water (with Cl-), so porewater Ca, Mg and Si reach typical paddy levels.
# --kaol=<wt fraction> puts kaolinite in the soil (spin-ups only; restarts inherit it). Without it secondary Al phases
#   must nucleate from a 1e-20 seed, so porewater Al builds up and stalls plagioclase dissolution for weeks.
# Output goes to $SCEPTER_OUT (default ./scepter_output/) under paddy/flood_<tag>.
import os,sys,shutil
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
nitro  = '--n' in flags
fallow = '--fallow' in flags
nfert  = float(opt.get('nfert',0))
zdry   = float(opt.get('zdry',0.15))
bas    = 'basalt' in opt
basfe3 = float(opt.get('basfe3',0))     # Fe(III) oxide in the basalt (as amfe3), wt fraction
bas_t  = float(opt.get('basalt',0))
rock   = opt.get('rock','')
basmgn = float(opt.get('basmgn',0))     # magnetite in the basalt, wt fraction
if rock:   # a characterised rock from data/dust_<rock>.in
    bas_min = {}
    for l in open(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),'data',f'dust_{rock}.in')).readlines()[1:]:
        w = l.split()
        if len(w) >= 2 and float(w[1]) > 0: bas_min[w[0]] = float(w[1])
else:      # default: data/dust_basalt.in without hornblende
    bas_min = {'fo':0.12,'ab':0.172,'an':0.258,'fa':0.05,'dp':0.189,'kfs':0.06}
bas_min = {k:v/sum(bas_min.values())*(1-basfe3-basmgn) for k,v in bas_min.items()}
if '--light' in flags: bas_min = {}   # controls/calibration: Si, Al and secondary minerals only (faster; not restart-compatible with dosed runs)
if basfe3 > 0: bas_min['amfe3'] = basfe3
if basmgn > 0: bas_min['mgn'] = basmgn
t_app  = 0.01                        # application period, yr
fdust  = bas_t*100/t_app             # t/ha -> g m-2, spread over t_app
sdry   = float(opt.get('sdry',0.6))
kaol   = float(opt.get('kaol',0))       # soil kaolinite, wt fraction
ref    = '--refsoil' in flags         # reference paddy soil: CEC 15, phytolith Si, river-type irrigation water
cec    = 15 if ref else 10              # cmol/kg
amsiwt = 0.01 if ref else 0            # amorphous (phytolith) silica, wt fraction
# urea [g m-2 yr-1] as a fraction of the OM input [g C m-2 yr-1]; applied at a constant rate (the run is one season)
urea_frac = nfert/10*60.056/28.014/0.31/omin

def awd_schedule(t_flood=0.08,t_drain=0.02,t_reflood=0.04,t_end=0.5,s_dry=sdry,ramp=0.003):
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

rundir = os.path.join(outdir,f'paddy/flood_{tag}')
os.makedirs(rundir,exist_ok=True)
if os.path.exists(os.path.join(rundir,'krxn_ext.in')): os.remove(os.path.join(rundir,'krxn_ext.in'))
if 'k' in opt: shutil.copy(opt['k'],os.path.join(rundir,'krxn_ext.in'))
with open(os.path.join(rundir,'floodwater.in'),'w') as f:
    f.write('transfer velocity [m/yr]\n%s\n' % opt.get('kfw','36.5'))

spinup.run_a_scepter_run(
    f'paddy/flood_{tag}',outdir,
    # ---- frame.in ----
    ztot=0.3,nz=30,ttot=ttot,temp=25,omrain=omin,zom=0.15,poro=0.5,
    moistsrf=0.5 if fallow else 1.0,  # water saturation 1 throughout: flooded
    zwater=zdry if awd else (0.5 if fallow else 1000),w=0,q=0.5,p=1e-5,nstep=10,rstrt=os.path.join(outdir,'paddy','flood_'+opt['rstrt']) if 'rstrt' in opt else 'self',
    restart='true' if 'rstrt' in opt else 'false',
    # ---- switches.in ----
    act_ON='false',cec_on='true' if fe2buf else 'false',display=1,report=1,season='true' if season else 'false',
    psd_bulk='false',psd_full='false',
    # ---- tracers ----
    sld_list=['inrt','g2','amfe3'] + (['sid','gr'] if fe2buf else []) + (['fes'] if ch4 else []) + (['urea'] if nfert > 0 else [])
        + ([m for m in bas_min if m != 'amfe3'] + [m for m in ['cc','ka','gb','amsi'] if m not in bas_min] if bas else []),
    aq_list=(['fe2','fe3','ca','mg','na','k'] if fe2buf else ['fe2','fe3','ca'])
        + (['so4','hs'] if ch4 else []) + (['no3'] if nitro else []) + (['si','al'] if bas else []) + (['cl'] if ref else []),
    gas_list=['pco2','po2'] + (['pch4'] if ch4 else []) + (['pnh3','pn2o'] if nitro else []),
    exrxn_list=rxns,
    # ---- boundary values ----
    pr_list=[('inrt',1-fe3wt-kaol-amsiwt),('amfe3',fe3wt)] + ([('ka',kaol)] if kaol > 0 else []) + ([('amsi',amsiwt)] if amsiwt > 0 else []),
    # rain inputs are free (primary) species with -Drain_input_as_primary:
    # free Fe3+ must be near Fe(OH)3 equilibrium, otherwise hydrolysis makes total Fe(III) huge
    # --refsoil: river-type irrigation water (Ca 0.5, Mg 0.2, Na 0.3, K 0.05, Si 0.2, SO4 0.2, Cl 1.15 mM), pH ~7.6 at
    # atmospheric pCO2; free Fe3+ and Al3+ lowered to near Fe(OH)3 and gibbsite equilibrium at that pH
    rain_list=[('ca',(5e-4 if ref else 2.05e-4) if ch4 else 5e-6),('fe3',1e-27 if ref else 1e-21),('fe2',1e-10)]
        + ([('si',2e-4 if ref else 1e-4),('al',1e-15 if ref else 1e-12)] if bas else [])
        + ([('mg',2e-4),('na',3e-4),('k',5e-5),('cl',1.15e-3)] if ref else [('mg',2e-6),('na',1e-5),('k',2e-6)] if fe2buf else [])
        + ([('so4',2e-4),('hs',1e-12)] if ch4 else []) + ([('no3',1e-6)] if nitro else []),  # irrigation water ~0.2 mM SO42-, balanced by Ca above
    atm_list=[('pco2',3.16e-4),('po2',0.21),('pnh3',1e-50),('pn2o',3.3e-7 if nitro else 1e-50),('pch4',1.9e-6)],
    # cec [cmol/kg], log KH\\X for Na, K, Ca, Mg, Al, beta; Fe2+ follows Mg (+0.158) in the code
    sld_varlist_cec=[('inrt',cec,5.9,4.8,10.47,10.786,16.47,3.4)] if fe2buf else [],sld_varlist_omrain=[('g2',1.0)] + ([('urea',urea_frac)] if nfert > 0 else []),
    srcfile_2ndslds='./data/2ndslds_paddy.in' if fe2buf else './data/2ndslds_def.in',
    use_local_storage=False,exename_src='scepter_paddy',show_runtime_res=False,
    # ---- seasonal forcing (only read when season is on) ----
    moist_temp=awd_schedule(),T_temp=[[0.0,0.99],[25,25]],q_temp=[[0.0,0.99],[0.5,0.5]],
    dust_temp=[[0.0,t_app,t_app+1e-4,0.99],[fdust,fdust,0,0]],
    # non-seasonal runs (fallow spin-ups) get no basalt
    sld_varlist_dust=list(bas_min.items()) if bas else [],zdust=0.15,
    )
