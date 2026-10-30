# Fig 6b data (v22): Si root-O2 feedback off (w 0), reference (w 0.2) and stronger (w 0.4) at the reference site,
# 10 yr at 0/10/20 t/ha/yr; each treatment against the control run with the same w. Usage (from /home/user/pd):
# python3 fig6b_si_v22.py <outdir> <json>
import sys,os,json; O=sys.argv[1]; out=sys.argv[2]; sys.argv=['x','--days=500','--alk=pda']
exec(open('paddy/ghg_budget.py').read().split('ctrl = season')[0])
AK='potentially durable ALK (t CO2/ha)'; K=lambda n: n*44/28*0.273
SETS={'0':('r22w0_ref_0','r22w0_ref_{d}'),'0.2':('r10_ref_0','r22_ref_{d}'),'0.4':('r22w4_ref_0','r22w4_ref_{d}')}
res={}
for w,(c,t) in SETS.items():
    yc=[season(f'{O}/flood_{c}_y{y}s')['CH4 (kg/ha)'] for y in range(1,11)]
    res[w]={'ctrl_ch4':yc}
    for d in (10,20):
        r=t.format(d=d); s=dict(ch4=0,n2o=0,pda=0,cc=0,cat=0); pct=[]
        for y in range(1,11):
            a={x:season(f'{O}/flood_{c}_y{y}{x}') for x in 'sf'}; b={x:season(f'{O}/flood_{r}_y{y}{x}') for x in 'sf'}
            dd=lambda k: sum(b[x][k]-a[x][k] for x in 'sf')
            s['ch4']+=dd('CH4 (kg/ha)')*27e-3; s['n2o']+=K(dd('N2O (kg N/ha)')); s['pda']+=dd(AK); s['cc']+=dd('calcite (t CO2/ha)')
            s['cat']+=sum(b[x]['cation release (t CO2/ha)'] for x in 'sf'); pct.append(100*(b['s']['CH4 (kg/ha)']/a['s']['CH4 (kg/ha)']-1))
        s['cc']=min(s['cc'],max(s['cat']-s['pda'],0)); s['net']=s['ch4']+s['n2o']-s['pda']-s['cc']; s['pct']=pct
        res[w][d]=s; print(w,d,{k:(round(v,1) if not isinstance(v,list) else [round(x) for x in v]) for k,v in s.items()})
json.dump(res,open(out,'w'),indent=1)
