#!/bin/bash
# v22 verification data (Fig 3): PHREEQC re-speciation of irr061 porewaters and element budgets of 8 season/fallow runs.
until grep -q "zv done" /home/user/runs/pan/statusZv; do sleep 30; done
V=/mnt/project-files/riceland-scepter/v22/verification; O=/home/user/runs/out/paddy/flood_r22v_irr061
cd /home/user/pd
TEMP_C=28.7 python3 -I $V/phq.py $V/phreeqc_respeciation.json ${O}_f0_y1s/prof:060:flood_ctrl ${O}_f20_y1s/prof:060:flood_b20 ${O}_f20_y1s/prof:100:flood_b20_late \
  ${O}_a20_y1s/prof:037:awd_drained ${O}_a20_y1s/prof:046:awd_reflooded ${O}_f20_y1f/prof:052:fallow_b20 ${O}_f0_y1f/prof:118:fallow_ctrl_end > $V/phq.log 2>&1
python3 -I $V/massbal2.py $V/mass_balance.json ${O}_{f0,f20,a0,a20}_y1{s,f} > $V/massbal.log 2>&1
cd $V && python3 -I fig3_verification.py
