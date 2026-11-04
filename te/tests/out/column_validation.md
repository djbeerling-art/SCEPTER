# TE column vs PHREEQC TRANSPORT

20 cells, 60 shifts, pure advection; exchanger 0.1 mol + Hfo (11.0 g) + gibbsite per cell; acid soil solution (pH 5) flushed with Ca-Mg-HCO3 water (pH 7.5), no metal inflow.
PHREEQC 0.6 s, Python column 28.0 s.

| metal | pool | max abs err / scale | scale (mol) |
|---|---|---|---|
| Ni | total | 1.76e-06 | 8.42e-05 |
| Ni | exchange | 2.83e-06 | 8.42e-05 |
| Ni | dissolved | 5.35e-05 | 5.76e-07 |
| Ni | effluent (last cell) | 1.84e-05 | 1.09e-06 |
| Ni | mass left after 60 shifts | PHREEQC 0.968, Python 0.968 of initial | |
| Cu | total | 1.82e-06 | 7.22e-05 |
| Cu | exchange | 3.70e-06 | 7.22e-05 |
| Cu | dissolved | 3.63e-05 | 4.96e-07 |
| Cu | effluent (last cell) | 1.91e-05 | 9.41e-07 |
| Cu | mass left after 60 shifts | PHREEQC 0.967, Python 0.967 of initial | |
| Pb | total | 5.98e-07 | 0.000277 |
| Pb | exchange | 4.85e-06 | 0.000277 |
| Pb | dissolved | 4.91e-05 | 4.36e-07 |
| Pb | effluent (last cell) | 2.59e-05 | 8.26e-07 |
| Pb | mass left after 60 shifts | PHREEQC 0.993, Python 0.993 of initial | |
| Cd | total | 1.31e-06 | 2.88e-06 |
| Cd | exchange | 3.07e-06 | 2.88e-06 |
| Cd | dissolved | 4.92e-05 | 1.44e-08 |
| Cd | effluent (last cell) | 1.84e-05 | 2.46e-08 |
| Cd | mass left after 60 shifts | PHREEQC 0.978, Python 0.978 of initial | |

Worst error 5.35e-05: PASS
