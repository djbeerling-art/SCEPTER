"""Mass-balance screening models for trace-element (TE) accumulation under ERW.

These are the total-concentration screens used in the literature to flag
regulatory risk. They are the comparators for the mechanistic SCEPTER-TE
results (same rock, soil, dose and years).

1. Dupla et al. (2023, Eur. J. Soil Sci. 74, e13343), Eqs. 1-6: three
   scenarios for total TE in the 0-20 cm layer after n annual applications.
   - 'minimal'  : 1 wt%/yr rock dissolution (Eqs. 1-4); released TE minus
                  annual leaching and crop export accumulates.
   - 'maximal'  : full dissolution in the year of application (Eq. 5).
   - 'inert'    : no dissolution; rock accumulates and dilutes/enriches (Eq. 6).
   Leaching fractions per year (Tyler 1978, pH 4.2, 1000 mm): Ni 0.166,
   Cu 0.049, Cr 0.026, Zn 0.098. Other metals are not given by Dupla
   (default 0, i.e. conservative); crop export as mg/kg soil per year (Table 1).
2. ERW-MAC (Olagaray et al., CDRXiv preprint https://cdrxiv.org/preprint/546),
   Eq. 1: all TE in the applied rock is released instantly and stays in the
   mixing layer (no leaching, no crop export, homogeneous mixing). ERW-MAC
   samples soil depth U(0.05, 0.3) m and bulk density N(1250, 250) truncated
   to 800-1700 kg/m3, basalt cumulative doses 25-125 t/ha:
       C_final = C_soil0 + C_rock * app_rate * n_app / (depth * BD).

All concentrations are mg/kg (total, aqua-regia-like). Doses in t/ha.
"""
import numpy as np

DUPLA_LEACH = {'Ni': 0.166, 'Cu': 0.049, 'Cr': 0.026, 'Zn': 0.098}
DEPTH, BD_SOIL0, BD_ROCK, R_WEATH = 0.20, 1300.0, 2000.0, 0.01
# Dupla et al. 2023 SI: GEOROC means (Table S1), Kabata-Pendias world soil
# means (start values), aqua-regia extractability used to adjust limits, and
# maize grain removal (Table S4, g/ha/yr).
DUPLA_LIP_MEAN = {'Cr': 165.3, 'Cu': 129.1, 'Ni': 79.8, 'Zn': 99.4}
DUPLA_SOIL0 = {'Cr': 59.5, 'Cu': 38.9, 'Ni': 29.0, 'Zn': 70.0}
DUPLA_EXTRACT = {'Ni': 0.884, 'Cr': 0.559, 'Cu': 0.906, 'Zn': 0.895}
DUPLA_CROP_G_HA = {'Ni': 2.30, 'Cr': 1.84, 'Cu': 17.26, 'Zn': 121.8}


def dupla(c_rock, c_soil0, dose_t_ha_yr, years, *, scenario='minimal', leach=0.0,
          crop_mg_kg_yr=0.0, depth=DEPTH, bd_soil0=BD_SOIL0, bd_rock=BD_ROCK,
          r_weath=R_WEATH):
    """Total TE (mg/kg) in the mixing layer after each year 1..years.
    Annual application at the start of each year. Returns an array."""
    vol = depth * 1.0                      # m3 soil per m2
    app = dose_t_ha_yr * 0.1               # kg rock / m2 / yr
    out = np.zeros(years)
    if scenario == 'inert':                # Eq. 6 (arithmetic-geometric sequence)
        f = app / (vol * bd_rock)          # volume fraction replaced per year
        for n in range(1, years + 1):
            out[n - 1] = c_rock + (c_soil0 - c_rock) * (1 - f) ** n
        return out
    if scenario == 'maximal':              # Eq. 5, fixed point of the yearly map
        keep = (1 - leach)
        c = c_soil0
        for n in range(1, years + 1):
            add = c_rock * app / (bd_soil0 * vol)
            c = c + add * keep - crop_mg_kg_yr
            out[n - 1] = c
        return out
    if scenario == 'minimal':
        # Eqs. 1-4 as published (reproduces Dupla's mean-LIP result: first Cu
        # limit after 6 yr, Ni after 11 yr at 40 t/ha/yr). Note Eq. 4.2 as
        # printed, Mass_i = Rate_w*Mass_{i-1} + Mass_app, so the rock that
        # weathers each year is ~one year's application, while the undissolved
        # rock (Eq. 2.1) accumulates all applications; losses (Eq. 4) are
        # applied once to the cumulative weathered TE.
        m_i, sum_m = app, 0.0
        for n in range(1, years + 1):
            if n > 1:
                m_i = r_weath * m_i + app
            sum_m += m_i
            v_rock = app / bd_rock * sum((1 - r_weath) ** i for i in range(1, n + 1))
            bd_n = (v_rock * bd_rock + (vol - v_rock) * bd_soil0) / vol
            nat = c_soil0 * (vol - v_rock) / vol
            unw = c_rock * v_rock / vol
            wth = c_rock * r_weath * sum_m / (bd_n * vol)
            out[n - 1] = nat + unw + wth * (1 - leach) - crop_mg_kg_yr
        return out
    # 'minimal_consistent': 1 %/yr of the whole undissolved stock weathers,
    # rock displaces soil, and leaching/crop losses act every year.
    mass_rock = 0.0                        # undissolved rock kg/m2
    weathered = 0.0                        # net TE from weathered rock, mg/m2
    for n in range(1, years + 1):
        mass_rock = mass_rock + app                      # add this year's rock
        mass_rock_after = mass_rock * (1 - r_weath)
        released = c_rock * r_weath * mass_rock          # mg/m2 this year
        weathered = (weathered + released) * (1 - leach) - crop_mg_kg_yr * vol * bd_soil0
        weathered = max(weathered, 0.0)
        v_rock = mass_rock_after / bd_rock
        bd_n = (v_rock * bd_rock + (vol - v_rock) * bd_soil0) / vol
        nat = c_soil0 * (vol - v_rock) * bd_soil0 / (vol * bd_n)
        unw = c_rock * mass_rock_after / (vol * bd_n)
        out[n - 1] = nat + unw + weathered / (vol * bd_n)
        mass_rock = mass_rock_after
    return out


def erw_mac(c_rock, c_soil0, dose_t_ha, n_app, *, depth=0.2, bd=1250.0):
    """ERW-MAC Eq. 1, instant full release, no losses: total TE (mg/kg)."""
    return c_soil0 + c_rock * dose_t_ha * 0.1 * n_app / (depth * bd)


def years_to_limit(series, limit):
    """First year (1-based) at which a yearly series exceeds limit, or None."""
    i = np.nonzero(np.asarray(series) > limit)[0]
    return int(i[0]) + 1 if len(i) else None


full_release = erw_mac   # generic alias
