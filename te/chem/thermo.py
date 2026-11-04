"""Read PHREEQC-format thermodynamic databases.

Only the parts the TE module needs are parsed: SOLUTION_MASTER_SPECIES,
SOLUTION_SPECIES, PHASES, EXCHANGE_MASTER_SPECIES/EXCHANGE_SPECIES and
SURFACE_MASTER_SPECIES/SURFACE_SPECIES, with log_k, delta_h, -analytic and
-gamma. Reading the same file PHREEQC reads keeps the constants identical
between the prototype, the PHREEQC benchmark and (later) the Fortran module.
"""
import math
import re
from dataclasses import dataclass, field

LN10 = math.log(10.0)
R_KJ = 8.314462618e-3  # kJ/mol/K
T25 = 298.15

_KEYWORDS = {
    'SOLUTION_MASTER_SPECIES', 'SOLUTION_SPECIES', 'PHASES',
    'EXCHANGE_MASTER_SPECIES', 'EXCHANGE_SPECIES',
    'SURFACE_MASTER_SPECIES', 'SURFACE_SPECIES',
}
_CHARGE = re.compile(r'([+-])(\d*\.?\d*)$')
_COEF = re.compile(r'^(\d+\.?\d*|\.\d+)([A-Za-z(].*)$')


def charge(name):
    """Charge from a PHREEQC species name (Ni+2, CrO4-2, Hfo_wOH2+, e-)."""
    m = _CHARGE.search(name)
    if not m:
        return 0.0
    mag = float(m.group(2)) if m.group(2) else 1.0
    return mag if m.group(1) == '+' else -mag


def parse_side(text):
    """'Ca+2 + 2 X-' -> [('Ca+2', 1.0), ('X-', 2.0)]."""
    toks = text.split()
    out, coef = [], None
    for t in toks:
        if t == '+':
            continue
        try:
            coef = float(t)
            continue
        except ValueError:
            pass
        m = _COEF.match(t)
        if m:
            c, name = float(m.group(1)), m.group(2)
        else:
            c, name = 1.0, t
        if coef is not None:
            c *= coef
            coef = None
        out.append((name, c))
    return out


@dataclass
class Reaction:
    """lhs = rhs, with log K at 25 C and its temperature dependence."""
    lhs: list
    rhs: list
    log_k: float = 0.0
    delta_h: float = 0.0          # kJ/mol
    analytic: tuple = None
    gamma: tuple = None           # (a, b) WATEQ Debye-Hueckel parameters
    defined: str = ''             # species (or phase) this reaction defines

    def logk(self, T):
        if self.analytic:
            a = list(self.analytic) + [0.0] * (6 - len(self.analytic))
            return (a[0] + a[1] * T + a[2] / T + a[3] * math.log10(T)
                    + a[4] / T ** 2 + a[5] * T ** 2)
        return self.log_k - self.delta_h / (R_KJ * LN10) * (1.0 / T - 1.0 / T25)


@dataclass
class Database:
    masters: dict = field(default_factory=dict)      # 'Cr(3)' -> 'Cr(OH)2+'
    species: dict = field(default_factory=dict)      # aqueous species -> Reaction
    phases: dict = field(default_factory=dict)       # phase name -> Reaction
    ex_masters: dict = field(default_factory=dict)   # 'X' -> 'X-'
    ex_species: dict = field(default_factory=dict)
    surf_masters: dict = field(default_factory=dict)  # 'Hfo_w' -> 'Hfo_wOH'
    surf_species: dict = field(default_factory=dict)


def _strip(line):
    return line.split('#', 1)[0].strip()


def _parse_dh(tokens):
    val = float(tokens[0])
    unit = tokens[1].lower() if len(tokens) > 1 else 'kj'
    if unit.startswith('kcal'):
        val *= 4.184
    elif unit.startswith('cal'):
        val *= 4.184e-3
    elif unit.startswith('j'):
        val *= 1e-3
    return val


def read(paths, db=None):
    """Parse one or more PHREEQC database/input files; later entries override."""
    db = db or Database()
    if isinstance(paths, str):
        paths = [paths]
    for path in paths:
        with open(path, encoding='latin-1') as f:
            lines = f.readlines()
        block, cur, phase_name = None, None, None
        for raw in lines:
            line = _strip(raw)
            if not line:
                continue
            word = line.split()[0]
            if word.upper() == word and word.isidentifier() and word[0].isalpha() \
                    and ('_' in word or word in ('PHASES', 'END', 'RATES')) \
                    and word not in db.species and '=' not in line:
                block = word if word in _KEYWORDS else None
                cur, phase_name = None, None
                continue
            if block is None:
                continue
            if block in ('SOLUTION_MASTER_SPECIES', 'EXCHANGE_MASTER_SPECIES',
                         'SURFACE_MASTER_SPECIES'):
                t = line.split()
                tgt = {'SOLUTION_MASTER_SPECIES': db.masters,
                       'EXCHANGE_MASTER_SPECIES': db.ex_masters,
                       'SURFACE_MASTER_SPECIES': db.surf_masters}[block]
                tgt[t[0]] = t[1]
                continue
            if '=' in line and not line.startswith('-'):
                l, r = line.split('=', 1)
                lhs, rhs = parse_side(l), parse_side(r)
                if block == 'PHASES':
                    cur = Reaction(lhs, rhs, defined=phase_name)
                    db.phases[phase_name] = cur
                else:
                    cur = Reaction(lhs, rhs, defined=rhs[0][0])
                    tgt = {'SOLUTION_SPECIES': db.species,
                           'EXCHANGE_SPECIES': db.ex_species,
                           'SURFACE_SPECIES': db.surf_species}[block]
                    tgt[cur.defined] = cur
                continue
            t = line.split()
            opt = t[0].lstrip('-').lower()
            if block == 'PHASES' and not line.startswith('-') and opt not in (
                    'log_k', 'logk', 'delta_h', 'deltah', 'analytic', 'analytical_expression',
                    'a_e', 'vm', 't_c', 'p_c', 'omega', 'add_logk', 'no_check'):
                phase_name, cur = t[0], None
                continue
            if cur is None:
                continue
            if opt in ('log_k', 'logk', 'log_k'):
                cur.log_k = float(t[1])
            elif opt in ('delta_h', 'deltah'):
                cur.delta_h = _parse_dh(t[1:])
            elif opt in ('analytic', 'analytical_expression', 'a_e', 'analytical'):
                cur.analytic = tuple(float(x) for x in t[1:7])
            elif opt == 'gamma':
                cur.gamma = (float(t[1]), float(t[2]) if len(t) > 2 else 0.0)
    return db
