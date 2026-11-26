"""
Cairo-gittret: alla resistanser ligger i Q + Q√3 + Q/π (forslag/bevis_cairo.md).
Samma residy- och Galoisargument som för snub square (forslag/bevis_formodan8_direkt.md): Cairo-gittret har exakt samma
determinant, alltså samma parametrisering X(s), Z_in(s); bara adjunkten byts (6 x 6).

Kontroller:
  C0  det L_Cairo = det M_snub som Laurentpolynom (parametriseringen D0 ärvs).
  C1  Z-ordningarna i adj L: för varje par (i, j) och varje n2 finns en post utan bidrag från Z = 0.
  C2  alla 36 g_kl har poler bara i {0} ∪ P (ordning <= 3), samma residy i s = 0, ingen polynomdel (härledd gradgräns).
  C3  alla R_ij(n) i fönstret |n1|, |n2| <= NW (alla i, j) får formen A - B/(2π): Σ_q a_q = 0, a_q imaginära,
      A ∈ Q(√3), B ∈ Q.
  C4  jämförelse mot runs/cairo/tabell.json (exakt värde från exact_snub.py cairo), alla 13 + 14 klasser.
  C5  Kirchhoffs lag med termen κ_x (g3 - g4), g3 - g4 = √3/18, exakt i alla noder vars grannar ligger i fönstret.
  C6  numerisk kontroll av några allmänna par mot direkt rutnätsintegration av Ekv. (1) (ca 10 siffror).
Körtid: ca 2 minuter.  Kör:  python3 bevis_cairo.py
"""
import json, os, sys, time, itertools
from fractions import Fraction
import sympy as sp
import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import exact_snub as es
import cairo_lattice as cl
import residue_tools as rt

NW = 5                                   # fönster för C3 och C5
ok = []


def check(name, cond):
    ok.append(bool(cond))
    print(('OK  ' if cond else 'FEL ') + name, flush=True)


t0 = time.time()
pos, deg, edges = cl.build()
Lm = cl.bloch_laplacian(deg, edges)
ADJ = Lm.adjugate().applyfunc(sp.expand)
check('C0: det L_Cairo = det M_snub (Laurentpolynom)', sp.expand(Lm.det() - es.DET) == 0)

ldict = lambda e: rt.laurent_dict(e, es.X, es.Z)
ADJd = {(i, j): ldict(ADJ[i, j]) for i in range(6) for j in range(6)}
DDZd = ldict(sp.diff(es.DET, es.Z))
M = rt.Machine(ADJd, DDZd, 2 * NW + 12)

# ---- C1
bad = []
for i, j in itertools.product(range(6), repeat=2):
    for n2 in range(-3, 4):
        try:
            M.entry(i, j, (0, n2))
        except ValueError:
            bad.append((i, j, n2))
check('C1: för alla 36 par (i, j) och alla n2 finns en post med Z-ordning >= -1 efter faktorn Z^n2 (Z = 0 ger inget bidrag)'
      + (f'; undantag {bad}' if bad else ''), not bad)

# ---- C2
okpf, bound, res0 = rt.partial_fraction_check(M)
check(f'C2: alla 36 g_kl = partialbråk med poler i {{0}} ∪ P, ordning <= 3, gemensam residy i 0 ({res0.p}), ingen polynomdel '
      f'(gradgräns {bound}, {bound + 1} punkter)', okpf)

# ---- C3
t1 = time.time()
RES = {}
nbad = 0
for i, j in itertools.product(range(6), repeat=2):
    for n in itertools.product(range(-NW, NW + 1), repeat=2):
        if (i, j) == (j, i) and n == (0, 0):
            RES[(i, j, n)] = (Fraction(0),) * 3
            continue
        r = M.resistance(i, j, n)
        if r is None:
            nbad += 1
            continue
        RES[(i, j, n)] = tuple(Fraction(int(x.p), int(x.q)) for x in r)
check(f'C3: R = A - B/(2π), A ∈ Q(√3), B ∈ Q, Σ a_q = 0, a_q imaginära för alla {len(RES) + nbad} par i fönstret '
      f'(alla i, j, |n| <= {NW}; {time.time() - t1:.0f} s)', nbad == 0)
print('    exempel:', {k: RES[k] for k in [(0, 3, (0, -1)), (0, 4, (0, 0)), (4, 5, (-1, -1))]}, flush=True)

# ---- C4
def parse(s):
    e = sp.expand(sp.sympify(s))
    c1 = e.coeff(sp.sqrt(3), 1)
    c2 = e.coeff(sp.pi, -1)
    c0 = sp.simplify(e - c1 * sp.sqrt(3) - c2 / sp.pi)
    return tuple(Fraction(int(sp.Rational(c).p), int(sp.Rational(c).q)) for c in (c0, c1, c2))


rows = json.load(open(os.path.join(os.path.dirname(HERE), 'runs', 'cairo', 'tabell.json'), encoding='utf-8'))['rows']
nok = 0
for src, (j, n), d2, mult, exact, digits in rows:
    nok += RES.get((src, j, tuple(n))) == parse(exact)
check(f'C4: residyformeln = exakta värdet i runs/cairo/tabell.json för alla {len(rows)} parklasser ({nok} av {len(rows)})', nok == len(rows))

# ---- C5
nb = {a: [] for a in range(6)}
for a, b, (n1, n2) in edges:
    nb[a].append((b, (n1, n2)))
    nb[b].append((a, (-n1, -n2)))
S3 = sp.sqrt(3)
g34 = (Fraction(0), Fraction(1, 18), Fraction(0))
kappa = {3: -2, 4: 4}
nk = nbad5 = 0
for src in (0, 4):
    for j in range(6):
        for m in itertools.product(range(-NW, NW + 1), repeat=2):
            inside = all(max(abs(m[0] + n[0]), abs(m[1] + n[1])) <= NW for _, n in nb[j])
            if not inside:
                continue
            lhs = [Fraction(0)] * 3
            Rx = RES[(src, j, m)]
            for jj, n in nb[j]:
                Ry = RES[(src, jj, (m[0] + n[0], m[1] + n[1]))]
                lhs = [lhs[t] + Ry[t] - Rx[t] for t in range(3)]
            delta = 2 if (j, m) == (src, (0, 0)) else 0
            want = [kappa[deg[j]] * g34[t] for t in range(3)]
            want[0] += delta
            nk += 1
            nbad5 += (lhs != want)
check(f'C5: Σ_(y~x)(R(y) - R(x)) = 2δ + κ_x(g3 - g4) exakt i {nk} noder (källa av grad 3 och 4)', nbad5 == 0 and nk > 100)

# ---- C6  (numeric_resistance i exact_snub.py är byggd för tabellens avstånd; här används direkt rutnätsintegration
# med Richardson, K = 400 och 800, som i numeric_check.py)
import numeric_check as nc
worst = 0.0
for (i, j, n) in [(0, 1, (2, -3)), (4, 5, (3, 1)), (1, 4, (-2, 2)), (0, 0, (4, -1)), (2, 3, (-3, -2)), (5, 0, (1, 4))]:
    nu = (4 * nc.cairo_grid_resistance(i, j, n, 800) - nc.cairo_grid_resistance(i, j, n, 400)) / 3
    c = RES[(i, j, n)]
    ex = float(c[0]) + float(c[1]) * 3**0.5 + float(c[2]) / 3.141592653589793
    worst = max(worst, abs(nu - ex))
check(f'C6: residyformeln = direkt rutnätsintegration av Ekv. (1) för 6 allmänna par  [max avvikelse {worst:.1e}]', worst < 1e-8)

print(f'\n{sum(ok)} godkända, {len(ok) - sum(ok)} fel  ({time.time() - t0:.0f} s)')
sys.exit(0 if all(ok) else 1)
