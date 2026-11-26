"""
Sandhögens höjd-1-sannolikhet P_1 och gradfördelningen i det likformiga spännträdet, exakt, för
  (a) Cairo-gittret, båda nodtyperna (grad 3 och grad 4), och
  (b) snub square, med exakt jämförelse mot Goodchilds P_1 och hans sex övergångsströmmar Y (GOODCHILD.md, Tabell II och III).
Resistanserna kommer från residy- och Galoismaskinen (residue_tools.py, forslag/bevis_formodan8_direkt.md och
forslag/bevis_cairo.md), alltså ur satsen själv; inget PSLQ.

Metod (som jobb_ust.py): för kanterna e_1..e_z vid en nod v, alla orienterade ut ur v, är
    Y(e, f) = (R(v,a) + R(v,b) - R(a,b)) / 2   (e = (v,a), f = (v,b)),    Y(e, e) = R(e),
    P(e in T för e i S, e inte i T annars) = det(M_S),  M_S = rader Y_e (e i S), I_e - Y_e (annars)   [Burton-Pemantle]
    P(grad = k) = Σ_{|S|=k} det(M_S).
Majumdar-Dhar: P_1 = det(I - Y') där Y' är Y utan en kant (z-1 kanter); Goodchild: P(grad 1) = z P_1.
Alla tal är polynom i √3 och p = 1/π med rationella koefficienter (p behandlas som en fri symbol; √3² reduceras till 3).

Kontroller (alla exakta, som polynom i (√3, p)):
  S1  Σ_f Y(e,f) = 1 för varje e (strömbevarande), Σ_k P(grad k) = 1, väntevärdet = Σ_e R(e).
  S2  P_1 = det(I - Y') är samma för alla z val av borttagen kant och P(grad 1) = z P_1  (Goodchilds korollarium).
  S3  Cairo: väntevärdena 2 - √3/9 (grad 3) och 2 + 2√3/9 (grad 4); genomsnitt över cellen = 2.
  S4  snub square: P_1 = Goodchilds Tabell II-uttryck exakt, och hans sex Y (Tabell III) finns i vår Y-matris vid en
      kantmärkning a, b, d, a', b' (hittas med sökning), exakt.
  S5  P_1 = τ(G\\v)/τ(G) på Cairo-torusar (determinanter av hela Laplacianen, N = 16, 32, Richardson) mot det exakta P_1.
  S6  Wilsons algoritm på Cairo-torusar mot de exakta gradsannolikheterna (valfritt: python3 sandhog_y.py sim).
Körtid: ca en minut utan simulering.  Kör:  python3 sandhog_y.py
"""
import itertools, os, sys, time
from fractions import Fraction
import sympy as sp
import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import exact_snub as es
import cairo_lattice as cl
import residue_tools as rt

s3, p = sp.symbols('s3 p')
ok = []


def check(name, cond):
    ok.append(bool(cond))
    print(('OK  ' if cond else 'FEL ') + name, flush=True)


def red(e):
    """polynom i s3, p: reducera s3^2 -> 3"""
    e = sp.expand(e)
    return sp.expand(sp.rem(sp.Poly(e, s3), sp.Poly(s3**2 - 3, s3)).as_expr())


def to_sym(t):
    return sp.Rational(t[0].numerator, t[0].denominator) + sp.Rational(t[1].numerator, t[1].denominator) * s3 \
        + sp.Rational(t[2].numerator, t[2].denominator) * p


def val(e, dps=25):
    return sp.N(e.subs({s3: sp.sqrt(3), p: 1 / sp.pi}), dps)


def build_machine(adj, det, nmax=14):
    X, Z = es.X, es.Z
    ADJd = {(i, j): rt.laurent_dict(adj[i, j], X, Z) for i in range(adj.shape[0]) for j in range(adj.shape[1])}
    return rt.Machine(ADJd, rt.laurent_dict(sp.diff(det, Z), X, Z), nmax)


def neighbours(edges, nv):
    nb = {i: [] for i in range(nv)}
    for i, j, n in edges:
        nb[i].append((j, tuple(n)))
        nb[j].append((i, (-n[0], -n[1])))
    return nb


class Lattice:
    def __init__(self, name, mach, nb, deg):
        self.name, self.mach, self.nb, self.deg = name, mach, nb, deg
        self._R = {}

    def R(self, ja, na, jb, nb_):
        """resistans mellan (ja, cell na) och (jb, cell nb_), exakt (polynom i s3, p)"""
        d = (nb_[0] - na[0], nb_[1] - na[1])
        if ja == jb and d == (0, 0):
            return sp.Integer(0)
        key = (ja, jb, d)
        if key not in self._R:
            self._R[key] = to_sym(self.mach.resistance(ja, jb, d))
        return self._R[key]

    def Ymatrix(self, v):
        nbr = self.nb[v]
        z = len(nbr)
        Rv = [self.R(v, (0, 0), j, n) for j, n in nbr]
        Y = sp.zeros(z, z)
        for x in range(z):
            for y in range(z):
                (ja, na), (jb, nb_) = nbr[x], nbr[y]
                Y[x, y] = red((Rv[x] + Rv[y] - self.R(ja, na, jb, nb_)) / 2) if x != y else Rv[x]
        return Y, Rv

    def analyse(self, v, tag):
        Y, Rv = self.Ymatrix(v)
        z = Y.shape[0]
        # S1
        rows = [red(sum(Y[e, f] for f in range(z))) for e in range(z)]
        probs = {}
        for k in range(z + 1):
            tot = 0
            for S in itertools.combinations(range(z), k):
                Mx = sp.zeros(z, z)
                for r in range(z):
                    for c in range(z):
                        Mx[r, c] = Y[r, c] if r in S else (int(r == c) - Y[r, c])
                tot += Mx.det(method='berkowitz')
            probs[k] = red(tot)
        mean = red(sum(k * v_ for k, v_ in probs.items()))
        check(f'S1 {tag}: Σ_f Y(e,f) = 1 för alla {z} kanter, Σ_k P(grad k) = 1, väntevärde = Σ_e R(e) = {sp.collect(red(sum(Rv)), s3)}',
              all(r == 1 for r in rows) and red(sum(probs.values())) == 1 and mean == red(sum(Rv)))
        # S2: P_1 = det(I - Y') för varje borttagen kant
        P1s = []
        for drop in range(z):
            idx = [e for e in range(z) if e != drop]
            Ysub = Y.extract(idx, idx)
            P1s.append(red((sp.eye(z - 1) - Ysub).det(method='berkowitz')))
        check(f'S2 {tag}: P_1 = det(I - Y\') är samma för alla {z} val av borttagen kant och P(grad 1) = {z}·P_1',
              all(q == P1s[0] for q in P1s) and red(z * P1s[0] - probs[1]) == 0)
        return Y, probs, P1s[0]


def show(tag, probs, P1):
    print(f'  {tag}:')
    print(f'    P_1 = {sp.collect(sp.expand(P1), p)}  = {val(P1)}')
    for k in sorted(probs):
        if k:
            print(f'    P(grad {k}) = {sp.collect(sp.expand(probs[k]), p)}  = {val(probs[k])}')


t0 = time.time()
RESULT = {}
# ------------------------------------------------------------------ Cairo
pos, deg, edges = cl.build()
Lc = cl.bloch_laplacian(deg, edges)
mc = build_machine(Lc.adjugate().applyfunc(sp.expand), es.DET)
cairo = Lattice('Cairo', mc, neighbours(edges, 6), deg)
print('Cairo-gittret', flush=True)
for v, name in ((0, 'grad 3'), (4, 'grad 4')):
    Y, probs, P1 = cairo.analyse(v, f'Cairo {name}')
    RESULT[name] = (probs, P1)
    show(f'Cairo, nod av {name}', probs, P1)
mean3 = red(sum(k * q for k, q in RESULT['grad 3'][0].items()))
mean4 = red(sum(k * q for k, q in RESULT['grad 4'][0].items()))
check('S3 Cairo: väntevärdet 2 - √3/9 (grad 3), 2 + 2√3/9 (grad 4); cellens medelvärde (4·grad3 + 2·grad4)/6 = 2',
      mean3 == red(2 - s3 / 9) and mean4 == red(2 + 2 * s3 / 9) and red((4 * mean3 + 2 * mean4) / 6) == 2)
# Cairo: z-entropin-kontroll mot Goodchild finns inte; däremot P(grad 1)/z = P_1 ovan (S2).

# ------------------------------------------------------------------ snub square
print('\nSnub square', flush=True)
ms = build_machine(es.ADJ, es.DET)
snub = Lattice('snub square', ms, neighbours(es.EDGES, 4), [5] * 4)
Ys, probs_s, P1_s = snub.analyse(0, 'snub square')
show('snub square', probs_s, P1_s)
P1_GC = sp.Rational(8953, 3888) - 281 * s3 / 216 + (44 * s3 / 27 - sp.Rational(11, 4)) * p + (sp.Rational(17, 36) - s3 / 3) * p**2
check('S4a snub square: P_1 = 8953/3888 - 281√3/216 + (44√3/27 - 11/4)/π + (17/36 - √3/3)/π² (Goodchild, Tabell II), exakt',
      red(P1_s - P1_GC) == 0)

# S4b: Goodchilds Tabell III. Kanter vid noden: ett triangel-triangel-d och fyra triangel-kvadrat-kanter.
GC_Y = {('a', 'b'): 1 - s3 / 2, ('a', 'd'): s3 / 9, ("a'", "b'"): (3 - s3) / 6, ('d', "a'"): sp.Rational(1, 2) - 2 * s3 / 9,
        ('a', "a'"): 5 * s3 / 9 - 1 + p / 2, ('a', "b'"): sp.Rational(1, 2) - s3 / 9 - p / 2}
nbr = snub.nb[0]
tt = [e for e, (j, n) in enumerate(nbr) if red(Ys[e, e] - 2 * s3 / 9) == 0]
ts = [e for e in range(5) if e not in tt]
found = []
if len(tt) == 1:
    d = tt[0]
    for perm in itertools.permutations(ts):
        lab = dict(zip(['a', 'b', "a'", "b'"], perm)); lab['d'] = d
        if all(red(Ys[lab[x], lab[y]] - v_) == 0 for (x, y), v_ in GC_Y.items()):
            found.append(lab)
check(f'S4b snub square: Goodchilds sex Y (Tabell III) finns exakt i vår Y-matris ({len(found)} märkningar a, b, a\', b\', d av kanterna hittade)',
      len(found) >= 1)
if found:
    lab = found[0]
    inv = {v_: k for k, v_ in lab.items()}
    print('    kanter vid noden (granne: (nod, cell)):', {inv[e]: nbr[e] for e in range(5)})
    print('    vår Y-matris (rad/kolumn i ordningen a, b, a\', b\', d):')
    order = ['a', 'b', "a'", "b'", 'd']
    for x in order:
        print('     ', x, [str(sp.collect(Ys[lab[x], lab[y]], p)) for y in order])
    print('    Y på diagonalen: R(TS) = 1/2 - √3/18 och R(TT) = 2√3/9 =', [str(Ys[lab[x], lab[x]]) for x in order])
    vals = sorted({str(sp.expand(Ys[x, y])) for x in range(5) for y in range(x + 1, 5)})
    print(f'    {len(vals)} olika Y-värden utanför diagonalen')

# ------------------------------------------------------------------ S5: oberoende kontroll på ändliga torusar
# P_1 = τ(G \ v)/τ(G) (Majumdar-Dhar: spännträd där v är löv fäst vid en given granne ↔ spännträd i G utan v), räknat med
# determinanter av hela Laplacianen på N x N celler, utan någon resistans eller formel; Richardson i 1/N^2.
import numpy as np


def torus_P1(lat, v, N):
    V = lat.deg.__len__() * N * N
    nv = len(lat.deg)
    idx = lambda j, n1, n2: ((n1 % N) * N + (n2 % N)) * nv + j
    A = np.zeros((V, V))
    for n1 in range(N):
        for n2 in range(N):
            for j in range(nv):
                for jj, (a, b) in lat.nb[j]:
                    A[idx(j, n1, n2), idx(jj, n1 + a, n2 + b)] += 1
    Lap = np.diag(A.sum(1)) - A
    v0 = idx(v, 0, 0)
    keep = [u for u in range(V) if u != v0]
    L2 = Lap[np.ix_(keep, keep)] - np.diag(A[keep, v0])
    ld = lambda M_: np.linalg.slogdet(M_[1:, 1:])[1]
    return float(np.exp(ld(L2) - ld(Lap)))


for v, name in ((0, 'grad 3'), (4, 'grad 4')):
    a16, a32 = torus_P1(cairo, v, 16), torus_P1(cairo, v, 32)
    extrap = (4 * a32 - a16) / 3
    exact = float(val(RESULT[name][1]))
    check(f'S5 Cairo {name}: P_1 = τ(G\\v)/τ(G) på tori N=16, 32 (Richardson) = {extrap:.7f}, exakt {exact:.7f}', abs(extrap - exact) < 2e-5)

print(f'\n{sum(ok)} godkända, {len(ok) - sum(ok)} fel  ({time.time() - t0:.0f} s)')

# ------------------------------------------------------------------ S6: Wilson-simulering (valfri)
if sys.argv[1:] == ['sim']:
    import random, multiprocessing as mpp
    NB = cairo.nb
    DEG = cairo.deg
    exact = {name: {k: float(val(q)) for k, q in RESULT[name][0].items()} for name in RESULT}

    def sim(args):
        N, trees, seed = args
        V = 6 * N * N
        idx = lambda j, n1, n2: ((n1 % N) * N + (n2 % N)) * 6 + j
        adj = [None] * V
        for n1 in range(N):
            for n2 in range(N):
                for j in range(6):
                    adj[idx(j, n1, n2)] = [idx(jj, n1 + a, n2 + b) for jj, (a, b) in NB[j]]
        rng = random.Random(seed)
        hist = {3: [0] * 5, 4: [0] * 6}
        for _ in range(trees):
            in_tree = [False] * V
            nxt = [-1] * V
            root = rng.randrange(V)
            in_tree[root] = True
            for start in range(V):
                u = start
                while not in_tree[u]:
                    nxt[u] = adj[u][int(rng.random() * len(adj[u]))]
                    u = nxt[u]
                u = start
                while not in_tree[u]:
                    in_tree[u] = True
                    u = nxt[u]
            dg = [0] * V
            for v_ in range(V):
                if v_ != root:
                    dg[v_] += 1
                    dg[nxt[v_]] += 1
            for v_ in range(V):
                hist[len(adj[v_])][dg[v_]] += 1
        return hist

    out = {}
    for N, trees in ((24, 4 * 250), (48, 4 * 60)):
        with mpp.Pool(4) as pool:
            hs = pool.map(sim, [(N, trees // 4, 1000 * N + s_) for s_ in range(4)])
        out[N] = {z_: [sum(h[z_][k] for h in hs) for k in range(len(hs[0][z_]))] for z_ in (3, 4)}
    print('\nS6 Wilson-simulering på Cairo-torusar (N x N celler); Richardson (4 v(48) - v(24))/3 mot exakt:')
    for z_, name in ((3, 'grad 3'), (4, 'grad 4')):
        for k in range(1, z_ + 1):
            v24 = out[24][z_][k] / sum(out[24][z_])
            v48 = out[48][z_][k] / sum(out[48][z_])
            print(f'    {name}, P(grad {k}): N=24 {v24:.5f}  N=48 {v48:.5f}  extrapolerat {(4 * v48 - v24) / 3:.5f}  exakt {exact[name][k]:.5f}')
sys.exit(0 if all(ok) else 1)
