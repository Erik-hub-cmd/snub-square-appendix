"""
Direkt bevis av Förmodan 8 för ALLA par (forslag/bevis_formodan8_direkt.md): residy- och Galoisargumentet från
diagonalkedjan (bevis_formodan8.py, del A) fungerar för varje R_1j(n1, n2), utan induktionen i del B.

Kontroller (exakt i Q(z), z = exp(iπ/12), om inget annat sägs):
  D0  parametriseringen: det M(X(s), Z_in(s)) = 0 som rationell funktion (antal punkter = härledd gradgräns + 1).
  D1  Z-ordningarna i adj M: inget bidrag från Z = 0 med valet G_j1 (n2 >= 1) resp. G_1j (n2 <= 0).
  D2  varje g_ij = r_ij (X'/X)/i har poler bara i de 8 primitiva 24:e enhetsrötterna och i 0 (residy 3/8):
      g_ij = (sin partialbråksutveckling) som rationella funktioner (antalet punkter = härledd gradgräns + 1, se D2).
  D3  R_1j(n) = A - B/(2π) med A = -(i/2) Σ_q ℓ_q a_q och B = Σ_q b_q; A ∈ Q(√3), B ∈ Q; jämförelse med HELA
      runs/tabell (1763 noder) exakt.
  D4  numerisk kontroll av 1D-integralen för några allmänna par.
Körtid: några minuter.  Kör:  .venv\\Scripts\\python kod\\bevis_formodan8_direkt.py
"""
import json, os, sys, time
import sympy as sp
import mpmath as mp
import flint

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'kod'))
from k24 import K, Z, I, S3, split_sqrt3
import exact_snub as es

ok = []


def check(name, cond):
    ok.append(bool(cond))
    print(('OK  ' if cond else 'FEL ') + name, flush=True)


def Q(a, b=1):
    return K(flint.fmpq(a, b))


EX = {1: 1, 5: 1, 7: 1, 11: 1, 13: -1, 17: -1, 19: -1, 23: -1}        # X(s)    = prod (s - z^k)^EX[k]
EZ = {1: 1, 5: 1, 19: 1, 23: 1, 7: -1, 11: -1, 13: -1, 17: -1}        # Z_in(s) = prod (s - z^k)^EZ[k]
ROOTS = [1, 5, 7, 11, 13, 17, 19, 23]
ZK = {k: Z**k for k in range(24)}


def ldict(expr):
    P = sp.Poly(sp.expand(expr * es.X**6 * es.Z**6), es.X, es.Z)
    return {(a - 6, b - 6): int(c) for (a, b), c in P.terms()}


ADJd = {(i, j): ldict(es.ADJ[i, j]) for i in range(4) for j in range(4)}
DDZd = ldict(sp.diff(es.DET, es.Z))

def width(d):
    """Grad (i s) för täljare och nämnare av L(X(s), Z_in(s)) för ett Laurentpolynom L(X, Z) = Σ c X^a Z^b.
    X = N_X/D_X och Z_in = N_Z/D_Z, där N och D har grad 4 i s. Multiplikation med N_X^p D_X^q N_Z^u D_Z^v,
    p = -min a, q = max a, u = -min b, v = max b, ger ett polynom i s av grad <= 4 (p+q+u+v); samma för nämnaren."""
    a = [k[0] for k in d]
    b = [k[1] for k in d]
    return 4 * ((max(a) - min(a)) + (max(b) - min(b)))


def eval_L(d, sv):
    """L(X(sv), Z_in(sv)) i K"""
    Xv = Zv = K(1)
    for k in ROOTS:
        Xv = Xv * (sv - ZK[k])**EX[k]
        Zv = Zv * (sv - ZK[k])**EZ[k]
    Xi, Zi = Xv.inv(), Zv.inv()
    return sum((K(c) * (Xv**a if a >= 0 else Xi**(-a)) * (Zv**b if b >= 0 else Zi**(-b)) for (a, b), c in d.items()), K(0))


# ---- D0: parametriseringen (Steg 1); ett rationellt uttryck av grad <= w som är noll i w + 1 punkter är identiskt noll
DETd = ldict(es.DET)
w0 = width(DETd)
check(f'D0: det M(X(s), Z_in(s)) = 0 som rationell funktion (gradgräns {w0}, kontrollerad i {w0 + 1} punkter)',
      all(eval_L(DETd, K(flint.fmpq(k + 2, 1))).is_zero() for k in range(w0 + 1)))

# ---- D1
zord = {ij: min(b for (_, b) in d) for ij, d in ADJd.items()}
check('D1: Z-ordning ≥ -2 för G_j1 (n2 ≥ 1), ≥ -1 för G_1j (n2 ≤ 0) och för diagonalen  ' +
      str({f'{j + 1}1': zord[(j, 0)] for j in range(4)}) + ' ' + str({f'1{j + 1}': zord[(0, j)] for j in range(4)}),
      all(zord[(j, 0)] >= -2 and zord[(0, j)] >= -1 and zord[(j, j)] >= -1 for j in range(4)))


# ---- Laurentserier i e = s - q över K:  (v, [c0, c1, ...]) = Σ c_k e^{v+k}
class L:
    __slots__ = ('v', 'c')

    def __init__(self, v, c):
        k = 0
        while k < len(c) and c[k].is_zero():
            k += 1
        self.v, self.c = v + k, c[k:]

    def __mul__(self, o):
        N = min(len(self.c), len(o.c))
        out = [K(0)] * N
        for i in range(N):
            a = self.c[i]
            if a.is_zero():
                continue
            for j in range(N - i):
                out[i + j] = out[i + j] + a * o.c[j]
        return L(self.v + o.v, out)

    def __add__(self, o):
        v = min(self.v, o.v)
        top = min(self.v + len(self.c), o.v + len(o.c))           # exakt kända ordningar < top
        out = [K(0)] * (top - v)
        for s_ in (self, o):
            for k, x in enumerate(s_.c):
                if s_.v + k < top:
                    out[s_.v + k - v] = out[s_.v + k - v] + x
        return L(v, out)

    def scale(self, x):
        return L(self.v, [y * x for y in self.c])

    def inv(self):
        assert self.c, 'ingen känd ledande term'
        a0i = self.c[0].inv()
        N = len(self.c)
        out = [a0i]
        for n in range(1, N):
            acc = K(0)
            for k in range(1, n + 1):
                acc = acc + self.c[k] * out[n - k]
            out.append(-acc * a0i)
        return L(-self.v, out)

    def coef(self, n):
        k = n - self.v
        assert k < len(self.c), 'för kort serie'
        return self.c[k] if k >= 0 else K(0)


def lin(c, N):                     # c + e
    return L(0, [c, K(1)] + [K(0)] * (N - 2)) if not c.is_zero() else L(1, [K(1)] + [K(0)] * (N - 1))


def lpow(x, n, cache):
    if n in cache:
        return cache[n]
    if n == 0:
        r = L(0, [K(1)] + [K(0)] * (len(x.c) - 1))
    elif n > 0:
        r = lpow(x, n - 1, cache) * x
    else:
        r = lpow(x, n + 1, cache) * cache['inv']
    cache[n] = r
    return r


def local(q, N):
    """Laurentdata vid s = q (q = 0 eller en primitiv rot): X, Z_in, X'/X och g_ij"""
    Xs = Zs = L(0, [K(1)] + [K(0)] * (N - 1))
    dl = L(0, [K(0)] * N)
    for k in ROOTS:
        f = lin(q - ZK[k], N)
        Xs = Xs * (f if EX[k] > 0 else f.inv())
        Zs = Zs * (f if EZ[k] > 0 else f.inv())
        dl = dl + f.inv().scale(K(EX[k]))
    cx, cz = {'inv': Xs.inv()}, {'inv': Zs.inv()}

    def ev(d):
        out = L(0, [K(0)] * N)
        for (a, b), c in d.items():
            out = out + (lpow(Xs, a, cx) * lpow(Zs, b, cz)).scale(K(c))
        return out
    den = (Zs * ev(DDZd)).inv()
    g = {ij: ev(ADJd[ij]) * den * dl.scale(I.inv()) for ij in [(0, 0), (1, 1), (2, 2), (3, 3)] +
         [(j, 0) for j in range(1, 4)] + [(0, j) for j in range(1, 4)]}
    return Xs, Zs, cx, cz, g


# ---- D2: partialbråk för varje g_ij, kontrollerade som rationella identiteter
NPF = 6
PF = {}
for q in [0] + ROOTS:
    qq = K(0) if q == 0 else ZK[q]
    _, _, _, _, g = local(qq, NPF + 4)
    for ij, ser in g.items():
        PF.setdefault(ij, []).append((qq, [ser.coef(-m) for m in range(1, 4)]))      # koefficienter för e^-1, e^-2, e^-3
okres0 = all(PF[ij][0][1][0] == Q(3, 8) and PF[ij][0][1][1].is_zero() for ij in PF)
check('D2: residyn i s = 0 är 3/8 (enkel pol) för alla g_ij som behövs', okres0)
okpf = True
# Gradgräns: g_ij = Num/Den med deg Num <= w_adj + w_ddz + 4 + 7 och deg Den <= w_adj + w_ddz + 4 + 8
# (täljare/nämnare av adj, av Z dZ det, faktorn 1/Z_in = D_Z/N_Z samt X'/X = Σ e_k/(s - z^k) = (täljare av grad 7)/(nämnare av grad 8)).
# Partialbråkssumman har nämnaren s^3 Π(s - z^k)^3 (grad 27) och täljare av grad <= 26. Täljaren i Num DenR - NumR Den
# har därför grad <= w_adj + w_ddz + 38.
w_ddz = width(DDZd)
bound = max(width(ADJd[ij]) for ij in PF) + w_ddz + 38
pts = [K(flint.fmpq(k + 2, 1)) for k in range(bound + 1)]
for ij in PF:
    for sv in pts:
        Xv = Zv = K(1)
        dl = K(0)
        for k in ROOTS:
            Xv = Xv * (sv - ZK[k])**EX[k]
            Zv = Zv * (sv - ZK[k])**EZ[k]
            dl = dl + K(EX[k]) / (sv - ZK[k])

        def evK(d):
            Xi, Zi = Xv.inv(), Zv.inv()
            return sum((K(c) * (Xv**a if a >= 0 else Xi**(-a)) * (Zv**b if b >= 0 else Zi**(-b)) for (a, b), c in d.items()), K(0))
        lhs = evK(ADJd[ij]) / (Zv * evK(DDZd)) * dl / I
        rhs = K(0)
        for qq, cs in PF[ij]:
            for m, c in enumerate(cs, 1):
                rhs = rhs + c / (sv - qq)**m
        okpf &= (lhs == rhs)
check('D2: g_ij = Σ_q Σ_m c_qm/(s-q)^m med q ∈ {0} ∪ {primitiva 24:e rötter}, m ≤ 3, och utan polynomdel (gradgräns ' + str(bound) + ', ' + str(bound + 1) + ' punkter, alla 10 g_ij)', okpf)

# ---- D3: hela tabellen
TAB = {}
for r in json.load(open(os.path.join(ROOT, 'runs', 'tabell', 'checkpoint.json'), encoding='utf-8'))['results'].values():
    for k, d in r.items():
        TAB[k] = d['coeffs']


def ell(k):
    """log(-z^k) = iπ ℓ (principalgren)"""
    m = (k + 12) % 24
    return flint.fmpq(m if m <= 12 else m - 24, 12)


t0 = time.time()
NMAX = 2 * 10 + 10
LOC = {k: local(ZK[k], NMAX) for k in ROOTS}
print(f'    (Laurentdata vid de 8 rötterna: {time.time() - t0:.0f} s)', flush=True)
nok = nbad = nA = nB = 0
nsum = nimag = 0
badlist = []
for key, cf in TAB.items():
    j, n1, n2 = map(int, key.split(':'))
    if n2 >= 1:
        ij, e1, e2 = (j, 0), n1, n2               # -2 Re(G_j1 X^n1 Z^n2)
    else:
        ij, e1, e2 = (0, j), -n1, -n2             # -2 Re(G_1j X^-n1 Z^-n2)
    A = K(0)
    B = K(0)
    asum = K(0)
    allimag = True
    for k in ROOTS:
        Xs, Zs, cx, cz, g = LOC[k]
        f = g[(0, 0)] + g[(j, j)] + (g[ij] * lpow(Xs, e1, cx) * lpow(Zs, e2, cz)).scale(K(-2))
        a = f.coef(-1)
        asum = asum + a
        allimag &= (a + a.conj()).is_zero()
        A = A + K(ell(k)) * a
        qi = ZK[k].inv()
        for m in range(1, -f.v):
            B = B + f.coef(-1 - m) * qi**m * Q((-1)**(m - 1), m)
    nsum += asum.is_zero()
    nimag += allimag
    A = A * I * Q(-1, 2)
    sA, rB = split_sqrt3(A), B.rational()
    nA += sA is not None
    nB += rB is not None
    if sA is None or rB is None:
        nbad += 1
        continue
    c0, c1, c2 = sA[0], sA[1], -rB / 2                   # R = c0 + c1 √3 + c2/π
    want = [flint.fmpq(sp.Rational(x).p, sp.Rational(x).q) for x in cf]
    if [c0, c1, c2] == want:
        nok += 1
    else:
        nbad += 1
        badlist.append(key)
print(f'    ({len(TAB)} noder på {time.time() - t0:.0f} s)')
check(f'D3: Σ_q a_q = 0 (dvs. f = O(s^-2) vid ∞; residyn i 0 är 3/8+3/8-2·3/8 = 0) för alla {len(TAB)} noder', nsum == len(TAB))
check(f'D3: varje residy a_q är rent imaginär (symmetrin s -> 1/conj s) för alla {len(TAB)} noder', nimag == len(TAB))
check(f'D3: A ∈ Q(√3) för alla {len(TAB)} noder', nA == len(TAB))
check(f'D3: B ∈ Q för alla {len(TAB)} noder', nB == len(TAB))
check(f'D3: R = A - B/(2π) ger EXAKT tabellvärdet för {nok} av {len(TAB)} noder' + (f'; avvikande: {badlist[:5]}' if badlist else ''),
      nok == len(TAB))

# ---- D4: numerisk kontroll av 1D-integralen (1/π) Re ∫_1^∞ f ds för några allmänna par
mp.mp.dps = 25
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
z = mp.expjpi(mp.mpf(1) / 12)
phi = (1 + mp.sqrt(5)) / 2


def gnum(ij, t):
    return sum(c.to_complex() * (1 / (t - (0 if qq.is_zero() else qq.to_complex())))**m
               for qq, cs in PF[ij] for m, c in enumerate(cs, 1))


def Xn(t):
    return mp.fprod((t - z**k)**EX[k] for k in ROOTS)


def Zn(t):
    return mp.fprod((t - z**k)**EZ[k] for k in ROOTS)


worst = 0
for key in ['1:0:0', '3:2:-3', '0:3:1', '2:-4:1', '1:5:5', '3:-2:-6', '0:7:0']:
    j, n1, n2 = map(int, key.split(':'))
    ij, e1, e2 = ((j, 0), n1, n2) if n2 >= 1 else ((0, j), -n1, -n2)
    fn = lambda t: (gnum((0, 0), t) + gnum((j, j), t) - 2 * gnum(ij, t) * Xn(t)**e1 * Zn(t)**e2).real
    v = mp.quad(fn, [1, 1.2, phi, 2, 4, 20, 200, mp.inf]) / mp.pi
    # referensvärde ur tabellens exakta koefficienter
    c = [sp.Rational(x) for x in TAB[key]]
    ref = mp.mpf(sp.N(c[0] + c[1] * sp.sqrt(3) + c[2] / sp.pi, 30))
    worst = max(worst, abs(v - ref))
check('D4: (1/π) Re ∫_1^∞ [g_11 + g_jj - 2 g X^a Z_in^b] ds = tabellvärdet för 7 allmänna par  [%s]' % mp.nstr(worst, 3), worst < 1e-18)

print('\nSUMMARY:', sum(ok), 'passed,', sum(not k for k in ok), 'failed')
