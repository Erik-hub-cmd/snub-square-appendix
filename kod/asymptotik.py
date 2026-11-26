"""
Asymptotik för resistansen på snub square-gittret (forslag/asymptotik.md):

    R_1j(n) = √3/36 + (3/(4π)) (γ + log(2√6 |x̃|)) - ((ε_1+ε_j)/2) sin 2φ/(48π |x̃|^2) + O(|x̃|^-3),
    x̃ = n + d_j  (cellkoordinater), φ = arg x̃, ε = (1, -1, 1, -1),

med d_j = nodernas harmoniska ("elektriska") lägen. I geometriska avstånd r: R = (3/(4π)) log r + C + O(1/r),
C = √3/36 + (3/(4π))(γ + log(6 - 2√3)).

Kontroller (exakt i Q(z), z = exp(iπ/12), om inget annat sägs):
  S1  alla g_kl har residyn 3/8 i s = 0, och h_kl = g_kl - 3/(8s) = O(s^-2) (summan av residyerna i rötterna är 0).
  S2  (log X)'(0) = 2i√6, (log Z_in)'(0) = -2√6; X(-s) X(s) = 1 och Z_in(-s) Z_in(s) = 1 (exponentmängderna).
  S3  K_j = (1/2π) ∫_0^∞ (h_11 + h_jj) ds = √3/36 för j = 1..4 (residyformeln med log(-s)).
  S4  harmoniska lägen d_j (varje nod i medelpunkten av sina grannar) och tensorn Σ e e^T = (4/3) I över cellens
      10 kanter; ger oberoende c = 1/(π·4/3) = 3/(4π).
  S5  1/r-termen: -(2/π) Re[h_g(0)/p] = (3/(4π)) d_j·n/|n|^2 exakt, för båda valen av post (n2 >= 1 och n2 <= 0).
  S6  alla 1 763 värden i runs/tabell: |R - formel| · |x̃|^2 (flyttal), med geometriska avstånd och med x̃.
  S7  numerisk 1D-integral (2) långt ut (|x̃| upp till ca 250): felet · |x̃|^2 mot -((ε_1+ε_j)/2) sin 2φ/(48π)
      (ε = +1 för nod 1, 3; -1 för nod 2, 4).
  S8  Taylorutveckling av Bloch-laplacianen i harmonisk gauge, M~ = U* M U, U = diag(e^{i d_k·θ}), θ = t(cos φ, sin φ):
      det M~ = 96 t^2 - 7 t^4 + O(t^6); för alla 16 poster G~_ji = (3/4)/t^2 + 0/t + A_0(φ) + O(t), och
      A_0 = konstant - ((ε_i+ε_j)/2) sin 2φ/96 (exakt, sympy).  Ger satsen med 1/r^2-termen (forslag/asymptotik.md).
  S9  hela runs/tabell med 1/r^2-termen: felet är O(|x̃|^-3).
  S10 tabellen "Table check in numbers" i forslag/asymptotik.md (halvöppna intervall i verkligt geometriskt avstånd).
  S11 Lemma F (Step 7 i beviset av Sats 2): gaussiska Fouriertransformen av e^{2iφ}, numeriskt.
  S12 enheterna: |x̃| i cellenheter, ρ = L|x̃| geometriskt, L^2 = 2+√3; sin 2φ = cos 2φ_E.
Enheter: |x̃| och φ i cellkoordinater (1/r^2-koefficient 1/(48π)); i geometriska enheter (kantlängd 1) är den (2+√3)/(48π).
Körtid: ca 1 minut.  Kör:  .venv\\Scripts\\python kod\\asymptotik.py
"""
import json, math, os, sys
sys.stdout.reconfigure(encoding="utf-8")          # utskrift till fil/pipe i Windows (cp1252) kraschar annars
from fractions import Fraction as F
import numpy as np
import sympy as sp

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'kod'))
import exact_snub as es
from k24 import K, Z, I, S3, split_sqrt3
from residue_tools import Machine, laurent_dict, ROOTS, ZK, EX, EZ, ell, Q

ok = []


def check(name, cond):
    ok.append(bool(cond))
    print(('OK  ' if cond else 'FEL ') + name, flush=True)


ADJd = {(i, j): laurent_dict(es.ADJ[i, j], es.X, es.Z) for i in range(4) for j in range(4)}
DDZd = laurent_dict(sp.diff(es.DET, es.Z), es.X, es.Z)
M = Machine(ADJd, DDZd, 8)
PF = {}                                           # PF[kl][q] = [c_q1, c_q2, c_q3]  (koefficienter för (s-q)^-m)
for q in [0] + ROOTS:
    Xs, Zs, _, _, g = M.local(K(0) if q == 0 else ZK[q], 8)
    if q == 0:
        X0, Z0, G0 = Xs, Zs, g
    for kl, ser in g.items():
        PF.setdefault(kl, {})[q] = [ser.coef(-m) for m in range(1, 4)]
    if q == 0:
        ord0 = all(ser.v >= -1 for ser in g.values())

# ---- S1
check('S1: residyn i s = 0 är 3/8 (enkel pol) för alla 16 g_kl', ord0 and all(PF[kl][0][0] == Q(3, 8) for kl in PF))
check('S1: Σ_q Res_q g_kl = 0 över rötterna, dvs. h_kl = g_kl - 3/(8s) = O(s^-2), för alla 16 g_kl',
      all(sum((PF[kl][q][0] for q in ROOTS), K(0)).is_zero() for kl in PF))

# ---- S2
S6 = (ZK[3] + ZK[21]) * S3                        # √2 = z^3 + z^-3
check("S2: (log X)'(0) = 2i√6 och (log Z_in)'(0) = -2√6", X0.coef(1) == 2 * I * S6 and Z0.coef(1) == -2 * S6)
odd = all(EX[(k + 12) % 24] == -EX[k] and EZ[(k + 12) % 24] == -EZ[k] for k in ROOTS)
check('S2: X(-s) = 1/X(s), Z_in(-s) = 1/Z_in(s) (k -> k+12 byter täljare och nämnare), så log X, log Z_in är udda', odd)

# ---- S3
def Kj(j):
    A, B = K(0), K(0)
    for q in ROOTS:
        c = [PF[(0, 0)][q][m] + PF[(j, j)][q][m] for m in range(3)]
        A = A + K(ell(q)) * c[0]
        for m in (2, 3):                          # Res_q[c (s-q)^-m log(s/q)] = c (-1)^m / ((m-1) q^(m-1))
            B = B + c[m - 1] * Q((-1) ** m, m - 1) * ZK[q].inv() ** (m - 1)
    a = split_sqrt3(A * I * Q(-1, 2))
    return (fr(a[0]), fr(a[1])), fr(B.rational())             # K_j = a0 + a1 √3 - B/(2π)


def fr(x):
    return F(int(x.p), int(x.q))
kv = [Kj(j) for j in range(4)]
check('S3: K_j = √3/36 för j = 1, 2, 3, 4  ' + str([f'{a[0]} + {a[1]}√3 - ({b})/(2π)' for a, b in kv]),
      all(a == (0, F(1, 36)) and b == 0 for a, b in kv))

# ---- S4
dsym = [sp.Matrix([0, 0])] + [sp.Matrix(sp.symbols(f'x{i} y{i}')) for i in (1, 2, 3)]
eqs = []
for v in range(4):
    S = sp.zeros(2, 1)
    for i, j, n in es.EDGES:
        if i == v:
            S += dsym[j] + sp.Matrix(n) - dsym[i]
        if j == v:
            S += dsym[i] - sp.Matrix(n) - dsym[j]
    eqs += list(S)
sol = sp.solve(eqs, list(sp.symbols('x1 y1 x2 y2 x3 y3')))
D = {0: (F(0), F(0))}
for i in (1, 2, 3):
    D[i] = (F(str(sol[sp.Symbol(f'x{i}')])), F(str(sol[sp.Symbol(f'y{i}')])))
T = [[F(0)] * 2 for _ in range(2)]
for i, j, n in es.EDGES:
    e = (D[j][0] + n[0] - D[i][0], D[j][1] + n[1] - D[i][1])
    for a in range(2):
        for b in range(2):
            T[a][b] += e[a] * e[b]
check(f'S4: harmoniska lägen d = {[(str(x), str(y)) for x, y in D.values()]}', D[1] == (F(1, 2), F(-1, 6)) and D[2] == (F(1, 3), F(-2, 3)) and D[3] == (F(-1, 6), F(-1, 2)))
check(f'S4: Σ_kanter e e^T = (4/3) I  {[[str(x) for x in r] for r in T]}  => c = 1/(π·4/3) = 3/(4π)', T == [[F(4, 3), 0], [0, F(4, 3)]])

# ---- S5:  h(0) = (√6/8)(u + iv) med rationella u, v;  -(2/π)Re[h/p] = -(u b - v a)/(8π|n|^2), p = 2√6(b - ia)
def uv(kl):
    w = G0[kl].coef(0) / (S6 * Q(1, 8))
    re = (w + w.conj()) * Q(1, 2)
    im = (w - w.conj()) * Q(1, 2) / I
    return fr(re.rational()), fr(im.rational())
s5 = True
for j in (1, 2, 3):
    u, v = uv((j, 0))                              # n2 >= 1:  (a, b) = (n1, n2)  =>  d·n = -(u n2 - v n1)/6
    s5 &= (v / 6, -u / 6) == D[j]
    u, v = uv((0, j))                              # n2 <= 0:  (a, b) = (-n1, -n2) =>  d·n = (u n2 - v n1)/6
    s5 &= (-v / 6, u / 6) == D[j]
    s5 &= G0[(j, j)].coef(0).is_zero()
check('S5: 1/r-termen -(2/π)Re[h_g(0)/p] är exakt (3/(4π)) d_j·n/|n|^2 för j = 2, 3, 4 och båda valen av post', s5 and G0[(0, 0)].coef(0).is_zero())

# ---- S6
g_ = 0.57721566490153286
Cn = math.sqrt(3) / 36 + 3 / (4 * math.pi) * (g_ + math.log(6 - 2 * math.sqrt(3)))
TAB = {}
for r in json.load(open(os.path.join(ROOT, 'runs', 'tabell', 'checkpoint.json'), encoding='utf-8'))['results'].values():
    TAB.update(r)


def asym(j, n1, n2):
    r = math.hypot(n1 + D[j][0], n2 + D[j][1])
    return math.sqrt(3) / 36 + 3 / (4 * math.pi) * (g_ + math.log(2 * math.sqrt(6) * r)), r
worstE = worstH = 0
for key, d in TAB.items():
    j, n1, n2 = map(int, key.split(':'))
    if (n1, n2) == (0, 0):
        continue
    Rv = float(d['float'])                       # 20 siffror (koefficienterna tar ut varandra kraftigt)
    A, r = asym(j, n1, n2)
    if r >= 6:
        worstH = max(worstH, abs(Rv - A) * r * r)
    if d['dist'] >= 10:
        worstE = max(worstE, abs(Rv - 3 / (4 * math.pi) * math.log(d['dist']) - Cn) * d['dist'])
print(f'    C = {Cn:.15f}')
check(f'S6: hela runs/tabell, geometriskt avstånd r >= 10: max |R - (3/(4π))log r - C|·r = {worstE:.4f} (felet är O(1/r))', worstE < 0.03)
check(f'S6: hela runs/tabell, elektriskt avstånd |x̃| >= 6: max |R - formel|·|x̃|^2 = {worstH:.4f} (felet är O(1/r^2))', worstH < 0.015)

# ---- S7: numerisk 1D-integral R = (1/π) Re ∫_1^∞ f ds (bevis_formodan8_direkt.md, (2))
z = np.exp(1j * np.pi / 12)
PFc = {kl: [(0 if q == 0 else z**q, [complex(c.to_complex()) for c in PF[kl][q]]) for q in [0] + ROOTS] for kl in PF}
xg, wg = np.polynomial.legendre.leggauss(40)


def gval(kl, s):
    return sum(c / (s - qq)**m for qq, cs in PFc[kl] for m, c in enumerate(cs, 1))


def Rnum(j, n1, n2):
    kl, a, b = ((j, 0), n1, n2) if n2 >= 1 else ((0, j), -n1, -n2)
    N = 60 + 8 * (abs(a) + abs(b))
    tt = np.concatenate([np.linspace(0, 12, N + 1), np.linspace(12, 40, 60)[1:]])     # t = log s
    tot = 0
    for t0, t1 in zip(tt[:-1], tt[1:]):
        s = np.exp((t1 - t0) / 2 * xg + (t1 + t0) / 2)
        Xv = np.prod([(s - z**k)**EX[k] for k in ROOTS], axis=0)
        Zv = np.prod([(s - z**k)**EZ[k] for k in ROOTS], axis=0).real
        f = gval((0, 0), s) + gval((j, j), s) - 2 * gval(kl, s) * Xv**a * Zv**b
        tot += np.sum(wg * (f * s).real) * (t1 - t0) / 2
    return tot / np.pi
w4 = max(abs(Rnum(*map(int, k.split(':'))) - float(TAB[k]['float'])) for k in ['1:5:5', '3:-2:-6', '0:7:0', '2:-4:1'])
check(f'S7: 1D-integralen reproducerar tabellvärden (fel {w4:.1e})', w4 < 1e-12)
EPS = {0: 1, 1: -1, 2: 1, 3: -1}
worst = 0
for (j, u, v) in [(0, -1, 1), (2, -1, 1), (1, -1, 1), (3, 1, 1), (0, 3, 1), (0, 2, 1), (2, 3, -1), (1, 3, 1), (0, 1, 0), (3, 2, 1)]:
    for k in (25, 80):
        n1, n2 = u * k, v * k
        A, r = asym(j, n1, n2)
        x = (n1 + float(D[j][0]), n2 + float(D[j][1]))
        pred = -(EPS[0] + EPS[j]) / 2 * (2 * x[0] * x[1] / r**2) / (48 * math.pi)
        worst = max(worst, abs((Rnum(j, n1, n2) - A) * r * r - pred))
check(f'S7: |x̃| = 25…250: (R - formel)·|x̃|^2 = -((ε_1+ε_j)/2) sin2φ/(48π) + O(1/r)  [största avvikelse {worst:.1e}]', worst < 2e-4)

# ---- S8: Taylorutveckling i harmonisk gauge
t, c, sn, ph = sp.symbols('t c s phi', real=True)
NT = 5


def ex(e):                                         # exp(i t e·ω), trunkerad
    u = sp.I * t * (sp.Rational(e[0]) * c + sp.Rational(e[1]) * sn)
    return sum(u**k / sp.factorial(k) for k in range(NT + 1))


def trunc(expr, n):
    P = sp.Poly(sp.expand(expr), t)
    return sum(cf * t**k for (k,), cf in P.terms() if k <= n)
Mt = 5 * sp.eye(4)
for i, j, m in es.EDGES:
    e = (D[j][0] + m[0] - D[i][0], D[j][1] + m[1] - D[i][1])
    Mt[i, j] -= ex(e)
    Mt[j, i] -= ex((-e[0], -e[1]))
Mt = Mt.applyfunc(lambda e: trunc(e, NT))
red = lambda e: sp.expand(sp.expand(e).subs(sn**2, 1 - c**2))
det = sp.expand(Mt.det(method='berkowitz'))
dc = [red(det.coeff(t, k)) for k in range(6)]
check(f'S8: det M~ = 96 t^2 - 7 t^4 + O(t^6)  (koefficienter t^0..t^5: {dc})', dc == [0, 0, 96, 0, -7, 0])
adj = Mt.adjugate()
s8 = True
for i in range(4):
    for j in range(4):
        a = sp.expand(adj[j, i])
        a0, a1, a2 = (red(a.coeff(t, k)) for k in range(3))
        A0 = sp.expand(a2 / 96 + 7 * a0 / 96**2)       # 1/det = (1/96t^2)(1 + 7t^2/96 + ...)
        A0 = sp.expand(sp.expand_trig(A0.subs({c: sp.cos(ph), sn: sp.sin(ph)})))
        F2 = sp.simplify(sp.fu(A0) - sp.integrate(A0, (ph, 0, 2 * sp.pi)) / (2 * sp.pi))
        want = -sp.Rational(EPS[i] + EPS[j], 2) * sp.sin(2 * ph) / 96
        s8 &= (a0 == 72) and (a1 == 0) and sp.simplify(F2 - want) == 0
check('S8: för alla 16 poster: A_-2 = (3/4)/t^2, A_-1 = 0 (=> ingen 1/r-term), A_0 - medelvärde = -((ε_i+ε_j)/2) sin2φ/96', s8)

# ---- S9: tabellen med 1/r^2-termen
worst3 = 0
for key, d in TAB.items():
    j, n1, n2 = map(int, key.split(':'))
    if (n1, n2) == (0, 0):
        continue
    A, r = asym(j, n1, n2)
    x = (n1 + float(D[j][0]), n2 + float(D[j][1]))
    A2 = A - (EPS[0] + EPS[j]) / 2 * (2 * x[0] * x[1] / r**2) / (48 * math.pi * r * r)
    if r >= 6:
        worst3 = max(worst3, abs(float(d['float']) - A2) * r**3)
check(f'S9: hela runs/tabell med 1/r^2-termen, |x̃| >= 6: max |R - formel|·|x̃|^3 = {worst3:.4f} (felet är O(1/r^3))', worst3 < 0.1)

# ---- S10: tabellen i forslag/asymptotik.md ("Table check in numbers"), med halvöppna intervall a <= r < b
TROWS = [((1.5, 2.5), 13, 6.8e-3, 0.016), ((7.5, 8.5), 54, 3.0e-3, 0.024), ((15.5, 16.5), 110, 1.4e-3, 0.023),
         ((23.5, 24.5), 43, 9.2e-4, 0.022), ((27, 29), 12, 7.3e-4, 0.020), ((10, math.inf), 1423, 2.2e-3, 0.023),
         ((0.5, math.inf), 1763, 2.3e-2, 0.028)]
ERR = [(d['dist'], float(d['float']) - 3 / (4 * math.pi) * math.log(d['dist']) - Cn) for d in TAB.values() if d['dist'] > 0.5]
s10 = True
print('    r in          par   max|e|    max|e|·r')
for (a, b), npar, me, mer in TROWS:
    sel = [(r, e) for r, e in ERR if a <= r < b]
    m1, m2 = max(abs(e) for _, e in sel), max(abs(e) * r for r, e in sel)
    print(f'    [{a}, {b})'.ljust(18) + f'{len(sel):5d}   {m1:.1e}   {m2:.3f}')
    s10 &= len(sel) == npar and float(f'{m1:.1e}') == me and round(m2, 3) == mer
check('S10: tabellen i forslag/asymptotik.md återskapas exakt (antal par och avrundade värden)', s10)

# ---- S11: Lemma F (Step 7) numeriskt: (1/4π²)∫ e^{-|θ|²/2} e^{2iφ} e^{ix·θ} dθ = -(e^{2iψ}/(π|x|²))(1 - e^{-q}(1+q)), q = |x|²/2
xr, wr = np.polynomial.legendre.leggauss(200)
rho = 7.5 * (xr + 1)
wr = 7.5 * wr
phi = 2 * np.pi * np.arange(512) / 512
s11 = 0
for xv in [(3.0, 1.7), (-0.6, 2.2), (5.0, -4.0)]:
    inner = np.exp(2j * phi)[None, :] * np.exp(1j * rho[:, None] * (xv[0] * np.cos(phi) + xv[1] * np.sin(phi))[None, :])
    val = np.sum(wr * rho * np.exp(-rho**2 / 2) * inner.mean(axis=1)) * 2 * np.pi / (4 * np.pi**2)
    X2 = xv[0]**2 + xv[1]**2
    q = X2 / 2
    pred = -np.exp(2j * math.atan2(xv[1], xv[0])) / (math.pi * X2) * (1 - math.exp(-q) * (1 + q))
    s11 = max(s11, abs(val - pred))
check(f'S11: Lemma F (gaussiska fallet) stämmer numeriskt (fel {s11:.1e})', s11 < 1e-12)

# ---- S12: enheterna (Units): L^2 = 2+√3, och sin 2φ (cellkoordinater) = cos 2φ_E (planet), a1 = (s, s), a2 = (s, -s)
sq = (1 + math.sqrt(3)) / 2
Lg = math.hypot(sq, sq)
s12 = abs(Lg**2 - (2 + math.sqrt(3))) < 1e-14
for xv in [(1.3, 0.4), (-2.0, 0.7), (0.2, -3.1)]:
    E = (sq * (xv[0] + xv[1]), sq * (xv[0] - xv[1]))
    s12 &= abs(math.sin(2 * math.atan2(xv[1], xv[0])) - math.cos(2 * math.atan2(E[1], E[0]))) < 1e-14
    s12 &= abs(math.hypot(*E) - Lg * math.hypot(*xv)) < 1e-12
check(f'S12: 1/r²-koefficienten är 1/(48π) = {1 / (48 * math.pi):.6f} i cellenheter och (2+√3)/(48π) = {(2 + math.sqrt(3)) / (48 * math.pi):.6f} geometriskt; sin 2φ = cos 2φ_E', s12)

print('\nSUMMARY:', sum(ok), 'passed,', sum(not k for k in ok), 'failed')
