"""
Isoradiell förklaring, steg 9.2: Kenyons formel för det isoradiella Z² med ALLMÄN vinkel Φ, exakt i t.

Parametrisering (steg 8): t = q = 1 + λ + √(λ(λ+2)), λ = (t−1)²/(2t), κ = √(λ(λ+2)) = (t²−1)/(2t),
    e^{iΦ} = i(t − i)/(t + i)   (dvs. cos Φ = 2t/(t²+1) = 1/(1+λ), sin Φ = (t²−1)/(t²+1) = κ/(1+λ)).
Tågspårsvinklar: α = {0, −Φ}, β = {π/2 − Φ, π/2}; alla poler ±e^{iγ} ligger i Q(i)(t).
Residysumman ger G = g₀(t) + g₁(t)/π + g₂(t)·Φ/π, där g₀, g₁, g₂ räknas exakt i Q(i)(t).

Skriptet kontrollerar (alla exakt):
  (K1) g₀, g₁, g₂ är reella, dvs. ligger i Q(t)                          [Lemma K, del (a)]
  (K2) logaritmens gren (heltalet k i γ + 2πk) är densamma för t i (0,1) och (1,∞)   [Lemma K, del (b)]
  (K3) t = 2 + √3 ger steg 6:s värden (λ = 1); t → 1 ger det vanliga Z² (G(e₁) = −1/4)
  (K4) ΔG = δ som identitet i t (alla tre delar), alltså för alla t > 0
Sedan, med A(t), B(t) från steg 8, för de 53 noderna med avstånd < 4 från nod 1:
  (R1) R = r₀ + r₁/π + r₂·Φ/π med r₀, r₁ invarianta och r₂ antiinvariant under τ: t ↦ 1/t
       (alltså R = a(λ) + b(λ)/π + c(λ)·κΦ/π med a, b, c ∈ Q(λ))                            [Sats, bevisets slutsteg]
  (R2) a, b, c i λ = 1/2, 2, 3 är exakt de rationella tal som PSLQ gav (runs/isoradiell_lambda/koefficienter_159.json)
  (R3) i λ = 1 ger a + c·√3/3 + b/π tabellvärdet i runs/tabell (om filen finns)
  (R4) nämnarna i a, b, c innehåller bara λ, λ + 1, λ + 2 (a, b, c ∈ Q[λ, 1/λ, 1/(1+λ), 1/(2+λ)])
Skriver de slutna formerna till runs/isoradiell_lambda/SLUTNA_FORMER.md.
Kör: python kod/isoradiell_steg9_kenyon_t.py
"""
import sys, os, json, math, cmath
sys.stdout.reconfigure(encoding="utf-8")
import sympy as sp
from fractions import Fraction as F
from functools import lru_cache

HERE = os.path.dirname(os.path.abspath(__file__))
t, X, Z, lam_s = sp.symbols('t X Z lambda')
K = sp.QQ_I.frac_field(t)
ONE, ZERO, I = K.one, K.zero, K.from_sympy(sp.I)
W = K.from_sympy(sp.I * (t - sp.I) / (t + sp.I))          # e^{iΦ}
WI = ONE / W
ok = []


def check(name, cond):
    ok.append(bool(cond)); print(('OK  ' if cond else 'FEL ') + name, flush=True)


def Phi_num(tv):
    return 2 * math.atan((tv - 1) / (tv + 1))


# vinkel = (a, ε): γ = aπ + εΦ;  e^{iγ} som element i K
ALPHA = {0: (F(0), 0), 1: (F(0), -1)}             # NE-spår u = x − y
BETA = {0: (F(1, 2), -1), 1: (F(1, 2), 0)}        # NW-spår w = x + y


def ekv(ang):
    a, e = ang
    r = {F(0): ONE, F(1, 2): I, F(1): -ONE, F(3, 2): -I}[a % 2]
    return r * (W if e == 1 else WI if e == -1 else ONE)


def steps(p1, p2):
    S = []
    for (a, b, tab) in ((p1[0] - p1[1], p2[0] - p2[1], ALPHA), (p1[0] + p1[1], p2[0] + p2[1], BETA)):
        if b > a: S += [(tab[k % 2], 1) for k in range(a, b)]
        else: S += [(tab[k % 2], -1) for k in range(b, a)]
    return S


def smul(a, b, n):
    out = [ZERO] * n
    for i in range(min(n, len(a))):
        if a[i] != ZERO:
            for j in range(min(n - i, len(b))): out[i + j] += a[i] * b[j]
    return out


def sinv(a, n):
    i0 = ONE / a[0]; out = [i0]
    for k in range(1, n):
        s = ZERO
        for j in range(1, min(k, len(a) - 1) + 1): s += a[j] * out[k - j]
        out.append(-(s * i0))
    return out


SAMPLES = [F(1, 50), F(1, 3), F(9, 10), F(11, 10), F(3), F(50)]


def branch(S, pole_ang, tv):
    """Heltalet k så att γ = aπ + εΦ + 2πk ligger i (θ₀ − π, θ₀ + π], θ₀ = arg(v₂ − v₁) i rombinbäddningen."""
    Ph = Phi_num(tv)
    disp = sum(s * cmath.exp(1j * (float(a) * math.pi + e * Ph)) for (a, e), s in S)
    th0 = cmath.phase(disp)
    g = float(pole_ang[0]) * math.pi + pole_ang[1] * Ph
    k = 0
    while g + 2 * math.pi * k > th0 + math.pi + 1e-9: k -= 1
    while g + 2 * math.pi * k <= th0 - math.pi + 1e-9: k += 1
    return k


BRANCH_OK = [True]


@lru_cache(maxsize=None)
def G(p1, p2):
    """(g0, g1, g2) i K med G(p1, p2) = g0 + g1/π + g2·Φ/π (Kenyon, residyer)."""
    if p1 == p2: return (ZERO, ZERO, ZERO)
    S = steps(p1, p2)
    fac = [(ekv(ang), s) for ang, s in S]
    poles = {}
    for ang, s in S:
        pa = ang if s == 1 else (ang[0] + 1, ang[1])
        key = (pa[0] % 2, pa[1])
        poles[key] = poles.get(key, 0) + 1
    s0 = s1 = s2 = ZERO                     # Σ Res = π·s0 + Φ·s2 + s1
    # Grenen: alla poler ligger i en öppen sektor med vinkel < π runt θ₀, så k ändras inte när t varierar; ett
    # gemensamt hopp för alla poler (θ₀ passerar ±π i atan2) ändrar inte G, eftersom Σ_p Res(E/z) = 1 − E(0) = 0.
    for i, tv in enumerate(SAMPLES):
        diff = {branch(S, pa, float(tv)) - branch(S, pa, 3.0) for pa in poles}
        if len(diff) != 1: BRANCH_OK[0] = False
    for pa, m in poles.items():
        k = branch(S, pa, 3.0)
        a_p, e_p = pa[0] + 2 * k, pa[1]
        e = ekv(pa)
        h = [ONE] + [ZERO] * (m - 1)
        for ek, s in fac:
            num, den = [e + ek, ONE], [e - ek, ONE]
            if s == -1: num, den = den, num
            if den[0] == ZERO: den = den[1:] + [ZERO]
            if num[0] == ZERO: raise ValueError('nollställe och pol möts')
            h = smul(h, smul(num, sinv(den, m), m), m)
        h = smul(h, sinv([e, ONE], m), m)
        alpha = h[m - 1]
        beta, pw, ei = ZERO, ONE, ONE / e
        for j in range(1, m):
            pw = pw * ei
            beta += h[m - 1 - j] * pw * K.from_sympy(sp.Rational((-1) ** (j + 1), j))
        s0 += I * K.from_sympy(sp.Rational(a_p.numerator, a_p.denominator)) * alpha
        s2 += I * K.from_sympy(sp.Integer(e_p)) * alpha
        s1 += beta
    q = K.from_sympy(sp.Rational(-1, 4))
    return (q * s0, q * s1, q * s2)


def conj(f):
    return K.from_sympy(K.to_sympy(f).subs(sp.I, -sp.I))


def tau(f):                                  # t ↦ 1/t
    return K.from_sympy(sp.cancel(K.to_sympy(f).subs(t, 1 / t)))


# ---------------- Kenyon-delen ----------------
PTS = [(0, 0), (1, 0), (2, 1), (1, 1), (3, 1), (2, 2), (3, 0)]
real = all(conj(g) == g for p in PTS for g in G((0, 0), p))
check('(K1) g₀, g₁, g₂ ∈ Q(t) (reella) för G((0,0), p), p bland 6 punkter', real)
s3 = sp.sqrt(3); T1 = 2 + s3
def Gval(p1, p2, tv):
    g0, g1, g2 = G(p1, p2); ph = 2 * sp.atan((tv - 1) / (tv + 1))
    return K.to_sympy(g0).subs(t, tv) + K.to_sympy(g1).subs(t, tv) / sp.pi + K.to_sympy(g2).subs(t, tv) * ph / sp.pi
v1 = sp.simplify(Gval((0, 0), (1, 0), T1) - (-sp.Rational(1, 6) - s3 / 12))
v2 = sp.simplify(sp.expand(Gval((0, 0), (2, 1), T1) - (sp.Rational(13, 12) + s3 / 3 - (3 + 2 * s3) / sp.pi)))
check('(K3) λ = 1: G((0,0),(1,0)) = −1/6 − √3/12 och G((0,0),(2,1)) = 13/12 + √3/3 − (3+2√3)/π (som steg 6)',
      v1 == 0 and v2 == 0)
def Gexpr(p1, p2):
    g0, g1, g2 = G(p1, p2)
    return K.to_sympy(g0) + K.to_sympy(g1) / sp.pi + K.to_sympy(g2) * 2 * sp.atan((t - 1) / (t + 1)) / sp.pi
lim1, lim2 = sp.limit(Gexpr((0, 0), (1, 0)), t, 1), sp.limit(Gexpr((0, 0), (1, 1)), t, 1)
check('(K3) t → 1 (vanligt Z²): G(e₁) → −1/4, G((0,0),(1,1)) → −1/π  (g₀, g₁, g₂ själva kan ha pol i t = 1)',
      sp.simplify(lim1 + sp.Rational(1, 4)) == 0 and sp.simplify(lim2 + 1 / sp.pi) == 0)

def cond(p, r):
    if p[1] != r[1]: return ONE                                   # vertikal kant
    x = min(p[0], r[0]); return [ONE / K.from_sympy(t), K.from_sympy(t)][(x + p[1]) % 2]   # h₀₀ = h₁₁ = 1/t, h₀₁ = h₁₀ = t
def lapG(y_):
    out = [ZERO, ZERO, ZERO]
    for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1)):
        w_ = (y_[0] + dx, y_[1] + dy); c = cond(y_, w_)
        gy, gw = G((0, 0), y_), G((0, 0), w_)
        for u in range(3): out[u] += c * (gy[u] - gw[u])
    return out
harm = True
for y_ in [(0, 0), (1, 0), (0, 1), (1, 1), (-1, 0), (2, 1), (-1, 2), (3, -2)]:
    lg = lapG(y_)
    harm = harm and lg[0] == (ONE if y_ == (0, 0) else ZERO) and lg[1] == ZERO and lg[2] == ZERO
check('(K4) ΔG = δ exakt som identitet i t (alla tre delar, 8 punkter); gäller alltså för alla t > 0, även t < 1', harm)

# ---------------- kombinationen med A(t), B(t) ----------------
src = open(os.path.join(HERE, 'isoradiell_steg8_allmant_lambda.py'), encoding='utf-8').read()
ns = {'sp': sp, 'X': X, 'Z': Z, 't': t}
exec(src[src.index('A = sp.Matrix'):src.index("\ncheck('")], ns)
A, B = ns['A'], ns['B']
D = 4 * (t + 1) ** 4 * (t ** 2 + 1) ** 2 / (t - 1) ** 4
Binv = (B.adjugate() / D)


def laurent(e):
    num, den = sp.fraction(sp.cancel(sp.together(e)))
    # den = c(t)·X^a·Z^b
    pd = sp.Poly(den, X, Z); assert len(pd.terms()) == 1
    (dx, dz), dc = pd.terms()[0]
    out = {}
    for (ex, ez), c in sp.Poly(sp.expand(num), X, Z).terms():
        out[(ex - dx, ez - dz)] = K.from_sympy(sp.cancel(c / dc))
    return out


AL = [[laurent(A[i, j]) for j in range(4)] for i in range(4)]
BL = [[laurent(Binv[i, j]) for j in range(4)] for i in range(4)]
PT = {0: (0, 0), 1: (0, 1), 2: (1, 0), 3: (1, 1)}
SIG = (0, 2, 3, 1); OFF = [(0, 0), (0, 0), (0, -1), (0, -1)]


def Gz(i, m, j):
    p1 = (PT[i][0] + 2 * m[0], PT[i][1] + 2 * m[1]); p2 = PT[j]
    return G(p1, p2) if p1 <= p2 else G(p2, p1)


def R_snub(a, b, n):
    w = [{} for _ in range(4)]
    w[a][(0, 0)] = w[a].get((0, 0), ZERO) + ONE
    w[b][(-n[0], -n[1])] = w[b].get((-n[0], -n[1]), ZERO) - ONE
    p = [{} for _ in range(4)]; r = [{} for _ in range(4)]
    for i in range(4):
        for k in range(4):
            for (ex, ez), c in AL[k][i].items():
                for (wx, wz), cw in w[k].items():
                    key = (-ex + wx, -ez + wz); p[i][key] = p[i].get(key, ZERO) + c * cw
            for (ex, ez), c in BL[i][k].items():
                for (wx, wz), cw in w[k].items():
                    key = (ex + wx, ez + wz); r[i][key] = r[i].get(key, ZERO) + c * cw
    tot = [ZERO, ZERO, ZERO]
    for i in range(4):
        for kk, cp in p[i].items():
            if cp == ZERO: continue
            for j in range(4):
                for ll, cr in r[j].items():
                    if cr == ZERO: continue
                    g = Gz(i, (ll[0] - kk[0], ll[1] - kk[1]), j)
                    c = cp * cr
                    for u in range(3): tot[u] += c * g[u]
    return tot


def to_lambda(f):
    """f ∈ Q(t) invariant under τ: t ↦ 1/t  ->  rationell funktion av λ.
    f = N·τ(D) / (D·τ(D)); både täljare och nämnare är då symmetriska Laurentpolynom, alltså polynom i
    y = t + 1/t = 2(λ + 1)."""
    if f == ZERO: return sp.Integer(0)
    y = sp.Symbol('y')
    N, D = sp.fraction(sp.cancel(K.to_sympy(f)))
    tD = sp.expand(D.subs(t, 1 / t))
    def in_y(P):
        P = sp.expand(P); c = {}
        for term in sp.Add.make_args(P):
            co, ex = term.as_coeff_exponent(t); c[int(ex)] = c.get(int(ex), 0) + co
        res = 0
        while any(v != 0 for v in c.values()):
            k = max(kk for kk, v in c.items() if v != 0); cc = c[k]
            if k < 0 or (k > 0 and c.get(-k, 0) != cc): raise ValueError('inte symmetrisk')
            res += cc * y ** k
            for j in range(k + 1):
                c[k - 2 * j] = c.get(k - 2 * j, 0) - cc * sp.binomial(k, j)
        return res
    e = in_y(N * tD) / in_y(D * tD)
    return sp.factor(sp.cancel(e.subs(y, 2 * (lam_s + 1))))


KAPPA = K.from_sympy((t ** 2 - 1) / (2 * t))


if __name__ == '__main__':

    coef = json.load(open(os.path.join(HERE, '..', 'runs', 'isoradiell_lambda', 'koefficienter_159.json'), encoding='utf-8'))
    nodes = sorted({k.split('|')[1] for k in coef}, key=lambda s: tuple(map(int, s.split(':'))))
    if len(sys.argv) > 1: nodes = nodes[:int(sys.argv[1])]
    tab = {}
    tp = os.path.join(HERE, '..', 'runs', 'tabell', 'checkpoint.json')
    if os.path.exists(tp):
        for tt in json.load(open(tp, encoding='utf-8'))['results'].values():
            tab.update(tt)
    nsym = nper = nl = n1 = 0
    rader = []
    for nd in nodes:
        j, n1_, n2_ = map(int, nd.split(':'))
        a_, b_ = SIG[0], SIG[j]; nn = (n1_ + OFF[j][0] - OFF[0][0], n2_ + OFF[j][1] - OFF[0][1])
        r0, r1, r2 = R_snub(a_, b_, nn)
        sym = tau(r0) == r0 and tau(r1) == r1 and tau(r2) == -r2
        nsym += sym
        c_t = r2 / KAPPA
        fa, fb, fc = (to_lambda(f) for f in (r0, r1, c_t))
        good = True
        for ls in ('1/2', '2', '3'):
            lv = sp.Rational(ls)
            got = [sp.Rational(x) for x in coef[f'{ls}|{nd}']]
            val = [sp.cancel(f.subs(lam_s, lv)) for f in (fa, fb, fc)]   # ingen nsimplify (gissar fel)
            good = good and got == val
        nper += good
        R1 = sp.radsimp(fa.subs(lam_s, 1) + fc.subs(lam_s, 1) * s3 / 3 + fb.subs(lam_s, 1) / sp.pi)
        t1 = 'saknas'
        if nd in tab:
            t1 = sp.simplify(R1 - sp.sympify(tab[nd]['expr'])) == 0
            n1 += t1
        nl += nd in tab
        print(f'{nd}: τ-symmetri {"ja" if sym else "NEJ"}, λ=1/2,2,3 {"ja" if good else "NEJ"}, λ=1 {t1}:  '
              f'a = {fa}, b = {fb}, c = {fc}', flush=True)
        rader.append((nd, fa, fb, fc, sym, good, t1))
    check(f'(K2) logaritmens gren densamma för t i (0,1) och (1,∞), sånär som på ett gemensamt 2π (prov t = {", ".join(map(str, SAMPLES))})', BRANCH_OK[0])
    check(f'(R1) r₀, r₁ invarianta, r₂ antiinvariant under t ↦ 1/t: {nsym} av {len(nodes)} noder', nsym == len(nodes))
    check(f'(R2) a(λ), b(λ), c(λ) i λ = 1/2, 2, 3 = PSLQ-koefficienterna: {nper} av {len(nodes)} noder (= {3 * nper} av {3 * len(nodes)} rader)',
          nper == len(nodes))
    if nl:
        check(f'(R3) λ = 1: a + c√3/3 + b/π = runs/tabell: {n1} av {nl} noder', n1 == nl)
    def okden(f):
        den = sp.fraction(sp.factor(f))[1]
        return all(sp.Poly(fac, lam_s).degree() == 0 or sp.expand(fac - lam_s) == 0 or sp.expand(fac - lam_s - 1) == 0
                   or sp.expand(fac - lam_s - 2) == 0 for fac, _ in sp.factor_list(den)[1]) if den.has(lam_s) else True
    nden = sum(all(okden(f) for f in r[1:4]) for r in rader)
    check(f'(R4) nämnarna i a, b, c innehåller bara λ, λ+1, λ+2: {nden} av {len(rader)} noder', nden == len(rader))
    L = ['# Slutna former för allmänt λ (steg 9.2)', '',
         'R(nod 1 → nod) = a(λ) + b(λ)/π + c(λ)·κ·Φ/π,  κ = √(λ(λ+2)),  Φ = arccos(1/(1+λ)) = arctan κ.',
         'Räknat exakt ur Kenyons formel för det isoradiella Z² och A(t), B(t) från steg 8 (kod/isoradiell_steg9_kenyon_t.py).',
         'Kolumnerna "λ = 1/2, 2, 3" och "λ = 1" säger om värdena stämmer exakt med PSLQ (KONTROLL_159.md) och runs/tabell.', '',
         '| nod j:n1:n2 | a(λ) | b(λ) | c(λ) | λ = 1/2, 2, 3 | λ = 1 |', '|---|---|---|---|---|---|']
    for nd, fa, fb, fc, sym, good, t1 in rader:
        L.append(f'| {nd} | ${sp.latex(fa)}$ | ${sp.latex(fb)}$ | ${sp.latex(fc)}$ | {"ja" if good else "NEJ"} | '
                 f'{"ja" if t1 is True else ("NEJ" if t1 is False else t1)} |')
    if len(sys.argv) == 1:
        open(os.path.join(HERE, '..', 'runs', 'isoradiell_lambda', 'SLUTNA_FORMER.md'), 'w', encoding='utf-8').write('\n'.join(L) + '\n')
    print(f'\n{sum(ok)} godkända, {len(ok) - sum(ok)} fel')
