"""Kontroller för den omskrivna artikel 1 (2026-10-21): påståenden som är nya i texten eller skrivs ut på nytt.

  (T1) artikel/ab_matrices.tex (Bilaga A) är exakt matriserna A(t), B(t) ur isoradiell_steg8_allmant_lambda.py,
       och M_snub(λ(t))·A = B·M_t, det A = det B = 4(t+1)^4(t^2+1)^2/(t-1)^4.
  (T2) Ordlistan mellan cellnumreringen (j, n) och Z²-modellen: (j, n) -> P_j + 2n är en grafisomorfi som skickar
       triangel–triangel-kanterna på diagonalerna.
  (T3) Z²_t är isoradiellt: rombvinklarna ur tågspårsvinklarna ger halvvinklarna π/4, π/4 ∓ Φ/2 och konduktanserna
       1, 1/t, t; vinkelsumman är 2π runt varje nod och varje yta (t = 2, 5/7, 2 + √3).
  (T4) Bilaga E: partialbråk, faktoriseringar och Euler-substitutionerna för allmänt λ, i den form texten skriver dem.
  (T5) Tabell 3 (a, b, c för de tolv klasserna): λ = 1 ger Tabell 1, och λ = 1/2, 2, 3 ger PSLQ-talen i
       runs/isoradiell_lambda/checkpoint.json (läses bara som data).
  (T6) Kenyons log det-formel (Thm 1.1 i arXiv-versionen) för Z²_t har derivatan 2R_TT(λ) och värdet 16G/π i λ = 0,
       alltså är den lika med m(λ) (sats F); dessutom numeriskt i fem punkter.
  (T7) Kenyons c·R = 2θ/π med Jacobi ger R_TT(λ) (anmärkningen i avsnitt 4).
Kör: .venv\\Scripts\\python kod\\artikel1_kontroll.py           (kontroll)
     .venv\\Scripts\\python kod\\artikel1_kontroll.py --skriv   (skriver om artikel/ab_matrices.tex)
"""
import sys, os, json
sys.stdout.reconfigure(encoding='utf-8')
import sympy as sp
import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, '..')
ok = []


def check(name, cond):
    ok.append(bool(cond)); print(('OK  ' if cond else 'FEL ') + name, flush=True)


def zero(e):
    return sp.simplify(sp.factor(sp.together(sp.expand(e)))) == 0


# ---------------------------------------------------------------- (T1) A, B
src = open(os.path.join(HERE, 'isoradiell_steg8_allmant_lambda.py'), encoding='utf-8').read()
ns = {'sp': sp}
exec(src[src.index('X, Z, t = sp.symbols'):src.index("check('M_snub")].replace('ok = []', ''), ns)
X, Z, t = ns['X'], ns['Z'], ns['t']
A, B, Ms, Mz = ns['A'], ns['B'], ns['Ms'], ns['Mz']


def laurent_tex(e):
    """e = Σ c_ij(t) X^i Z^j  ->  LaTeX med faktoriserade koefficienter, monom sist."""
    num, den = sp.fraction(sp.cancel(sp.together(e)))
    pd = sp.Poly(den, X, Z)
    (dx, dz), dc = pd.terms()[0]
    assert len(pd.terms()) == 1
    terms = {}
    for (ex, ez), c in sp.Poly(sp.expand(num), X, Z).terms():
        terms[(ex - dx, ez - dz)] = sp.factor(c / dc)
    out = ''
    order = sorted(terms, key=lambda m: (m[0] != 0 or m[1] != 0, m))
    for m in order:
        c = terms[m]
        mon = ''
        for v, k in (('X', m[0]), ('Z', m[1])):
            if k == 1: mon += v
            elif k != 0: mon += f'{v}^{{{k}}}'
        if mon and c == 1: s = mon
        elif mon and c == -1: s = '-' + mon
        else:
            s = sp.latex(c)
            if mon:
                s = (s + '\\,' + mon) if not isinstance(c, sp.Add) else f'\\left({s}\\right){mon}'
        if out and not s.startswith('-'): out += '+'
        out += s
    return out or '0'


def ab_fragment():
    lines = ['% Genererad av kod/artikel1_kontroll.py --skriv (ändra inte för hand).']
    for nm, M in (('A', A), ('B', B)):
        lines.append('\\begin{align*}')
        rows = []
        for i in range(4):
            for j in range(4):
                rows.append(f'{nm}_{{{i + 1}{j + 1}}}&={laurent_tex(M[i, j])}')
        lines.append(',\\\\\n'.join(rows) + '.')
        lines.append('\\end{align*}')
    return '\n'.join(lines) + '\n'


FRAG = os.path.join(ROOT, 'artikel', 'ab_matrices.tex')
if '--skriv' in sys.argv:
    open(FRAG, 'w', encoding='utf-8').write(ab_fragment())
    print('skrev', FRAG)
check('artikel/ab_matrices.tex = A(t), B(t) ur steg 8 (genererad text, ordagrant)',
      os.path.exists(FRAG) and open(FRAG, encoding='utf-8').read() == ab_fragment())
check('M_snub·A = B·M_t identiskt i X, Z, t', all(zero(e) for e in (Ms * A - B * Mz)))
Dd = 4*(t + 1)**4*(t**2 + 1)**2/(t - 1)**4
check('det A = det B = 4(t+1)^4(t^2+1)^2/(t-1)^4', zero(A.det() - Dd) and zero(B.det() - Dd))
lam_t = (t - 1)**2/(2*t)
check('det M_snub(λ(t)) = det M_t (identiskt i X, Z, t)', zero(Ms.det() - Mz.det()))
one = sp.ones(4, 1)
check('A(1,1)·1 = (t+1)^2/(t-1)·1 och 1^T B(1,1) = 2(t^2+1)/(t-1)·1^T',
      all(zero(e) for e in (A.subs({X: 1, Z: 1})*one - (t + 1)**2/(t - 1)*one)) and
      all(zero(e) for e in (one.T*B.subs({X: 1, Z: 1}) - 2*(t**2 + 1)/(t - 1)*one.T)))
# nämnarna i A, B, adj A, adj B innehåller bara t, t-1, t+1, t^2+1
okden = True
for M in (A, B, A.adjugate(), B.adjugate()):
    for e in M:
        den = sp.fraction(sp.factor(sp.cancel(e)))[1]
        for f, _ in sp.factor_list(sp.Poly(den, X, Z, t).as_expr())[1]:
            if f.free_symbols - {t} and f not in (X, Z):
                okden = False
            if f.free_symbols == {t} and sp.factor(f) not in (t, t - 1, t + 1, t**2 + 1):
                okden = False
check('elementen i A, B, adj A, adj B är Laurentpolynom i X, Z med nämnare bara ur t, t±1, t²+1', okden)

# ---------------------------------------------------------------- (T2) ordlistan
EDGES = [(1, 2, (-1, 0)), (1, 2, (0, 0)), (1, 4, (0, 0)), (1, 4, (0, 1)), (2, 3, (0, 0)),
         (2, 3, (0, 1)), (3, 4, (0, 0)), (3, 4, (1, 0)), (1, 3, (0, 1)), (2, 4, (1, 0))]
P = {1: (0, 0), 2: (1, 0), 3: (1, -1), 4: (0, -1)}
img = lambda j, n: (P[j][0] + 2*n[0], P[j][1] + 2*n[1])
diag = lambda p, q: abs(p[0] - q[0]) == 1 and abs(p[1] - q[1]) == 1
model = set()
for x in range(-2, 2):
    for y in range(-2, 2):
        model.add(frozenset({(x, y), (x + 1, y)})); model.add(frozenset({(x, y), (x, y + 1)}))
        if x % 2 == 0 and y % 2 == 0: model.add(frozenset({(x, y), (x + 1, y + 1)}))
        if x % 2 == 1 and y % 2 == 1: model.add(frozenset({(x + 1, y), (x, y + 1)}))
good = True
for (i, j, n) in EDGES:
    e = frozenset({img(i, (0, 0)), img(j, n)})
    good &= e in model
    good &= diag(*tuple(e)) == ((i, j) in ((1, 3), (2, 4)))
# antal kanter per 2x2-block i modellen = 10 och grad 5
deg = {}
for e in model:
    for v in e: deg[v] = deg.get(v, 0) + 1
check('(j, n) -> P_j + 2n skickar de 10 kanterna på kanter i Z²-modellen; TT-kanterna blir diagonalerna', good)
check('Z²-modellen: alla inre noder har grad 5', all(deg[(x, y)] == 5 for x in range(-1, 1) for y in range(-1, 1)))

# ---------------------------------------------------------------- (T3) isoradiellt Z²_t
for tv in (sp.Integer(2), sp.Rational(5, 7), 2 + sp.sqrt(3)):
    Phi = 2*sp.atan(tv) - sp.pi/2
    al = {0: 0, 1: -Phi}; be = {0: sp.pi/2 - Phi, 1: sp.pi/2}
    good = True
    for x in range(2):
        for y in range(2):
            par = (x + y) % 2
            # horisontell kant (x,y)-(x+1,y): rombkanter e^{iα_par}, e^{iβ_par}
            th_h = (be[par] - al[par])/2
            # vertikal kant (x,y)-(x,y+1): -e^{iα_{1-par}} = e^{i(α+π)}, e^{iβ_par}
            th_v = ((al[1 - par] + sp.pi) - be[par])/2
            c_h = 1/tv if par == 0 else tv
            good &= sp.simplify(sp.tan(th_h) - c_h) == 0 and sp.simplify(th_v - sp.pi/4) == 0
    # vinkelsumma runt en nod: Σ 2θ = 2π; runt en yta: Σ (π - 2θ) = 2π
    thE, thO = (be[0] - al[0])/2, (be[1] - al[1])/2
    good &= sp.simplify(2*(thE + thO + sp.pi/4 + sp.pi/4) - 2*sp.pi) == 0
    good &= sp.simplify((sp.pi - 2*thE) + (sp.pi - 2*thO) + 2*(sp.pi/2) - 2*sp.pi) == 0
    good &= all(0 < float(b_ - a_) < float(sp.pi) for a_ in al.values() for b_ in be.values())
    check(f'Z²_t isoradiellt, t = {tv}: halvvinklar och konduktanser 1/t, t, 1; vinkelsummor 2π', good)

# ---------------------------------------------------------------- (T4) Bilaga E
lam = sp.Symbol('lambda', positive=True)
a, b = sp.symbols('a b', real=True)
sg, tt, tau = sp.symbols('sigma t tau', positive=True)
sys.path.insert(0, HERE)
import snub_viktad as sv
Dl = a**2 + b**2 - (4*lam**2 + 8*lam + 2)*a*b - 4*(lam + 2)**2*(a + b) + 4*(3*lam + 4)*(lam + 2)
NTT = 4*(3*lam + 5 - (lam + 2)*(a + b) - (lam + 1)*a*b)
NTS = (a**2 - (lam**2 + 3*lam + 1)*a*b - (lam**2 + 8*lam + 11)*a - (lam**2 + 4*lam + 5)*b
       + 3*lam**2 + 15*lam + 16)
check('E: medelvärdesbildade täljare N_TT, N_TS (kant (1,3;0,1) resp. (1,2;0,0)) och det M_λ = 4D_λ',
      zero(sv.to_ab(sv.numerator(sv.TT_EDGE))/4 - NTT.subs(lam, sv.lam)) and
      zero(sv.to_ab(sv.numerator(sv.TS_EDGE))/4 - NTS.subs(lam, sv.lam)) and
      zero(sp.expand(sv.bloch().det()) - 4*Dl.subs(lam, sv.lam).subs({a: (sv.X + 1/sv.X)/2, b: (sv.Z + 1/sv.Z)/2})))
L_ = sg**2 - 2
Aa = 4/(tt**2 - L_) - 1
rho = 4*tt/(tt**2 - L_)
Lam = tt**2 - L_
good = True
pp, qq = 2*L_**2 + 4*L_ + 1, 2*(L_ + 2)**2
for s_ in (1, -1):
    bs = pp*Aa + qq + s_*2*(L_ + 1)*sg*rho
    bo = pp*Aa + qq - s_*2*(L_ + 1)*sg*rho
    Ls = (2*sg**2 - 1)*tt + s_*sg*(2*sg**2 - 3)
    good &= zero(Dl.subs(lam, L_).subs({a: Aa, b: bs}))
    good &= zero(bs + 1 - 4*(sg**2 + s_*sg*tt - 1)**2/Lam) and zero(bs - 1 - 2*(tt + s_*sg)*Ls/Lam)
    cTT = NTT.subs(lam, L_).subs({a: Aa, b: bs})/(bo - bs)
    cTS = NTS.subs(lam, L_).subs({a: Aa, b: bs})/(bo - bs)
    # c_± = 2(t±σ)^2 |σ²±σt-1|/(σ t Λ)  och  c^TS_± = (t±σ)|σ²±σt-1| P_±/(2σtΛ), P_± = (λ+4)t ± λσ
    Ps = (L_ + 4)*tt + s_*L_*sg
    good &= zero(cTT - s_*2*(tt + s_*sg)**2*(sg**2 + s_*sg*tt - 1)/(sg*tt*Lam))
    good &= zero(cTS - s_*(tt + s_*sg)*(sg**2 + s_*sg*tt - 1)*Ps/(2*sg*tt*Lam))
check('E: rötterna b_±, faktoriseringarna av b_± ± 1, koefficienterna c_± (TT och TS)', good)
check('E: ρ² = (a+1)(λa+λ+4), 1-a = 2(t²-σ²)/Λ, 1+a = 4/Λ, diskriminanten',
      zero(rho**2 - (Aa + 1)*(L_*Aa + L_ + 4)) and zero(1 - Aa - 2*(tt**2 - sg**2)/Lam) and zero(1 + Aa - 4/Lam) and
      zero(((pp*a + qq)**2 - (a**2 - 2*qq*a + 4*(3*L_ + 4)*(L_ + 2))) - 4*(L_ + 1)**2*(L_ + 2)*(a + 1)*(L_*a + L_ + 4)))
# Euler: τ² = L_+/(t-σ) och τ² = L_-/(t+σ)
kk = 2*sg**2 - 1
Q = tau**4 + (4*(sg**2 - 1)**2 - 2)*tau**2 + 1
good = True
for s_ in (1, -1):
    tv = sg*(tau**2 + 2*sg**2 - 3)/(tau**2 - kk) if s_ > 0 else sg*(tau**2 + 2*sg**2 - 3)/(kk - tau**2)
    good &= zero((tv - s_*sg)*tau**2 - ((2*sg**2 - 1)*tv + s_*sg*(2*sg**2 - 3)))
    # TT-termen: (2/σ)(t±σ)/sqrt((t∓σ)L_±) dt/Λ = (2/σ)(t±σ)/((t∓σ)τ) dt/Λ
    gTT = (2/sg)*(tv + s_*sg)/((tv - s_*sg)*tau)/(tv**2 - L_)*sp.diff(tv, tau)
    gTS = (1/(2*sg))*((L_ + 4)*tv + s_*L_*sg)/((tv - s_*sg)*tau)/(tv**2 - L_)*sp.diff(tv, tau)
    sgn = -1 if s_ > 0 else 1        # τ går ∞ -> √k (term +) resp. 1 -> √k (term -)
    good &= zero(sgn*gTT - 4*(tau**2 - 1)/Q) and zero(sgn*gTS - (sg**2*tau**2 + 3*sg**2 - 4)/Q)
check('E: Euler-substitutionerna ger (1/π)·4(τ²-1)/Q och (1/π)(σ²τ²+3σ²-4)/Q för båda termerna', good)
check('E: Q/τ² = (τ+1/τ)² + 4κ² = (τ-1/τ)² + 4(λ+1)²; σ²τ²+3σ²-4 = -λ(τ²-1) + 2(λ+1)(τ²+1)',
      zero(Q/tau**2 - (tau + 1/tau)**2 - 4*L_*(L_ + 2)) and zero(Q/tau**2 - (tau - 1/tau)**2 - 4*(L_ + 1)**2) and
      zero(sg**2*tau**2 + 3*sg**2 - 4 + L_*(tau**2 - 1) - 2*(L_ + 1)*(tau**2 + 1)))
check('E: för t > σ är σ²+σt-1 > 0, σ²-σt-1 < 0 och (2σ²-1)t - σ(2σ²-3) > 2σ > 0 (prov i 200 punkter)',
      all((s2**2 + s2*x - 1 > 0) and (s2**2 - s2*x - 1 < 0) and ((2*s2**2 - 1)*x - s2*(2*s2**2 - 3) > 2*s2 - 1e-12)
          for s2 in [mp.sqrt(2 + l) for l in (0.01, 0.3, 1, 3, 30)] for x in [s2*(1 + k/8) for k in range(1, 41)]))

# ---------------------------------------------------------------- (T5) Tabell 3
L = sp.Symbol('lambda')
T3 = {  # klass: (j, n1, n2) i 0-baserad numrering som i runs/  ->  (a, b, c)
    (1, 0, 0): (sp.Rational(1, 2), 0, -1/(2*(L + 2))),
    (2, 0, 1): (0, 0, 2/(L*(L + 2))),
    (2, 0, 0): (0, 1, 1/(L*(L + 2))),
    (2, -1, 0): (-L, 0, 2*(L + 1)**2/(L*(L + 2))),
    (1, -1, 1): (-sp.Rational(1, 2), 0, (3*L + 8)/(2*L*(L + 2))),
    (0, -1, 0): (L + 2, -1, -(2*L**2 + 6*L + 3)/(L*(L + 2))),
    (1, -1, -1): (-(2*L**2 + 4*L + 1)/2, 2*(L + 2), (4*L**2 + 12*L + 7)/(2*(L + 2))),
    (2, 1, 1): (-2*(L + 2), (5*L + 8)/L, (4*L**3 + 12*L**2 + 5*L - 8)/(L**2*(L + 2))),
    (0, -1, -1): (L, -4*(L + 2)/L, -2*(L**2 - 2*L - 2)/L**2),
    (0, -1, 1): (0, -2*(L + 1)/L, 2*(L**2 + 3*L + 1)/(L**2*(L + 2))),
    (2, 1, 0): (-2*(L + 2), 2*(L + 1)/L, 2*(2*L**3 + 7*L**2 + 6*L - 1)/(L**2*(L + 2))),
    (1, -2, 0): ((6*L**2 + 20*L + 17)/2, -6*(L + 2), -(12*L**3 + 52*L**2 + 69*L + 24)/(2*L*(L + 2))),
}
r3 = sp.sqrt(3)
T1 = {(1, 0, 0): 1/sp.Integer(2) - r3/18, (2, 0, 1): 2*r3/9, (2, 0, 0): r3/9 + 1/sp.pi, (2, -1, 0): 8*r3/9 - 1,
      (1, -1, 1): 11*r3/18 - sp.Rational(1, 2), (0, -1, 0): 3 - 11*r3/9 - 1/sp.pi,
      (1, -1, -1): 23*r3/18 - sp.Rational(7, 2) + 6/sp.pi, (2, 1, 1): 13*r3/9 - 6 + 13/sp.pi,
      (0, -1, -1): 1 + 2*r3 - 12/sp.pi, (0, -1, 1): 10*r3/9 - 4/sp.pi, (2, 1, 0): 28*r3/9 - 6 + 4/sp.pi,
      (1, -2, 0): sp.Rational(43, 2) - 157*r3/18 - 18/sp.pi}
check('Tabell 3 i λ = 1 (κΦ/π = √3/3) ger Tabell 1, alla 12 rader',
      all(sp.simplify(sp.sympify(T3[k][0]).subs(L, 1) + sp.sympify(T3[k][1]).subs(L, 1)/sp.pi
                      + sp.sympify(T3[k][2]).subs(L, 1)*r3/3 - T1[k]) == 0 for k in T3))
cp = os.path.join(ROOT, 'runs', 'isoradiell_lambda', 'checkpoint.json')
if os.path.exists(cp):
    mp.mp.dps = 50
    data = json.load(open(cp, encoding='utf-8'))['results']
    n_ok = n_tot = 0
    for key, v in data.items():          # nyckel 'λ|dps|j:n1:n2'
        ls, _, node = key.split('|')
        k = tuple(int(x) for x in node.split(':'))
        if k not in T3:
            continue
        lv = sp.Rational(ls)
        val = mp.mpf(v['value'])
        if True:
            kap = mp.sqrt(lv*(lv + 2)); Ph = mp.acos(1/(1 + mp.mpf(lv)))
            f = (mp.mpf(sp.sympify(T3[k][0]).subs(L, lv)) + mp.mpf(sp.sympify(T3[k][1]).subs(L, lv))/mp.pi
                 + mp.mpf(sp.sympify(T3[k][2]).subs(L, lv))*kap*Ph/mp.pi)
            n_tot += 1; n_ok += abs(f - val) < mp.mpf(10)**-40
    check(f'Tabell 3 mot de numeriska värdena för λ = 1/2, 2, 3 i runs/isoradiell_lambda ({n_ok} av {n_tot})',
          n_tot > 0 and n_ok == n_tot)
else:
    print('(T5b hoppas över: runs/isoradiell_lambda/checkpoint.json saknas)')

# ---------------------------------------------------------------- (T6) Kenyons log det
mp.mp.dps = 40
Lob = lambda x: mp.clsin(2, 2*x)/2


def kenyon_m(lv):
    Ph = mp.acos(1/(1 + lv))
    f = lambda th: (2/mp.pi)*(th*mp.log(mp.tan(th)) + Lob(th) + Lob(mp.pi/2 - th))
    return 4*f(mp.pi/4) + 2*f(mp.pi/4 + Ph/2) + 2*f(mp.pi/4 - Ph/2)


def RTT(lv):
    kap = mp.sqrt(lv*(lv + 2)); return 2*mp.atan(kap)/(mp.pi*kap)


check('Kenyons formel för Z²_t: dK/dλ = 2R_TT(λ) i λ = 0.3, 1, 4 och K(0+) = 16G/π',
      all(abs(mp.diff(kenyon_m, mp.mpf(lv)) - 2*RTT(mp.mpf(lv))) < mp.mpf(10)**-25 for lv in (0.3, 1, 4)) and
      abs(kenyon_m(mp.mpf(10)**-30) - 16*mp.catalan/mp.pi) < mp.mpf(10)**-12)
check('Kenyons formel = Clausenformen i sats F (λ = 1/4, 1, 3, 10)',
      all(abs(kenyon_m(mp.mpf(lv)) - (8*mp.catalan/mp.pi + 4/mp.pi*(mp.acos(1/(1 + mp.mpf(lv)))*mp.log(1 + lv + mp.sqrt(lv*(lv + 2)))
          + mp.clsin(2, mp.pi/2 + mp.acos(1/(1 + mp.mpf(lv)))) + mp.clsin(2, mp.pi/2 - mp.acos(1/(1 + mp.mpf(lv))))))) < mp.mpf(10)**-30
          for lv in (0.25, 1, 3, 10)))

# ---------------------------------------------------------------- (T7) R_TT via Kenyon c·R = 2θ/π och Jacobi
# m(t) = Σ_e (Kenyon) ; dm/dt = Σ_e R_e dc_e/dt = 2(2θ_-/π)(1/c_-)(-1/t²) + 2(2θ_+/π)(1/c_+)(1)
good = True
for lv in (mp.mpf('0.2'), mp.mpf(1), mp.mpf(7)):
    kap = mp.sqrt(lv*(lv + 2)); tv = 1 + lv + kap; Ph = mp.atan(kap)
    thp, thm = mp.pi/4 + Ph/2, mp.pi/4 - Ph/2
    dm_dt = 2*(2*thm/mp.pi)*tv*(-1/tv**2) + 2*(2*thp/mp.pi)/tv
    dt_dl = tv/kap
    good &= abs(dm_dt*dt_dl - 2*RTT(lv)) < mp.mpf(10)**-30 and abs(4*Ph/(mp.pi*tv)*dt_dl - 2*RTT(lv)) < mp.mpf(10)**-30
check('Anmärkning: Jacobi + Kenyons c·R = 2θ/π på Z²_t ger 2R_TT = 4Φ/(πκ)', good)

# ---------------------------------------------------------------- (T8) gradgränsen för Cairo (granskning 2, punkt 4.2)
# Gradgränsen w = 4·(bredd i X + bredd i Z) i residue_tools.width() gäller när varje exponentlåda innehåller 0.
try:
    import exact_snub as es, cairo_lattice as cl, residue_tools as rt
    pos_, deg_, edges_ = cl.build()
    Lc = cl.bloch_laplacian(deg_, edges_)
    ADJc = Lc.adjugate().applyfunc(sp.expand)
    ld = lambda e: rt.laurent_dict(e, es.X, es.Z)
    boxes = [ld(ADJc[i, j]) for i in range(6) for j in range(6)] + [ld(sp.diff(es.DET, es.Z))]
    inne = lambda d: (min(a for a, _ in d) <= 0 <= max(a for a, _ in d)) and (min(b for _, b in d) <= 0 <= max(b for _, b in d))
    check('Cairo: alla 36 exponentlådor i adj L och lådan för ∂_Z det innehåller 0 (gradgränsen 86 i C2 är giltig)',
          all(inne(d) for d in boxes) and max(rt.width(d) for d in boxes[:-1]) + rt.width(boxes[-1]) + 38 == 86)
except ImportError as e:
    print('(T8 hoppas över: ' + str(e) + ')')

print(f'\n{sum(ok)} godkända, {len(ok) - sum(ok)} fel')
