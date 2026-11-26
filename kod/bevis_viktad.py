"""
Kontroll av härledningen för snub square-gittret med två konduktanser (kandidat A4):
konduktans c på triangel–triangel-kanterna, d = 1 på triangel–kvadrat-kanterna, λ = c/d > 0.

Resultat (κ = sqrt(λ(λ+2)), Φ = arccos(1/(1+λ)) = arctan κ):
    R_TT(λ) = 2 arctan(κ) / (π κ)
    R_TS(λ) = 1/2 - λ arctan(κ) / (2 π κ)
    m(λ)    = 8G/π + (4/π) [ Φ log(1+λ+κ) + Cl2(π/2+Φ) + Cl2(π/2-Φ) ]      (m = ∫∫ log det M)

Varje steg kontrolleras exakt med sympy för allmänt λ (via σ = sqrt(λ+2)), eller numeriskt till
minst 50 siffror. Körtid: under en minut. Kör:  .venv\\Scripts\\python kod\\bevis_viktad.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sympy as sp
import mpmath as mp

import snub_viktad as sv

mp.mp.dps = 70
TOL = mp.mpf(10)**-55
ok = []


def check(name, cond):
    ok.append(bool(cond))
    print(('OK  ' if cond else 'FEL ') + name, flush=True)


def zero(e):
    return sp.simplify(sp.factor(sp.together(sp.expand(e)))) == 0


X, Z, lam, a, b = sv.X, sv.Z, sv.lam, sv.a, sv.b
t, sg, tau, w = sp.symbols('t sigma tau w', positive=True)

# ---- Steg 1: determinanten ----------------------------------------------------------------
p = 2*lam**2 + 4*lam + 1
q = 2*(lam + 2)**2
r = 4*(3*lam + 4)*(lam + 2)
D = b**2 - 2*b*(p*a + q) + a**2 - 2*q*a + r       # = a²+b²-(4λ²+8λ+2)ab-4(λ+2)²(a+b)+4(3λ+4)(λ+2)
det = sp.expand(sv.bloch().det())
check('det M = 4 D_λ(a, b) exakt (utan medelvärdesbildning)',
      zero(det - 4*D.subs({a: (X + 1/X)/2, b: (Z + 1/Z)/2})))
D1 = a**2 + b**2 - 14*a*b - 36*(a + b) + 84
check('λ = 1 ger artikelns D(a, b) (Lemma 2)', zero(D.subs(lam, 1) - D1))
check('D_λ(1, 1) = 0 (kurvan går genom (1,1))', zero(D.subs({a: 1, b: 1})))
check('D_λ(-1, b) = (b - 4λ - 7)² (tangent till a = -1)', zero(D.subs(a, -1) - (b - 4*lam - 7)**2))
check('D_λ(1, b) = (b - 1)(b - 8λ² - 24λ - 17)', zero(D.subs(a, 1) - (b - 1)*(b - 8*lam**2 - 24*lam - 17)))
check('diskriminant: (pa+q)² - (a² - 2qa + r) = 4(λ+1)²(λ+2)(a+1)(λa+λ+4)',
      zero((p*a + q)**2 - (a**2 - 2*q*a + r) - 4*(lam + 1)**2*(lam + 2)*(a + 1)*(lam*a + lam + 4)))

# ---- Steg 2: symmetriserade täljare --------------------------------------------------------
NTT = sv.to_ab(sv.numerator(sv.TT_EDGE)) / 4          # integrand = NTT / D
NTS = sv.to_ab(sv.numerator(sv.TS_EDGE)) / 4
check('täljare TT: 4(3λ+5 - (λ+2)(a+b) - (λ+1)ab)', zero(NTT - 4*(3*lam + 5 - (lam + 2)*(a + b) - (lam + 1)*a*b)))
check('λ = 1 ger artikelns 4(8 - 3a - 3b - 2ab) (Lemma 3)', zero(NTT.subs(lam, 1) - 4*(8 - 3*a - 3*b - 2*a*b)))
check('täljarna har grad 1 i b (inget polynomdel vid division med D)',
      sp.degree(NTT, b) == 1 and sp.degree(NTS, b) == 1)

# ---- Steg 3: rationell parameter t, σ = sqrt(λ+2) -----------------------------------------
L_ = sg**2 - 2                                        # λ uttryckt i σ
sub = lambda e: e.subs(lam, L_)
A = 4/(t**2 - L_) - 1
rho = 4*t/(t**2 - L_)
check('ρ² = (a+1)(λa+λ+4) med a = 4/(t²-λ) - 1, ρ = 4t/(t²-λ)', zero(rho**2 - (A + 1)*(L_*A + L_ + 4)))
k = 2*sg**2 - 1
fac = {}
for s_ in (1, -1):
    bs = sub(p)*A + sub(q) + s_*2*(L_ + 1)*sg*rho
    bo = sub(p)*A + sub(q) - s_*2*(L_ + 1)*sg*rho
    check(f'D_λ(a, b_{"+" if s_ > 0 else "-"}) = 0', zero(sub(D).subs({a: A, b: bs})))
    sq = (sg**2 + s_*sg*t - 1)
    Lpm = k*t + s_*sg*(2*sg**2 - 3)
    check(f'b{"+" if s_ > 0 else "-"} + 1 = 4(σ² {"+" if s_ > 0 else "-"} σt - 1)²/(t²-λ)',
          zero(bs + 1 - 4*sq**2/(t**2 - L_)))
    check(f'b{"+" if s_ > 0 else "-"} - 1 = 2(t {"+" if s_ > 0 else "-"} σ)((2σ²-1)t {"+" if s_ > 0 else "-"} σ(2σ²-3))/(t²-λ)',
          zero(bs - 1 - 2*(t + s_*sg)*Lpm/(t**2 - L_)))
    for nm, N in (('TT', NTT), ('TS', NTS)):
        fac[(nm, s_)] = sp.factor(sp.cancel(sub(N).subs({a: A, b: bs}) / (bo - bs)))
check('1 - a = 2(t-σ)(t+σ)/(t²-λ), 1 + a = 4/(t²-λ)',
      zero(1 - A - 2*(t - sg)*(t + sg)/(t**2 - L_)) and zero(1 + A - 4/(t**2 - L_)))

# ---- Steg 4: Euler-substitution; integranden blir rationell i τ -----------------------------
# Inre integralen: (1/π)∫ c/(β-b) db/sqrt(1-b²) = c/sqrt(β²-1). Med faktoriseringarna ovan och
# |da|/sqrt(1-a²) = 2√2 t dt/((t²-λ) sqrt(t²-σ²)) blir termen för tecknet ±:
#   c± (t²-λ) / (2√2 (σ² ± σt - 1) sqrt((t±σ) L±)) * 2√2 t / ((t²-λ) sqrt(t²-σ²)) dt
# och τ² = L±/(t∓σ) gör sqrt((t∓σ) L±) = (t∓σ) τ.
Q = tau**4 + (4*(sg**2 - 1)**2 - 2)*tau**2 + 1
res = {}
for nm in ('TT', 'TS'):
    tot = 0
    for s_ in (1, -1):
        c = fac[(nm, s_)]
        sq = (sg**2 + s_*sg*t - 1)
        # integrand i t (utan rottecken): c * t / (sq * sqrt((t+s_σ)L) * sqrt(t²-σ²))
        #   = c t / (sq (t + s_σ) * (t - s_σ) τ)   eftersom sqrt((t+s_σ)L)·sqrt((t-s_σ)(t+s_σ)) = (t+s_σ)(t-s_σ)τ
        tt = sg*(tau**2 + 2*sg**2 - 3)/(tau**2 - k) if s_ > 0 else sg*(tau**2 + 2*sg**2 - 3)/(k - tau**2)
        check(f'τ-substitution ({nm}, {"+" if s_ > 0 else "-"}): (t ∓ σ)τ² = L±',
              zero((tt - s_*sg)*tau**2 - (k*tt + s_*sg*(2*sg**2 - 3))))
        # sqrt((b+1)) = 2|σ² ± σt - 1|/sqrt(t²-λ); för t > σ är σ² + σt - 1 > 0 och σ² - σt - 1 < 0
        g = c.subs(t, tt) * tt / (s_*sq.subs(t, tt) * (tt + s_*sg) * (tt - s_*sg) * tau) * sp.diff(tt, tau)
        g = sp.factor(sp.cancel(g))
        # term +: t går σ→∞ när τ går ∞→√k (byt tecken); term -: τ går 1→√k
        tot += (-g if s_ > 0 else g)
        res[(nm, s_)] = (-g if s_ > 0 else g)
    check(f'{nm}: de två termerna ger samma rationella funktion av τ (ihopfogas till ∫_1^∞)',
          zero(res[(nm, 1)] - res[(nm, -1)]))
gTT, gTS = res[('TT', 1)], res[('TS', 1)]
check('TT: integrand (1/π)·4(τ²-1)/Q(τ), Q = τ⁴ + (4(σ²-1)²-2)τ² + 1', zero(gTT - 4*(tau**2 - 1)/Q))
check('TS: integrand (1/π)(σ²τ² + 3σ² - 4)/Q(τ)', zero(gTS - (sg**2*tau**2 + 3*sg**2 - 4)/Q))
check('λ = 1: Q = τ⁴ + 14τ² + 1 (artikeln, Steg 4)', zero(Q.subs(sg, sp.sqrt(3)) - (tau**4 + 14*tau**2 + 1)))

# ---- Steg 5: de elementära integralerna ---------------------------------------------------
kap2 = L_*(L_ + 2)                                         # κ² = λ(λ+2) = σ²(σ²-2)
check('w = τ + 1/τ: Q/τ² = w² + 4λ(λ+2)', zero(Q/tau**2 - ((tau + 1/tau)**2 + 4*kap2)))
check('v = τ - 1/τ: Q/τ² = v² + 4(λ+1)²', zero(Q/tau**2 - ((tau - 1/tau)**2 + 4*(L_ + 1)**2)))
check('σ²τ² + 3σ² - 4 = -λ(τ² - 1) + 2(λ+1)(τ² + 1)',
      zero(sg**2*tau**2 + 3*sg**2 - 4 - (-L_*(tau**2 - 1) + 2*(L_ + 1)*(tau**2 + 1))))
# ∫_1^∞ (τ²-1)/Q dτ = ∫_2^∞ dw/(w² + 4κ²) = arctan(κ)/(2κ) ;  ∫_1^∞ (τ²+1)/Q dτ = ∫_0^∞ dv/(v² + 4(λ+1)²) = π/(4(λ+1))
kp = sp.Symbol('kappa', positive=True)
I1 = sp.integrate(1/(w**2 + 4*kp**2), (w, 2, sp.oo))
check('∫_2^∞ dw/(w² + 4κ²) = arctan(κ)/(2κ)', zero(sp.simplify(I1 - sp.atan(kp)/(2*kp))))


def R_TT(lv):
    kv = mp.sqrt(lv*(lv + 2))
    return 2*mp.atan(kv)/(mp.pi*kv)


def R_TS(lv):
    kv = mp.sqrt(lv*(lv + 2))
    return mp.mpf(1)/2 - lv*mp.atan(kv)/(2*mp.pi*kv)


def m_closed(lv):
    kv = mp.sqrt(lv*(lv + 2))
    F = mp.acos(1/(1 + lv))
    cl2 = lambda x: mp.clsin(2, x)
    return 8*mp.catalan/mp.pi + 4/mp.pi*(F*mp.log(1 + lv + kv) + cl2(mp.pi/2 + F) + cl2(mp.pi/2 - F))


# kontroll av integralerna i τ numeriskt för några λ
for lv in (mp.mpf(1)/4, mp.mpf(1), mp.mpf(3)):
    s2 = lv + 2
    Qf = lambda x: x**4 + (4*(s2 - 1)**2 - 2)*x**2 + 1
    vTT = mp.quad(lambda x: 4*(x**2 - 1)/Qf(x), [1, 3, mp.inf]) / mp.pi
    vTS = mp.quad(lambda x: (s2*x**2 + 3*s2 - 4)/Qf(x), [1, 3, mp.inf]) / mp.pi
    check(f'λ = {mp.nstr(lv, 3)}: τ-integralerna = slutna formlerna (R_TT, R_TS)',
          abs(vTT - R_TT(lv)) < TOL and abs(vTS - R_TS(lv)) < TOL)

# ---- Steg 6: kontroller av formlerna -------------------------------------------------------
check('λ = 1: R_TT = 2√3/9', abs(R_TT(mp.mpf(1)) - 2*mp.sqrt(3)/9) < TOL)
check('λ = 1: R_TS = 1/2 - √3/18', abs(R_TS(mp.mpf(1)) - (mp.mpf(1)/2 - mp.sqrt(3)/18)) < TOL)
check('Foster 2λR_TT + 8R_TS = 4 (exakt i λ)',
      zero(2*lam*2*sp.atan(sp.sqrt(lam*(lam + 2)))/(sp.pi*sp.sqrt(lam*(lam + 2)))
           + 8*(sp.Rational(1, 2) - lam*sp.atan(sp.sqrt(lam*(lam + 2)))/(2*sp.pi*sp.sqrt(lam*(lam + 2)))) - 4))
check('λ → 0: R_TT → 2/π (diagonal i en kvadrat i kvadratgittret), R_TS → 1/2 (kvadratgittrets kant)',
      abs(R_TT(mp.mpf(10)**-80) - 2/mp.pi) < TOL and abs(R_TS(mp.mpf(10)**-80) - mp.mpf(1)/2) < TOL)
check('λ → ∞: λ·R_TT → 1, R_TS → 1/4 (kvadratgittret med dubbla kanter)',
      abs(mp.mpf(10)**60*R_TT(mp.mpf(10)**60) - 1) < 1e-50 and abs(R_TS(mp.mpf(10)**60) - mp.mpf(1)/4) < 1e-50)
check('λ = √2 - 1: R_TT = 1/2 (Φ = π/4)', abs(R_TT(mp.sqrt(2) - 1) - mp.mpf(1)/2) < TOL)

# ---- Steg 7: spännträdskonstanten m(λ) = ∫∫ log det M -------------------------------------
# d m/dλ = Σ_{TT-kanter per cell} R_e = 2 R_TT (Jacobis formel: d log det M = tr(M⁻¹ dM)).
F_ = sp.Symbol('Phi', positive=True)
# med λ = sec Φ - 1: κ = tan Φ, dλ = sec Φ tan Φ dΦ; d/dΦ[Cl2(π/2+Φ) + Cl2(π/2-Φ)] = -log tan(π/4+Φ/2)
dm_dPhi = 4/sp.pi*(sp.log(1/sp.cos(F_) + sp.tan(F_)) + F_/sp.cos(F_) - sp.log(sp.tan(sp.pi/4 + F_/2)))
check('sec Φ + tan Φ = tan(π/4 + Φ/2)', abs(mp.sec(0.7) + mp.tan(0.7) - mp.tan(mp.pi/4 + 0.35)) < TOL)
check('dm/dλ = 2 R_TT symboliskt (via Φ)',
      zero(sp.simplify(dm_dPhi.subs(sp.log(sp.tan(sp.pi/4 + F_/2)), sp.log(1/sp.cos(F_) + sp.tan(F_)))
                       / (sp.sin(F_)/sp.cos(F_)**2) - 2*2*F_/(sp.pi*sp.tan(F_)))))
for lv in (mp.mpf(1)/2, mp.mpf(2)):
    check(f'dm/dλ = 2R_TT numeriskt vid λ = {mp.nstr(lv, 3)}', abs(mp.diff(m_closed, lv) - 2*R_TT(lv)) < mp.mpf(10)**-40)
check('λ = 0: m = 16G/π (kvadratgittret, 4 noder per cell)', abs(m_closed(mp.mpf(0)) - 16*mp.catalan/mp.pi) < TOL)
check('λ = 1: m = (40/3)G/π + (4/3)log(2+√3) (bevisat tidigare)',
      abs(m_closed(mp.mpf(1)) - (40*mp.catalan/(3*mp.pi) + mp.mpf(4)/3*mp.log(2 + mp.sqrt(3)))) < TOL)
check('λ → ∞: m - 2 log λ → 8G/π + 2 log 2 (kvadratgittret med dubbla kanter, 2 noder per cell)',
      abs(m_closed(mp.mpf(10)**30) - 2*mp.log(mp.mpf(10)**30) - (8*mp.catalan/mp.pi + 2*mp.log(2))) < mp.mpf(10)**-25)

print()
print(f'{sum(ok)} godkända, {len(ok) - sum(ok)} fel')
