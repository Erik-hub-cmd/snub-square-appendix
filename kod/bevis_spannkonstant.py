"""
Kontroll av beviset för  m(P) = (40/3) G/pi + (4/3) log(2+sqrt3)  (se forslag/bevis_spannkonstant.md).

Varje steg kontrolleras exakt (sympy) eller numeriskt till minst 50 siffror (mpmath).
Körtid: några sekunder. Kör:  .venv\\Scripts\\python kod\\bevis_spannkonstant.py
"""
import sympy as sp
import mpmath as mp

mp.mp.dps = 70
TOL = mp.mpf(10)**-55
ok = []


def check(name, cond):
    ok.append(bool(cond))
    print(('OK  ' if cond else 'FEL ') + name)


s = sp.symbols('s')
r2, r3, r6 = sp.sqrt(2), sp.sqrt(3), sp.sqrt(6)

# ---- Steg 1: parametrisering ------------------------------------------------------------
q = s - 1/s
mu = (1 + q**2) / (1 - q**2)
m = r3 * mu
a = sp.cancel((5 - m**2) / (m**2 - 1))
sint = sp.cancel(2 * r6 * (mu + 1) * q / (3 * mu**2 - 1))
b = sp.cancel((11 * m**2 + 17 + 16 * r3 * m) / (m**2 - 1))
w = sp.cancel(2 * r6 * (r3 * m + 2) * (mu + 1) * (s + 1/s) / (m**2 - 1))     # = sqrt(b^2 - 1)

Xd = s**8 - s**4 + 1
check('a = (s^8-12s^6+23s^4-12s^2+1)/(s^8-s^4+1)', sp.simplify(a - (s**8 - 12*s**6 + 23*s**4 - 12*s**2 + 1) / Xd) == 0)
check('a^2 + sin^2 = 1', sp.simplify(a**2 + sint**2 - 1) == 0)
check('b = 7a+18+4sqrt3*rho*sign(m), rho = 4|m|/(m^2-1)', sp.simplify(b - (7*a + 18 + 4*r3*4*m/(m**2 - 1))) == 0)
check('w^2 = b^2 - 1', sp.simplify(w**2 - (b**2 - 1)) == 0)
# D(a,b) = 0 (b är en rot till determinanten)
check('D(a,b) = 0', sp.simplify(a**2 + b**2 - 14*a*b - 36*(a + b) + 84) == 0)

# ---- Steg 2: faktorisering i 24:e enhetsrötter --------------------------------------------
z = (r6 + r2)/4 + sp.I*(r6 - r2)/4          # exp(i pi/12)
check('z^12 = -1', sp.simplify(sp.expand(z**12) + 1) == 0)


def prod(ks):
    e = sp.Integer(1)
    for k in ks:
        e = sp.expand(e * (s - sp.expand(z**k)))
    return e


def is_zero(expr):
    p = sp.Poly(sp.expand(expr), s)
    return all(sp.simplify(sp.radsimp(c)) == 0 for c in p.all_coeffs())


X = a + sp.I*sint
Z = b + w
eX = {1: 1, 5: 1, 7: 1, 11: 1, 13: -1, 17: -1, 19: -1, 23: -1}
eZ = {7: 1, 11: 1, 13: 1, 17: 1, 1: -1, 5: -1, 19: -1, 23: -1}
Xn, Xd2 = sp.fraction(sp.together(X))
Zn, Zd = sp.fraction(sp.together(Z))
check('X = prod_{k=1,5,7,11} (s-z^k)/(s-z^-k)',
      is_zero(Xn * prod([13, 17, 19, 23]) - Xd2 * prod([1, 5, 7, 11])))
check('Z = (s-z^7)(s-z^11)(s-z^13)(s-z^17)/((s-z)(s-z^5)(s-z^19)(s-z^23))',
      is_zero(Zn * prod([1, 5, 19, 23]) - Zd * prod([7, 11, 13, 17])))

# ---- numeriska funktioner ----------------------------------------------------------------
fa = sp.lambdify(s, a, 'mpmath')
fs = sp.lambdify(s, sint, 'mpmath')
fb = sp.lambdify(s, b, 'mpmath')
fZ = sp.lambdify(s, Z, 'mpmath')
phi = (1 + mp.sqrt(5)) / 2
# |Z| = exp(arccosh b) och b = b_+ på (1,phi), b_- på (phi,oo)
S3 = mp.sqrt(3)
err = 0
for sv in [mp.mpf('1.1'), mp.mpf('1.4'), mp.mpf('1.6'), mp.mpf('1.7'), mp.mpf(3), mp.mpf(20)]:
    av = fa(sv)
    rho = mp.sqrt((av + 1)*(av + 5))
    bt = 7*av + 18 + (1 if sv < phi else -1) * 4*S3*rho
    err = max(err, abs(fb(sv) - bt), abs(mp.log(abs(fZ(sv))) - mp.acosh(bt)))
check('b = b_+ för 1<s<phi, b = b_- för s>phi, log|Z| = arccosh b', err < TOL)

# ---- Steg 3: m(P) = (1/pi) int_1^oo log|Z| d arg X ------------------------------------------
def fA(t):
    av = mp.cos(t)
    rho = mp.sqrt((av + 1)*(av + 5))
    return mp.acosh(7*av + 18 + 4*S3*rho) + mp.acosh(7*av + 18 - 4*S3*rho)


mA = mp.quad(fA, [0, mp.pi/4, mp.pi/2, 3*mp.pi/4, mp.pi]) / mp.pi
dargX = sp.lambdify(s, sp.simplify(sp.diff(sp.atan2(sint, a), s)), 'mpmath')
mB = mp.quad(lambda t: mp.log(abs(fZ(t))) * dargX(t), [1, mp.mpf('1.3'), phi, 3, 10, mp.inf]) / mp.pi
check('m(P) = (1/pi) int_1^oo log|Z(s)| d arg X(s)  [%s]' % mp.nstr(abs(mA - mB), 3), abs(mA - mB) < TOL)


# ---- Steg 4: varje term int_1^oo eta(s-al, s-be) -----------------------------------------
def D(x):
    return mp.im(mp.polylog(2, x)) + mp.arg(1 - x) * mp.log(abs(x))


def cl2(th):
    return mp.clsin(2, th)


def zeta(k):
    return mp.expjpi(mp.mpf(k)/12)


worst = 0
worst_k = 0
for ka in eX:
    for kb in eZ:
        if ka == kb:
            continue
        al, be = zeta(ka), zeta(kb)
        f = lambda t: mp.log(abs(t - al)) * mp.im(1/(t - be)) - mp.log(abs(t - be)) * mp.im(1/(t - al))
        direct = mp.quad(f, [1, 2, 10, mp.inf])
        u0 = (1 - al)/(be - al)
        lemma = -D(u0) + mp.log(abs(be - al)) * (mp.arg(1 - al) - mp.arg(1 - be))
        worst = max(worst, abs(direct - lemma))
        # Kummer/Zagier: D(u0) = (1/2)[Cl2(th_a) - Cl2(th_b) + Cl2(th_b - th_a)]
        ta, tb = mp.pi*ka/12, mp.pi*kb/12
        worst_k = max(worst_k, abs(D(u0) - (cl2(ta) - cl2(tb) + cl2(tb - ta))/2))
check('int_1^oo eta(s-al,s-be) = -D(u0) + log|be-al| (arg(1-al)-arg(1-be)), alla 56 par', worst < TOL)
check('D(u0) = (Cl2(th_a) - Cl2(th_b) + Cl2(th_b-th_a))/2, alla 56 par', worst_k < TOL)

# ---- Steg 5: summering (exakt) -------------------------------------------------------------
from collections import Counter
W = Counter()
Lc = Counter()
for ka, ea in eX.items():
    for kb, eb in eZ.items():
        W[(kb - ka) % 24] += ea*eb
        if ka != kb:
            n = abs(kb - ka)//2
            n = min(n, 12 - n)                      # sin(n pi/12) = sin((12-n) pi/12)
            Lc[n] += ea*eb*sp.Rational(ka - kb, 24)
check('sum e_a = sum e_b = 0', sum(eX.values()) == 0 and sum(eZ.values()) == 0)
Wd = {k: v for k, v in W.items() if v}
check('vikter för Cl2(n pi/12): {2:4, 6:8, 10:4, 14:-4, 18:-8, 22:-4}', Wd == {2: 4, 6: 8, 10: 4, 14: -4, 18: -8, 22: -4})
Ld = {k: v for k, v in Lc.items() if v}
check('log-vikter: pi*[(4/3)log(2sin(pi/12)) - (4/3)log(2sin(5pi/12))]', Ld == {1: sp.Rational(4, 3), 5: sp.Rational(-4, 3)})
check('sin(pi/12)/sin(5pi/12) = 2 - sqrt3', sp.simplify(sp.sin(sp.pi/12)/sp.sin(5*sp.pi/12) - (2 - r3)) == 0)

# ---- Steg 6: Cl2(pi/6) + Cl2(5pi/6) = (4/3) G -------------------------------------------
# 2 sin(n pi/2) cos(n pi/3) för n = 1..12 ska vara chi_{-4}(n) * (1 om 3 ∤ n, -2 om 3 | n)
chk = all(sp.simplify(sp.sin(n*sp.pi/6) + sp.sin(5*n*sp.pi/6)
                      - (0 if n % 2 == 0 else (1 if n % 4 == 1 else -1)) * (1 if n % 3 else -2)) == 0
          for n in range(1, 13))
check('sin(n pi/6)+sin(5n pi/6) = chi_-4(n)*(1 eller -2), period 12', chk)
G = mp.catalan
check('Cl2(pi/6) + Cl2(5pi/6) = (4/3) G', abs(cl2(mp.pi/6) + cl2(5*mp.pi/6) - 4*G/3) < TOL)

# ---- Slutresultat --------------------------------------------------------------------------
C = -(8*(cl2(mp.pi/6) + cl2(5*mp.pi/6)) + 16*G)/2 - mp.mpf(4)/3*mp.pi*mp.log(2 + S3)
target = 40*G/(3*mp.pi) + mp.mpf(4)/3*mp.log(2 + S3)
check('m(P) = -C/pi = (40/3)G/pi + (4/3)log(2+sqrt3)  [%s]' % mp.nstr(abs(-C/mp.pi - mA), 3),
      abs(-C/mp.pi - target) < TOL and abs(mA - target) < TOL)

print()
print(f'{sum(ok)} godkända, {len(ok) - sum(ok)} fel')
