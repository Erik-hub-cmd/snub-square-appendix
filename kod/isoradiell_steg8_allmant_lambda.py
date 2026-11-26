"""
Isoradiell förklaring, steg 8: samma konstruktion för allmänt λ > 0.

Parametrisera med t = q = 1 + λ + √(λ(λ+2)) > 1, alltså λ = (t − 1)²/(2t) (och cos Φ = 1/(1 + λ)).
Då är alla vikter rationella i t, och det finns A(X, Z), B(X, Z) med koefficienter i Q(t) och
    M_snub(λ) · A = B · M_Z²(q),   det A = det B = 4 (t + 1)⁴ (t² + 1)² / (t − 1)⁴   (konstant i X, Z).
A och B har samma stöd (monom) som i steg 5 och fås ur ett linjärt ekvationssystem, normerat med A[1,1] = 1;
lösningen är då entydig. Alltså M_snub⁻¹ = A · M_Z²⁻¹ · B⁻¹ för varje λ.
Tillägg 2026-10-19 (granskning 4): entydigheten kontrolleras (full rang i t = 2 och t = 5/7, alltså också över Q(t)), och
konstanterna A(1,1)·1 = α·1, 1ᵀB(1,1) = β·1ᵀ som används i identifieringen med Kenyons G (anteckningen, steg 10.1).
För M_snub⁻¹ = A·M_Z²⁻¹·B⁻¹ räcker att det B är konstant.
Kör: python kod/isoradiell_steg8_allmant_lambda.py
"""
import sys
sys.stdout.reconfigure(encoding="utf-8")
import sympy as sp

X, Z, t = sp.symbols('X Z t')
ok = []


def check(name, cond):
    ok.append(bool(cond)); print(('OK  ' if cond else 'FEL ') + name, flush=True)


def lap(edges):
    idx = {(i, j): 2*i + j for i in range(2) for j in range(2)}
    M = sp.zeros(4, 4)
    for p, r, w in edges:
        a, b = idx[(p[0] % 2, p[1] % 2)], idx[(r[0] % 2, r[1] % 2)]
        sx = (r[0] - r[0] % 2)//2 - (p[0] - p[0] % 2)//2; sy = (r[1] - r[1] % 2)//2 - (p[1] - p[1] % 2)//2
        ph = X**sx * Z**sy
        M[a, a] += w; M[b, b] += w; M[a, b] -= w*ph; M[b, a] -= w/ph
    return M


def z2e(h, v):
    E = []
    for i in range(2):
        for j in range(2):
            E.append(((i, j), (i + 1, j), h[2*i + j])); E.append(((i, j), (i, j + 1), v[2*i + j]))
    return E


lam = (t - 1)**2/(2*t)
Ms = lap(z2e([1]*4, [1]*4) + [((0, 0), (1, 1), lam), ((2, 1), (1, 2), lam)])     # snub square, vikt λ
Mz = lap(z2e([1/t, t, t, 1/t], [1]*4))                                          # isoradiellt Z², q = t
A = sp.Matrix([
    [(X*t**5 + 3*X*t**4 + 4*X*t**3 + 4*X*t**2 + 3*X*t + X - t + 1)/(X*t**2*(t - 1)), -(t + 1)*(t**2 + 1)/(Z*t*(t - 1)), -(t + 1)*(X*t**2 + X + t**4 + 1)/(X*t**2*(t - 1)), -1/(X*Z*t)],
    [Z*(t + 1)*(t**2 + 1)/(t*(t - 1)), 1, Z/t, 0],
    [-(t + 1)*(X*t**2 + X + 2)/(t*(t - 1)), -X/Z, (-X*t + X + 2*t**3 + 6*t**2 + 6*t + 2)/(t*(t - 1)), -2*(t + 1)/(Z*(t - 1))],
    [Z, 0, 2*Z*(t + 1)/(t - 1), t]])
B = sp.Matrix([
    [(-2*X*t**3 + 2*X*t**2 + t**5 + 7*t**4 + 8*t**3 + 8*t**2 + 7*t + 1)/(2*t**2*(t - 1)), (3*t**3 - t**2 + 5*t + 1)/(2*t*(t - 1)), -(4*X*t**3 + 4*X*t**2 + 3*t**3 - t**2 + 5*t + 1)/(2*X*t**2*(t - 1)), t],
    [-(t + 1)*(t**2 + 1)/(t*(t - 1)), 1, -1/(X*t), 0],
    [-(t + 1)*(2*X + t**2 + 1)/(t*(t - 1)), 1, (2*X*t**3 + 6*X*t**2 + 6*X*t + 2*X - t + 1)/(X*t*(t - 1)), 2*(t + 1)/(t - 1)],
    [-(2*X*t**2 + t**4 - 1)/(2*t**2), (t - 1)**2/(2*t), -(4*X*t**3 + 4*X*t**2 + t**3 - 3*t**2 + 3*t - 1)/(2*X*t**2*(t - 1)), t]])

check('M_snub · A = B · M_Z² (exakt, alla t)', sp.simplify(sp.expand(sp.together(Ms*A - B*Mz))) == sp.zeros(4, 4))
D = 4*(t + 1)**4*(t**2 + 1)**2/(t - 1)**4
for name, M in (('A', A), ('B', B)):
    check(f'det {name} = 4(t+1)⁴(t²+1)²/(t−1)⁴ (konstant)', sp.cancel(M.det() - D) == 0)
check('t = 2 + √3 ger λ = 1', sp.simplify(lam.subs(t, 2 + sp.sqrt(3)) - 1) == 0)
u = (t - 1)/(t + 1)                                      # u = tan(Φ/2) när t = tan(π/4 + Φ/2)
check('cos Φ = (1 − u²)/(1 + u²) = 1/(1 + λ)', sp.cancel((1 - u**2)/(1 + u**2) - 1/(1 + lam)) == 0)
# --- Granskning 4, fynd 4: "entydig lösning" kontrolleras ---
# Okända: en koefficient per monom i varje element av A och B (samma monom som i A, B ovan), normering A[1,1] = 1.
# Om ekvationssystemet M_snub·A = B·M_Z² har en entydig lösning i en rationell punkt t₀, så har koefficientmatrisen full
# rang där, alltså också över Q(t) (rangen kan bara sjunka vid specialisering), och lösningen över Q(t) är entydig.
def monom(e):
    num, den = sp.fraction(sp.cancel(sp.together(e)))
    pd = sp.Poly(den, X, Z); (dx, dz), _ = pd.terms()[0]
    return [(ex - dx, ez - dz) for (ex, ez), _ in sp.Poly(sp.expand(num), X, Z).terms()]


def entydig(t0):
    unk, Au, Bu = [], sp.zeros(4, 4), sp.zeros(4, 4)
    for Mx, U, nm in ((A, Au, 'a'), (B, Bu, 'b')):
        for i in range(4):
            for j in range(4):
                for (ex, ez) in monom(Mx[i, j]):
                    s = sp.Symbol(f'{nm}{i}{j}_{ex}_{ez}'.replace('-', 'm')); unk.append(s)
                    U[i, j] += s * X**ex * Z**ez
    eqs = []
    E = (Ms * Au - Bu * Mz).subs(t, t0)
    for e in E:
        num = sp.fraction(sp.together(sp.expand(e)))[0]
        if num != 0: eqs += [c for c in sp.Poly(sp.expand(num), X, Z).coeffs()]
    eqs.append(Au[1, 1] - 1)
    Mlin, rhs = sp.linear_eq_to_matrix(eqs, unk)
    r = Mlin.rank()
    sol = sp.linsolve((Mlin, rhs), unk)
    ok_sol = len(sol) == 1 and all(not v.free_symbols for v in next(iter(sol)))
    if ok_sol:
        vals = dict(zip(unk, next(iter(sol))))
        same = all(sp.simplify(Au.subs(vals)[k] - A.subs(t, t0)[k]) == 0 and sp.simplify(Bu.subs(vals)[k] - B.subs(t, t0)[k]) == 0
                   for k in range(16))
    return len(unk), r, ok_sol and same


for t0 in (sp.Integer(2), sp.Rational(5, 7)):
    n, r, good = entydig(t0)
    check(f'entydig lösning med samma monom och A[1,1] = 1 i t = {t0} ({n} okända, rang {r}), lika med A, B ovan', good and r == n)

# --- Granskning 4, avsnitt 2: konstanterna i (X, Z) = (1, 1) (används i identifieringen med Kenyons G) ---
one = sp.ones(4, 1)
A1, B1 = A.subs({X: 1, Z: 1}), B.subs({X: 1, Z: 1})
check('A(1,1)·1 = ((t+1)²/(t−1))·1', sp.simplify(A1 * one - (t + 1)**2/(t - 1) * one) == sp.zeros(4, 1))
check('1ᵀ·B(1,1) = (2(t²+1)/(t−1))·1ᵀ', sp.simplify(one.T * B1 - 2*(t**2 + 1)/(t - 1) * one.T) == sp.zeros(1, 4))

print(f'\n{sum(ok)} godkända, {len(ok) - sum(ok)} fel')
