"""
Uppföljning 2026-10-18, kontroll 2: varför har snub square genus 0 för alla λ?

Påstående (exakt): snub square med konduktans λ på triangel–triangel-kanterna och 1 på övriga har SAMMA
spektralkurva (samma Bloch-determinant, identiskt i X, Z) som kvadratgittret Z² med perioden 2×2, vertikala
konduktanser 1 och horisontella omväxlande 1/q och q, där q = 1 + λ + κ = tan(π/4 + Φ/2), cos Φ = 1/(1+λ).
Det viktade Z² är isoradiellt (Kenyons kritiska vikter), så genus 0 följer av den kända satsen, och Kenyons
c·R = 2θ/π på Z² ger via Jacobis formel R_TT(λ) = 2Φ/(πκ) = 2 arctan κ/(πκ).

Snub square ritas här som Z² med diagonaler (de två triangel–triangel-kanterna i varje cell):
cell = 2×2 rutor av Z², diagonal (0,0)–(1,1) och (2,1)–(1,2). Kör:  python kod/uppfoljning_isoradiell.py
"""
import sys
sys.stdout.reconfigure(encoding="utf-8")
import sympy as sp

X, Z = sp.symbols('X Z')
ok = []


def check(name, cond):
    ok.append(bool(cond))
    print(('OK  ' if cond else 'FEL ') + name, flush=True)


def lap(edges):
    """Bloch-Laplacian på cellen {0,1}², kanter (p, q, vikt) med p, q i Z²; X, Z = förflyttning med 2 steg."""
    idx = {(i, j): 2 * i + j for i in range(2) for j in range(2)}
    M = sp.zeros(4, 4)
    for p, q, w in edges:
        a, b = idx[(p[0] % 2, p[1] % 2)], idx[(q[0] % 2, q[1] % 2)]
        sx = (q[0] - q[0] % 2) // 2 - (p[0] - p[0] % 2) // 2
        sy = (q[1] - q[1] % 2) // 2 - (p[1] - p[1] % 2) // 2
        ph = X**sx * Z**sy
        M[a, a] += w; M[b, b] += w; M[a, b] -= w * ph; M[b, a] -= w / ph
    return M


def z2_edges(h, v):
    E = []
    for i in range(2):
        for j in range(2):
            E.append(((i, j), (i + 1, j), h[2 * i + j]))
            E.append(((i, j), (i, j + 1), v[2 * i + j]))
    return E


def snub_edges(lam):
    return z2_edges([1] * 4, [1] * 4) + [((0, 0), (1, 1), lam), ((2, 1), (1, 2), lam)]


lam = sp.Symbol('lambda', positive=True)
q = sp.Symbol('q', positive=True)

# 1. grafen är snub square: alla noder har grad 5, ytor 3.3.4.3.4 (4 trianglar, 2 kvadrater per cell)
deg = {}
for p, r, w in snub_edges(1):
    for s in (p, r):
        deg[(s[0] % 2, s[1] % 2)] = deg.get((s[0] % 2, s[1] % 2), 0) + 1
check('Z² + diagonaler: fyra noder per cell, alla av grad 5', sorted(deg.values()) == [5, 5, 5, 5])

# 2. determinanten är artikelns 4·D_λ(a, b), a = (X+1/X)/2, b = (Z+1/Z)/2
a, b = (X + 1 / X) / 2, (Z + 1 / Z) / 2
D = a**2 + b**2 - (4*lam**2 + 8*lam + 2)*a*b - 4*(lam + 2)**2*(a + b) + 4*(3*lam + 4)*(lam + 2)
dsnub = lap(snub_edges(lam)).det()
check('det M_snub(λ) = 4·D_λ(a, b) (samma kurva som i viktad_avsnitt.tex)', sp.simplify(sp.expand(dsnub - 4*D)) == 0)

# 3. huvudidentiteten
L_of_q = (q + 1/q)/2 - 1                          # q = 1+λ+κ  <=>  λ = (q + 1/q)/2 − 1
dz2 = lap(z2_edges([1/q, q, q, 1/q], [1, 1, 1, 1])).det()
check('det M_snub(λ) = det Δ_Z²(h = 1/q, q, q, 1/q; v = 1) identiskt i X, Z',
      sp.simplify(sp.expand(lap(snub_edges(L_of_q)).det() - dz2)) == 0)
kap = sp.sqrt(lam*(lam + 2))
check('q = 1 + λ + κ uppfyller λ = (q + 1/q)/2 − 1',
      sp.simplify(L_of_q.subs(q, 1 + lam + kap) - lam) == 0)

# 4. de viktade Z² är isoradiella: c = tan θ med θ_v = π/4, θ_h = π/4 ∓ Φ/2;
#    vinkelsumma π vid varje nod och i varje ruta (rombvinklarna 2θ fyller 2π)
Phi = sp.Symbol('Phi', positive=True)
th_p, th_q, th_v = sp.pi/4 - Phi/2, sp.pi/4 + Phi/2, sp.pi/4
u = sp.Symbol('u', positive=True)                # u = tan(Φ/2)
tq = (1 + u)/(1 - u)                                # tan(π/4 + Φ/2), additionsformeln
secpt = (1 + u**2)/(1 - u**2) + 2*u/(1 - u**2)      # sec Φ + tan Φ uttryckt i u
check('tan(π/4 + Φ/2) = sec Φ + tan Φ = q och tan(π/4 − Φ/2) = (1 − u)/(1 + u) = 1/q',
      sp.simplify(tq - secpt) == 0 and sp.simplify(tq*(1 - u)/(1 + u) - 1) == 0
      and abs(sp.N(sp.tan(th_q).subs(Phi, 0.7) - (1/sp.cos(0.7) + sp.tan(0.7)))) < 1e-12)
theta = {}
for i in range(2):
    for j in range(2):
        theta[('h', i, j)] = th_p if (i + j) % 2 == 0 else th_q    # h00 = h11 = 1/q, h01 = h10 = q
        theta[('v', i, j)] = th_v
node = all(sp.simplify(theta[('h', i, j)] + theta[('h', (i - 1) % 2, j)] + theta[('v', i, j)]
                       + theta[('v', i, (j - 1) % 2)] - sp.pi) == 0 for i in range(2) for j in range(2))
face = all(sp.simplify(theta[('h', i, j)] + theta[('h', i, (j + 1) % 2)] + theta[('v', i, j)]
                       + theta[('v', (i + 1) % 2, j)] - sp.pi) == 0 for i in range(2) for j in range(2))
check('isoradiellt: Σθ = π vid varje nod och i varje ruta, alla θ i (0, π/2) för 0 < Φ < π/2', node and face)
# tågspår: NE-vinklar {0, −Φ}, NW-vinklar {π/2 − Φ, π/2}; θ = (β − α)/2 för kanten där spåren korsas
al = {sp.Rational(1, 2): 0, sp.Rational(3, 2): -Phi}
be = {sp.Rational(1, 2): sp.pi/2 - Phi, sp.Rational(3, 2): sp.pi/2}
tt = True
for i in range(2):
    for j in range(2):
        for kind, (u, w) in (('h', (i - j + sp.Rational(1, 2), i + j + sp.Rational(1, 2))),
                             ('v', (i - j - sp.Rational(1, 2), i + j + sp.Rational(1, 2)))):
            th = (be[w % 2] - al[u % 2]) / 2
            tt = tt and sp.simplify(th - theta[(kind, i, j)]) == 0
check('rombinbäddning finns: tågspårsvinklar α = {0, −Φ}, β = {π/2 − Φ, π/2} ger alla θ', tt)

# 5. Jacobi + Kenyon: dm/dq = Σ R_e dc_e/dq med c_e R_e = 2θ_e/π  =>  2 R_TT = dm/dλ
Rp, Rq = (2*th_p/sp.pi)*q, (2*th_q/sp.pi)/q
dm_dq = sp.simplify(2*Rp*(-1/q**2) + 2*Rq)
check('Kenyon på Z²: dm/dq = 4Φ/(πq)', sp.simplify(dm_dq - 4*Phi/(sp.pi*q)) == 0)
qq = 1 + lam + kap
RTT = sp.simplify((dm_dq.subs(q, qq) * sp.diff(qq, lam) / 2).subs(Phi, sp.atan(kap)))
check('R_TT(λ) = (1/2)·dm/dλ = 2 arctan κ/(πκ) (Sats i viktad_avsnitt.tex)',
      sp.simplify(RTT - 2*sp.atan(kap)/(sp.pi*kap)) == 0)

# 6. snub square själv har ingen periodisk isoradiell inbäddning (för något λ):
#    φ_e = halva medelpunktsvinkeln, Σφ = π i varje yta. Per cell: 2 kvadrater (bara TS-kanter, 8 st),
#    4 trianglar (varje TT-kant i två trianglar, varje TS-kant i en). Summera över cellens ytor:
sTT, sTS = sp.symbols('sTT sTS')                 # summan av φ över cellens 2 TT- resp. 8 TS-kanter
sol = sp.solve([sTS - 2*sp.pi, 2*sTT + sTS - 4*sp.pi], [sTT, sTS])
check('ytsummorna tvingar Σ_TT φ = π, alltså medel φ_TT = π/2 (urartat, θ_TT = 0); omöjligt om varje φ < π/2',
      sol[sTT] == sp.pi)

# 7. numerisk kontroll av Kenyon c·R = 2θ/π för en horisontell kant i det viktade Z² (λ = 1, q = 2 + √3)
try:
    import numpy as np
    qv = 2 + np.sqrt(3); N = 600
    t = (np.arange(N) + 0.5) * 2*np.pi/N
    T1, T2 = np.meshgrid(t, t, indexing='ij')
    Mq = lap(z2_edges([1/q, q, q, 1/q], [1]*4)).subs(q, qv)
    fs = [[sp.lambdify((X, Z), e, 'numpy') for e in row] for row in Mq.tolist()]
    Xs, Zs = np.exp(1j*T1), np.exp(1j*T2)
    one = np.ones_like(Xs)
    Ms = np.moveaxis(np.array([[e * one for e in row] for row in [[g(Xs, Zs) for g in r] for r in fs]], dtype=complex), [0, 1], [-2, -1])
    G = np.linalg.inv(Ms)
    # kant (0,0)–(1,0): index 0 och 2, samma cell
    R = np.real(G[..., 0, 0] + G[..., 2, 2] - G[..., 0, 2] - G[..., 2, 0]).mean()
    check(f'numeriskt: (1/q)·R_h00 = {R/qv:.6f} ≈ 2θ_p/π = 1/6 = {1/6:.6f} (λ = 1, Φ = π/3)', abs(R/qv - 1/6) < 1e-3)
except ImportError:
    print('(numpy saknas; numerisk kontroll 7 hoppas över)')

print(f'\n{sum(ok)} godkända, {len(ok) - sum(ok)} fel')
