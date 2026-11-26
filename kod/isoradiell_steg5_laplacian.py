"""
Isoradiell förklaring, steg 5: Laplacian-nivå (λ = 1, q = 2 + √3).

Lokala matriser A(X, Z), B(X, Z) (Laurentpolynom, koefficienter i Q(√3)) med
    M_snub · A = B · M_Z²,   det A = det B = 576 q²  (konstanter).
Alltså M_snub = B · M_Z² · A⁻¹ med A, B unimodulära, och Greenfunktionerna uppfyller M_snub⁻¹ = A · M_Z²⁻¹ · B⁻¹:
varje resistans på snub square är en ändlig linjärkombination av potentialkärnans värden på det isoradiella Z².
A fås ur G i steg 4 (A = S⁻¹ G S, S = byte från cellvektor till transfertillstånd); B = M_snub A adj(M_Z²)/det.
Kör: python kod/isoradiell_steg5_laplacian.py
"""
import sys
sys.stdout.reconfigure(encoding="utf-8")
import sympy as sp
from sympy import sqrt

X, Z = sp.symbols('X Z')
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


q = 2 + sqrt(3)
Ms = lap(z2e([1]*4, [1]*4) + [((0, 0), (1, 1), 1), ((2, 1), (1, 2), 1)])     # snub square, λ = 1
Mz = lap(z2e([1/q, q, q, 1/q], [1]*4))                                        # isoradiellt Z²
A = sp.Matrix([
    [(24*sqrt(3)*X - 7 + 4*sqrt(3))/X, -4*sqrt(3)/Z, 12 + (-8*sqrt(3)*X - 14*sqrt(3))/X, (-2 + sqrt(3))/(X*Z)],
    [4*sqrt(3)*Z, 1, Z*(2 - sqrt(3)), 0],
    [-4*sqrt(3)*X - 4*sqrt(3) + 6, -X/Z, -2*X + sqrt(3)*X + 12*sqrt(3), -2*sqrt(3)/Z],
    [Z, 0, 2*sqrt(3)*Z, sqrt(3) + 2]])
B = sp.Matrix([
    [-X + 20*sqrt(3), 1 + 4*sqrt(3), (-2*sqrt(3)*X - 7*sqrt(3) + 10)/X, sqrt(3) + 2],
    [-4*sqrt(3), 1, (-2 + sqrt(3))/X, 0],
    [-4*sqrt(3)*X + 6*X - 4*sqrt(3), 1, (12*sqrt(3)*X - 2 + sqrt(3))/X, 2*sqrt(3)],
    [-X - 4*sqrt(3), 1, (-2*sqrt(3)*X - 2 + sqrt(3))/X, sqrt(3) + 2]])

check('M_snub · A = B · M_Z² (exakt)', all(sp.simplify(sp.expand(sp.together(e)*X**3*Z**3)) == 0 for e in (Ms*A - B*Mz)))
for name, M in (('A', A), ('B', B)):
    d = sp.expand(sp.cancel(sp.together(M.det()))*X**4*Z**4 - 576*q**2*X**4*Z**4)
    check(f'det {name} = 576 q² (konstant)', sp.simplify(d) == 0)
print(f'\n{sum(ok)} godkända, {len(ok) - sum(ok)} fel')
