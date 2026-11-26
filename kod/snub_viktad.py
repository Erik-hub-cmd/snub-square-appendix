"""
Kandidat A4: snub square-gittret med två konduktanser.
c på triangel–triangel-kanterna (2 per cell), d på triangel–kvadrat-kanterna (8 per cell), λ = c/d.
Vi sätter d = 1 (resistanser för allmänt d fås genom R(c, d) = R(λ, 1)/d).

Samma enhetscell och kantlista som artikeln (Sektion 2) och kod/exact_snub.py.
"""
import sympy as sp

X, Z = sp.symbols('X Z')
lam = sp.Symbol('lambda', positive=True)
a, b = sp.symbols('a b', real=True)
t1, t2 = sp.symbols('theta1 theta2', real=True)

# (i, j, (n1, n2), klass) ; klass 'ts' = triangel–kvadrat (vikt 1), 'tt' = triangel–triangel (vikt λ)
EDGES = [(0, 1, (-1, 0), 'ts'), (0, 1, (0, 0), 'ts'), (0, 3, (0, 0), 'ts'), (0, 3, (0, 1), 'ts'),
         (1, 2, (0, 0), 'ts'), (1, 2, (0, 1), 'ts'), (2, 3, (0, 0), 'ts'), (2, 3, (1, 0), 'ts'),
         (0, 2, (0, 1), 'tt'), (1, 3, (1, 0), 'tt')]
TT_EDGE = (0, 2, (0, 1))       # R_TT
TS_EDGE = (0, 1, (0, 0))       # R_TS


def weight(cls, lv=lam):
    return lv if cls == 'tt' else 1


def bloch(lv=lam):
    """viktad Bloch-laplacian M(X, Z) (Laurentpolynom); varje nod har grad 4 + λ"""
    M = (4 + lv) * sp.eye(4)
    for i, j, (n1, n2), cls in EDGES:
        w = weight(cls, lv)
        ph = X**n1 * Z**n2
        M[i, j] -= w * ph
        M[j, i] -= w / ph
    return M


def to_ab(expr):
    """medelvärde över teckenbyten av θ1, θ2; uttryck i a = cos θ1, b = cos θ2"""
    f = sp.expand(expr.subs({X: sp.exp(sp.I * t1), Z: sp.exp(sp.I * t2)}))
    g = (f + f.subs(t1, -t1) + f.subs(t2, -t2) + f.subs({t1: -t1, t2: -t2}, simultaneous=True)) / 4
    g = sp.expand_trig(sp.expand(sp.expand_complex(sp.expand(g.rewrite(sp.cos)))))
    g = sp.expand(g.subs({sp.sin(t1)**2: 1 - sp.cos(t1)**2, sp.sin(t2)**2: 1 - sp.cos(t2)**2}))
    return sp.expand(sp.simplify(g.subs({sp.cos(t1): a, sp.cos(t2): b})))


def spectral_poly(lv=lam):
    """P(X, Z) = X^2 Z^2 det M som polynom"""
    return sp.expand(sp.cancel(X**2 * Z**2 * bloch(lv).det()))


def numerator(edge, lv=lam):
    i, j, n = edge
    M = bloch(lv)
    A = M.adjugate()
    ph = X**n[0] * Z**n[1]
    return sp.expand(A[i, i] + A[j, j] - A[i, j] / ph - A[j, i] * ph)


if __name__ == '__main__':
    M = bloch()
    det = sp.expand(M.det())
    Dab = sp.factor(to_ab(det))
    print('det M(a, b; λ) =', Dab)
    D1 = a**2 + b**2 - 14*a*b - 36*(a + b) + 84
    print('λ = 1 ger 4 D(a, b):', sp.expand(Dab.subs(lam, 1) - 4 * D1) == 0)
    print('P(X, Z) =', sp.factor(spectral_poly()))
