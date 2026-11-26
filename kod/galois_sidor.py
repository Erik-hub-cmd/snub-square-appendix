"""Sidpolynomen till det M för snub square (allmänt λ) och isoradiellt Z²_t. Kräver isoradiell_steg8_allmant_lambda.py."""
src = open('isoradiell_steg8_allmant_lambda.py', encoding='utf-8').read().split("check('M_snub")[0]
exec(src)
import sympy as sp
L = sp.symbols('lambda', positive=True)
def lapL(edges):
    return lap(edges)
Ms_l = lap(z2e([1]*4, [1]*4) + [((0, 0), (1, 1), L), ((2, 1), (1, 2), L)])
def sidor(M, name, var):
    P = sp.cancel(M.det())
    num, den = sp.fraction(sp.together(P))
    pd = sp.Poly(den, X, Z); (dx, dz), _ = pd.terms()[0]
    terms = {(ex - dx, ez - dz): c / pd.coeffs()[0] for (ex, ez), c in sp.Poly(sp.expand(num), X, Z).terms()}
    pts = list(terms)
    print(f'== {name}: Newtonpolygonens hörn/stöd: x+z in', sorted({p[0]+p[1] for p in pts}), ' x-z in', sorted({p[0]-p[1] for p in pts}))
    for lin, sgn in (((1, 1), 1), ((1, 1), -1), ((1, -1), 1), ((1, -1), -1)):
        val = max(sgn*(lin[0]*p[0] + lin[1]*p[1]) for p in pts)
        side = sorted([p for p in pts if sgn*(lin[0]*p[0] + lin[1]*p[1]) == val])
        w = sp.Symbol('w')
        poly = sum(sp.factor(terms[p]) * w**k for k, p in enumerate(side))
        print(f'  sida {sgn}*({lin}): punkter {side}')
        print(f'     sidpolynom = {sp.factor(poly)}')
        print(f'     diskriminant = {sp.factor(sp.discriminant(sp.Poly(sp.expand(poly), w)))}' if len(side) == 3 else '')
sidor(Ms_l, 'snub square (vikt λ)', L)
sidor(Mz, 'isoradiellt Z²_t', t)
