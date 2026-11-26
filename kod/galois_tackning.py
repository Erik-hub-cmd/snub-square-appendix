"""Galois-hindret på övertorusar: sidpolynomen till det M för snub square (λ = 1) på T_Λ har inga rationella rötter."""
import sympy as sp
from fractions import Fraction
from torus_sok import snub_karta, edges_of, tails
from galois_lemma import detpoly, sidpolynom, w
import torus_sok
def snub_cond(B, lam):
    rot, volt = snub_karta(B)
    red, idx = torus_sok.reducer(B)
    reps = sorted({red((x, y))[0] for x in range(-20, 21) for y in range(-20, 21)})
    e = 0; cond = {}
    for (x, y) in reps:
        dirs = [(1, 0), (0, 1)] + ([(1, 1)] if x % 2 == 0 and y % 2 == 0 else []) + ([(-1, 1)] if x % 2 == 0 and y % 2 == 1 else [])
        for dl in dirs:
            cond[e] = lam if dl[0] * dl[1] != 0 else sp.Integer(1); e += 1
    return rot, volt, cond
for namn, Bc in (('egen', ((1, 0), (0, 1))), ('2x1', ((2, 0), (0, 1))), ('diagonal <a,b>', ((1, -1), (1, 1))), ('3x1', ((3, 0), (0, 1))), ('2x2', ((2, 0), (0, 2)))):
    B = tuple((2 * x, 2 * y) for x, y in Bc)
    rot, volt, cond = snub_cond(B, sp.Integer(1))
    P = detpoly(rot, volt, cond, B)
    print(f'== {namn}:')
    for nn, g, f in sidpolynom(P):
        fl = sp.factor_list(f, w)[1]
        rat = [fa for fa, m in fl if sp.Poly(fa, w).degree() == 1]
        print(f'   normal {nn}, längd {g}: faktorer {[(fa, m) for fa, m in fl]}  rationella rötter: {len(rat)}')

# Gitternät (det enda minimala nätet med denna Newtonpolygon) på samma torusar, slumpvisa rationella konduktanser
import random
from galois_lemma import z2_karta
rng = random.Random(11)
for namn, Bc in (('egen', ((1, 0), (0, 1))), ('2x1', ((2, 0), (0, 1))), ('diagonal <a,b>', ((1, -1), (1, 1))), ('3x1', ((3, 0), (0, 1))), ('2x2', ((2, 0), (0, 2)))):
    B = tuple((2 * x, 2 * y) for x, y in Bc)
    rot, volt = z2_karta(B)
    cond = {e: sp.Rational(rng.randint(1, 9), rng.randint(1, 9)) for e in edges_of(rot)}
    P = detpoly(rot, volt, cond, B)
    lin = all(sp.Poly(fa, w).degree() <= 1 for nn, g, f in sidpolynom(P) for fa, m in sp.factor_list(f, w)[1])
    print(f'gitter {namn}: {len(edges_of(rot))} kanter, alla sidpolynom i rationella linjära faktorer: {lin}')
