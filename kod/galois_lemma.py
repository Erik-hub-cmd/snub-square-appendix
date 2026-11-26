"""Kontroll av lemmat bakom Galois-hindret:
 (1) Z² (cell 2x2) med 8 godtyckliga konduktanser: varje sidpolynom till det M faktoriserar i linjära faktorer över Q(konduktanserna).
 (2) (borttagen; se lemma4_test.py och lemma4_pl.py)
 (3) Snub square med λ = 1: sidpolynomet w² − 14w + 1 är irreducibelt över Q (rötter (2±√3)²)."""
import random, itertools
from fractions import Fraction
import sympy as sp
from torus_sok import snub_karta, edges_of, tails, ta_bort, drag, genus
X, Z, w = sp.symbols('X Z w')

def detpoly(rot, volt, cond, per):
    vs = sorted(rot); ix = {v: i for i, v in enumerate(vs)}; t = tails(rot)
    n = len(vs); M = sp.zeros(n, n)
    for v, L in rot.items():
        for d in L:
            c = cond[d >> 1]; (vx, vy) = volt[d]
            # spänning i Λ-koordinater: volt = s*P0 + r*P1
            (a, b), (cc, dd) = per; D = a * dd - b * cc
            s = Fraction(vx * dd - vy * cc, D); r = Fraction(a * vy - b * vx, D)
            assert s.denominator == 1 and r.denominator == 1
            M[ix[v], ix[v]] += c
            M[ix[v], ix[t[d ^ 1]]] -= c * X**int(s) * Z**int(r)
    return sp.cancel(M.det())

def sidpolynom(P):
    num, den = sp.fraction(sp.together(P))
    pd = sp.Poly(den, X, Z); (dx, dz), _ = pd.terms()[0]
    terms = {(ex - dx, ez - dz): c / pd.coeffs()[0] for (ex, ez), c in sp.Poly(sp.expand(num), X, Z).terms()}
    pts = list(terms); ut = []
    # sidornas normaler: hitta konvexa höljets sidor genom att pröva normaler
    normals = set()
    for p, q in itertools.combinations(pts, 2):
        d = (q[0] - p[0], q[1] - p[1]); nrm = (d[1], -d[0])
        for s in (1, -1):
            nn = (s * nrm[0], s * nrm[1]); g = sp.gcd(nn[0], nn[1]); nn = (nn[0] // g, nn[1] // g)
            val = max(nn[0] * x + nn[1] * y for x, y in pts)
            side = [p_ for p_ in pts if nn[0] * p_[0] + nn[1] * p_[1] == val]
            if len(side) >= 2: normals.add(nn)
    for nn in sorted(normals):
        val = max(nn[0] * x + nn[1] * y for x, y in pts)
        side = sorted([p_ for p_ in pts if nn[0] * p_[0] + nn[1] * p_[1] == val])
        # parametrisera längs sidan med primitivt steg
        d = (side[-1][0] - side[0][0], side[-1][1] - side[0][1]); g = abs(sp.gcd(d[0], d[1])); st = (d[0] // g, d[1] // g)
        poly = sum(terms.get((side[0][0] + k * st[0], side[0][1] + k * st[1]), 0) * w**k for k in range(g + 1))
        ut.append((nn, g, sp.factor(poly)))
    return ut

import torus_sok
def z2_karta(B):
    rot, volt = snub_karta(B)
    # bygg om utan diagonaler: snub_karta numrerar kanter i ordning; diagonalerna är de tredje riktningarna
    red, idx = torus_sok.reducer(B)
    reps = sorted({red((x, y))[0] for x in range(-20, 21) for y in range(-20, 21)})
    e = 0; tt = []
    for (x, y) in reps:
        dirs = [(1, 0), (0, 1)] + ([(1, 1)] if x % 2 == 0 and y % 2 == 0 else []) + ([(-1, 1)] if x % 2 == 0 and y % 2 == 1 else [])
        for dl in dirs:
            if dl[0] * dl[1] != 0: tt.append(e)
            e += 1
    r = ta_bort(rot, volt, None, tt)
    return r[0], r[1]

def main():
  pass
if __name__ == "__main__":
    ok = True
    # (1) Z² med symboliska konduktanser
    rot, volt = snub_karta(((2, 0), (0, 2)))
    t = tails(rot)
    diag = [e for e in edges_of(rot) if abs(volt[2*e][0]) + abs(volt[2*e][1]) not in (0, 2) or False]
    # TT-kanterna är de med Z²-förskjutning (1,1) eller (-1,1): identifiera via positioner
    rot, volt = z2_karta(((2, 0), (0, 2)))
    cs = sp.symbols('c0:8', positive=True)
    E = edges_of(rot); cond = {e: cs[i] for i, e in enumerate(E)}
    P = detpoly(rot, volt, cond, ((2, 0), (0, 2)))
    print('(1) Z², 8 symboliska konduktanser:')
    for nn, g, f in sidpolynom(P):
        lin = all(sp.Poly(fa, w).degree() <= 1 for fa, m in sp.factor_list(f, w)[1])
        ok &= lin
        print(f'    normal {nn}, längd {g}: {f}   linjära faktorer: {lin}')

    # (2) borttagen 2026-11-25: Z² har inga noder av grad 3 och inga trianglar, så vandringarna gjorde 0 drag och
    #     prövade bara startgittret. Lemma G.5 kontrolleras av lemma4_test.py och lemma4_pl.py.
    # (3)
    f = w**2 - 14*w + 1
    print('(3) snub λ=1, sidpolynom', f, 'irreducibelt över Q:', sp.Poly(f, w).is_irreducible, ' rötter', sp.solve(f, w))
    print('ALLT OK' if ok else 'FEL')
