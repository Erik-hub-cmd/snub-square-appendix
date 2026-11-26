"""Sidpolynom för de fyra Galois-hindrade arkimediska näten med en vikt per kantklass (symboliskt)."""
import json, math, sys
sys.path.insert(0, '.')
from natgalois import *
z, w = sp.symbols('z w')
def symdet(nV, E, W):
    M = sp.zeros(nV, nV)
    for (i, j, n, _), c in zip(E, W):
        x = z**n[0] * w**n[1]
        if i == j:
            M[i, i] += c * (2 - x - 1/x); continue
        M[i, i] += c; M[j, j] += c; M[i, j] -= c * x; M[j, i] -= c / x
    D = sp.expand(sp.cancel(M.det(method='berkowitz')))
    P = {}
    for term in sp.Add.make_args(D):
        c, m = term.as_coeff_Mul()
        # Laurent: använd as_powers_dict
        pd = term.as_powers_dict(); a = int(pd.get(z, 0)); b = int(pd.get(w, 0))
        coef = term / (z**a * w**b)
        P[(a, b)] = sp.expand(P.get((a, b), 0) + coef)
    return {k: v for k, v in P.items() if v != 0}
d = json.load(open('uniform12.json'))
it = {t['name']: from_json(t) for t in d if t['summary']['k'] == 1}
s = 1 / math.sqrt(2); L = 1 + math.sqrt(2); pos = [s, 1j * s, -s, -1j * s]
E48 = geom_edges(pos, L, 1j * L)
E48 = [(i, j, n, 0 if n == (0, 0) else 1) for (i, j, n, _) in E48]
it['4.8^2'] = (pos, L, 1j * L, E48)
a, b, c = sp.symbols('a b c', positive=True)
if __name__ == "__main__":
  for name in ['3^2.4.3.4', '3^4.6', '4.8^2', '4.6.12']:
      pos, v1, v2, E = it[name]
      cls = sorted(set(e[3] for e in E)); sym = [1, a, b, c][:len(cls)]
      W = [sym[cls.index(e[3])] for e in E]
      from collections import Counter
      print('==', name, 'kantklasser', dict(Counter(e[3] for e in E)), '-> vikter', sym)
      P = symdet(len(pos), E, W)
      H = hull(list(P))
      for A, B in list(zip(H, H[1:] + H[:1]))[:len(H)//2]:
          dd = (B[0]-A[0], B[1]-A[1]); g = math.gcd(abs(dd[0]), abs(dd[1])); st = (dd[0]//g, dd[1]//g)
          poly = sum(P.get((A[0]+k*st[0], A[1]+k*st[1]), 0) * U**k for k in range(g+1))
          f = sp.factor(poly)
          disc = sp.factor(sp.discriminant(poly, U)) if g == 2 else None
          print('   sida', st, 'längd', g, ':', f, '   diskriminant:', disc)
