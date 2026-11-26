"""Separat kontroll: det Δ med sympy (Berkowitz, symboliskt i z, w) i stället för modulär evaluering + CRT."""
import json, sys, math
sys.path.insert(0, '.')
from natgalois import from_json, hull, U, det_poly
import sympy as sp
from symb_klass import symdet
d = json.load(open('uniform12.json'))
for idx in [int(x) for x in sys.argv[1].split(',')]:
    t = d[idx]; pos, v1, v2, E = from_json(t)
    P = symdet(len(pos), E, [1] * len(E))
    Pi = {k: int(v) for k, v in P.items()}
    Q = det_poly(len(pos), E, [1] * len(E))
    H = hull(list(P)); sid = []
    for A, B in list(zip(H, H[1:] + H[:1]))[:len(H) // 2]:
        dd = (B[0] - A[0], B[1] - A[1]); g = math.gcd(abs(dd[0]), abs(dd[1])); st = (dd[0] // g, dd[1] // g)
        sid.append(str(sp.factor(sum(P.get((A[0] + k * st[0], A[1] + k * st[1]), 0) * U ** k for k in range(g + 1)))))
    print(f"[{idx}] {t['name']}: sympy == CRT: {Pi == Q};  sidor: {sid}", flush=True)
