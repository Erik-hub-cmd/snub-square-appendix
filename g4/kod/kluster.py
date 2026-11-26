"""Klusterförmodan (GALOIS_UTVIDGNING.md, Förmodan K): om N(det Δ) = Z_G så är sidpolynomet längs s en produkt
över kluster C (sammanhängande komponenter av korsningsrelationen bland spåren ∥ s) av faktorer Q_C av grad
Σ_{T∈C} k_T (T har klass k_T·s), med rotprodukt ±∏_{T∈C} μ_T^{±1}."""
import json, sys, random, math, itertools
sys.path.insert(0, '.')
from natgalois import *
from fractions import Fraction

def cluster_check(pos, v1, v2, E, W):
    T = tracks(pos, v1, v2, E); mi = minimal_info(E, T); mu = monodromy(T, W)
    P = det_poly(len(pos), E, W); H, S = sides(P)
    NZ = sorted((2 * x, 2 * y) for x, y in H) == zonogon(mi['h'])
    out = []
    for sd in S[:len(S) // 2]:
        s = sd['dir']
        par = [i for i, h in enumerate(mi['h']) if h != (0, 0) and h[0] * s[1] - h[1] * s[0] == 0]
        k = {i: abs(mi['h'][i][0] // s[0]) if s[0] else abs(mi['h'][i][1] // s[1]) for i in par}
        # kluster
        comp = {i: i for i in par}
        def f(i):
            while comp[i] != i: i = comp[i]
            return i
        for i, j in itertools.combinations(par, 2):
            if mi['C'][i][j] > 0: comp[f(i)] = f(j)
        for i in par:
            if mi['C'][i][i] > 0: pass
        cl = {}
        for i in par: cl.setdefault(f(i), []).append(i)
        clusters = []
        for c in cl.values():
            D = sum(k[i] for i in c); m = Fraction(1)
            for i in c: m *= mu[i]
            clusters.append((D, m, c))
        c0, fl = sp.factor_list(sd['poly'].as_expr(), U)
        facs = []
        for fa, mult in fl:
            p = sp.Poly(fa, U); dg = p.degree(); co = p.all_coeffs()
            r = Fraction(int(co[-1]), int(co[0])) * (-1) ** dg
            facs += [(dg, r)] * mult
        ok = match(facs, clusters)
        out.append({'dir': s, 'len': sd['len'], 'deg_sum': sum(c[0] for c in clusters), 'clusters': [(c[0], str(c[1]), len(c[2])) for c in clusters],
                    'factors': [(d, str(r)) for d, r in facs], 'match': ok,
                    'irr_cluster': any(d >= 2 for d, _ in facs)})
    hs = mi['h']; Hh = all(h != (0, 0) for h in hs)
    for i, j in itertools.combinations(range(len(hs)), 2):
        wd = abs(hs[i][0] * hs[j][1] - hs[i][1] * hs[j][0])
        if wd and mi['C'][i][j] != wd: Hh = False
    return {'NZ': NZ, 'H': Hh, 'minimal': mi['minimal'], 'sides': out}

def match(facs, clusters):
    n = len(clusters)
    if sum(d for d, _ in facs) != sum(c[0] for c in clusters): return False
    def rec(i, deg, prod):
        if i == len(facs):
            for (D, m, _), dd, pp in zip(clusters, deg, prod):
                if dd != D: return False
                if not any(pp == sg * m ** e for sg in (1, -1) for e in (1, -1)): return False
            return True
        d, r = facs[i]
        for j in range(n):
            if deg[j] + d <= clusters[j][0]:
                deg[j] += d; prod[j] *= r
                if rec(i + 1, deg, prod): return True
                deg[j] -= d; prod[j] /= r
        return False
    return rec(0, [0] * n, [Fraction(1)] * n)

if __name__ == '__main__':
    files = sys.argv[1].split(','); mode = sys.argv[2]
    random.seed(11)
    items = []
    for fn in files:
        for t in json.load(open(fn)):
            items.append((t.get('name') or ' ; '.join(t['types']), t['summary']['k']) + from_json(t))
    if 'uniform12.json' in files:
        s = 1 / math.sqrt(2); L = 1 + math.sqrt(2); pos = [s, 1j * s, -s, -1j * s]
        items.append(('4.8^2', 1, pos, L, 1j * L, geom_edges(pos, L, 1j * L)))
    stat = {}
    for idx, (nm, kk, pos, v1, v2, E) in enumerate(items):
        W = [1] * len(E) if mode == 'enhet' else [random.randint(1, 97) for _ in E]
        r = cluster_check(pos, v1, v2, E, W)
        allm = all(s['match'] for s in r['sides'])
        key = ('N=Z' if r['NZ'] else 'N≠Z', '(H)' if r['H'] else 'ej (H)', 'minimal' if r['minimal'] else 'ej minimal', 'kluster stämmer' if allm else 'kluster STÄMMER EJ')
        stat[key] = stat.get(key, 0) + 1
        flag = '' if allm or not (r['NZ'] and r['H']) else '   <== AVVIKELSE'
        print(f"[{idx}] {nm} (k={kk}, V={len(pos)}): {key}{flag}")
        for s in r['sides']:
            print(f"     sida {s['dir']} längd {s['len']}: kluster (grad, ∏μ, #spår) {s['clusters']}  faktorer (grad, rotprodukt) {s['factors']}  {'ok' if s['match'] else 'NEJ'}")
    print('SAMMANFATTNING', stat)
