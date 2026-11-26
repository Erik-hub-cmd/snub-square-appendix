"""Galois-testet för periodiska resistornät (GALOIS_UTVIDGNING.md).

För ett nät G på sin egen torus (cellbas v1, v2) räknas:
  * tågspåren (zig-zag-stigarna: sväng max vänster, max höger, omväxlande), deras klasser h_i i Z^2,
    korsningstalen, och om G är MINIMALT (alla h_i != 0, inga självkorsningar, och spår i, j korsas exakt |h_i ^ h_j| gånger);
  * det Δ_G(z, w) exakt (evaluering i rötter av enheten mod stora primtal + CRT), Newtonpolygonen och sidpolynomen;
  * faktorisering av sidpolynomen över Q, och Galois-hindret:
      H0 (egen torus):   något sidpolynom har en irreducibel faktor av grad >= 2;
      H∞ (alla torusar): någon rot ρ har en konjugerad rot med annat belopp  (=> ρ^m ∉ Q för alla m >= 1).
Indata: lista av noder (komplexa lägen), v1, v2 (komplexa), kanter (i, j, (n1, n2), vikt): nod i i cell 0 -- nod j i cell n.
"""
import cmath, math, random, itertools
from functools import reduce
import flint
import sympy as sp

ZETA = cmath.exp(1j * math.pi / 6)


def z12(c):
    return sum(ck * ZETA ** k for k, ck in enumerate(c))


def from_json(t):
    pos = [z12(v[0]) for v in t['verts']]
    v1, v2 = z12(t['v1']), z12(t['v2'])
    E = [(e[0], e[1], tuple(e[2]), e[3]) for e in t['edges']]
    return pos, v1, v2, E


def geom_edges(pos, v1, v2, R=2, tol=1e-7):
    """alla par på avstånd 1 (en representant per kant)"""
    E = []
    n = len(pos)
    for i in range(n):
        for j in range(n):
            for a in range(-R, R + 1):
                for b in range(-R, R + 1):
                    if (j, a, b) <= (i, 0, 0) and not (i == j and (a, b) == (0, 0)):
                        pass
                    d = pos[j] + a * v1 + b * v2 - pos[i]
                    if abs(abs(d) - 1) < tol:
                        key = (i, j, a, b)
                        rk = (j, i, -a, -b)
                        if rk < key or (rk == key):
                            continue
                        E.append((i, j, (a, b), 0))
    return E


def check_geom(pos, v1, v2, E):
    for (i, j, n, _) in E:
        d = pos[j] + n[0] * v1 + n[1] * v2 - pos[i]
        assert abs(abs(d) - 1) < 1e-6, (i, j, n, abs(d))


# ---------------------------------------------------------------- tågspår
def darts(pos, v1, v2, E):
    out = {i: [] for i in range(len(pos))}
    for k, (i, j, n, _) in enumerate(E):
        d = pos[j] + n[0] * v1 + n[1] * v2 - pos[i]
        out[i].append((cmath.phase(d), (i, j, n, k)))
        out[j].append((cmath.phase(-d), (j, i, (-n[0], -n[1]), k)))
    rot = {}
    for i in out:
        lst = [x[1] for x in sorted(out[i])]
        rot[i] = lst
    idx = {}
    for i, lst in rot.items():
        for p, d in enumerate(lst):
            idx[(d, )] = (i, p)
    return rot


def tracks(pos, v1, v2, E):
    rot = darts(pos, v1, v2, E)
    where = {}
    for i, lst in rot.items():
        for p, d in enumerate(lst):
            where[d] = (i, p)

    def rev(d):
        i, j, n, k = d
        return (j, i, (-n[0], -n[1]), k)

    seen = set(); T = []
    for i, lst in rot.items():
        for d0 in lst:
            for par0 in (0, 1):
                if (d0, par0) in seen: continue
                d, par = d0, par0; off = (0, 0); seq = []
                while (d, par) not in seen:
                    seen.add((d, par)); seq.append((d[3], off, par))
                    off = (off[0] + d[2][0], off[1] + d[2][1])
                    r = rev(d); v, p = where[r]; L = rot[v]
                    d = L[(p + 1) % len(L)] if par == 0 else L[(p - 1) % len(L)]
                    par ^= 1
                assert (d, par) == (d0, par0)
                T.append({'edges': seq, 'h': off, 'start': (d0, par0)})
    # varje oorienterat spår förekommer två gånger (två orienteringar); behåll en
    uniq = []; used = set()
    for t in T:
        key = frozenset(x[0] for x in t['edges']), len(t['edges'])
        es = sorted(x[0] for x in t['edges'])
        sig = (tuple(es), )
        # hitta par via gemensam kantmultimängd och motsatt klass
        uniq.append(t)
    # para ihop: spår t och s är samma oorienterade om de har samma kantmultimängd och h_s = -h_t
    res = []; taken = [False] * len(uniq)
    for a in range(len(uniq)):
        if taken[a]: continue
        ea = sorted(x[0] for x in uniq[a]['edges'])
        for b in range(a + 1, len(uniq)):
            if taken[b]: continue
            if sorted(x[0] for x in uniq[b]['edges']) == ea and uniq[b]['h'] == (-uniq[a]['h'][0], -uniq[a]['h'][1]):
                taken[a] = taken[b] = True; res.append(uniq[a]); break
        else:
            raise RuntimeError('oparat spår')
    return res


def minimal_info(E, T):
    nE = len(E)
    on = [[] for _ in range(nE)]
    for ti, t in enumerate(T):
        for (k, off, *_r) in t['edges']:
            on[k].append(ti)
    assert all(len(x) == 2 for x in on), 'varje kant ska korsas av två spårpassager'
    m = len(T)
    C = [[0] * m for _ in range(m)]
    for x in on:
        a, b = x
        if a == b: C[a][a] += 1
        else: C[a][b] += 1; C[b][a] += 1
    hs = [t['h'] for t in T]
    wedge = lambda p, q: p[0] * q[1] - p[1] * q[0]
    reasons = []
    for i in range(m):
        if hs[i] == (0, 0): reasons.append(f'spår {i} kontraktibelt')
        if C[i][i]: reasons.append(f'spår {i} korsar sig själv {C[i][i]} ggr')
        for j in range(i + 1, m):
            if C[i][j] != abs(wedge(hs[i], hs[j])):
                reasons.append(f'spår {i},{j}: {C[i][j]} korsningar, |h∧h|={abs(wedge(hs[i], hs[j]))}')
    return {'minimal': not reasons, 'reasons': reasons, 'h': hs, 'C': C,
            'lengths': [len(t['edges']) for t in T]}


# ---------------------------------------------------------------- det Δ
def _primes(N, count, start=2 ** 61):
    ps = []; k = (start // N) + 1
    while len(ps) < count:
        p = k * N + 1
        if sp.isprime(p): ps.append(p)
        k += 1
    return ps


def _root(p, N):
    g = sp.primitive_root(p)
    return pow(g, (p - 1) // N, p)


def det_poly(nV, E, W):
    """det Δ som dict {(a,b): heltal}; W = lista av heltalsvikter per kant"""
    D1 = sum(abs(e[2][0]) for e in E); D2 = sum(abs(e[2][1]) for e in E)
    N1, N2 = 2 * D1 + 1, 2 * D2 + 1
    # Hadamard-gräns
    rows = [0] * nV
    for (i, j, n, _), c in zip(E, W):
        rows[i] += 2 * abs(c); rows[j] += 2 * abs(c)
    bound = reduce(lambda a, b: a * b, [max(r, 1) for r in rows], 1) * N1 * N2
    N = N1 * N2 // math.gcd(N1, N2)
    need = bound.bit_length() // 60 + 2
    ps = _primes(N, need)
    residues = []
    for p in ps:
        om1 = pow(_root(p, N), N // N1, p); om2 = pow(_root(p, N), N // N2, p)
        vals = [[0] * N2 for _ in range(N1)]
        for a in range(N1):
            za = pow(om1, a, p)
            for b in range(N2):
                wb = pow(om2, b, p)
                M = [[0] * nV for _ in range(nV)]
                for (i, j, n, _), c in zip(E, W):
                    x = pow(za, n[0] % N1, p) * pow(wb, n[1] % N2, p) % p
                    xi = pow(x, p - 2, p)
                    M[i][i] += c; M[j][j] += c
                    M[i][j] -= c * x; M[j][i] -= c * xi
                d = flint.nmod_mat([[v % p for v in r] for r in M], p).det()
                # multiplicera med z^D1 w^D2 så att exponenterna blir 0..2D
                vals[a][b] = int(d) * pow(za, D1, p) * pow(wb, D2, p) % p
        # invers DFT
        i1 = pow(om1, p - 2, p); i2 = pow(om2, p - 2, p)
        tmp = [[sum(vals[a][b] * pow(i1, a * k, p) for a in range(N1)) % p for b in range(N2)] for k in range(N1)]
        co = [[sum(tmp[k][b] * pow(i2, b * l, p) for b in range(N2)) * pow(N1 * N2, p - 2, p) % p for l in range(N2)]
              for k in range(N1)]
        residues.append(co)
    out = {}
    Mtot = 1
    for p in ps: Mtot *= p
    for k in range(N1):
        for l in range(N2):
            x = 0
            for p, co in zip(ps, residues):
                Mp = Mtot // p
                x = (x + co[k][l] * Mp * pow(Mp, -1, p)) % Mtot
            if x > Mtot // 2: x -= Mtot
            if x: out[(k - D1, l - D2)] = x
    return out


def hull(pts):
    pts = sorted(set(pts))
    if len(pts) <= 2: return pts
    cr = lambda o, a, b: (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])
    lo, up = [], []
    for p in pts:
        while len(lo) >= 2 and cr(lo[-2], lo[-1], p) <= 0: lo.pop()
        lo.append(p)
    for p in reversed(pts):
        while len(up) >= 2 and cr(up[-2], up[-1], p) <= 0: up.pop()
        up.append(p)
    return lo[:-1] + up[:-1]


U = sp.Symbol('u')


def sides(P):
    H = hull(list(P))
    out = []
    for a, b in zip(H, H[1:] + H[:1]):
        d = (b[0] - a[0], b[1] - a[1]); g = math.gcd(abs(d[0]), abs(d[1])); s = (d[0] // g, d[1] // g)
        poly = sum(P.get((a[0] + k * s[0], a[1] + k * s[1]), 0) * U ** k for k in range(g + 1))
        out.append({'from': a, 'dir': s, 'len': g, 'poly': sp.Poly(poly, U)})
    return H, out


def analyse_side(poly):
    c, fl = sp.factor_list(poly.as_expr(), U)
    facs = []
    H0 = False; Hinf = False; reell = True; signs = set()
    for f, m in fl:
        fp = sp.Poly(f, U); dg = fp.degree()
        rts = [complex(r) for r in sp.Poly(f, U).nroots(n=30)] if dg >= 1 else []
        for r in rts:
            if abs(r.imag) > 1e-12 * max(1, abs(r)): reell = False
            else: signs.add(1 if r.real > 0 else -1)
        stab = None
        if dg >= 2:
            H0 = True
            ab = [abs(r) for r in rts]
            stab = (max(ab) - min(ab)) > 1e-9 * max(ab)
            if stab: Hinf = True
        facs.append({'f': str(f), 'deg': dg, 'mult': m, 'stabil': stab, 'rotter': rts})
    return {'factors': facs, 'H0': H0, 'Hinf': Hinf, 'reell': reell, 'tecken': sorted(signs)}


def zonogon(hs):
    """hörnen i Σ[-h/2, h/2] (dubblade koordinater: 2*hörn)"""
    vs = set()
    for sg in itertools.product((1, -1), repeat=len(hs)):
        vs.add((sum(s * h[0] for s, h in zip(sg, hs)), sum(s * h[1] for s, h in zip(sg, hs))))
    return sorted(hull(list(vs)))


def analyse(name, pos, v1, v2, E, W=None, verbose=True):
    if W is None: W = [1] * len(E)
    check_geom(pos, v1, v2, E)
    T = tracks(pos, v1, v2, E)
    mi = minimal_info(E, T)
    P = det_poly(len(pos), E, W)
    H, S = sides(P)
    res = {'name': name, 'V': len(pos), 'E': len(E), 'min': mi, 'hull': H, 'sides': []}
    for s in S:
        a = analyse_side(s['poly'])
        a.update(s); res['sides'].append(a)
    res['H0'] = any(s['H0'] for s in res['sides'])
    res['Hinf'] = any(s['Hinf'] for s in res['sides'])
    res['split'] = not res['H0']
    res['N_eq_Z'] = sorted((2 * x, 2 * y) for x, y in H) == zonogon(mi['h'])
    return res


def report(r):
    mi = r['min']
    L = [f"== {r['name']}: V={r['V']} E={r['E']}  spår={len(mi['h'])}  minimalt={'JA' if mi['minimal'] else 'nej'}"]
    from collections import Counter
    L.append(f"   spårklasser: {dict(Counter(tuple(sorted([h, (-h[0], -h[1])])[1]) for h in mi['h']))}")
    if not mi['minimal']:
        L.append('   icke-minimalt: ' + '; '.join(mi['reasons'][:4]) + (' …' if len(mi['reasons']) > 4 else ''))
    L.append(f"   Newtonpolygon (hörn): {r['hull']}   N = zonogon(spår): {r['N_eq_Z']}")
    for s in r['sides'][: (len(r['sides']) + 1) // 2]:
        fs = ', '.join(f"({f['f']})" + (f"^{f['mult']}" if f['mult'] > 1 else '') for f in s['factors'] if f['deg'] >= 1)
        tag = ''
        if s['H0']: tag = '   <-- IRREDUCIBEL FAKTOR' + (' (stabil på alla övertorusar)' if s['Hinf'] else ' (stabilitet oklar)')
        L.append(f"   sida riktning {s['dir']} längd {s['len']}: {fs}{tag}")
    L.append(f"   GALOIS-HINDER: egen torus {'JA' if r['H0'] else 'nej'}, alla torusar {'JA' if r['Hinf'] else 'nej'}")
    return '\n'.join(L)


def monodromy(T, W):
    """μ_T = ∏_{par=0} c / ∏_{par=1} c (exakt bråk)"""
    from fractions import Fraction
    out = []
    for t in T:
        m = Fraction(1)
        for (k, off, par) in t['edges']:
            m = m * W[k] if par == 0 else m / W[k]
        out.append(m)
    return out
