"""Minimala resistornät från geodetiska tågspår på R^2/Z^2.
Tågspår i: linjen b_i x - a_i y = o_i (mod 1), klass h_i=(a_i,b_i) primitiv.
Område = golvvektor k_i = floor(b_i x - a_i y - o_i) modulo k -> k + A l, A rader (b_i,-a_i).
Svarta områden (sum k jämn) = noder, korsningar = kanter."""
import itertools, random, math
from fractions import Fraction as Fr
import sympy as sp

def bygg(H, O):
    m = len(H); A = [(b, -a) for (a, b) in H]
    assert sum(a for a,b in H) % 2 == 0 and sum(b for a,b in H) % 2 == 0
    # två icke-parallella rader for kanonisering
    i1, i2 = next((i, j) for i, j in itertools.combinations(range(m), 2) if A[i][0]*A[j][1]-A[i][1]*A[j][0] != 0)
    B = ((A[i1][0], A[i1][1]), (A[i2][0], A[i2][1])); d = B[0][0]*B[1][1]-B[0][1]*B[1][0]
    def kanon(k):
        v = (k[i1], k[i2])
        # l = floor(B^{-1} v)
        l0 = math.floor(Fr(B[1][1]*v[0]-B[0][1]*v[1], d)); l1 = math.floor(Fr(-B[1][0]*v[0]+B[0][0]*v[1], d))
        # B^{-1} v: B = [[A_i1],[A_i2]] verkar som k_i = A_i . l
        return tuple(k[i]-(A[i][0]*l0+A[i][1]*l1) for i in range(m)), (l0, l1)
    noder = {}; kanter = []
    for i, j in itertools.combinations(range(m), 2):
        det = A[i][0]*A[j][1]-A[i][1]*A[j][0]
        if det == 0: continue
        R = 4*(abs(A[i][0])+abs(A[i][1])+abs(A[j][0])+abs(A[j][1]))+4
        for p in range(-R, R):
            for q in range(-R, R):
                # A_i.(x,y) = o_i+p, A_j.(x,y)=o_j+q
                ri, rj = O[i]+p, O[j]+q
                x = (ri*A[j][1]-A[i][1]*rj)/det; y = (A[i][0]*rj-A[j][0]*ri)/det
                if not (0 <= x < 1 and 0 <= y < 1): continue
                k = [math.floor(A[l][0]*x+A[l][1]*y-O[l]) for l in range(m)]
                k[i] = p; k[j] = q
                k1 = list(k); k2 = list(k)
                if sum(k) % 2 == 0: k2[i] -= 1; k2[j] -= 1
                else: k1[i] -= 1; k2[j] -= 1
                c1, l1 = kanon(tuple(k1)); c2, l2 = kanon(tuple(k2))
                for c in (c1, c2): noder.setdefault(c, len(noder))
                kanter.append((noder[c1], noder[c2], (l2[0]-l1[0], l2[1]-l1[1])))
    return len(noder), kanter

X, Z, w = sp.symbols('X Z w')
def blochdet(n, kanter, cond):
    M = sp.zeros(n, n)
    for (u, v, (lx, ly)), c in zip(kanter, cond):
        M[u, u] += c; M[v, v] += c
        M[u, v] -= c*X**lx*Z**ly; M[v, u] -= c*X**(-lx)*Z**(-ly)
    P = sp.cancel(M.det(method='berkowitz'))
    num, den = sp.fraction(sp.together(P))
    pd = sp.Poly(den, X, Z); (dx, dz), _ = pd.terms()[0]
    return {(ex-dx, ez-dz): c/pd.coeffs()[0] for (ex, ez), c in sp.Poly(sp.expand(num), X, Z).terms()}

def hull(pts):
    pts = sorted(set(pts))
    def cr(o, a, b): return (a[0]-o[0])*(b[1]-o[1])-(a[1]-o[1])*(b[0]-o[0])
    lo, up = [], []
    for p in pts:
        while len(lo) >= 2 and cr(lo[-2], lo[-1], p) <= 0: lo.pop()
        lo.append(p)
    for p in reversed(pts):
        while len(up) >= 2 and cr(up[-2], up[-1], p) <= 0: up.pop()
        up.append(p)
    return lo[:-1]+up[:-1]

def zonogon(H):
    # Σ [-h/2, h/2]; hörn = 1/2 Σ s_i h_i
    pts = []
    for ang in [2*math.pi*(t+0.5)/720 for t in range(720)]:
        n = (math.cos(ang), math.sin(ang))
        s = [1 if n[0]*a+n[1]*b > 0 else -1 for a, b in H]
        pts.append((sum(si*a for si, (a, b) in zip(s, H))//2, sum(si*b for si, (a, b) in zip(s, H))//2))
    return hull(pts)

def sidpolynom(terms):
    pts = list(terms); hv = hull(pts); ut = []
    for p, q in zip(hv, hv[1:]+hv[:1]):
        dd = (q[0]-p[0], q[1]-p[1]); g = math.gcd(abs(dd[0]), abs(dd[1])); st = (dd[0]//g, dd[1]//g)
        poly = sum(terms.get((p[0]+k*st[0], p[1]+k*st[1]), 0)*w**k for k in range(g+1))
        ut.append((st, g, sp.factor(poly)))
    return ut
