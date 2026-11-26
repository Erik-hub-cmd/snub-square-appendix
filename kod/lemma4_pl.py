"""Separat skrivet test av Lemma 4 (Lemma G.5 i artikeln) för ICKE-geodetiska minimala nät.

Tågspåren är slumpvisa styckvis linjära (PL) slutna kurvor på R^2/Z^2 med givna primitiva klasser h_i,
med stora utböjningar. Konfigurationen behålls bara om den är minimal i GK:s mening (inga självkorsningar,
spår i och j korsas exakt |h_i ^ h_j| gånger, alla ytor är skivor). Medialgrafen, ytorna och schackfärgningen
räknas kombinatoriskt (rotationssystem + ytspårning med gitterförskjutningar), utan att anta något om
Y-Delta-sammanhang. Sedan: Bloch-determinanten exakt (sympy), N(det) mot zonogonen, sidpolynomens faktorer
över Q, och rötterna mot tågspårsmonodromierna (alternerande konduktansprodukter längs varje spår) i nätet själv.

Körning: python3 lemma4_pl.py <frö> <antal> [maxkanter]
"""
import sys, math, random, itertools
import sympy as sp
from fractions import Fraction as Fr

def wedge(a, b): return a[0]*b[1]-a[1]*b[0]

def kurva(h, rng, N, brus):
    L = math.hypot(*h); perp = (-h[1]/L, h[0]/L)
    o = (rng.random(), rng.random())
    pts = []
    for k in range(N):
        s = k/N + rng.uniform(-0.3, 0.3)/N
        a = rng.uniform(-brus, brus); b = rng.uniform(-brus, brus)
        pts.append((o[0]+s*h[0]+a*perp[0]+b*0.2*h[0]/L, o[1]+s*h[1]+a*perp[1]+b*0.2*h[1]/L))
    return pts  # punkt N = pts[0] + h

def segment(pts, h, k):
    N = len(pts); p = pts[k]
    q = pts[k+1] if k+1 < N else (pts[0][0]+h[0], pts[0][1]+h[1])
    return p, q

def skarning(p, q, r, s):
    d1 = (q[0]-p[0], q[1]-p[1]); d2 = (s[0]-r[0], s[1]-r[1])
    den = wedge(d1, d2)
    if abs(den) < 1e-14: return None
    w_ = (r[0]-p[0], r[1]-p[1])
    t = wedge(w_, d2)/den; u = wedge(w_, d1)/den
    if 0 <= t < 1 and 0 <= u < 1: return t, u
    return None

def bygg(H, rng, N=6, brus=0.25):
    m = len(H)
    C = [kurva(h, rng, N, brus) for h in H]
    # korsningar: (i,ki,ti,j,kj,tj,v) där punkten = seg_i(ti) = seg_j(tj)+v
    kors = []
    for i in range(m):
        for j in range(i, m):
            for ki in range(N):
                p, q = segment(C[i], H[i], ki)
                for kj in range(N):
                    r0, s0 = segment(C[j], H[j], kj)
                    for vx in range(-4, 5):
                        for vy in range(-4, 5):
                            if i == j and (vx, vy) == (0, 0): continue
                            r = (r0[0]+vx, r0[1]+vy); s = (s0[0]+vx, s0[1]+vy)
                            st = skarning(p, q, r, s)
                            if st is None: continue
                            if i == j: return None  # självkorsning (också mellan olika lyft)
                            kors.append((i, ki, st[0], j, kj, st[1], (vx, vy)))
    # minimalitet: antal korsningar per par
    MINIMAL = '--icke-minimal' not in sys.argv
    extra = 0
    for i, j in itertools.combinations(range(m), 2):
        d = sum(1 for c in kors if c[0] == i and c[3] == j) - abs(wedge(H[i], H[j]))
        extra += d
        if MINIMAL and d != 0: return None
    if not MINIMAL and extra == 0: return None
    nx = len(kors)
    # förekomster längs varje spår: (param, korsning x, delta = lyftförskjutning relativt pos(x), riktning)
    pos = []; forek = [[] for _ in range(m)]
    for x, (i, ki, ti, j, kj, tj, v) in enumerate(kors):
        p, q = segment(C[i], H[i], ki); r, s = segment(C[j], H[j], kj)
        pos.append((p[0]+ti*(q[0]-p[0]), p[1]+ti*(q[1]-p[1])))
        forek[i].append((ki+ti, x, (0, 0), (q[0]-p[0], q[1]-p[1])))
        forek[j].append((kj+tj, x, (-v[0], -v[1]), (s[0]-r[0], s[1]-r[1])))
    # halvkanter: vid x, fyra utgående: (x, spår, fram/bak). Granne och förskjutning.
    ut = {}  # (x, idx) -> (y, idx_tillbaka, offset) ; offset = (lyft av y) - ... i "pos"-koordinater
    rot = {x: [] for x in range(nx)}
    for k in range(m):
        F = sorted(forek[k]); n_ = len(F)
        for a in range(n_):
            pa, x, da, dirx = F[a]; pb, y, db, diry = F[(a+1) % n_]
            wrap = H[k] if a+1 == n_ else (0, 0)
            # punkt x i lyftet = pos(x)+da ; y = pos(y)+db+wrap. offset för y relativt x:s lyft:
            off = (db[0]+wrap[0]-da[0], db[1]+wrap[1]-da[1])
            ut[(x, k, +1)] = (y, (k, -1), off)
            ut[(y, k, -1)] = (x, (k, +1), (-off[0], -off[1]))
            rot[x].append(((k, +1), math.atan2(dirx[1], dirx[0])))
            rot[y].append(((k, -1), math.atan2(-diry[1], -diry[0])))
    ordn = {}
    for x in range(nx):
        L = sorted(rot[x], key=lambda e: e[1]); assert len(L) == 4
        ordn[x] = [e[0] for e in L]  # moturs
    # ytspårning: halvkant (x,e). nästa i ytan: gå till y, ta halvkanten före tillbakakanten moturs (=medurs nästa)
    sett = {}; ytor = []
    for x in range(nx):
        for e in ordn[x]:
            if (x, e) in sett: continue
            f = len(ytor); hörn = []; cx, ce, o = x, e, (0, 0)
            while (cx, ce) not in sett:
                sett[(cx, ce)] = f
                y, eb, off = ut[(cx, ce[0], ce[1])]
                o = (o[0]+off[0], o[1]+off[1])
                L = ordn[y]; i_ = L.index(eb); ne = L[(i_-1) % 4]
                hörn.append((y, o, eb, ne))  # hörnet vid y mellan eb och ne, y ligger på pos(y)+o
                cx, ce = y, ne
            if o != (0, 0): return None  # ytan är inte en skiva (icke-kontraktibel rand)
            ytor.append(hörn)
    # Euler: V - E + F = 0 på torusen
    if nx - 2*nx + len(ytor) != 0: return None
    # färgning: två ytor över en medialkant har olika färg
    ytaav = sett
    adj = {f: set() for f in range(len(ytor))}
    for (x, e), f in sett.items():
        y, eb, off = ut[(x, e[0], e[1])]
        g = sett[(y, eb)]
        adj[f].add(g); adj[g].add(f)
    farg = {0: 0}; stack = [0]
    while stack:
        f = stack.pop()
        for g in adj[f]:
            if g not in farg: farg[g] = 1-farg[f]; stack.append(g)
            elif farg[g] == farg[f]: return None
    svarta = [f for f in range(len(ytor)) if farg[f] == 0]
    nod = {f: a for a, f in enumerate(svarta)}
    # kanter: vid varje korsning x, de två svarta hörnen
    hornvid = {x: [] for x in range(nx)}
    for f, H_ in enumerate(ytor):
        for (y, o, eb, ne) in H_:
            hornvid[y].append((f, o))
    kanter = []; kant_vid = {}
    for x in range(nx):
        sv = [(f, o) for f, o in hornvid[x] if farg[f] == 0]
        assert len(sv) == 2, sv
        (f1, o1), (f2, o2) = sv
        kant_vid[x] = len(kanter)
        kanter.append((nod[f1], nod[f2], (o1[0]-o2[0], o1[1]-o2[1])))
    # korsningsordning längs varje spår (för monodromin)
    sparordn = [[x for _, x, _, _ in sorted(forek[k])] for k in range(m)]
    # icke-geodetisk kombinatorik: på ett rakt spår k är korsningarna med j och l två aritmetiska följder,
    # så antalet l-korsningar mellan två på varandra följande j-korsningar varierar högst 1.
    partner = {x: (kors[x][0], kors[x][3]) for x in range(nx)}
    global ICKEGEO
    ICKEGEO = False
    for k in range(m):
        seq = [(set(partner[x]) - {k}).pop() for x in sparordn[k]]
        for j in set(seq):
            for l in set(seq) - {j}:
                idx = [t for t, q in enumerate(seq) if q == j]
                L = len(seq); cnt = []
                for a in range(len(idx)):
                    b0 = idx[a]; b1 = idx[(a+1) % len(idx)] + (L if a+1 == len(idx) else 0)
                    cnt.append(sum(1 for t in range(b0+1, b1) if seq[t % L] == l))
                if max(cnt) - min(cnt) >= 2: ICKEGEO = True
    return len(svarta), kanter, kant_vid, sparordn

X, Z, w = sp.symbols('X Z w')
def blochdet(n, kanter, cond):
    M = sp.zeros(n, n)
    for (u, v, (lx, ly)), c in zip(kanter, cond):
        M[u, u] += c; M[v, v] += c
        M[u, v] -= c*X**lx*Z**ly; M[v, u] -= c*X**(-lx)*Z**(-ly)
    P = sp.together(M.det(method='berkowitz'))
    num, den = sp.fraction(P)
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
    pts = []
    for t in range(720):
        ang = 2*math.pi*(t+0.5)/720; n = (math.cos(ang), math.sin(ang))
        s = [1 if n[0]*a+n[1]*b > 0 else -1 for a, b in H]
        pts.append((sum(si*a for si, (a, b) in zip(s, H))//2, sum(si*b for si, (a, b) in zip(s, H))//2))
    return hull(pts)

def prim(d):
    g = math.gcd(abs(d[0]), abs(d[1])); return (d[0]//g, d[1]//g), g

def main():
    seed = int(sys.argv[1]); antal = int(sys.argv[2]); maxk = int(sys.argv[3]) if len(sys.argv) > 3 else 14
    rng = random.Random(seed)
    DIRS = [(1,0),(0,1),(1,1),(1,-1),(2,1),(1,2),(2,-1),(1,-2),(3,1),(1,3),(3,2),(2,3)]
    ok = True; gjort = 0; forsok = 0; ickegeo = 0
    while gjort < antal:
        forsok += 1
        m = rng.randint(2, 6)
        H = [rng.choice(DIRS) for _ in range(m)]
        if sum(a for a, b in H) % 2 or sum(b for a, b in H) % 2: continue
        if len({prim(h)[0] if (h[0] > 0 or (h[0] == 0 and h[1] > 0)) else prim((-h[0], -h[1]))[0] for h in H}) < 2: continue
        E = sum(abs(wedge(H[i], H[j])) for i in range(m) for j in range(i+1, m))
        if E > maxk: continue
        res = bygg(H, rng, N=rng.randint(3, 7), brus=rng.choice([0.1, 0.2, 0.3, 0.45]))
        if res is None: continue
        n, K, kant_vid, sparordn = res
        if '--bara-ickegeo' in sys.argv and not ICKEGEO: continue
        if n > 9: continue
        cond = [sp.Rational(rng.randint(1, 12), rng.randint(1, 12)) for _ in K]
        T = blochdet(n, K, cond)
        N = hull(list(T)); Zg = zonogon(H)
        likaN = sorted(N) == sorted(Zg)
        # monodromier: för spår k, alternerande produkter längs korsningsordningen
        mono = {}
        for k in range(m):
            xs = sparordn[k]
            a = sp.prod([cond[kant_vid[x]] for x in xs[0::2]]); b = sp.prod([cond[kant_vid[x]] for x in xs[1::2]])
            mono[k] = a/b
        sidok = True; monook = True; info = []
        for p, q in zip(N, N[1:]+N[:1]):
            st, g = prim((q[0]-p[0], q[1]-p[1]))
            poly = sp.Poly(sum(T.get((p[0]+k*st[0], p[1]+k*st[1]), 0)*w**k for k in range(g+1)), w)
            fl = sp.factor_list(poly.as_expr(), w)[1]
            lin = all(sp.Poly(f, w).degree() <= 1 for f, _ in fl)
            sidok &= lin
            par = [k for k in range(m) if wedge(H[k], st) == 0]
            if len(par) != g: sidok = False
            rotter = sorted([abs(r) for r in sp.roots(poly, w, multiple=True)])
            # jämn/udda-märkningen längs varje spår är godtycklig: pröva m_k eller 1/m_k per spår
            mm = any(sorted([abs(mono[k]) if e else abs(1/mono[k]) for k, e in zip(par, val)]) == rotter
                     for val in itertools.product([0, 1], repeat=len(par)))
            monook &= mm
            info.append((st, g, lin, mm))
        geo = all(len(xs) for xs in sparordn)
        bra = likaN and sidok and len(K) == E
        ok &= bra
        if '--icke-minimal' in sys.argv: E = len(K)
        gjort += 1
        print(f'H={H} noder={n} kanter={len(K)} (pred {E}) N=Z:{likaN} sidor linjära:{sidok} rötter=monodromi(i nätet självt):{monook} icke-geodetisk kombinatorik:{ICKEGEO}', flush=True)
        if not bra: print('   AVVIKELSE', N, Zg, info)
    print('ALLT OK' if ok else 'AVVIKELSE', f'({gjort} nät, {forsok} försök)')

if __name__ == '__main__':
    main()
