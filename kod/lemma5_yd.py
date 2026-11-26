"""Separat skriven kontroll av Lemma 5 (Anm. G.6 i artikeln).

Bygger hexagonala gittret på en k x k-torus med slumpvisa rationella konduktanser, gör slumpvisa Y->Delta-drag
(rotationssystemet uppdateras kombinatoriskt), och kontrollerar
  (1) att den kanoniska monodromin mu_T (konduktanser vid L->R-korsningar delat med R->L-korsningar) är oförändrad
      för varje tågspår efter varje drag, och att antalet korsningar per spår är jämnt;
  (2) att varje sidpolynom till det Delta(X,Z) är konst * prod_T (u + eps*mu_T^sigma) över spåren T parallella med
      sidan, med sigma, eps gemensamma för sidan.
Determinanten räknas exakt genom evaluering + interpolation (python-flint).
Körning: python3 lemma5_yd.py <frö> <k> <antal drag> <antal nät>
"""
import sys, random, math
from fractions import Fraction as Fr
import flint

class Nat:
    def __init__(s):
        s.he = []      # halvkant: [nod, motstående halvkant, förskjutning (dx,dy), kant-id]
        s.rot = {}     # nod -> cyklisk lista av halvkanter (moturs)
        s.c = {}       # kant-id -> konduktans
        s.nid = 0
    def kant(s, u, v, sh, c):
        e = s.nid; s.nid += 1
        a = len(s.he); s.he.append([u, a+1, sh, e]); s.he.append([v, a, (-sh[0], -sh[1]), e])
        s.c[e] = c
        return a, a+1

def hexnat(k, rng):
    N = Nat()
    e1 = (1.0, 0.0); e2 = (0.5, math.sqrt(3)/2)
    B = lambda i, j: ('B', i % k, j % k); W = lambda i, j: ('W', i % k, j % k)
    pos = {}
    tmp = {}
    for i in range(k):
        for j in range(k):
            for (di, dj) in [(0, 0), (-1, 0), (0, -1)]:
                ii, jj = i+di, j+dj
                sh = (ii//k - 0, jj//k - 0)  # floor-division ger förskjutning
                sh = (math.floor(ii/k), math.floor(jj/k))
                c = Fr(rng.randint(1, 9), rng.randint(1, 9))
                a, b = N.kant(B(i, j), W(ii, jj), sh, c)
                vec = ((e1[0]+e2[0])/3 + di*e1[0] + dj*e2[0], (e1[1]+e2[1])/3 + di*e1[1] + dj*e2[1])
                tmp.setdefault(B(i, j), []).append((math.atan2(vec[1], vec[0]), a))
                tmp.setdefault(W(ii, jj), []).append((math.atan2(-vec[1], -vec[0]), b))
    for v, L in tmp.items():
        N.rot[v] = [h for _, h in sorted(L)]
    return N

def nasta(N, h, d):
    v = N.he[h][0]; R = N.rot[v]; i = R.index(h)
    return R[(i + d) % len(R)]

def spar(N):
    """Orienterade tågspår som listor av (kant, typ, förskjutning); typ +1 = L->R."""
    sett = set(); ut = []
    for h0 in range(len(N.he)):
        for d0 in (1, -1):
            if (h0, d0) in sett or N.he[h0][3] not in N.c: continue
            st = (h0, d0); seq = []; cls = [0, 0]
            while st not in sett:
                sett.add(st); h, d = st
                h2 = nasta(N, h, d)
                seq.append((N.he[h2][3], d)); cls[0] += N.he[h2][2][0]; cls[1] += N.he[h2][2][1]
                st = (N.he[h2][1], -d)
            if st != (h0, d0): raise RuntimeError('spåret sluts inte i startläget')
            ut.append((seq, tuple(cls)))
    # varje geometriskt spår förekommer två gånger (båda orienteringarna): ta bort dubbletter
    sedda = {}; res = []
    for seq, cls in ut:
        kanter = tuple(sorted(e for e, _ in seq))
        # samma kantmultimängd + motsatt klass = samma spår baklänges
        k1 = (kanter, cls); k2 = (kanter, (-cls[0], -cls[1]))
        if k2 in sedda and sedda[k2] > 0:
            sedda[k2] -= 1; continue
        sedda[k1] = sedda.get(k1, 0) + 1
        res.append((seq, cls))
    return res

def mu(N, seq):
    m = Fr(1)
    for e, t in seq:
        m = m * N.c[e] if t == 1 else m / N.c[e]
    return m

def y_till_delta(N, o):
    R = N.rot[o]
    assert len(R) == 3
    legs = []
    for h in R:  # moturs: a, b, c
        ho = N.he[h][1]
        legs.append((N.he[ho][0], ho, N.he[h][2], N.c[N.he[h][3]]))
    (a, ha, sa, x), (b, hb, sb, y), (cn, hc, sc, z) = legs
    if len({a, b, cn}) < 3 or o in (a, b, cn):
        return False   # håll det enkelt: bara Y med tre olika grannar
    S = x + y + z
    sub = lambda p, q: (p[0]-q[0], p[1]-q[1])
    ab = N.kant(a, b, sub(sb, sa), x*y/S)
    bc = N.kant(b, cn, sub(sc, sb), y*z/S)
    ca = N.kant(cn, a, sub(sa, sc), z*x/S)
    # vid a: ersätt (a->o) med [a->b, a->c] moturs; osv.
    def ers(v, hold, h1, h2):
        L = N.rot[v]; i = L.index(hold); N.rot[v] = L[:i] + [h1, h2] + L[i+1:]
    ers(a, ha, ab[0], ca[1]); ers(b, hb, bc[0], ab[1]); ers(cn, hc, ca[0], bc[1])
    for h in R:
        N.c.pop(N.he[h][3], None)
    del N.rot[o]
    return True

def bloch(N, X, Z):
    nod = sorted(N.rot); idx = {v: i for i, v in enumerate(nod)}; n = len(nod)
    M = [[flint.fmpq(0)]*n for _ in range(n)]
    for h in range(0, len(N.he), 2):
        u, h2, sh, e = N.he[h]
        if e not in N.c: continue
        v = N.he[h2][0]; c = flint.fmpq(N.c[e].numerator, N.c[e].denominator)
        w = X**sh[0] * Z**sh[1]
        i, j = idx[u], idx[v]
        M[i][i] += c; M[j][j] += c; M[i][j] -= c*w; M[j][i] -= c/w
    return flint.fmpq_mat(M).det()

def interp1(xs, ys):
    # Newton -> koefficientlista
    n = len(xs); a = list(ys)
    for j in range(1, n):
        for i in range(n-1, j-1, -1):
            a[i] = (a[i]-a[i-1])/(xs[i]-xs[i-j])
    p = flint.fmpq_poly([0])
    for i in range(n-1, -1, -1):
        p = p*flint.fmpq_poly([-xs[i], 1]) + a[i]
    return [p[i] for i in range(n)]

def detpoly(N):
    Dx = sum(abs(N.he[h][2][0]) for h in range(0, len(N.he), 2) if N.he[h][3] in N.c)
    Dz = sum(abs(N.he[h][2][1]) for h in range(0, len(N.he), 2) if N.he[h][3] in N.c)
    xs = [flint.fmpq(i+2) for i in range(2*Dx+1)]; zs = [flint.fmpq(i+2) for i in range(2*Dz+1)]
    tab = [[bloch(N, X, Z) * X**Dx * Z**Dz for Z in zs] for X in xs]
    rad = [interp1(zs, tab[i]) for i in range(len(xs))]
    P = {}
    for b in range(len(zs)):
        col = interp1(xs, [rad[i][b] for i in range(len(xs))])
        for a, v in enumerate(col):
            if v != 0: P[(a-Dx, b-Dz)] = v
    return P

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
    return lo[:-1] + up[:-1]

def kontroll_sidor(N, P, S):
    H = hull(list(P)); ok = True; rapport = []
    for i in range(len(H)):
        p, q = H[i], H[(i+1) % len(H)]
        g = math.gcd(q[0]-p[0], q[1]-p[1]); s = ((q[0]-p[0])//g, (q[1]-p[1])//g)
        coef = [P.get((p[0]+j*s[0], p[1]+j*s[1]), flint.fmpq(0)) for j in range(g+1)]
        par = [mu(N, seq) for seq, cls in S if cls == s or cls == (-s[0], -s[1])]
        if len(par) != g: ok = False; rapport.append(('fel antal spår', s, g, len(par))); continue
        lyckad = None
        for sig in (1, -1):
            for eps in (1, -1):
                Q = flint.fmpq_poly([1])
                for m in par:
                    r = m**sig
                    Q = Q * flint.fmpq_poly([flint.fmpq(r.numerator, r.denominator)*eps, 1])
                Pc = flint.fmpq_poly(coef)
                if Pc * Q[g] == Q * Pc[g]: lyckad = (sig, eps)
        if lyckad is None: ok = False; rapport.append(('sidpolynom ≠ produkt', s))
        else: rapport.append((s, g, lyckad))
    return ok, rapport

def main():
    fro, k, drag, antal = map(int, sys.argv[1:5])
    rng = random.Random(fro); fel = 0; ndrag = 0
    for t in range(antal):
        N = hexnat(k, rng)
        S0 = spar(N)
        for _ in range(drag):
            kand = [v for v in N.rot if len(N.rot[v]) == 3]
            rng.shuffle(kand)
            mu_fore = sorted(mu(N, seq) for seq, _ in spar(N))
            cls_fore = sorted(max(cls, (-cls[0], -cls[1])) for _, cls in spar(N))
            gjort = False
            for o in kand:
                if y_till_delta(N, o): gjort = True; break
            if not gjort: break
            ndrag += 1
            S1 = spar(N)
            if any(len(seq) % 2 for seq, _ in S1): fel += 1; print('udda spår')
            if sorted(mu(N, seq) for seq, _ in S1) != mu_fore or sorted(max(c, (-c[0], -c[1])) for _, c in S1) != cls_fore:
                fel += 1; print('mu ändrades', t)
        P = detpoly(N)
        ok, rap = kontroll_sidor(N, P, spar(N))
        print(t, 'noder', len(N.rot), 'kanter', len(N.c), 'spår', len(spar(N)), 'sidor ok' if ok else 'FEL', rap)
        fel += (not ok)
    print('drag', ndrag, 'fel', fel)

main()
