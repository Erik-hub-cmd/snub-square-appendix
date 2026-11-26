"""
Separat skriven kontroll (2026-10-18; arbetsanteckningen GRANSKNING2.md publiceras inte).

Skriptet importerar INGEN projektkod (inte exact_snub, residue_tools, k24, cairo_lattice ...). Allt byggs från grunden:
  * snub square-gittret ur geometrin: enhetskvadrater centrerade i (0,0) och (L/2, L/2), L = sqrt(2+sqrt3), vridna ±alpha;
    kanter = alla nodpar på avstånd 1. Sedan speglas y -> -y så att numrering och cellkoordinater blir projektets
    (de harmoniska lägena blir exakt d_j i forslag/asymptotik.md).
  * Cairo-gittret som dualgrafen (snub-ytor -> noder, gemensam kant -> kant).
  * Resistanser direkt ur Ekv. (1): θ2-integralen med numeriska residyer i Z (det och adj som Laurentpolynom via FFT,
    rötter numeriskt; för θ1 < 0,05 mpmath med 40 siffror eftersom rötterna nästan sammanfaller vid Z = 1),
    θ1-integralen med Gauss-Legendre på ett nät som förfinas mot 0 och π. Ingen X(s), Z_in(s), inga partialbråk.

Kontroller:
  K1  geometrin: alpha = 15°, 10 kanter per cell, grad 5; Cairo 4 noder av grad 3, 2 av grad 4.
  K2  det M_snub = 4D(cos θ1, cos θ2) (sympy, exakt) och det L_Cairo = det M_snub (slumpade θ).
  K3  Z-ordningarna i adj M = kontroll D1, och alla exponentlådor innehåller 0 (villkoret för width() i
      bevis_formodan8_direkt.py), w_adj = 20.
  K4  harmoniska lägen = d_j; Σ e e^T = (4/3) I (snub) och (3/4) I (Cairo).
  K5  P3: θ1 = arg X(s) växer 0 -> π på (1, φ) och -π -> 0 på (φ, ∞); Z_in(s) = b+- resp. b--roten.
  K6  kantresistanser: snub 1/2 - √3/18, 2√3/9; Cairo 1 - 2√3/9, 1/2 + √3/18.
  K7  40 slumpvalda värden i runs/tabell/checkpoint.json (läses bara som data).
  K8  sandhög/spännträd: P1 och gradfördelning för snub square (mot Goodchild och sandhog_goodchild.md) och för båda
      Cairo-nodtyperna, ur egna numeriska resistanser.
  K9  Cairo-stickprov: R(4-4) = 4√3/9 och PSLQ-identifiering i Q + Q√3 + Q/π (stickprov, inget bevis).
  K10 asymptotiken vid |n| ≈ 20, 40, 80 (64 par per radie, alla nodtyper, åtta riktningar):
      konstanten C, felet efter Sats 2 är O(|x̃|^-3), sin 2φ-amplituden 1/(48π), ingen 1/r²-term för blandade par,
      ingen cos 4φ-term, och felet med geometriskt avstånd är O(1/r).

Körtid: ca 7 minuter (3 processer, låg prioritet; 431 s vid första körningen 2026-10-18, 14 godkända).  Kör:  .venv\\Scripts\\python kod\\granskning2_kontroll.py
"""
import itertools, json, os, random, sys, time
sys.stdout.reconfigure(encoding="utf-8")          # utskrift till fil/pipe i Windows (cp1252) kraschar annars
import numpy as np
import mpmath as mp

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
S3 = np.sqrt(3)
L = np.sqrt(2 + S3)
GAMMA = 0.5772156649015329
C_ASYM = S3 / 36 + 3 / (4 * np.pi) * (GAMMA + np.log(6 - 2 * S3))
EPS = [1, -1, 1, -1]


# ================================================================== gittren
def square_vertices(c, ang):
    return [np.array(c) + np.array([np.cos(ang + np.pi / 4 + k * np.pi / 2), np.sin(ang + np.pi / 4 + k * np.pi / 2)]) / np.sqrt(2)
            for k in range(4)]


def mismatch(alpha):
    v0 = square_vertices((0, 0), alpha)
    v1 = square_vertices((L / 2, L / 2), -alpha)
    return min(np.linalg.norm(a - b) for a in v0 for b in v1)


def build_snub():
    """(alpha, V, E) i projektets konvention (efter spegling y -> -y)"""
    from scipy.optimize import minimize_scalar
    al = np.linspace(-np.pi / 4, np.pi / 4, 20001)
    a0 = al[np.argmin([mismatch(a) for a in al])]
    alpha = minimize_scalar(mismatch, bracket=(a0 - 1e-4, a0, a0 + 1e-4), tol=1e-14).x
    assert mismatch(alpha) < 1e-9
    red = []
    for p in square_vertices((0, 0), alpha) + square_vertices((L / 2, L / 2), -alpha):
        q = p - L * np.floor(p / L + 1e-12)
        if not any(np.linalg.norm(((q - r) / L + 0.5) % 1 - 0.5) < 1e-9 for r in red):
            red.append(q)
    V = np.array(red)
    E = []
    for i, j in itertools.product(range(len(V)), repeat=2):
        for m1, m2 in itertools.product(range(-2, 3), repeat=2):
            if (i, j, m1, m2) <= (j, i, -m1, -m2) and not (i == j and m1 == m2 == 0):
                if abs(np.linalg.norm(V[j] + L * np.array([m1, m2]) - V[i]) - 1) < 1e-9:
                    E.append((i, j, (m1, m2)))
    # spegling y -> -y: då blir numreringen och cellkoordinaterna projektets
    V = np.array([[x, -y] for x, y in V])
    E = [(i, j, (m1, -m2)) for i, j, (m1, m2) in E]
    # positioner måste stämma med kantvektorerna efter speglingen
    for i, j, m in E:
        assert abs(np.linalg.norm(V[j] + L * np.array(m) - V[i]) - 1) < 1e-9
    return alpha, V, E


def neighbours(E, v):
    out = []
    for i, j, m in E:
        if i == v:
            out.append((j, m))
        if j == v:
            out.append((i, (-m[0], -m[1])))
    return out


def build_cairo(V, E):
    """dualgrafen: (centroider, grader, kanter)"""
    nb = {v: neighbours(E, v) for v in range(4)}
    pos = lambda v, m: V[v] + L * np.array(m)
    ang = lambda d: np.arctan2(d[1], d[0])
    faces = {}
    for v in range(4):
        for b in nb[v]:
            cyc, prev, cur = [(v, (0, 0))], (v, (0, 0)), b
            while cur != (v, (0, 0)):
                cyc.append(cur)
                cand = [(w, (cur[1][0] + mm[0], cur[1][1] + mm[1])) for w, mm in nb[cur[0]]]
                pc = pos(*cur)
                ain = ang(pos(*prev) - pc)
                prev, cur = cur, min((c for c in cand if c != prev), key=lambda c: (ain - ang(pos(*c) - pc)) % (2 * np.pi))
                assert len(cyc) < 10
            cen = np.mean([pos(*c) for c in cyc], axis=0)
            cr = cen - L * np.floor(cen / L + 1e-12)
            faces.setdefault(tuple(np.round(cr, 9)), (len(cyc), cr))
    keys = list(faces)
    cen = np.array([faces[k][1] for k in keys])
    deg = [faces[k][0] for k in keys]
    CE = []
    for i, j, m in E:
        a, b = V[i], V[j] + L * np.array(m)
        mid, t = (a + b) / 2, b - a
        nrm = np.array([-t[1], t[0]])
        sides = []
        for sgn in (1, -1):
            best = None
            for f in range(len(cen)):
                for m1, m2 in itertools.product(range(-2, 3), repeat=2):
                    c = cen[f] + L * np.array([m1, m2])
                    if sgn * np.dot(c - mid, nrm) > 1e-9 and abs(np.dot(c - mid, t)) < 0.5:
                        dd = np.linalg.norm(c - mid)
                        if best is None or dd < best[0]:
                            best = (dd, f, (m1, m2))
            sides.append(best[1:])
        (f1, c1), (f2, c2) = sides
        CE.append((f1, f2, (c2[0] - c1[0], c2[1] - c1[1])))
    cdeg = [sum((e[0] == v) + (e[1] == v) for e in CE) for v in range(len(cen))]
    assert cdeg == deg
    return cen, deg, CE


def bloch(nv, E, deg, t1, t2):
    M = np.diag(np.array(deg, dtype=complex))
    for i, j, (m1, m2) in E:
        ph = np.exp(1j * (t1 * m1 + t2 * m2))
        M[i, j] -= ph
        M[j, i] -= np.conj(ph)
    return M


def harmonic(nv, E, V):
    """harmoniska lägen (geometriska enheter, h_0 = V_0) och Σ e e^T i cellenheter"""
    A = np.zeros((nv, nv))
    b = np.zeros((nv, 2))
    for i, j, m in E:
        mv = L * np.array(m, float)
        A[i, i] -= 1; A[i, j] += 1; b[i] -= mv
        A[j, j] -= 1; A[j, i] += 1; b[j] += mv
    h = np.vstack([V[0], np.linalg.solve(A[1:, 1:], b[1:] - np.outer(A[1:, 0], V[0]))])
    T = np.zeros((2, 2))
    for i, j, m in E:
        e = (h[j] + L * np.array(m) - h[i]) / L
        T += np.outer(e, e)
    return h, T


# ================================================================== resistanser ur Ekv. (1)
NF = 16          # FFT-punkter på |Z| = 1 (Laurentgraden i Z är <= 4)
OFF = 4          # Z^OFF · (Laurentpolynom) blir ett vanligt polynom


def _sum_terms(inner, pv_num, pv_den, Pq, Q, terms, sh, t1, expj, zero):
    """Σ residyer i de inre rötterna (samlad täljare) + eventuell pol i Z = 0"""
    tot = zero
    for z in inner:
        num = zero
        for (r, s, m1, m2, w) in terms:
            num += w * expj(m1 * t1) * z**(m2 - OFF + sh - 1) * pv_num(Q(r, s), z)
        tot += num / pv_den(z)
    for (r, s, m1, m2, w) in terms:
        ex = m2 - OFF + sh - 1
        if ex < 0:                                   # koefficienten för Z^(-ex-1) i Q/Pq
            Qs, ser = Q(r, s), []
            for i in range(-ex):
                acc = Qs[i] if i < len(Qs) else 0
                for j in range(1, min(i, len(Pq) - 1) + 1):
                    acc -= Pq[j] * ser[i - j]
                ser.append(acc / Pq[0])
            tot += w * expj(m1 * t1) * ser[-ex - 1]
    return tot


def _terms(k, l, n1, n2):
    """F = G_kk + G_ll - 2 Re(G_lk e^{iθn}); använd G_kl e^{-iθn} om n2 < 0 (hermitesk)"""
    a, b, e1, e2 = (l, k, n1, n2) if n2 >= 0 else (k, l, -n1, -n2)
    return ((k, k, 0, 0, 1), (l, l, 0, 0, 1), (a, b, e1, e2, -2))


def theta2_np(nv, E, deg, t1, pairs):
    zs = np.exp(2j * np.pi * np.arange(NF) / NF)
    Ms = np.array([bloch(nv, E, deg, t1, np.angle(z)) for z in zs])
    dets = np.linalg.det(Ms)
    adjs = np.linalg.inv(Ms) * dets[:, None, None]
    ks = np.arange(NF)
    ks = np.where(ks < NF // 2, ks, ks - NF)

    def poly(c):
        out = np.zeros(NF, dtype=complex)
        for k, v in zip(ks, c):
            if 0 <= k + OFF < NF:
                out[k + OFF] += v
        return out
    P = poly(np.fft.fft(dets) / NF)
    P[np.abs(P) < 1e-11 * np.abs(P).max()] = 0
    P = np.trim_zeros(P, 'b')
    lo = int(np.argmax(np.abs(P) > 0))
    Pq = P[lo:]
    ca = np.fft.fft(adjs, axis=0) / NF
    cache = {}
    Q = lambda r, s: cache.setdefault((r, s), poly(ca[:, r, s]))
    roots = np.roots(Pq[::-1])
    inner = roots[np.abs(roots) < 1]
    dP = np.polyder(Pq[::-1])
    return np.array([_sum_terms(inner, lambda q, z: np.polyval(q[::-1], z), lambda z: np.polyval(dP, z), Pq, Q,
                                _terms(*p), OFF - lo, t1, lambda x: np.exp(1j * x), 0j).real for p in pairs])


def theta2_mp(nv, E, deg, t1, pairs, dps=40):
    with mp.workdps(dps):
        t1 = mp.mpf(t1)
        X = mp.expj(t1)
        zs = [mp.expjpi(mp.mpf(2 * k) / NF) for k in range(NF)]
        dets, adjs = [], []
        for z in zs:
            M = mp.matrix(nv, nv)
            for v in range(nv):
                M[v, v] = deg[v]
            for i, j, (m1, m2) in E:
                ph = X**m1 * z**m2
                M[i, j] -= ph
                M[j, i] -= mp.conj(ph)
            d = mp.det(M)
            dets.append(d)
            adjs.append(mp.inverse(M) * d)
        four = lambda vals, k: mp.fsum(vals[q] * zs[q]**(-k) for q in range(NF)) / NF
        P = [four(dets, k - OFF) for k in range(NF)]
        mx = max(abs(c) for c in P)
        P = [c if abs(c) > mp.mpf(10)**(-dps + 8) * mx else mp.mpf(0) for c in P]
        while P[-1] == 0:
            P.pop()
        lo = next(i for i, c in enumerate(P) if c != 0)
        Pq = P[lo:]
        inner = [r for r in mp.polyroots(Pq[::-1], maxsteps=200, extraprec=4 * dps) if abs(r) < 1]
        dP = [k * Pq[k] for k in range(1, len(Pq))]
        pv = lambda cs, z: mp.fsum(c * z**i for i, c in enumerate(cs))
        cache = {}
        Q = lambda r, s: cache[(r, s)] if (r, s) in cache else cache.setdefault(
            (r, s), [four([A[r, s] for A in adjs], k - OFF) for k in range(NF)])
        return np.array([float(mp.re(_sum_terms(inner, pv, lambda z: pv(dP, z), Pq, Q, _terms(*p), OFF - lo, t1,
                                                 mp.expj, mp.mpf(0)))) for p in pairs])


def theta1_nodes(nmax):
    """Gauss-Legendre (24 punkter) på [0, π], geometriskt förfinat mot 0 och π, panelbredd <= 2/nmax"""
    xg, wg = np.polynomial.legendre.leggauss(24)
    br = [0.0] + [np.pi / 2 * 2.0**-k for k in range(40, 0, -1)] + [np.pi / 2] \
        + [np.pi - np.pi / 2 * 2.0**-k for k in range(1, 41)] + [np.pi]
    br = sorted(set(br))
    fine = []
    for a, b in zip(br[:-1], br[1:]):
        m = max(1, int(np.ceil((b - a) * max(nmax, 4) / 2.0)))
        fine += list(np.linspace(a, b, m + 1)[:-1])
    fine.append(np.pi)
    t, w = [], []
    for a, b in zip(fine[:-1], fine[1:]):
        t += list((b - a) / 2 * xg + (a + b) / 2)
        w += list((b - a) / 2 * wg)
    return np.array(t), np.array(w)


def resist(nv, E, deg, pairs):
    """R((k,0), (l,n)) för pairs = [(k, l, n1, n2)];  R = (1/π) ∫_0^π h(θ1) dθ1 (jämnhet)"""
    t, w = theta1_nodes(max(abs(p[2]) + abs(p[3]) for p in pairs) + 1)
    acc = np.zeros(len(pairs))
    for ti, wi in zip(t, w):
        acc += wi * (theta2_mp if ti < 0.05 else theta2_np)(nv, E, deg, ti, pairs)
    return acc / np.pi


# ================================================================== deluppgifter (körs parallellt)
def task_tabell():
    path = os.path.join(ROOT, 'runs', 'tabell', 'checkpoint.json')
    if not os.path.exists(path):
        return None
    TAB = {}
    for r in json.load(open(path, encoding='utf-8'))['results'].values():
        TAB.update(r)
    random.seed(7)
    pick = random.sample([k for k in sorted(TAB) if k.split(':')[1:] != ['0', '0']], 40)
    _, V, E = build_snub()
    pairs = [tuple([0] + list(map(int, k.split(':')))) for k in pick]
    R = resist(4, E, [5] * 4, pairs)
    dR = max(abs(r - float(TAB[k]['float'])) for k, r in zip(pick, R))
    dd = max(abs(np.linalg.norm(V[p[1]] + L * np.array(p[2:]) - V[0]) - TAB[k]['dist']) for k, p in zip(pick, pairs))
    return dict(dR=dR, dd=dd, maxdist=max(TAB[k]['dist'] for k in pick))


def ust(nv, E, deg, v):
    nb = neighbours(E, v)
    z = len(nb)
    pairs = [(v, j, m[0], m[1]) for j, m in nb]
    pairs += [(ja, jb, mb[0] - ma[0], mb[1] - ma[1]) for (ja, ma), (jb, mb) in itertools.combinations(nb, 2)]
    R = resist(nv, E, deg, pairs)
    Rv, Rab = R[:z], iter(R[z:])
    Y = np.diag(Rv)
    for x, y in itertools.combinations(range(z), 2):
        Y[x, y] = Y[y, x] = (Rv[x] + Rv[y] - next(Rab)) / 2
    P = [sum(np.linalg.det(np.array([Y[r] if r in S else np.eye(z)[r] - Y[r] for r in range(z)]))
             for S in itertools.combinations(range(z), k)) for k in range(z + 1)]
    P1 = [np.linalg.det(np.eye(z - 1) - np.delete(np.delete(Y, d, 0), d, 1)) for d in range(z)]
    return dict(Rv=list(Rv), Y=Y.tolist(), P=P, P1=P1)


def task_sand():
    _, V, E = build_snub()
    cen, deg, CE = build_cairo(V, E)
    v3, v4 = deg.index(3), deg.index(4)
    return dict(snub=ust(4, E, [5] * 4, 0), c3=ust(6, CE, deg, v3), c4=ust(6, CE, deg, v4))


def task_cairo():
    _, V, E = build_snub()
    cen, deg, CE = build_cairo(V, E)
    d4 = [v for v in range(6) if deg[v] == 4]
    pairs = [(d4[0], d4[1], m1, m2) for m1 in (-1, 0, 1) for m2 in (-1, 0, 1)] + [(d4[0], d4[0], 1, 0), (d4[0], d4[0], 1, 1)]
    R = resist(6, CE, deg, pairs)
    ids = [mp.identify(r, ['sqrt(3)', '1/pi'], tol=1e-12) for r in R]
    return dict(R=list(R), ids=ids)


def task_asym(Rc):
    _, V, E = build_snub()
    h, _ = harmonic(4, E, V)
    pairs = []
    for k in (0, 1):
        for l in range(4):
            for ph in np.radians([3, 17, 31, 52, 74, 96, 118, 140]):
                pairs.append((k, l, int(round(Rc * np.cos(ph))), int(round(Rc * np.sin(ph)))))
    R = resist(4, E, [5] * 4, pairs)
    out = []
    for (k, l, n1, n2), r in zip(pairs, R):
        geo = np.linalg.norm(V[l] + L * np.array([n1, n2]) - V[k])
        x = (h[l] + L * np.array([n1, n2]) - h[k]) / L
        out.append(dict(k=k, l=l, R=float(r), geo=float(geo), x=float(x[0]), y=float(x[1])))
    return Rc, out


def run_task(t):
    try:
        import psutil
        psutil.Process().nice(psutil.BELOW_NORMAL_PRIORITY_CLASS if os.name == 'nt' else 10)
    except Exception:
        pass
    if t[0] == 'asym':
        return t, task_asym(t[1])
    return t, {'tabell': task_tabell, 'sand': task_sand, 'cairo': task_cairo}[t[0]]()


# ================================================================== huvudprogram
ok = []


def check(name, cond):
    ok.append(bool(cond))
    print(('OK  ' if cond else 'FEL ') + name, flush=True)


def main():
    import multiprocessing as mpp
    import sympy as sp
    t0 = time.time()
    tasks = [('asym', 80), ('asym', 40), ('asym', 20), ('tabell',), ('sand',), ('cairo',)]
    pool = mpp.Pool(3)
    async_res = pool.map_async(run_task, tasks)

    # ---- K1-K5 (snabba, i huvudprocessen)
    alpha, V, E = build_snub()
    cen, deg, CE = build_cairo(V, E)
    degs = [sum((e[0] == v) + (e[1] == v) for e in E) for v in range(4)]
    check(f'K1: alpha = {np.degrees(alpha):.9f}°, {len(E)} kanter per cell, grader {degs}; Cairo-grader {sorted(deg)}',
          abs(np.degrees(abs(alpha)) - 15) < 1e-8 and len(E) == 10 and degs == [5] * 4 and sorted(deg) == [3, 3, 3, 3, 4, 4])

    X, Z = sp.symbols('X Z')
    M = sp.diag(*[5] * 4)
    for i, j, (m1, m2) in E:
        M[i, j] -= X**m1 * Z**m2
        M[j, i] -= X**(-m1) * Z**(-m2)
    det = sp.expand(M.det())
    a, b = (X + 1 / X) / 2, (Z + 1 / Z) / 2
    D = a**2 + b**2 - 14 * a * b - 36 * (a + b) + 84
    rng = np.random.default_rng(1)
    dc = max(abs(np.linalg.det(bloch(6, CE, deg, *th)) - np.linalg.det(bloch(4, E, [5] * 4, *th)))
             for th in rng.uniform(-np.pi, np.pi, (20, 2)))
    check(f'K2: det M_snub = 4D(cos θ1, cos θ2) exakt (sympy), det L_Cairo = det M_snub (20 slumpade θ, max skillnad {dc:.1e})',
          sp.simplify(det - 4 * D) == 0 and dc < 1e-9)

    A = M.adjugate()

    def box(e):
        P = sp.Poly(sp.expand(e * X**8 * Z**8), X, Z)
        ex = [(p - 8, q - 8) for (p, q), _ in P.terms()]
        return (min(x for x, _ in ex), max(x for x, _ in ex)), (min(y for _, y in ex), max(y for _, y in ex))
    need = [(j, j) for j in range(4)] + [(j, 0) for j in (1, 2, 3)] + [(0, j) for j in (1, 2, 3)]
    bx = {ij: box(A[ij]) for ij in need}
    zord = {f'{i + 1}{j + 1}': bx[(i, j)][1][0] for i, j in need if i != j or i == 0}
    D1 = {'11': -1, '21': -1, '31': -2, '41': -2, '12': -1, '13': 0, '14': -1}
    wadj = max(4 * ((x[1] - x[0]) + (y[1] - y[0])) for x, y in bx.values())
    check(f'K3: Z-ordningar {zord} = kontroll D1; alla tio exponentlådor innehåller 0; w_adj = {wadj}',
          zord == D1 and all(x[0] <= 0 <= x[1] and y[0] <= 0 <= y[1] for x, y in bx.values()) and wadj == 20)

    h, T = harmonic(4, E, V)
    d = (h - h[0]) / L
    dwant = np.array([[0, 0], [1 / 2, -1 / 6], [1 / 3, -2 / 3], [-1 / 6, -1 / 2]])
    _, Tc = harmonic(6, CE, cen)
    check(f'K4: harmoniska lägen = d_j i asymptotik.md; Σ e e^T = {np.round(T[0, 0], 12)} I (snub), {np.round(Tc[0, 0], 12)} I (Cairo)',
          np.abs(d - dwant).max() < 1e-12 and np.abs(T - 4 / 3 * np.eye(2)).max() < 1e-12 and np.abs(Tc - 0.75 * np.eye(2)).max() < 1e-12)

    z = np.exp(1j * np.pi / 12)
    phi = (1 + 5**0.5) / 2
    EX = {1: 1, 5: 1, 7: 1, 11: 1, 13: -1, 17: -1, 19: -1, 23: -1}
    EZ = {1: 1, 5: 1, 19: 1, 23: 1, 7: -1, 11: -1, 13: -1, 17: -1}
    Xs = lambda s: np.prod([(s - z**k)**e for k, e in EX.items()], axis=0)
    Zs = lambda s: np.prod([(s - z**k)**e for k, e in EZ.items()], axis=0)
    p5 = True
    for s, lo_hi, sign in ((np.linspace(1 + 1e-9, phi - 1e-9, 200001), (0, np.pi), +1),
                           (np.geomspace(phi + 1e-9, 1e6, 200001), (-np.pi, 0), -1)):
        th = np.unwrap(np.angle(Xs(s)))
        A_ = np.cos(th)
        bb = 7 * A_ + 18 + sign * 4 * S3 * np.sqrt((A_ + 1) * (A_ + 5))
        p5 &= bool(np.all(np.diff(th) > 0)) and abs(th[0] - lo_hi[0]) < 1e-3 and abs(th[-1] - lo_hi[1]) < 1e-3
        p5 &= np.abs(Zs(s).real - (bb - np.sqrt(bb**2 - 1))).max() < 1e-8 and np.abs(Zs(s).imag).max() < 1e-12
    check('K5: P3: θ1 = arg X(s) växer 0 -> π på (1, φ) och -π -> 0 på (φ, ∞); Z_in(s) = b+-roten resp. b--roten', p5)

    # ---- resultat från de parallella delarna
    res = dict(async_res.get())
    pool.close()
    sand, cai, tab = res[('sand',)], res[('cairo',)], res[('tabell',)]

    p = 1 / np.pi
    Rs, R3, R4 = sand['snub']['Rv'], sand['c3']['Rv'], sand['c4']['Rv']
    check(f'K6: kantresistanser snub {sorted(set(np.round(Rs, 13)))}, Cairo {sorted(set(np.round(R3 + R4, 13)))}',
          np.allclose(sorted(Rs), sorted([2 * S3 / 9] + [1 / 2 - S3 / 18] * 4), atol=1e-13)
          and np.allclose(sorted(R3), sorted([1 - 2 * S3 / 9] + [1 / 2 + S3 / 18] * 2), atol=1e-13)
          and np.allclose(R4, [1 / 2 + S3 / 18] * 4, atol=1e-13))

    if tab is None:
        print('    (K7 hoppas över: runs/tabell/checkpoint.json saknas)')
    else:
        check(f'K7: 40 slumpvalda tabellvärden (avstånd upp till {tab["maxdist"]:.1f}): max |R_egen - R_tabell| = {tab["dR"]:.1e}, '
              f'avstånden på {tab["dd"]:.0e}', tab['dR'] < 1e-12 and tab['dd'] < 1e-9)

    GC = 8953 / 3888 - 281 * S3 / 216 + (44 * S3 / 27 - 11 / 4) * p + (17 / 36 - S3 / 3) * p**2
    Psnub = [5 * GC,
             (55 * S3 / 9 - 179 / 18) * p**2 + (977 / 18 - 853 * S3 / 27) * p + 30245 * S3 / 972 - 51893 / 972,
             (44 / 3 - 25 * S3 / 3) * p**2 + (401 * S3 / 9 - 232 / 3) * p - 8125 * S3 / 162 + 56447 / 648,
             (5 * S3 - 161 / 18) * p**2 + (281 / 6 - 727 * S3 / 27) * p + 10685 * S3 / 324 - 55511 / 972,
             (67 / 36 - 10 * S3 / 9) * p**2 + (157 * S3 / 27 - 361 / 36) * p - 14455 * S3 / 1944 + 50057 / 3888]
    GCY = sorted([1 - S3 / 2, S3 / 9, (3 - S3) / 6, 1 / 2 - 2 * S3 / 9, 5 * S3 / 9 - 1 + p / 2, 1 / 2 - S3 / 9 - p / 2])
    Ys = np.array(sand['snub']['Y'])
    ys = sorted(set(np.round(Ys[np.triu_indices(5, 1)], 12)))
    sn = sand['snub']
    check(f'K8: snub square: P1 = {sn["P1"][0]:.15f} (Goodchild {GC:.15f}) för alla 5 val, P(grad 1..5) = sandhog_goodchild.md, '
          f'de 6 Y-värdena = Goodchilds Tabell III',
          max(abs(x - GC) for x in sn['P1']) < 1e-13 and np.allclose(sn['P'][1:], Psnub, atol=1e-13) and abs(sn['P'][0]) < 1e-13
          and len(ys) == 6 and np.allclose(ys, GCY, atol=1e-11))
    P3c = [S3 / 3 - 2 / 9, 13 / 9 - 5 * S3 / 9, 2 * S3 / 9 - 2 / 9]
    P4c = [-2 * p**3 + (2 - 8 * S3 / 9) * p**2 + (4 * S3 / 9 - 10 / 27) * p - 4 * S3 / 243 + 2 / 27,
           6 * p**3 + (8 * S3 / 3 - 9) * p**2 + (46 / 9 - 22 * S3 / 9) * p + 40 * S3 / 81 - 13 / 27,
           -6 * p**3 + (12 - 8 * S3 / 3) * p**2 + (32 * S3 / 9 - 82 / 9) * p - 94 * S3 / 81 + 74 / 27,
           2 * p**3 + (8 * S3 / 9 - 5) * p**2 + (118 / 27 - 14 * S3 / 9) * p + 166 * S3 / 243 - 4 / 3]
    c3, c4 = sand['c3'], sand['c4']
    check(f'K8: Cairo: P1 = {c3["P1"][0]:.13f} (grad 3; √3/9 - 2/27) och {c4["P1"][0]:.13f} (grad 4), gradfördelningarna = '
          f'sandhog_goodchild.md', np.allclose(c3['P1'], S3 / 9 - 2 / 27, atol=1e-13) and np.allclose(c4['P1'], P4c[0] / 4, atol=1e-13)
          and np.allclose(c3['P'][1:], P3c, atol=1e-13) and np.allclose(c4['P'][1:], P4c, atol=1e-13))

    print('    Cairo-stickprov (grad 4 - grad 4):', [f'{r:.12f} = {i}' for r, i in zip(cai['R'], cai['ids'])])
    check('K9: Cairo: 4√3/9 förekommer (bevis_cairo.md, R(4-4, d² = 3)), och alla 11 stickprov identifieras i Q + Q√3 + Q/π',
          any(abs(r - 4 * S3 / 9) < 1e-13 for r in cai['R']) and all(i is not None for i in cai['ids']))

    # ---- K10 asymptotiken
    print(f'    C = {C_ASYM:.15f}')
    summary = {}
    for Rc in (20, 40, 80):
        rows = res[('asym', Rc)][1]
        cest, e3, egeo, fit = [], [], [], {}
        for r in rows:
            rx = np.hypot(r['x'], r['y'])
            ph = np.arctan2(r['y'], r['x'])
            sats1 = S3 / 36 + 3 / (4 * np.pi) * (GAMMA + np.log(2 * np.sqrt(6) * rx))
            t2 = -(EPS[r['k']] + EPS[r['l']]) / 2 * np.sin(2 * ph) / (48 * np.pi * rx**2)
            cest.append(r['R'] - 3 / (4 * np.pi) * np.log(L * rx) - t2)
            e3.append(abs(r['R'] - sats1 - t2) * rx**3)
            egeo.append(abs(r['R'] - 3 / (4 * np.pi) * np.log(r['geo']) - C_ASYM) * r['geo'])
            same = EPS[r['k']] == EPS[r['l']]
            fit.setdefault((r['k'], same), ([], []))
            fit[(r['k'], same)][0].append([1, np.cos(2 * ph), np.sin(2 * ph), np.cos(4 * ph), np.sin(4 * ph)])
            fit[(r['k'], same)][1].append((r['R'] - sats1) * rx**2)
        coef = {key: np.linalg.lstsq(np.array(Xm), np.array(y), rcond=None)[0] for key, (Xm, y) in fit.items()}
        summary[Rc] = (min(cest), max(cest), max(e3), max(egeo), coef)
        print(f'    |n| ≈ {Rc}: C-skattning [{min(cest):.12f}, {max(cest):.12f}], max |R - Sats 2|·|x̃|³ = {max(e3):.5f}, '
              f'max |R - (3/4π)log r - C|·r = {max(egeo):.4f}; sin2φ-koefficient: nod 1 {coef[(0, True)][2]:+.5f}, '
              f'nod 2 {coef[(1, True)][2]:+.5f} (väntat ∓1/(48π) = ∓{1 / (48 * np.pi):.5f})')
    lo80, hi80 = summary[80][0], summary[80][1]
    check(f'K10: konstanten C ligger i C-skattningen vid |n| ≈ 80 ([{lo80:.11f}, {hi80:.11f}]), bredd < 3e-8',
          lo80 <= C_ASYM <= hi80 and hi80 - lo80 < 3e-8)
    check('K10: felet efter Sats 2 är O(|x̃|^-3): max |fel|·|x̃|³ < 0,005 vid alla tre radier, och spridningen i C-skattningen '
          'minskar med en faktor ≥ 6 per fördubbling',
          all(summary[r][2] < 0.005 for r in summary)
          and all((summary[r][1] - summary[r][0]) / (summary[2 * r][1] - summary[2 * r][0]) > 6 for r in (20, 40)))
    c80 = summary[80][4]
    amp = 1 / (48 * np.pi)
    check('K10: vid |n| ≈ 80: sin 2φ-amplituden = ∓1/(48π) inom 1 % för samma paritet (tecken - för nod 1, + för nod 2), '
          'blandade par och cos 2φ, cos 4φ, sin 4φ < 1e-4',
          abs(c80[(0, True)][2] + amp) < 0.01 * amp and abs(c80[(1, True)][2] - amp) < 0.01 * amp
          and all(abs(c80[(k, False)][i]) < 1e-4 for k in (0, 1) for i in range(5))
          and all(abs(c80[(k, True)][i]) < 1e-4 for k in (0, 1) for i in (0, 1, 3, 4)))
    check('K10: med geometriskt avstånd är felet O(1/r): max |fel|·r är 0,020–0,023 vid alla tre radier',
          all(0.020 < summary[r][3] < 0.023 for r in summary))

    print(f'\nSUMMARY: {sum(ok)} godkända, {len(ok) - sum(ok)} fel  ({time.time() - t0:.0f} s)')
    sys.exit(0 if all(ok) else 1)


if __name__ == '__main__':
    main()
