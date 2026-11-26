"""
Isoradiell förklaring, steg 1 (forslag/isoradiell_forklaring.md): når elektriska drag Z² från snub square?

Söker över alla följder av Δ→Y, Y→Δ (en bana av trianglar/noder per drag, periodiskt), med seriekoppling,
parallellkoppling och borttagning av hängande noder direkt efter varje drag. Rapporterar antalet grafer som nås
och om någon har formen av Z² (4 noder per cell, alla av grad 4, 8 kanter). Kör: python kod/isoradiell_steg1_elektrisk.py
"""
import sys, copy, itertools
sys.stdout.reconfigure(encoding="utf-8")
import numpy as np

DJUP = 7


def snub(lam):
    """snub square som Z² med diagonaler, cell 2×2; kanter (u, v, förflyttning, konduktans)"""
    E = []
    idx = lambda i, j: 2 * (i % 2) + (j % 2)
    sh = lambda i, j: ((i - (i % 2)) // 2, (j - (j % 2)) // 2)

    def add(p, q, w):
        s0, s1 = sh(*p), sh(*q)
        E.append((idx(*p), idx(*q), (s1[0] - s0[0], s1[1] - s0[1]), w))
    for i in range(2):
        for j in range(2):
            add((i, j), (i + 1, j), 1.0); add((i, j), (i, j + 1), 1.0)
    add((0, 0), (1, 1), lam); add((2, 1), (1, 2), lam)
    return {'V': set(range(4)), 'E': E}


def bloch(G, X, Z):
    V = sorted(G['V']); ix = {v: k for k, v in enumerate(V)}
    M = np.zeros((len(V), len(V)), complex)
    for u, v, (sx, sy), w in G['E']:
        ph = X**sx * Z**sy; a, b = ix[u], ix[v]
        M[a, a] += w; M[b, b] += w; M[a, b] -= w * ph; M[b, a] -= w / ph
    return M


def nbrs(G, u):
    out = []
    for k, (a, b, s, w) in enumerate(G['E']):
        if a == u: out.append((b, s, w, k))
        if b == u: out.append((a, (-s[0], -s[1]), w, k))
    return out


def sub(s, t):
    return (s[0] - t[0], s[1] - t[1])


def simplify(G):
    """serie, parallell, hängande noder, öglor utan förflyttning; returnerar faktorn f med det_ny = f·det_gammal"""
    f, changed = 1.0, True
    while changed:
        changed = False
        d = {}
        for a, b, s, w in G['E']:
            if a == b and s == (0, 0):
                continue
            key = (a, b, s) if (a, s) <= (b, (-s[0], -s[1])) else (b, a, (-s[0], -s[1]))
            if a == b: key = (a, a, max(s, (-s[0], -s[1])))
            if key in d: d[key] += w; changed = True
            else: d[key] = w
        G['E'] = [(k[0], k[1], k[2], w) for k, w in d.items()]
        for u in list(G['V']):
            nb = nbrs(G, u)
            if len(nb) == 2 and all(x[0] != u for x in nb):
                (a, s1, w1, k1), (b, s2, w2, k2) = nb
                G['E'] = [e for i, e in enumerate(G['E']) if i not in (k1, k2)] + [(a, b, sub(s2, s1), w1 * w2 / (w1 + w2))]
                G['V'].remove(u); f /= (w1 + w2); changed = True; break
            if len(nb) == 1:
                (a, s1, w1, k1), = nb
                G['E'] = [e for i, e in enumerate(G['E']) if i != k1]; G['V'].remove(u); f /= w1; changed = True; break
    return f


def triangles(G):
    T = set()
    for u in G['V']:
        nb = [x for x in nbrs(G, u) if x[0] != u]
        for (v, s1, w1, k1), (x, s2, w2, k2) in itertools.combinations(nb, 2):
            for (y, s3, w3, k3) in nbrs(G, v):
                if y == x and s3 == sub(s2, s1) and k3 not in (k1, k2):
                    T.add(tuple(sorted((k1, k2, k3))))
    return sorted(T)


def DY(G, tri):
    E = G['E']; a, b, s, wab = E[tri[0]]
    rest = [E[tri[1]], E[tri[2]]]

    def rel(e, base):
        x, y, t, w = e
        if x == base: return (y, t, w)
        if y == base: return (x, (-t[0], -t[1]), w)
    r = [z for z in (rel(e, a) for e in rest) if z and not (z[0] == b and z[1] == s)]
    c, sc, wac = r[0]
    wbc = [e for e in rest if not (rel(e, a) and rel(e, a)[0] == c and rel(e, a)[1] == sc)][0][3]
    S = wab * wac + wab * wbc + wac * wbc
    n = max(G['V']) + 1
    newE = [e for i, e in enumerate(E) if i not in tri] + [(n, a, (0, 0), S / wbc), (n, b, s, S / wac), (n, c, sc, S / wab)]
    return {'V': G['V'] | {n}, 'E': newE}, S / wbc + S / wac + S / wab


def YD(G, u):
    nb = nbrs(G, u)
    if len(nb) != 3 or any(x[0] == u for x in nb): return None
    S = sum(x[2] for x in nb)
    newE = [e for i, e in enumerate(G['E']) if i not in [x[3] for x in nb]]
    for (a, s1, w1, _), (b, s2, w2, _) in itertools.combinations(nb, 2):
        newE.append((a, b, sub(s2, s1), w1 * w2 / S))
    return {'V': G['V'] - {u}, 'E': newE}, 1 / S


def z2like(G):
    return (len(G['V']) == 4 and len(G['E']) == 8 and not any(a == b for a, b, s, w in G['E'])
            and all(len(nbrs(G, u)) == 4 for u in G['V']))


PTS = [(np.exp(1j * a), np.exp(1j * b)) for a, b in [(0.3, 1.1), (2.0, -0.7), (1.3, 0.4)]]
G0 = snub(0.7)
d0 = [np.linalg.det(bloch(G0, *p)) for p in PTS]
seen, hits, faktor_ok = set(), [], True


def dfs(G, F, depth):
    global faktor_ok
    key = (len(G['V']), len(G['E']), tuple(sorted(round(w, 9) for *_, w in G['E'])))
    seen.add(key)
    # kontroll av buffringen: det(G) = F·det(snub) i varje steg
    faktor_ok = faktor_ok and all(abs(np.linalg.det(bloch(G, *p)) - F * d) < 1e-8 * abs(F * d) for p, d in zip(PTS, d0))
    if z2like(G): hits.append(G)
    if depth == 0 or len(G['V']) > 12: return
    for tri in triangles(G):
        H, f = DY(copy.deepcopy(G), tri); f2 = simplify(H); dfs(H, F * f * f2, depth - 1)
    for u in sorted(G['V']):
        r = YD(copy.deepcopy(G), u)
        if r: H, f = r; f2 = simplify(H); dfs(H, F * f * f2, depth - 1)


dfs(G0, 1.0, DJUP)
print(f'OK  Bloch-determinanten följer med rätt faktor i varje drag: {faktor_ok}')
print(f'{"OK " if not hits else "OBS"} olika grafer nådda upp till djup {DJUP}: {len(seen)}; med Z²-form: {len(hits)}')
