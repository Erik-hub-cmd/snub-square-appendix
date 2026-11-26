"""Kombinatorisk sökning efter drag (Y→Δ, Δ→Y, serie, parallell, döda grenar, kontraktibla öglor) som för ett
periodiskt nät på sin egen torus till ett minimalt nät. Inbäddningen hålls som rotationssystem (darts med cellförskjutning).
Används som NUMERISKT stöd för omvändningen (Förmodan O i GALOIS_UTVIDGNING.md); en misslyckad sökning bevisar inget.
"""
import cmath, itertools, random, sys
from collections import deque


class Net:
    def __init__(self):
        self.dart = {}   # id -> [tail, head, (a,b)]
        self.rev = {}
        self.rot = {}    # vertex -> list of dart ids (ccw)
        self.nid = 0; self.nv = 0

    @staticmethod
    def from_geom(pos, v1, v2, E):
        G = Net(); G.nv = len(pos)
        tmp = {i: [] for i in range(len(pos))}
        for (i, j, n, _) in E:
            a, b = G._new_edge(i, j, tuple(n))
            d = pos[j] + n[0] * v1 + n[1] * v2 - pos[i]
            tmp[i].append((cmath.phase(d), a)); tmp[j].append((cmath.phase(-d), b))
        G.rot = {i: [x[1] for x in sorted(tmp[i])] for i in tmp}
        return G

    def _new_edge(self, i, j, n):
        a, b = self.nid, self.nid + 1; self.nid += 2
        self.dart[a] = [i, j, n]; self.dart[b] = [j, i, (-n[0], -n[1])]
        self.rev[a] = b; self.rev[b] = a
        return a, b

    def copy(self):
        G = Net(); G.dart = {k: list(v) for k, v in self.dart.items()}; G.rev = dict(self.rev)
        G.rot = {k: list(v) for k, v in self.rot.items()}; G.nid = self.nid; G.nv = self.nv
        return G

    def E(self): return len(self.dart) // 2

    def _del_edge(self, a):
        b = self.rev[a]
        for d in (a, b):
            self.rot[self.dart[d][0]].remove(d)
        for d in (a, b):
            del self.dart[d]; del self.rev[d]

    def _replace(self, v, old, new_list):
        L = self.rot[v]; p = L.index(old); self.rot[v] = L[:p] + new_list + L[p + 1:]

    # --------------------------------------------------------------- reduktioner (minskar E)
    def reduce_all(self):
        changed = True
        while changed:
            changed = False
            for v in list(self.rot):
                L = self.rot.get(v)
                if L is None: continue
                if len(L) == 0:
                    del self.rot[v]; changed = True; continue
                # kontraktibel ögla
                for d in L:
                    t, h, n = self.dart[d]
                    if h == v and n == (0, 0):
                        self._del_edge(d); changed = True; break
                if changed: break
                L = self.rot[v]
                if len(L) == 1:
                    d = L[0]; t, h, n = self.dart[d]
                    if h != v:
                        self._del_edge(d); del self.rot[v]; changed = True; break
                if len(L) == 2:
                    a, b = L
                    if self.rev[a] != b and self.dart[a][1] != v and self.dart[b][1] != v:
                        w1, oa = self.dart[a][1], self.dart[a][2]; w2, ob = self.dart[b][1], self.dart[b][2]
                        ra, rb = self.rev[a], self.rev[b]
                        x, y = self._new_edge(w1, w2, (ob[0] - oa[0], ob[1] - oa[1]))
                        self._replace(w1, ra, [x]); self._replace(w2, rb, [y])
                        for d in (a, b, ra, rb):
                            del self.dart[d]; del self.rev[d]
                        del self.rot[v]; changed = True; break
                # parallella kanter (samma ände, samma förskjutning)
                seen = {}
                for d in self.rot[v]:
                    t, h, n = self.dart[d]
                    key = (h, n)
                    if key in seen and self.rev[d] != seen[key]:
                        self._del_edge(d); changed = True; break
                    seen[key] = d
                if changed: break
        return self

    # --------------------------------------------------------------- Y-Δ
    def faces(self):
        seen = set(); F = []
        for d0 in self.dart:
            if d0 in seen: continue
            f = []; d = d0
            while d not in seen:
                seen.add(d); f.append(d)
                r = self.rev[d]; v = self.dart[r][0]; L = self.rot[v]
                d = L[(L.index(r) - 1) % len(L)]
            F.append(f)
        return F

    def delta_moves(self):
        out = []
        for f in self.faces():
            if len(f) == 3:
                off = [sum(self.dart[d][2][k] for d in f) for k in (0, 1)]
                es = {min(d, self.rev[d]) for d in f}
                if off == [0, 0] and len(es) == 3: out.append(('DY', tuple(f)))
        for v, L in self.rot.items():
            if len(L) == 3 and all(self.dart[d][1] != v for d in L): out.append(('YD', v))
        return out

    def apply(self, mv):
        G = self.copy()
        if mv[0] == 'DY':
            d1, d2, d3 = mv[1]
            vs = [G.dart[d][0] for d in (d1, d2, d3)]
            o1 = G.dart[d1][2]; o2 = G.dart[d2][2]
            offs = [(0, 0), o1, (o1[0] + o2[0], o1[1] + o2[1])]
            c = G.nv; G.nv += 1; G.rot[c] = []
            legs = []
            ds = [d1, d2, d3]
            for i in range(3):
                x, y = G._new_edge(c, vs[i], offs[i]); legs.append(x)
                di = ds[i]; dprev = G.rev[ds[i - 1]]
                L = G.rot[vs[i]]; p = L.index(di)
                assert L[(p + 1) % len(L)] == dprev
                if p + 1 < len(L):
                    G.rot[vs[i]] = L[:p] + [y] + L[p + 2:]
                else:
                    G.rot[vs[i]] = [y] + L[1:p]
            for d in ds:
                r = G.rev[d]
                del G.dart[d]; del G.dart[r]; del G.rev[d]; del G.rev[r]
            G.rot[c] = legs
        else:
            c = mv[1]; a = G.rot[c]
            ws = [G.dart[d][1] for d in a]; os_ = [G.dart[d][2] for d in a]
            new = {}
            for i in range(3):
                j = (i + 1) % 3
                x, y = G._new_edge(ws[i], ws[j], (os_[j][0] - os_[i][0], os_[j][1] - os_[i][1]))
                new[(i, j)] = x; new[(j, i)] = y
            for i in range(3):
                G._replace(ws[i], G.rev[a[i]], [new[(i, (i + 1) % 3)], new[(i, (i - 1) % 3)]])
            for d in a:
                r = G.rev[d]; del G.dart[d]; del G.dart[r]; del G.rev[d]; del G.rev[r]
            del G.rot[c]
        return G

    # --------------------------------------------------------------- tågspår och minimalitet
    def tracks(self):
        where = {}
        for v, L in self.rot.items():
            for p, d in enumerate(L): where[d] = (v, p)
        seen = set(); T = []
        for d0 in self.dart:
            for par0 in (0, 1):
                if (d0, par0) in seen: continue
                d, par = d0, par0; off = (0, 0); seq = []
                while (d, par) not in seen:
                    seen.add((d, par)); seq.append(min(d, self.rev[d]))
                    n = self.dart[d][2]; off = (off[0] + n[0], off[1] + n[1])
                    r = self.rev[d]; v, p = where[r]; L = self.rot[v]
                    d = L[(p + 1) % len(L)] if par == 0 else L[(p - 1) % len(L)]
                    par ^= 1
                T.append((sorted(seq), off))
        res = []; used = [False] * len(T)
        for a in range(len(T)):
            if used[a]: continue
            for b in range(a + 1, len(T)):
                if not used[b] and T[b][0] == T[a][0] and T[b][1] == (-T[a][1][0], -T[a][1][1]):
                    used[a] = used[b] = True; res.append(T[a]); break
        return res

    def is_minimal(self):
        T = self.tracks()
        on = {}
        for ti, (es, h) in enumerate(T):
            for e in es: on.setdefault(e, []).append(ti)
        m = len(T); C = [[0] * m for _ in range(m)]
        for x in on.values():
            if len(x) != 2: return False
            a, b = x
            C[a][b] += 1
            if a != b: C[b][a] += 1
        hs = [h for _, h in T]
        for i in range(m):
            if hs[i] == (0, 0) or C[i][i]: return False
            for j in range(i + 1, m):
                if C[i][j] != abs(hs[i][0] * hs[j][1] - hs[i][1] * hs[j][0]): return False
        return True

    def key(self):
        """invariant (WL-förfining, inte kanonisk form) för att undvika upprepningar"""
        lab = {v: (len(L),) for v, L in self.rot.items()}
        for _ in range(3):
            lab = {v: hash((lab[v], tuple(sorted(lab[self.dart[d][1]] for d in L)))) for v, L in self.rot.items()}
        fl = sorted(len(f) for f in self.faces())
        T = sorted(tuple(sorted([h, (-h[0], -h[1])])[1]) + (len(es),) for es, h in self.tracks())
        return (self.E(), tuple(sorted(lab.values())), tuple(fl), tuple(T))


def search(G0, maxdepth=4, maxstates=20000, rounds=30, log=print):
    G = G0.copy().reduce_all()
    if G.is_minimal(): return G, 'minimalt efter reduktioner', []
    path = []
    for rnd in range(rounds):
        E0 = G.E()
        Q = deque([(G, 0, [])]); seen = {G.key()}; best = None; n = 0
        while Q and n < maxstates:
            H, dep, p = Q.popleft(); n += 1
            if dep >= maxdepth: continue
            for mv in H.delta_moves():
                K = H.apply(mv).reduce_all()
                if K.is_minimal(): return K, f'MINIMALT nått (runda {rnd + 1})', path + p + [mv]
                if K.E() < E0:
                    best = (K, p + [mv]); break
                k = K.key()
                if k not in seen:
                    seen.add(k); Q.append((K, dep + 1, p + [mv]))
            if best: break
        if not best:
            return G, f'fastnade vid E={G.E()} (djup {maxdepth}, {n} tillstånd)', path
        G, pp = best; path += pp
        log(f'   runda {rnd + 1}: E {E0} -> {G.E()}')
    return G, f'slut på rundor vid E={G.E()}', path


def random_walks(G0, walks=2000, length=40, seed=0):
    """slumpvandringar av Y-Δ-drag med reduktioner; börjar om från bästa nät hittills"""
    rng = random.Random(seed)
    best = G0.copy().reduce_all()
    if best.is_minimal(): return best, 'minimalt efter reduktioner'
    for w in range(walks):
        G = best
        for s in range(length):
            mv = G.delta_moves()
            if not mv: break
            G = G.apply(rng.choice(mv)).reduce_all()
            if G.is_minimal(): return G, f'MINIMALT nått (slumpvandring {w + 1})'
            if G.E() < best.E():
                best = G; break
    return best, f'ej minimalt efter {walks} vandringar (lägst E={best.E()})'
