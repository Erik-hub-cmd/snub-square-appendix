"""Separat skriven kontroll av datorsökningen i Anm. 3.5 (arbetsanteckningen ekvivalens_teori.md, avsnitt 9.3, publiceras inte).

Skrivet utan att läsa kod/storre_torus*.py eller kod/isoradiell_steg1_elektrisk.py.

Modell. Snub square = Z^2 med diagonalen (2i,2j)-(2i+1,2j+1) och antidiagonalen (2i+2,2j+1)-(2i+1,2j+2).
Övertorusen T_Λ ges av två generatorer för Λ i Z^2-koordinater (Λ ⊂ 2Z^2). Nätet lagras som en
kombinatorisk karta på torusen: rotationssystem (darter moturs runt varje nod) och en spänning
(voltage) i Λ på varje dart. Dart d och d^1 är de två riktningarna av samma kant.

Drag (alla periodiska med perioden Λ, alla elektriskt giltiga med positiva konduktanser):
  Y->Δ   nod av grad 3 (utan ögla),
  Δ->Y   yta av grad 3 med tre olika kanter och tre olika lyfta hörn,
  serie  nod av grad 2 med två olika kanter,
  parallell  två kanter med samma lyfta ändpunkter (oavsett om de begränsar en yta),
  död gren   en delgraf som hänger fast i en enda nod och har ändligt lyft (inkl. kontraktibla öglor).
Inga kanter delas. Sökningen är uttömmande (bredden först, alla nåbara kartor upp till isomorfi
och gauge), om inte --max-tillstand sätter ett tak.

Kör:  .venv\\Scripts\\python kod\\granskning5\\torus_sok.py                         (alla fall utom 2x2 och 4x1, ~5 min)
      .venv\\Scripts\\python kod\\granskning5\\torus_sok.py --det --invariant      (dessutom kontrollerna)
      .venv\\Scripts\\python kod\\granskning5\\torus_sok.py --fall 2x2 --kompakt   (de stora fallen, en i taget; timmar)
Kör med PYTHONIOENCODING=utf-8 om konsolen visar å/ä/ö fel.
"""
import argparse
import cmath
import hashlib
from array import array
import math
import random
import sys
import time
from fractions import Fraction
from collections import deque

import numpy as np


def add(p, q):
    return (p[0] + q[0], p[1] + q[1])


def sub(p, q):
    return (p[0] - q[0], p[1] - q[1])


def neg(p):
    return (-p[0], -p[1])


# ---------------------------------------------------------------- start: snub square på T_Λ

def reducer(B):
    (a, b), (c, d) = B  # generatorer (a,b), (c,d)
    det = a * d - b * c
    assert det != 0

    def red(p):
        x, y = p
        # p = s*(a,b) + t*(c,d)
        s = Fraction(x * d - y * c, det)
        t = Fraction(a * y - b * x, det)
        fs, ft = math.floor(s), math.floor(t)
        lam = (fs * a + ft * c, fs * b + ft * d)
        return sub(p, lam), lam
    return red, abs(det)


def snub_karta(B):
    red, idx = reducer(B)
    reps = set()
    for x in range(-20, 21):
        for y in range(-20, 21):
            reps.add(red((x, y))[0])
    reps = sorted(reps)
    assert len(reps) == idx, (len(reps), idx)
    vid = {r: i for i, r in enumerate(reps)}
    darts_at = {i: [] for i in range(len(reps))}  # (vinkel, dart)
    volt = {}
    e = 0
    for r in reps:
        x, y = r
        dirs = [(1, 0), (0, 1)]
        if x % 2 == 0 and y % 2 == 0:
            dirs.append((1, 1))
        if x % 2 == 0 and y % 2 == 1:
            dirs.append((-1, 1))
        for dlt in dirs:
            q = add(r, dlt)
            qr, lam = red(q)
            d0, d1 = 2 * e, 2 * e + 1
            volt[d0] = lam
            volt[d1] = neg(lam)
            darts_at[vid[r]].append((math.atan2(dlt[1], dlt[0]), d0))
            darts_at[vid[qr]].append((math.atan2(-dlt[1], -dlt[0]), d1))
            e += 1
    rot = {v: tuple(d for _, d in sorted(L)) for v, L in darts_at.items()}
    return rot, volt


# ---------------------------------------------------------------- kartverktyg

def tails(rot):
    t = {}
    for v, L in rot.items():
        for d in L:
            t[d] = v
    return t


def sigma(rot):
    s = {}
    for v, L in rot.items():
        n = len(L)
        for i, d in enumerate(L):
            s[d] = L[(i + 1) % n]
    return s


def faces(rot):
    s = sigma(rot)
    seen = set()
    F = []
    for d in s:
        if d in seen:
            continue
        orb = []
        x = d
        while x not in seen:
            seen.add(x)
            orb.append(x)
            x = s[x ^ 1]
        F.append(orb)
    return F


def genus(rot):
    V = len(rot)
    D = sum(len(L) for L in rot.values())
    E = D // 2
    F = len(faces(rot))
    chi = V - E + F
    return (2 - chi) // 2, chi


def ny_dart(volt):
    m = max(volt) if volt else -1
    return (m | 1) + 1  # nästa jämna id


def edges_of(rot):
    return sorted({d >> 1 for L in rot.values() for d in L})


# ---------------------------------------------------------------- kanonisk form

def kanon(rot, volt):
    t = tails(rot)
    s = sigma(rot)
    best = None
    darts = list(s)
    # starta bara i darter med minsta (isomorfi-invarianta) etikett: grader i ändarna och ytlängder
    flen = {}
    for f in faces(rot):
        for x in f:
            flen[x] = len(f)
    inv = {d: (len(rot[t[d]]), len(rot[t[d ^ 1]]), flen[d], flen[d ^ 1]) for d in darts}
    mi = min(inv.values())
    for d0 in [d for d in darts if inv[d] == mi]:
        lab = {d0: 0}
        order = [d0]
        pos = {t[d0]: (0, 0)}
        i = 0
        code = []
        while i < len(order):
            x = order[i]
            for y in (s[x], x ^ 1):
                if y not in lab:
                    lab[y] = len(order)
                    order.append(y)
                    if t[y] not in pos:
                        # y = x^1, nås från tail(x) via x
                        pos[t[y]] = add(pos[t[x]], volt[x])
            i += 1
        if len(order) != len(darts):
            raise RuntimeError("osammanhängande karta")
        for x in order:
            res = sub(add(pos[t[x]], volt[x]), pos[t[x ^ 1]])
            code.append((lab[s[x]], lab[x ^ 1], res[0], res[1]))
        code = tuple(code)
        if best is None or code < best:
            best = code
    return best


# ---------------------------------------------------------------- dragen

def y_till_delta(rot, volt, cond, v):
    d = rot[v]
    t = tails(rot)
    r = [x ^ 1 for x in d]
    if any(t[x] == v for x in r):
        return None  # ögla vid v
    w = [volt[x] for x in d]
    base = ny_dart(volt)
    X = [base + 2 * i for i in range(3)]
    nv = dict(volt)
    for x in d:
        del nv[x], nv[x ^ 1]
    for i in range(3):
        j = (i + 1) % 3
        nv[X[i]] = sub(w[j], w[i])
        nv[X[i] ^ 1] = neg(nv[X[i]])
    repl = {r[i]: (X[i], X[(i - 1) % 3] ^ 1) for i in range(3)}
    nr = {}
    for u, L in rot.items():
        if u == v:
            continue
        out = []
        for x in L:
            if x in repl:
                out.extend(repl[x])
            else:
                out.append(x)
        nr[u] = tuple(out)
    nc = None
    if cond is not None:
        nc = dict(cond)
        a = [cond[x >> 1] for x in d]
        S = sum(a)
        for x in d:
            del nc[x >> 1]
        for i in range(3):
            j = (i + 1) % 3
            nc[X[i] >> 1] = a[i] * a[j] / S
    return nr, nv, nc


def delta_till_y(rot, volt, cond, f):
    f1, f2, f3 = f
    if len({f1 >> 1, f2 >> 1, f3 >> 1}) < 3:
        return None
    t = tails(rot)
    p = [(0, 0), volt[f1], add(volt[f1], volt[f2])]
    if add(p[2], volt[f3]) != (0, 0):
        return None  # ytan är inte en kontraktibel triangel
    lifted = {(t[f1], p[0]), (t[f2], p[1]), (t[f3], p[2])}
    if len(lifted) < 3:
        return None
    base = ny_dart(volt)
    G = [base + 2 * i for i in range(3)]
    v = max(rot) + 1
    nv = dict(volt)
    for x in f:
        del nv[x], nv[x ^ 1]
    for i in range(3):
        nv[G[i]] = p[i]
        nv[G[i] ^ 1] = neg(p[i])
    fl = [f1, f2, f3]
    repl = {}
    drop = set()
    for i in range(3):
        repl[fl[(i - 1) % 3] ^ 1] = G[i] ^ 1
        drop.add(fl[i])
    nr = {}
    for u, L in rot.items():
        out = []
        for x in L:
            if x in drop:
                continue
            out.append(repl.get(x, x))
        nr[u] = tuple(out)
    nr[v] = (G[0], G[2], G[1])
    nc = None
    if cond is not None:
        nc = dict(cond)
        g = [cond[x >> 1] for x in fl]  # g[i] = kant t_i t_{i+1}
        P = g[0] * g[1] + g[1] * g[2] + g[2] * g[0]
        for x in fl:
            del nc[x >> 1]
        for i in range(3):
            opp = g[(i + 1) % 3]  # kanten t_{i+1} t_{i+2} ligger mitt emot t_i
            nc[G[i] >> 1] = P / opp
    return nr, nv, nc


def serie(rot, volt, cond, v):
    d1, d2 = rot[v]
    if (d1 >> 1) == (d2 >> 1):
        return None
    r1, r2 = d1 ^ 1, d2 ^ 1
    X = ny_dart(volt)
    nv = dict(volt)
    for x in (d1, d2):
        del nv[x], nv[x ^ 1]
    nv[X] = sub(volt[d2], volt[d1])
    nv[X ^ 1] = neg(nv[X])
    nr = {}
    for u, L in rot.items():
        if u == v:
            continue
        nr[u] = tuple(X if x == r1 else (X ^ 1 if x == r2 else x) for x in L)
    nc = None
    if cond is not None:
        nc = dict(cond)
        c1, c2 = cond[d1 >> 1], cond[d2 >> 1]
        del nc[d1 >> 1], nc[d2 >> 1]
        nc[X >> 1] = c1 * c2 / (c1 + c2)
    return nr, nv, nc


def ta_bort(rot, volt, cond, kanter, extra_cond=None):
    ks = set(kanter)
    nr = {}
    for u, L in rot.items():
        out = tuple(x for x in L if (x >> 1) not in ks)
        if out:
            nr[u] = out
    nv = {d: w for d, w in volt.items() if (d >> 1) not in ks}
    nc = None
    if cond is not None:
        nc = {e: c for e, c in cond.items() if e not in ks}
        if extra_cond:
            for e, c in extra_cond.items():
                nc[e] += c
    return nr, nv, nc


def parallella(rot, volt):
    t = tails(rot)
    grupper = {}
    for e in edges_of(rot):
        d = 2 * e
        a, b, w = t[d], t[d ^ 1], volt[d]
        if a == b:
            if w == (0, 0):
                continue
            w = max(w, neg(w))
            key = (a, a, w)
        else:
            if a > b:
                a, b, w = b, a, neg(w)
            key = (a, b, w)
        grupper.setdefault(key, []).append(e)
    return [g for g in grupper.values() if len(g) >= 2]


def doda_grenar(rot, volt):
    """Alla döda grenar: (nod v, kantmängd H) där H hänger fast i resten bara genom v och har ändligt lyft."""
    t = tails(rot)
    E = edges_of(rot)
    res = []
    for e in E:  # kontraktibla öglor
        d = 2 * e
        if t[d] == t[d ^ 1] and volt[d] == (0, 0):
            res.append((t[d], [e]))
    adj = {u: [] for u in rot}
    for e in E:
        d = 2 * e
        adj[t[d]].append((t[d ^ 1], d))
        adj[t[d ^ 1]].append((t[d], d ^ 1))
    for v in rot:
        sett = {v}
        for start in rot:
            if start in sett:
                continue
            comp = {start}
            st = [start]
            while st:
                u = st.pop()
                for w, _ in adj[u]:
                    if w != v and w not in comp:
                        comp.add(w)
                        st.append(w)
            sett |= comp
            if len(comp) + 1 == len(rot):
                # hela resten; då är "resten av grafen" bara v (plus öglor vid v)
                rest_kanter = [e for e in E if t[2 * e] == v and t[2 * e + 1] == v]
                if not rest_kanter:
                    continue
            H = [e for e in E if t[2 * e] in comp or t[2 * e + 1] in comp]
            # ändligt lyft: konsistenta positioner från v
            pos = {v: (0, 0)}
            ok = True
            st = [v]
            Hs = set(H)
            while st and ok:
                u = st.pop()
                for w, d in adj[u]:
                    if (d >> 1) not in Hs:
                        continue
                    q = add(pos[u], volt[d])
                    if w in pos:
                        if pos[w] != q:
                            ok = False
                            break
                    else:
                        pos[w] = q
                        st.append(w)
            if ok:
                res.append((v, H))
    return res


def drag(rot, volt, cond=None):
    """Alla drag från en karta: lista av (namn, (rot, volt, cond))."""
    ut = []
    for v, L in rot.items():
        if len(L) == 3:
            r = y_till_delta(rot, volt, cond, v)
            if r:
                ut.append(("Y->D", r))
        elif len(L) == 2:
            r = serie(rot, volt, cond, v)
            if r:
                ut.append(("serie", r))
    for f in faces(rot):
        if len(f) == 3:
            r = delta_till_y(rot, volt, cond, f)
            if r:
                ut.append(("D->Y", r))
    for g in parallella(rot, volt):
        for e2 in g[1:]:
            extra = {g[0]: cond[e2]} if cond is not None else None
            ut.append(("parallell", ta_bort(rot, volt, cond, [e2], extra)))
    for v, H in doda_grenar(rot, volt):
        ut.append(("dod gren", ta_bort(rot, volt, cond, H)))
    return ut


# ---------------------------------------------------------------- tågspår (zig-zag-vägar) för rapporten

def tagspar(rot, volt):
    """Zig-zag-vägar: omväxlande σ och σ^{-1}. Returnerar lista av (längd, klass i Z^2)."""
    s = sigma(rot)
    si = {b: a for a, b in s.items()}
    seen = set()
    ut = []
    for d0 in s:
        for p0 in (0, 1):
            if (d0, p0) in seen:
                continue
            x, p, w, n = d0, p0, (0, 0), 0
            while (x, p) not in seen:
                seen.add((x, p))
                w = add(w, volt[x])
                n += 1
                x = s[x ^ 1] if p == 0 else si[x ^ 1]
                p ^= 1
            if (x, p) == (d0, p0):
                ut.append((n, w))
    # varje väg hittas två gånger (en per riktning); behåll en per par ±
    uniq = {}
    for n, w in ut:
        key = (n, max(w, neg(w)))
        uniq[key] = uniq.get(key, 0) + 1
    return sorted((k[0], k[1], c // 2) for k, c in uniq.items())


def familjekorsningar(rot, volt):
    """Delar kanterna (= korsningar av tågspår) efter familj: a = spår med klass parallell med (1,-1),
    b = parallell med (1,1). Returnerar (antal a-a, antal b-b, antal a-b, antal övriga)."""
    s = sigma(rot)
    si = {b: a for a, b in s.items()}
    seen = set()
    pass_av = {}
    sid = 0
    for d0 in s:
        for p0 in (0, 1):
            if (d0, p0) in seen:
                continue
            orb = []
            x, p, w = d0, p0, (0, 0)
            while (x, p) not in seen:
                seen.add((x, p))
                orb.append(x)
                w = add(w, volt[x])
                x = s[x ^ 1] if p == 0 else si[x ^ 1]
                p ^= 1
            fam = "a" if (w[0] == -w[1] and w != (0, 0)) else ("b" if (w[0] == w[1] and w != (0, 0)) else "?")
            for x in orb:
                pass_av.setdefault(x >> 1, []).append((sid, fam))
            sid += 1
    # varje oorienterat spår förekommer som två banor (en per riktning), så varje kant har 4 pass
    cnt = {"aa": 0, "bb": 0, "ab": 0, "?": 0}
    for e, L in pass_av.items():
        fams = sorted(f for _, f in L)
        if fams == ["a"] * 4:
            cnt["aa"] += 1
        elif fams == ["b"] * 4:
            cnt["bb"] += 1
        elif fams == ["a", "a", "b", "b"]:
            cnt["ab"] += 1
        else:
            cnt["?"] += 1
    return cnt["aa"], cnt["bb"], cnt["ab"], cnt["?"]


# ---------------------------------------------------------------- sökningen

def packa(code):
    return array("h", [x for row in code for x in row]).tobytes()


def fran_kod(buf):
    """Bygger upp kartan ur en (packad) kanonisk kod: dart i har σ-etikett, α-etikett och spänning
    (i den gauge där alla valda lyft ligger i origo är spänningen just residualen)."""
    a = array("h")
    a.frombytes(buf)
    rows = [tuple(a[i:i + 4]) for i in range(0, len(a), 4)]
    nid = {}
    k = 0
    for i, (sl, al, rx, ry) in enumerate(rows):
        if i < al:
            nid[i] = 2 * k
            nid[al] = 2 * k + 1
            k += 1
    volt = {nid[i]: (r[2], r[3]) for i, r in enumerate(rows)}
    rot = {}
    seen = set()
    v = 0
    for i in range(len(rows)):
        if i in seen:
            continue
        cyc = []
        x = i
        while x not in seen:
            seen.add(x)
            cyc.append(nid[x])
            x = rows[x][0]
        rot[v] = tuple(cyc)
        v += 1
    return rot, volt


def sok_kompakt(B, max_tillstand=None, logg=None):
    """Som sok(), men lagrar bara 16-bytes hashvärden av de kanoniska koderna och packade koder på stacken
    (djupet först). Räcker för 40 kanter. Ingen invariant- eller exempelutskrift."""
    rot, volt = snub_karta(B)
    E0 = len(edges_of(rot))
    c0 = kanon(rot, volt)
    seen = {hashlib.blake2b(packa(c0), digest_size=16).digest()}
    stack = [packa(c0)]
    perE = {E0: 1}
    antal_drag = {}
    best = E0
    t0 = time.time()
    nxt = 50000
    while stack:
        rot, volt = fran_kod(stack.pop())
        for namn, (nr, nv, _) in drag(rot, volt):
            gg, _ = genus(nr)
            if gg != 1:
                raise RuntimeError(f"draget {namn} gav genus {gg}")
            p = packa(kanon(nr, nv))
            h = hashlib.blake2b(p, digest_size=16).digest()
            if h in seen:
                continue
            seen.add(h)
            E = len(edges_of(nr))
            perE[E] = perE.get(E, 0) + 1
            antal_drag[namn] = antal_drag.get(namn, 0) + 1
            best = min(best, E)
            stack.append(p)
        if logg and len(seen) >= nxt:
            nxt += 50000
            print(f"    ... {len(seen)} kartor, stack {len(stack)}, lägst {best}, per E {dict(sorted(perE.items()))},"
                  f" {time.time()-t0:.0f} s", file=logg, flush=True)
        if max_tillstand and len(seen) >= max_tillstand:
            break
    uttomd = not stack
    return dict(E0=E0, best=best, tillstand=len(seen), uttomd=uttomd, perE=dict(sorted(perE.items())),
                antal_drag=antal_drag, tid=time.time() - t0)


def sok(B, max_tillstand=None, logg=None, invariant=False):
    rot, volt = snub_karta(B)
    g, chi = genus(rot)
    assert g == 1 and chi == 0
    E0 = len(edges_of(rot))
    inv = {"ab": set(), "min_aa": None, "min_bb": None, "fel": 0}

    def kolla(r, v):
        aa, bb, ab, ovr = familjekorsningar(r, v)
        inv["ab"].add(ab)
        inv["min_aa"] = aa if inv["min_aa"] is None else min(inv["min_aa"], aa)
        inv["min_bb"] = bb if inv["min_bb"] is None else min(inv["min_bb"], bb)
        if ovr or aa + bb + ab != len(edges_of(r)):
            inv["fel"] += 1
    if invariant:
        kolla(rot, volt)
    start = kanon(rot, volt)
    seen = {start: len(edges_of(rot))}
    q = deque([(rot, volt, 0)])
    best = E0
    best_state = (rot, volt)
    best_depth = 0
    exempel = {E0: (rot, volt)}
    antal_drag = {}
    t0 = time.time()
    while q:
        rot, volt, dep = q.popleft()
        for namn, (nr, nv, _) in drag(rot, volt):
            gg, chi = genus(nr)
            if gg != 1:
                raise RuntimeError(f"draget {namn} gav genus {gg}")
            k = kanon(nr, nv)
            if k in seen:
                continue
            E = len(edges_of(nr))
            seen[k] = E
            if E not in exempel:
                exempel[E] = (nr, nv)
            if invariant:
                kolla(nr, nv)
            antal_drag[namn] = antal_drag.get(namn, 0) + 1
            if E < best:
                best, best_state, best_depth = E, (nr, nv), dep + 1
            q.append((nr, nv, dep + 1))
            if max_tillstand and len(seen) >= max_tillstand:
                q.clear()
                break
        if logg and len(seen) % 20000 < 5:
            print(f"    ... {len(seen)} tillstånd, kö {len(q)}, bäst {best}, {time.time()-t0:.0f} s", file=logg, flush=True)
    uttomd = not (max_tillstand and len(seen) >= max_tillstand)
    perE = {}
    for E in seen.values():
        perE[E] = perE.get(E, 0) + 1
    return dict(E0=E0, best=best, best_depth=best_depth, tillstand=len(seen), uttomd=uttomd,
                perE=dict(sorted(perE.items())), best_state=best_state, antal_drag=antal_drag,
                tid=time.time() - t0, inv=inv if invariant else None, exempel=exempel)


def forutsagt(Bcell):
    """Satsen i 9.2: [L:Λ](8 + 1/N_a + 1/N_b), Λ i cellbasen."""
    (p, q), (r, s) = Bcell
    idx = abs(p * s - q * r)

    def N(v):
        for n in range(1, 1000):
            x, y = n * v[0], n * v[1]
            # n*v ∈ Λ ?
            a = Fraction(x * s - y * r, p * s - q * r)
            b = Fraction(p * y - q * x, p * s - q * r)
            if a.denominator == 1 and b.denominator == 1:
                return n
    Na, Nb = N((1, -1)), N((1, 1))
    return Fraction(idx) * (8 + Fraction(1, Na) + Fraction(1, Nb)), idx, Na, Nb


# ---------------------------------------------------------------- determinantkontroll

def bloch_det(rot, volt, cond, theta):
    vs = sorted(rot)
    ix = {v: i for i, v in enumerate(vs)}
    t = tails(rot)
    n = len(vs)
    M = np.zeros((n, n), dtype=complex)
    for v, L in rot.items():
        for d in L:
            c = cond[d >> 1]
            w = volt[d]
            chi = cmath.exp(1j * (theta[0] * w[0] + theta[1] * w[1]))
            M[ix[v], ix[v]] += c
            M[ix[v], ix[t[d ^ 1]]] -= c * chi
    return np.linalg.det(M)


def detkontroll(B, vandringar=20, steg=40, seed=1):
    """Korta slumpvandringar med konduktanser. Kvoten mellan Bloch-determinanterna före och efter ett drag
    ska vara densamma i alla fyra slumpvalda punkter θ (en konstant pivot)."""
    rng = random.Random(seed)
    worst = 0.0
    typer = {}
    for _ in range(vandringar):
        rot, volt = snub_karta(B)
        cond = {e: rng.uniform(0.5, 2.0) for e in edges_of(rot)}
        thetas = [(rng.uniform(0, 1), rng.uniform(0, 1)) for _ in range(4)]
        vals = [bloch_det(rot, volt, cond, th) for th in thetas]
        for _ in range(steg):
            M = drag(rot, volt, cond)
            if not M:
                break
            # föredra förenklingar när de finns, så att de också prövas
            red = [m for m in M if m[0] not in ("Y->D", "D->Y")]
            namn, (nr, nv, nc) = rng.choice(red if red and rng.random() < 0.7 else M)
            nvals = [bloch_det(nr, nv, nc, th) for th in thetas]
            kv = [a / b for a, b in zip(vals, nvals)]
            spread = max(abs(k - kv[0]) for k in kv) / abs(kv[0])
            worst = max(worst, spread)
            typer[namn] = typer.get(namn, 0) + 1
            rot, volt, cond, vals = nr, nv, nc, nvals
    return worst, typer


# ---------------------------------------------------------------- fallen

FALL = [
    # namn, Λ i cellbasen (cell = 2x2 i Z^2)
    ("egen 1x1", ((1, 0), (0, 1))),
    ("diagonal index 2, <a,b>", ((1, -1), (1, 1))),
    ("2x1", ((2, 0), (0, 1))),
    ("1x2", ((1, 0), (0, 2))),
    ("snett index 3, <(3,0),(1,1)>", ((3, 0), (1, 1))),
    ("3x1", ((3, 0), (0, 1))),
    ("2x2", ((2, 0), (0, 2))),
    ("4x1", ((4, 0), (0, 1))),
    ("diagonal index 4, <2a,b>", ((2, -2), (1, 1))),
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fall", nargs="*", help="välj fall efter namn (delsträng)")
    ap.add_argument("--max-tillstand", type=int, default=None)
    ap.add_argument("--det", action="store_true", help="kör determinantkontrollen")
    ap.add_argument("--invariant", action="store_true",
                    help="räkna tågspårsfamiljernas korsningar i varje besökt karta (bevisets invariant)")
    ap.add_argument("--kompakt", action="store_true",
                    help="minnessnål uttömmande sökning (djupet först, hashvärden) för de stora fallen")
    a = ap.parse_args()
    try:  # låg prioritet: datorn används samtidigt
        import psutil
        psutil.Process().nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
    except Exception:
        pass
    fall = [f for f in FALL if f[0] not in ("2x2", "4x1")]  # de två stora: välj med --fall 2x2 4x1 --kompakt
    if a.fall:
        fall = [f for f in FALL if any(s in f[0] for s in a.fall)]
    for namn, Bc in fall:
        B = tuple((2 * x, 2 * y) for x, y in Bc)
        pred, idx, Na, Nb = forutsagt(Bc)
        print(f"== {namn}: Λ (cellbas) = {Bc}, index {idx}, N_a = {Na}, N_b = {Nb}, satsen: {pred} kanter, minimalt nät: {8*idx}")
        if a.det:
            worst, typer = detkontroll(B)
            print(f"   determinantkontroll (20 slumpvandringar à 40 drag): största relativa avvikelse {worst:.1e}, drag {typer}")
        if a.kompakt:
            r = sok_kompakt(B, a.max_tillstand, logg=sys.stdout)
            print(f"   start {r['E0']} kanter; {r['tillstand']} olika kartor "
                  f"({'uttömmande' if r['uttomd'] else 'AVBRUTEN vid taket'}), {r['tid']:.0f} s")
            print(f"   kartor per antal kanter: {r['perE']}")
            print(f"   drag som gav en ny karta: {r['antal_drag']}")
            print(f"   lägsta antal kanter: {r['best']}; stämmer med satsen: {'JA' if r['best'] == pred else 'NEJ'}",
                  flush=True)
            continue
        r = sok(B, a.max_tillstand, logg=sys.stdout, invariant=a.invariant)
        if r["inv"]:
            iv = r["inv"]
            print(f"   invariant i alla besökta kartor: a-b-korsningar {sorted(iv['ab'])} (satsen: {8*idx}),"
                  f" minst {iv['min_aa']} a-a- och {iv['min_bb']} b-b-korsningar (satsen: {idx//Na} och {idx//Nb}),"
                  f" kartor där kanterna inte delas upp så: {iv['fel']}")
        rot, volt = r["best_state"]
        deg = sorted(len(L) for L in rot.values())
        print(f"   start {r['E0']} kanter; {r['tillstand']} olika kartor ({'uttömmande' if r['uttomd'] else 'AVBRUTEN vid taket'}),"
              f" {r['tid']:.0f} s")
        print(f"   kartor per antal kanter: {r['perE']}")
        print(f"   drag som gav en ny karta: {r['antal_drag']}")
        print(f"   lägsta antal kanter: {r['best']} (första gången på djup {r['best_depth']});"
              f" {len(rot)} noder, grader {deg}")
        print("   tågspår (längd, klass i Z^2, antal) i den första kartan med E kanter:")
        for E in sorted(r["exempel"], reverse=True):
            print(f"      E = {E}: {tagspar(*r['exempel'][E])}")
        print(f"   stämmer med satsen: {'JA' if r['best'] == pred else 'NEJ'}", flush=True)


if __name__ == "__main__":
    main()
