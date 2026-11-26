"""Separat kontroll av Sats A/B: egna, separat skrivna verktyg (används av granskning4_k2 och granskning4_k3).
- grund(M, a, b, n): R(a@cell 0, b@cell n) = (2π)^-2 ∫∫ w* M^-1 w, w = e_a − X^{-n1} Z^{-n2} e_b (inre integral med residyer,
  yttre med tanh-sinh), M = Bloch-matris (sympy, symboler X, Z, ev. t).
- kenyon(v0, v1, t): Kenyons Green-funktion på Z²_t som NUMERISK konturintegral −(1/(8π² i)) ∮_C g_{v1}(z) log z dz
  längs en kontur av Kenyons typ (Sats 4.2/7.1 i arXiv math-ph/0202018), helt utan residyer.
"""
import sympy as sp
import mpmath as mp
from functools import lru_cache

X, Z, t = sp.symbols('X Z t')
TYPES = [(0, 0), (0, 1), (1, 0), (1, 1)]


# ---------------- grundintegralen ----------------
def _laurent_in_Z(expr):
    d = {}
    for term in sp.Add.make_args(sp.expand(expr)):
        c, p = term.as_coeff_exponent(Z)
        d[int(p)] = d.get(int(p), 0) + c
    return {k: v for k, v in d.items() if v != 0}


class Grund:
    def __init__(self, M, tval=None):
        self.M = M; self.tval = tval
        self.D = _laurent_in_Z(sp.expand(sp.cancel(M.det() * X**4 * Z**4)) / (X**4 * Z**4))
        self.adj = M.adjugate()
        self.fD = {k: sp.lambdify((X, t), v, 'mpmath') for k, v in self.D.items()}

    def R(self, a, b, n, dps=30):
        ph = X**n[0] * Z**n[1]
        N = sp.expand(self.adj[a, a] + self.adj[b, b] - self.adj[a, b] / ph - self.adj[b, a] * ph)
        Nd = _laurent_in_Z(N)
        fN = {k: sp.lambdify((X, t), v, 'mpmath') for k, v in Nd.items()}
        lowD, lowN = min(self.D), min(Nd)
        degD = max(self.D) - lowD; degN = max(Nd) - lowN
        e = lowN - lowD - 1                                   # F/Z = Ñ(Z) Z^e / D̃(Z)
        pal = sorted(self.D) == [-2, -1, 0, 1, 2] and sp.simplify(self.D[2] - 1) == 0 and sp.simplify(self.D[-2] - 1) == 0 and sp.simplify(self.D[1] - self.D[-1]) == 0
        mp.mp.dps = 2 * dps + 10
        tv = mp.mpf(0) if self.tval is None else self.tval

        def inner(x):
            x = min(max(x, mp.mpf(10) ** (-(dps + 3))), 2 * mp.pi - mp.mpf(10) ** (-(dps + 3)))
            Xv = mp.expj(x)
            P = [mp.mpc(self.fD.get(k + lowD, lambda *_: 0)(Xv, tv)) for k in range(degD + 1)]   # stigande potenser
            Q = [mp.mpc(fN.get(k + lowN, lambda *_: 0)(Xv, tv)) for k in range(degN + 1)]
            if pal:                                           # D̃ palindromisk av grad 4: s = Z + 1/Z löser s² + c1 s + c0 − 2 = 0
                c0, c1 = P[2], P[1]
                sq = mp.sqrt(c1**2 - 4 * (c0 - 2)); roots = []
                for sv in ((-c1 + sq) / 2, (-c1 - sq) / 2):
                    rr = mp.sqrt(sv**2 - 4); roots += [(sv + rr) / 2, (sv - rr) / 2]
            else:
                roots = mp.polyroots(P[::-1], maxsteps=3000, extraprec=2 * dps + 40)
            dP = [k * P[k] for k in range(1, len(P))]
            tot = mp.mpc(0)
            for r in roots:
                if abs(r) < 1:
                    tot += mp.polyval(Q[::-1], r) * r**e / mp.polyval(dP[::-1], r)
            if e < 0:                                         # pol i 0 av ordning −e: koefficient för Z^{−e−1} i Ñ/D̃
                k = -e - 1; h = []
                for i in range(k + 1):
                    s = (Q[i] if i < len(Q) else 0) - mp.fsum(P[r_] * h[i - r_] for r_ in range(1, min(i, len(P) - 1) + 1))
                    h.append(s / P[0])
                tot += h[k]
            return tot.real

        v = mp.quad(inner, [0, mp.mpf(1) / 2, 2, mp.mpf(4.5), 2 * mp.pi - mp.mpf(1) / 2, 2 * mp.pi]) / (2 * mp.pi)   # undvik X = −1 som brytpunkt
        mp.mp.dps = dps
        return +v


# ---------------- Kenyon, numerisk konturintegral ----------------
def Phi_of(tv):
    return 2 * mp.atan(tv) - mp.pi / 2


def alpha(k, Ph): return mp.mpf(0) if k % 2 == 0 else -Ph            # spår u = x − y, korsas från k till k+1
def beta(k, Ph): return mp.pi / 2 - Ph if k % 2 == 0 else mp.pi / 2   # spår w = x + y


def steps(p0, p1, Ph):
    S = []
    for (a, b, f) in ((p0[0] - p0[1], p1[0] - p1[1], alpha), (p0[0] + p0[1], p1[0] + p1[1], beta)):
        if b > a: S += [mp.expj(f(k, Ph)) for k in range(a, b)]
        else: S += [-mp.expj(f(k, Ph)) for k in range(b, a)]
    return S


def pos(p, Ph):
    return mp.fsum(steps((0, 0), p, Ph)) if p != (0, 0) else mp.mpc(0)


def cond_geom(p, r, Ph):
    """tan(halva rombvinkeln i primalnoden p) för kanten p–r, ur rombinbäddningen."""
    S = steps(p, r, Ph)                       # två steg: p -> dual -> r
    assert len(S) == 2
    # rombens kanter i p är S0 (p -> dual 1) och S1 (p -> dual 2 är parallell med dual 1 -> r):
    ang = abs(mp.arg(S[0] / S[1]))
    return mp.tan(ang / 2)


def kenyon(p0, p1, tv, dps=30, eps=0.15):
    """G(p0, p1) = −(1/(8π² i)) ∮_C g_{p1}(z) log z dz, g_{p1}(z) = (1/z) Π (z + s)/(z − s), s = stegvektorer p0 -> p1."""
    mp.mp.dps = dps + 10
    if p0 == p1: return mp.mpf(0)
    Ph = Phi_of(mp.mpf(tv))
    S = steps(p0, p1, Ph)
    disp = mp.fsum(S); th0 = mp.arg(disp)
    # alla poler måste ligga i (θ0 − π + 2ε, θ0 + π − 2ε)
    for s in S:
        d = mp.arg(s / mp.expj(th0))
        assert abs(d) < mp.pi - 2 * eps, ('pol nära snittet', p0, p1)

    def g(z):
        v = 1 / z
        for s in S: v *= (z + s) / (z - s)
        return v

    def lg(z):                                 # log z med arg i (θ0 − π, θ0 + π)
        return mp.log(abs(z)) + 1j * (th0 + mp.arg(z / mp.expj(th0)))

    Rb, rb = mp.mpf(3), mp.mpf(1) / 3
    a1, a2 = th0 - mp.pi + eps, th0 + mp.pi - eps
    f_out = lambda th: g(Rb * mp.expj(th)) * lg(Rb * mp.expj(th)) * 1j * Rb * mp.expj(th)
    f_in = lambda th: g(rb * mp.expj(th)) * lg(rb * mp.expj(th)) * 1j * rb * mp.expj(th)
    f_rad2 = lambda r: g(r * mp.expj(a2)) * lg(r * mp.expj(a2)) * mp.expj(a2)
    f_rad1 = lambda r: g(r * mp.expj(a1)) * lg(r * mp.expj(a1)) * mp.expj(a1)
    I = (mp.quad(f_out, mp.linspace(a1, a2, 9)) - mp.quad(f_rad2, [rb, mp.mpf(0.8), 1, mp.mpf(1.25), Rb])
         - mp.quad(f_in, mp.linspace(a1, a2, 9)) + mp.quad(f_rad1, [rb, mp.mpf(0.8), 1, mp.mpf(1.25), Rb]))
    G = -I / (8 * mp.pi**2 * 1j)
    mp.mp.dps = dps
    return G


# ---------------- egna Bloch-matriser (samma som i k1_identitet.py) ----------------
def _typ(p):
    return TYPES.index((p[0] % 2, p[1] % 2)), (p[0] // 2, p[1] // 2)


def _has_diag(p, r):
    if abs(p[0] - r[0]) != 1 or abs(p[1] - r[1]) != 1: return False
    c = (min(p[0], r[0]) % 2, min(p[1], r[1]) % 2); main = (p[0] - r[0]) == (p[1] - r[1])
    return (c == (0, 0) and main) or (c == (1, 1) and not main)


def neighbours(p, lam, h):
    out = []
    for d in [(1, 0), (-1, 0), (0, 1), (0, -1)]:
        r = (p[0] + d[0], p[1] + d[1])
        out.append((r, h(min(p[0], r[0]), p[1]) if d[1] == 0 else 1))
    if lam is not None:
        for d in [(1, 1), (-1, -1), (1, -1), (-1, 1)]:
            r = (p[0] + d[0], p[1] + d[1])
            if _has_diag(p, r): out.append((r, lam))
    return out


def bloch(lam, h):
    M = sp.zeros(4, 4)
    for a, p in enumerate(TYPES):
        for r, w in neighbours(p, lam, h):
            b, c = _typ(r)
            M[a, a] += w; M[a, b] -= w * X**c[0] * Z**c[1]
    return M


def hz(tv):
    return lambda x, y: 1 / tv if (x + y) % 2 == 0 else tv
