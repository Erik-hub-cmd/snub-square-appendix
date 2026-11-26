"""
Residy- och Galoismaskineriet från forslag/bevis_formodan8_direkt.md, för ett godtyckligt gitter vars spektralkurva
har samma parametrisering X(s), Z_in(s) (24:e enhetsrötter): snub square och Cairo-gittret (samma determinant).

Maskin(ADJd, DDZd, nmax) tar adjunkten (Laurentpolynom i X, Z, som ordböcker {(a, b): heltalskoefficient}) och
Z dZ det L, och ger exakt R_ij(n) = c0 + c1 √3 + c2/π som tre rationella tal:
    R = (1/π) Re ∫_1^∞ f ds = A - B/(2π),  f = g_ii + g_jj - 2 g Φ,  A = -(i/2) Σ_q ℓ_q a_q,  B = Σ_q b_q.
Används av bevis_cairo.py och sandhog_y.py.  Exakt aritmetik i Q(z), z = exp(iπ/12), via k24.py.
"""
import os, sys
import flint

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from k24 import K, Z, I, S3, split_sqrt3

EX = {1: 1, 5: 1, 7: 1, 11: 1, 13: -1, 17: -1, 19: -1, 23: -1}        # X(s)    = prod (s - z^k)^EX[k]
EZ = {1: 1, 5: 1, 19: 1, 23: 1, 7: -1, 11: -1, 13: -1, 17: -1}        # Z_in(s) = prod (s - z^k)^EZ[k]
ROOTS = [1, 5, 7, 11, 13, 17, 19, 23]
ZK = {k: Z**k for k in range(24)}


def Q(a, b=1):
    return K(flint.fmpq(a, b))


# ---- Laurentserier i e = s - q över K:  (v, [c0, c1, ...]) = Σ c_k e^{v+k}
class L:
    __slots__ = ('v', 'c')

    def __init__(self, v, c):
        k = 0
        while k < len(c) and c[k].is_zero():
            k += 1
        self.v, self.c = v + k, c[k:]

    def __mul__(self, o):
        N = min(len(self.c), len(o.c))
        out = [K(0)] * N
        for i in range(N):
            a = self.c[i]
            if a.is_zero():
                continue
            for j in range(N - i):
                out[i + j] = out[i + j] + a * o.c[j]
        return L(self.v + o.v, out)

    def __add__(self, o):
        v = min(self.v, o.v)
        top = min(self.v + len(self.c), o.v + len(o.c))           # exakt kända ordningar < top
        out = [K(0)] * (top - v)
        for s_ in (self, o):
            for k, x in enumerate(s_.c):
                if s_.v + k < top:
                    out[s_.v + k - v] = out[s_.v + k - v] + x
        return L(v, out)

    def scale(self, x):
        return L(self.v, [y * x for y in self.c])

    def inv(self):
        assert self.c, 'ingen känd ledande term'
        a0i = self.c[0].inv()
        N = len(self.c)
        out = [a0i]
        for n in range(1, N):
            acc = K(0)
            for k in range(1, n + 1):
                acc = acc + self.c[k] * out[n - k]
            out.append(-acc * a0i)
        return L(-self.v, out)

    def coef(self, n):
        k = n - self.v
        assert k < len(self.c), 'för kort serie'
        return self.c[k] if k >= 0 else K(0)


def lin(c, N):                     # c + e
    return L(0, [c, K(1)] + [K(0)] * (N - 2)) if not c.is_zero() else L(1, [K(1)] + [K(0)] * (N - 1))


def lpow(x, n, cache):
    if n in cache:
        return cache[n]
    if n == 0:
        r = L(0, [K(1)] + [K(0)] * (len(x.c) - 1))
    elif n > 0:
        r = lpow(x, n - 1, cache) * x
    else:
        r = lpow(x, n + 1, cache) * cache['inv']
    cache[n] = r
    return r


def ell(k):
    """log(-z^k) = iπ ℓ (principalgren)"""
    m = (k + 12) % 24
    return flint.fmpq(m if m <= 12 else m - 24, 12)


def laurent_dict(expr, X, Zs, shift=8):
    import sympy as sp
    P = sp.Poly(sp.expand(expr * X**shift * Zs**shift), X, Zs)
    return {(a - shift, b - shift): int(c) for (a, b), c in P.terms()}


def width(d):
    """Grad (i s) för täljare och nämnare av L(X(s), Z_in(s)) för ett Laurentpolynom L = Σ c X^a Z^b:
    4 · (bredden av X-exponenterna + bredden av Z-exponenterna); se bevis_formodan8_direkt.md, (P1)."""
    a = [k[0] for k in d]
    b = [k[1] for k in d]
    return 4 * ((max(a) - min(a)) + (max(b) - min(b)))


class Machine:
    def __init__(self, ADJd, DDZd, nmax):
        self.ADJd, self.DDZd = ADJd, DDZd
        self.keys = sorted(ADJd)
        self.zord = {ij: min(b for (_, b) in d) for ij, d in ADJd.items()}
        self.nmax = nmax
        self.LOC = {k: self.local(ZK[k], nmax) for k in ROOTS}

    def local(self, q, N):
        """Laurentdata vid s = q (q = 0 eller en primitiv rot): X, Z_in och g_kl = r_kl (X'/X)/i för alla poster kl"""
        Xs = Zs = L(0, [K(1)] + [K(0)] * (N - 1))
        dl = L(0, [K(0)] * N)
        for k in ROOTS:
            f = lin(q - ZK[k], N)
            Xs = Xs * (f if EX[k] > 0 else f.inv())
            Zs = Zs * (f if EZ[k] > 0 else f.inv())
            dl = dl + f.inv().scale(K(EX[k]))
        cx, cz = {'inv': Xs.inv()}, {'inv': Zs.inv()}

        def ev(d):
            out = L(0, [K(0)] * N)
            for (a, b), c in d.items():
                out = out + (lpow(Xs, a, cx) * lpow(Zs, b, cz)).scale(K(c))
            return out
        den = (Zs * ev(self.DDZd)).inv()
        g = {ij: ev(self.ADJd[ij]) * den * dl.scale(I.inv()) for ij in self.keys}
        return Xs, Zs, cx, cz, g

    def entry(self, i, j, n):
        """vilken post av G och vilka exponenter som används (inget bidrag från Z = 0): G_ji X^n1 Z^n2 om n2 >= 1,
        G_ij X^-n1 Z^-n2 om n2 <= -1, och för n2 = 0 den av G_ji, G_ij vars Z-ordning är >= -1"""
        n1, n2 = n
        if n2 >= 1:
            cands = [((j, i), n1, n2)]
        elif n2 <= -1:
            cands = [((i, j), -n1, -n2)]
        else:
            cands = [((j, i), n1, 0), ((i, j), -n1, 0)]
        for ij, e1, e2 in cands:
            if self.zord[ij] + e2 >= -1 and self.zord[(i, i)] >= -1 and self.zord[(j, j)] >= -1:
                return ij, e1, e2
        raise ValueError(f'ingen post med Z-ordning >= -1 för paret {(i, j, n)}')

    def resistance(self, i, j, n):
        """(c0, c1, c2) med R_ij(n) = c0 + c1 √3 + c2/π, plus kontrollerna Σ a_q = 0 och a_q imaginär, eller None"""
        ij, e1, e2 = self.entry(i, j, n)
        A = K(0)
        B = K(0)
        asum = K(0)
        allimag = True
        for k in ROOTS:
            Xs, Zs, cx, cz, g = self.LOC[k]
            f = g[(i, i)] + g[(j, j)] + (g[ij] * lpow(Xs, e1, cx) * lpow(Zs, e2, cz)).scale(K(-2))
            a = f.coef(-1)
            asum = asum + a
            allimag &= (a + a.conj()).is_zero()
            A = A + K(ell(k)) * a
            qi = ZK[k].inv()
            for m in range(1, -f.v):
                B = B + f.coef(-1 - m) * qi**m * Q((-1)**(m - 1), m)
        A = A * I * Q(-1, 2)
        sA, rB = split_sqrt3(A), B.rational()
        ok = asum.is_zero() and allimag and sA is not None and rB is not None
        if not ok:
            return None
        return sA[0], sA[1], -rB / 2


def partial_fraction_check(M, pts_from=2):
    """Kontroll (D2) för alla poster g_kl: g_kl = Σ_q Σ_{m<=3} c_qm/(s-q)^m med q i {0} ∪ P, utan polynomdel, med
    samma residy c i s = 0 för alla poster.  Gradgräns: täljare och nämnare av g_kl har grad <= w_adj + w_ddz + 11
    resp. + 12, partialbråkssumman har nämnare av grad 27 och täljare av grad <= 26, så skillnadens täljare har
    grad <= w_adj + w_ddz + 38; identiteten kontrolleras i en punkt mer än så.  Returnerar (ok, gräns, residyn i 0)."""
    PF = {}
    for q in [0] + ROOTS:
        qq = K(0) if q == 0 else ZK[q]
        _, _, _, _, g = M.local(qq, 10)
        for ij, ser in g.items():
            if ser.v < -3:
                return False, None, None
            PF.setdefault(ij, []).append((qq, [ser.coef(-m) for m in range(1, 4)]))
    res0 = {ij: (PF[ij][0][1][0], PF[ij][0][1][1].is_zero() and PF[ij][0][1][2].is_zero()) for ij in PF}
    vals = {str(v[0].p) for v in res0.values()}
    ok0 = len(vals) == 1 and all(v[1] for v in res0.values())
    w_ddz = width(M.DDZd)
    bound = max(width(M.ADJd[ij]) for ij in PF) + w_ddz + 38
    okpf = True
    for ij in PF:
        for k in range(bound + 1):
            sv = K(flint.fmpq(k + pts_from, 1))
            Xv = Zv = K(1)
            dl = K(0)
            for r in ROOTS:
                Xv = Xv * (sv - ZK[r])**EX[r]
                Zv = Zv * (sv - ZK[r])**EZ[r]
                dl = dl + K(EX[r]) / (sv - ZK[r])
            Xi, Zi = Xv.inv(), Zv.inv()

            def evK(d):
                return sum((K(c) * (Xv**a if a >= 0 else Xi**(-a)) * (Zv**b if b >= 0 else Zi**(-b))
                            for (a, b), c in d.items()), K(0))
            lhs = evK(M.ADJd[ij]) / (Zv * evK(M.DDZd)) * dl / I
            rhs = K(0)
            for qq, cs in PF[ij]:
                for m, c in enumerate(cs, 1):
                    rhs = rhs + c / (sv - qq)**m
            okpf &= (lhs == rhs)
            if not okpf:
                return False, bound, res0
    return ok0 and okpf, bound, res0[next(iter(res0))][0]
