"""
Exact arithmetic in the cyclotomic field Q(z), z = exp(i pi/12) (primitive 24th root of unity),
and truncated power series over it.  Used by bevis_formodan8.py.
Elements are polynomials in z of degree < 8 modulo Phi_24(z) = z^8 - z^4 + 1 (python-flint fmpq_poly).
"""
import flint
import mpmath as mp

PHI = flint.fmpq_poly([1, 0, 0, 0, -1, 0, 0, 0, 1])


class K:
    __slots__ = ('p',)

    def __init__(self, p):
        if isinstance(p, K):
            p = p.p
        if not isinstance(p, flint.fmpq_poly):
            p = flint.fmpq_poly(p if isinstance(p, list) else [p])
        self.p = p % PHI

    def __add__(self, o): return K(self.p + _k(o).p)
    __radd__ = __add__
    def __sub__(self, o): return K(self.p - _k(o).p)
    def __rsub__(self, o): return K(_k(o).p - self.p)
    def __neg__(self): return K(-self.p)
    def __mul__(self, o): return K(self.p * _k(o).p)
    __rmul__ = __mul__

    def inv(self):
        g, a, b = self.p.xgcd(PHI)          # a*p + b*PHI = g (constant, since PHI is irreducible)
        assert g.degree() == 0 and not self.p.is_zero()
        return K(a * (1 / g.coeffs()[0]))

    def __truediv__(self, o): return self * _k(o).inv()
    def __rtruediv__(self, o): return _k(o) * self.inv()

    def __pow__(self, n):
        if n < 0:
            return self.inv() ** (-n)
        r, b = K(1), self
        while n:
            if n & 1:
                r = r * b
            b, n = b * b, n >> 1
        return r

    def __eq__(self, o): return (self.p - _k(o).p).is_zero()
    def is_zero(self): return self.p.is_zero()

    def galois(self, j):
        """the automorphism z -> z^j (j coprime to 24)"""
        zj = Z ** j
        r = K(0)
        for c in reversed(self.p.coeffs()):
            r = r * zj + K(c)
        return r

    def conj(self): return self.galois(23)

    def to_complex(self):
        zz = mp.expjpi(mp.mpf(1) / 12)
        return sum((mp.mpf(int(c.p)) / int(c.q)) * zz**k for k, c in enumerate(self.p.coeffs()))

    def rational(self):
        """the rational number if self lies in Q, else None"""
        c = self.p.coeffs()
        if len(c) <= 1:
            return c[0] if c else flint.fmpq(0)
        return None

    def __repr__(self): return f'K({self.p})'


def _k(o):
    return o if isinstance(o, K) else K(o)


Z = K([0, 1])
I = Z ** 6                       # i
S3 = Z ** 2 - Z ** 10            # sqrt(3) = 2 cos(pi/6)


def split_sqrt3(x):
    """(a, b) rational with x = a + b sqrt3, or None if x is not in Q(sqrt3)"""
    if not (x.galois(13) == x and x.galois(23) == x):
        return None
    b = (x - x.galois(5)) / (2 * S3)       # galois 5 maps sqrt3 -> -sqrt3
    a = (x + x.galois(5)) / 2
    ra, rb = a.rational(), b.rational()
    assert ra is not None and rb is not None
    return ra, rb


# ---------------------------------------------------------------- power series in e over K
def s_add(a, b): return [x + y for x, y in zip(a, b)]
def s_scale(a, c): return [x * c for x in a]


def s_mul(a, b, N):
    out = [K(0)] * N
    for i, x in enumerate(a[:N]):
        if x.is_zero():
            continue
        for j, y in enumerate(b[:N - i]):
            out[i + j] = out[i + j] + x * y
    return out


def s_pow(a, n, N):
    r, b = [K(1)] + [K(0)] * (N - 1), a
    while n:
        if n & 1:
            r = s_mul(r, b, N)
        b, n = s_mul(b, b, N), n >> 1
    return r


def s_inv_linear(c, N):
    """series of 1/(c + e)"""
    ci = c.inv()
    out, t = [], ci
    for m in range(N):
        out.append(t)
        t = -t * ci
    return out


def s_linear(c, N):
    """series of c + e"""
    return [c, K(1)] + [K(0)] * (N - 2)
