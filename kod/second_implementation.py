"""
Second, separately written implementation for "Exact two-point resistances on the snub square lattice" (Appendix H.2).

It uses neither the substitutions of Steps 2-4 nor residues, and does not import exact_snub.py:
  1. Bloch matrix M(X, Z) from the edge list; integrand of Eq. (1) as a Laurent polynomial N(X, Z) / det M.
  2. Average over the sign changes (theta1, theta2) -> (+-theta1, +-theta2); X^k Z^l -> T_|k|(a) T_|l|(b),
     a = cos theta1, b = cos theta2, T_k = Chebyshev polynomial.
  3. For fixed a the b-integral is done exactly:  N/det = q(b) + (r1 b + r0)/(b^2 - B1 b + B0), with
        (1/pi) int_{-1}^{1} b^k / sqrt(1-b^2) db = binom(k, k/2) / 2^k  (k even, else 0),
        (1/pi) int_{-1}^{1} db / (sqrt(1-b^2) (beta - b)) = 1 / sqrt(beta^2 - 1)   (beta > 1),
     the roots beta_+- of the denominator being real and > 1 for a < 1.
  4. The theta1 integral is done numerically (tanh-sinh, mpmath).
The result is compared with the exact entries of Table 1.
Requires: sympy, mpmath.   Run:  python3 second_implementation.py [digits]     (default 55; takes a few minutes)
"""
import sys
import sympy as sp
import mpmath as mp

X, Z = sp.symbols('X Z')
EDGES = [(0, 1, (-1, 0)), (0, 1, (0, 0)), (0, 3, (0, 0)), (0, 3, (0, 1)), (1, 2, (0, 0)),
         (1, 2, (0, 1)), (2, 3, (0, 0)), (2, 3, (1, 0)), (0, 2, (0, 1)), (1, 3, (1, 0))]
s3, pi = sp.sqrt(3), sp.pi

# Table 1 of the paper: (i, j, (n1, n2)) -> exact value
TABLE = [((0, 2, (0, 1)), 2*s3/9),
         ((0, 1, (0, 0)), sp.Rational(1, 2) - s3/18),
         ((0, 2, (0, 0)), s3/9 + 1/pi),
         ((0, 2, (-1, 0)), 8*s3/9 - 1),
         ((0, 1, (-1, 1)), 11*s3/18 - sp.Rational(1, 2)),
         ((0, 0, (-1, 0)), 3 - 11*s3/9 - 1/pi),
         ((0, 1, (-1, -1)), 23*s3/18 - sp.Rational(7, 2) + 6/pi),
         ((0, 2, (1, 1)), 13*s3/9 - 6 + 13/pi),
         ((0, 0, (-1, -1)), 1 + 2*s3 - 12/pi),
         ((0, 0, (-1, 1)), 10*s3/9 - 4/pi),
         ((0, 2, (1, 0)), 28*s3/9 - 6 + 4/pi),
         ((0, 1, (-2, 0)), sp.Rational(43, 2) - 157*s3/18 - 18/pi)]

M = 5 * sp.eye(4)
for i, j, (n1, n2) in EDGES:
    M[i, j] -= X**n1 * Z**n2
    M[j, i] -= X**-n1 * Z**-n2
ADJ = M.adjugate().applyfunc(sp.expand)
DET = sp.expand(M.det())


def laurent_terms(expr):
    out = {}
    for term in sp.Add.make_args(sp.expand(expr)):
        c, kx = term.as_coeff_exponent(X)
        c, kz = c.as_coeff_exponent(Z)
        out[(int(kx), int(kz))] = out.get((int(kx), int(kz)), 0) + c
    return out


def chebyshev_form(expr):
    """sign-averaged Laurent polynomial -> polynomial in a, b via X^k Z^l -> T_|k|(a) T_|l|(b)"""
    a, b = sp.symbols('a b')
    terms = laurent_terms(expr)
    avg = {}                                             # average over the four sign changes
    for (k, l), c in terms.items():
        for kk, ll in ((k, l), (-k, l), (k, -l), (-k, -l)):
            avg[(kk, ll)] = avg.get((kk, ll), 0) + sp.Rational(1, 4) * c
    poly = sum(c * sp.chebyshevt(abs(k), a) * sp.chebyshevt(abs(l), b) for (k, l), c in avg.items())
    return sp.Poly(sp.expand(poly), b), a, b


def prepare(i, j, n):
    ph = X**n[0] * Z**n[1]
    num = ADJ[i, i] + ADJ[j, j] - ADJ[i, j] / ph - ADJ[j, i] * ph
    Pn, a, b = chebyshev_form(num)
    Pd, _, _ = chebyshev_form(DET)
    assert Pd.degree() == 2
    # coefficients as functions of a, ready for mpmath
    fn = [sp.lambdify(a, Pn.coeff_monomial(b**k), 'mpmath') for k in range(Pn.degree() + 1)]
    fd = [sp.lambdify(a, Pd.coeff_monomial(b**k), 'mpmath') for k in range(3)]
    return fn, fd


def b_mean(fn, fd, av):
    """(1/pi) int_{-1}^{1} N(av, b) / det(av, b) db / sqrt(1 - b^2), exactly in b"""
    num = [f(av) for f in fn]                            # ascending coefficients in b
    d0, d1, d2 = [f(av) for f in fd]
    # polynomial division num / det
    rem = num[:]
    quo = [mp.mpf(0)] * max(len(num) - 2, 1)
    for k in range(len(num) - 1, 1, -1):
        c = rem[k] / d2
        quo[k - 2] = c
        rem[k] -= c * d2; rem[k - 1] -= c * d1; rem[k - 2] -= c * d0
    r0, r1 = rem[0], (rem[1] if len(rem) > 1 else 0)
    total = mp.fsum(c * mp.binomial(k, k // 2) / mp.mpf(2)**k for k, c in enumerate(quo) if k % 2 == 0)
    disc = mp.sqrt((d1/d2)**2 - 4*d0/d2)
    bp, bm = (-d1/d2 + disc) / 2, (-d1/d2 - disc) / 2    # roots of det in b
    # (r1 b + r0) / (d2 (b - bp)(b - bm)) = cp/(bp - b) + cm/(bm - b)
    cp = -(r1*bp + r0) / (d2*(bp - bm))
    cm = -(r1*bm + r0) / (d2*(bm - bp))
    for beta, c in ((bp, cp), (bm, cm)):
        if abs(beta) <= 1:
            raise ValueError("denominator root inside [-1, 1]")
        total += c * mp.sign(beta) / mp.sqrt(beta**2 - 1)
    return total


def resistance(i, j, n, dps):
    fn, fd = prepare(i, j, n)
    old = mp.mp.dps
    mp.mp.dps = 2*dps + 20                               # strong cancellation near theta1 = 0
    try:
        eps = mp.mpf(10)**-(dps + 10)                     # tanh-sinh nodes closer to the ends than eps carry no weight;
        # at theta1 = pi the two roots of det in b coincide (D(-1, b) = (b - 11)^2), so stay off the endpoint
        f = lambda t: b_mean(fn, fd, mp.cos(min(max(t, eps), mp.pi - eps)))
        val = mp.quad(f, [0, mp.pi/16, mp.pi/4, mp.pi/2, mp.pi]) / mp.pi
    finally:
        mp.mp.dps = old
    return val


if __name__ == '__main__':
    digits = int(sys.argv[1]) if len(sys.argv) > 1 else 55
    worst = 10**9
    for (i, j, n), ex in TABLE:
        mp.mp.dps = digits + 10
        val = resistance(i, j, n, digits)
        mp.mp.dps = digits + 10
        err = abs(val - mp.mpf(sp.N(ex, digits + 10)))
        ok = int(-mp.log10(err)) if err > 0 else digits + 10
        worst = min(worst, ok)
        print(f"R_{i+1}{j+1}{n}  exact = {str(ex):34s} numeric = {mp.nstr(val, 20)}   agreement: {ok} digits", flush=True)
    print(f"worst agreement: {worst} digits")
    assert worst >= digits, f"agreement of {worst} digits is below the required {digits}"
