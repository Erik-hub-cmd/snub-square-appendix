"""
Exact two-point resistances on the snub square lattice (3^2.4.3.4), unit resistors.
Ancillary file for "Exact two-point resistances on the snub square lattice".

Pipeline (Section 5 of the paper):
  1. Bloch Laplacian M(theta) from the explicit edge list (Section 2).
  2. Integrand of Eq. (1), averaged over (theta1,theta2) -> (+-theta1,+-theta2),
     written as N(a,b) / (4 D(a,b)),  a = cos(theta1), b = cos(theta2).
  3. Division by D in b; polynomial part -> rational number.
  4. Remainder -> partial fractions in b, Eq. (3), substitution a=(5-m^2)/(m^2-1),
     m = sqrt3*mu, Euler substitutions -> integrals of rational functions of tau.
  5. Exact integration with SymPy; numerical cross-check of Eq. (1) with mpmath.
Requires: sympy, mpmath.   Run:  python3 exact_snub.py
      python3 exact_snub.py cairo   -> the same reduction for the Cairo lattice (needs cairo_lattice.py):
                                       det is identical, only the numerator changes
"""
import itertools
import sympy as sp
import mpmath as mp

X, Z = sp.symbols('X Z')                 # X = e^{i theta1}, Z = e^{i theta2}
t1, t2 = sp.symbols('theta1 theta2', real=True)
a, b = sp.symbols('a b', real=True)
m, mu, tau = sp.symbols('m mu tau', positive=True)
s3 = sp.sqrt(3)

# vertices 1..4 -> 0..3 ; edge (i, j, (n1, n2)) joins (i, cell 0) to (j, cell n)
EDGES = [(0, 1, (-1, 0)), (0, 1, (0, 0)), (0, 3, (0, 0)), (0, 3, (0, 1)), (1, 2, (0, 0)),
         (1, 2, (0, 1)), (2, 3, (0, 0)), (2, 3, (1, 0)),      # triangle-square edges
         (0, 2, (0, 1)), (1, 3, (1, 0))]                       # triangle-triangle edges

M = 5 * sp.eye(4)
for i, j, (n1, n2) in EDGES:
    ph = X**n1 * Z**n2
    M[i, j] -= ph
    M[j, i] -= 1 / ph
DET = sp.expand(M.det())
ADJ = M.adjugate().applyfunc(sp.expand)
D = a**2 + b**2 - 14*a*b - 36*(a + b) + 84

def to_ab(expr):
    """average over sign changes of theta1, theta2; express in a = cos theta1, b = cos theta2"""
    f = sp.expand(expr.subs({X: sp.exp(sp.I*t1), Z: sp.exp(sp.I*t2)}))
    g = (f + f.subs(t1, -t1) + f.subs(t2, -t2) + f.subs({t1: -t1, t2: -t2}, simultaneous=True)) / 4
    g = sp.expand_trig(sp.expand(sp.expand_complex(sp.expand(g.rewrite(sp.cos)))))
    g = sp.expand(g.subs({sp.sin(t1)**2: 1 - sp.cos(t1)**2, sp.sin(t2)**2: 1 - sp.cos(t2)**2}))
    return sp.expand(sp.simplify(g.subs({sp.cos(t1): a, sp.cos(t2): b})))

assert sp.expand(to_ab(DET) - 4*D) == 0, "Lemma 2 failed"

def numerator(i, j, n):
    ph = X**n[0] * Z**n[1]
    return sp.expand(ADJ[i, i] + ADJ[j, j] - ADJ[i, j] / ph - ADJ[j, i] * ph)

def arcsine_mean(poly, var):
    """(1/pi) int_{-1}^{1} var^k / sqrt(1 - var^2) = binom(k, k/2) / 2^k  (k even)"""
    P = sp.Poly(sp.expand(poly), var)
    return sp.expand(sum(c * sp.binomial(k, k // 2) / sp.Integer(2)**k
                         for (k,), c in P.terms() if k % 2 == 0))

def exact_resistance(i, j, n, N=None):
    """exact value of Eq. (1); N (optional) = symmetrised numerator as a polynomial in a, b, for any
    Bloch Laplacian with determinant 4 D(a, b) (e.g. the Cairo lattice, see cairo_lattice.py)"""
    if N is None:
        N = to_ab(numerator(i, j, n))                   # integrand = N / (4 D)
    q, r = sp.div(sp.Poly(sp.expand(N / 4), b), sp.Poly(D, b))
    total = arcsine_mean(arcsine_mean(q.as_expr(), b), a)
    r1, r0 = r.coeff_monomial(b), r.coeff_monomial(1)
    am, rho = (5 - m**2) / (m**2 - 1), 4*m / (m**2 - 1)
    for sg in (+1, -1):
        bpm = (11*m**2 + 17 + sg*16*s3*m) / (m**2 - 1)            # roots b_+-
        c = -(r1.subs(a, am)*bpm + r0.subs(a, am)) / (sg*8*s3*rho)  # c/(b_+- - b) form
        # (1/pi) int da/sqrt(1-a^2) * c/sqrt(b^2-1), radicals as in Steps 2-3 of the proof
        f_m = c * (2*sp.sqrt(2)*m/(m**2 - 1)) / (2*sp.sqrt(2)*(s3*m + sg*2)/(m**2 - 1)) / (m + sg*s3)
        f_mu = f_m.subs(m, s3*mu)               # dm = sqrt3 dmu ; remaining radical = sqrt3*sqrt(...)
        if sg == 1:
            mut = (tau**2 + 3) / (tau**2 - 5); rad = (mut - 1)*tau; lo, hi = sp.oo, sp.sqrt(5)
        else:
            mut = (tau**2 + 3) / (5 - tau**2); rad = (mut + 1)*tau; lo, hi = 1, sp.sqrt(5)
        g = sp.cancel(sp.together(f_mu.subs(mu, mut) / rad * sp.diff(mut, tau)))
        F = sp.integrate(sp.apart(sp.factor(g, extension=[s3]), tau, extension=[s3]), tau)
        val = sp.limit(F, tau, hi) - sp.limit(F, tau, lo)
        total += val / sp.pi
    total = total.subs(sp.atan(7*s3/3), sp.atan(s3/6) + sp.pi/3)   # arctan identity, Sec. 7
    raw = sp.simplify(total)
    out = sp.expand(raw)                     # exact derivation only; no numerical recognition
    if out.atoms(sp.Function):
        raise ArithmeticError("SymPy left transcendental functions in the result (%s); "
                              "the value is not derived exactly, refusing to guess it numerically" % out)
    return out

def _roots(P):
    """roots of sum_k P[k] z^k: mpmath polyroots, falling back to companion-matrix eigenvalues + Newton"""
    try:
        return mp.polyroots(P[::-1], maxsteps=400, extraprec=3 * mp.mp.dps)
    except mp.libmp.libhyper.NoConvergence:
        n = len(P) - 1
        C = mp.matrix(n, n)
        for k in range(1, n):
            C[k, k - 1] = 1
        for k in range(n):
            C[k, n - 1] = -P[k] / P[n]
        out = []
        cf, dcf = P[::-1], [k * P[k] for k in range(1, len(P))][::-1]
        for z in mp.eig(C, left=False, right=False):
            for _ in range(80):
                step = mp.polyval(cf, z) / mp.polyval(dcf, z)
                z -= step
                if abs(step) < mp.mpf(10) ** (-mp.mp.dps + 5):
                    break
            out.append(z)
        return out

def numeric_resistance(i, j, n, dps=30, Nn=None):
    """Eq. (1) numerically: theta2-integral by residues, theta1-integral by tanh-sinh.
    Nn (optional) = Laurent numerator in X, Z for another lattice with the same determinant."""
    mp.mp.dps = dps + 20
    if Nn is None:
        Nn = numerator(i, j, n)
    def coeffs(expr):
        d = {}
        for term in sp.Add.make_args(sp.expand(expr)):
            c, p = term.as_coeff_exponent(Z); d[int(p)] = d.get(int(p), 0) + c
        return d
    Dc, Nc = coeffs(DET), coeffs(Nn)
    low = -min(Dc)
    Pf = [sp.lambdify(X, Dc.get(k - low, 0), 'mpmath') for k in range(max(Dc) + low + 1)]
    Qf = [sp.lambdify(X, Nc.get(k - low, 0), 'mpmath') for k in range(max(Nc) + low + 1)]
    def inner(x):
        x = max(x, mp.mpf(10)**(-(dps // 2 + 3)))
        Xv = mp.expj(x)
        P = [mp.mpc(f(Xv)) for f in Pf]; Q = [mp.mpc(f(Xv)) for f in Qf]
        dP = [k*P[k] for k in range(1, len(P))]
        tot = mp.polyval(Q[::-1], 0) / P[0]
        for z in _roots(P):
            if abs(z) < 1:
                tot += mp.polyval(Q[::-1], z) / (mp.polyval(dP[::-1], z) * z)
        return tot.real
    v = mp.quad(inner, [0, mp.pi/8, mp.pi/2, mp.pi]) / mp.pi
    mp.mp.dps = dps
    return +v

def cairo_table(dps=40):
    """exact + numerical resistances on the Cairo lattice for all pairs at distance < 3, and g3 - g4"""
    import cairo_lattice as cl
    pos, deg, edges = cl.build()
    L = cl.bloch_laplacian(deg, edges)
    assert sp.expand(L.det() - DET) == 0, "Cairo determinant differs"
    print("Cairo: det L = det M (snub square) exactly")
    for src in (4, 0):
        for y, d2, orb in cl.targets(pos, src):
            Nn = cl.numerator(src, y[0], y[1])
            ex = exact_resistance(src, y[0], y[1], N=cl.laurent_to_ab(Nn, a, b))
            nu = numeric_resistance(src, y[0], y[1], dps=dps, Nn=Nn)
            err = abs(nu - mp.mpf(str(sp.N(ex, dps + 10))))
            print(f"R_{src+1}{y[0]+1}{y[1]}  deg {deg[src]}-{deg[y[0]]}  d^2 = {str(d2):10s} x{len(orb)}  "
                  f"exact = {str(ex):32s} = {sp.N(ex, 25)}   |exact - numeric| = {mp.nstr(err, 3)}")
            assert err < mp.mpf(10)**-28, "exact and numerical Cairo values differ"
    Nd = cl.diag_difference(0, 4)
    g34 = exact_resistance(0, 0, (0, 0), N=cl.laurent_to_ab(Nd, a, b))
    print("g3 - g4 =", g34)
    assert sp.simplify(g34 - s3/18) == 0, "g3 - g4 differs from sqrt(3)/18"


if __name__ == '__main__':
    import sys
    if sys.argv[1:] == ['cairo']:
        cairo_table()
        sys.exit()
    print("Lemma 2: det M = 4 D(a,b) with D =", D)
    Rtt = exact_resistance(0, 2, (0, 1))
    DPS = 40                                 # numerical precision of the cross-check; the paper claims >= 30 digits
    TOL = mp.mpf(10)**-31
    Rtt_num = numeric_resistance(0, 2, (0, 1), dps=DPS)
    print("Theorem 1: R_tt =", Rtt, " numeric:", Rtt_num)
    assert abs(Rtt_num - mp.mpf(str(sp.N(Rtt, 50)))) < TOL, "Theorem 1 differs from the numerical value"
    rows = [((0, 1, (0, 0)), 'triangle-square edge'), ((0, 2, (0, 0)), 'square diagonal'),
            ((0, 2, (-1, 0)), 'rhombus long diagonal'), ((0, 1, (-1, 1)), ''), ((0, 0, (-1, 0)), 'translation'),
            ((0, 1, (-1, -1)), ''), ((0, 2, (1, 1)), ''), ((0, 0, (-1, -1)), 'translation'),
            ((0, 0, (-1, 1)), 'translation'), ((0, 2, (1, 0)), ''), ((0, 1, (-2, 0)), '')]
    for (i, j, n), name in rows:
        ex = exact_resistance(i, j, n)
        nu = numeric_resistance(i, j, n, dps=DPS)
        assert abs(nu - mp.mpf(str(sp.N(ex, 50)))) < TOL, "exact and numerical values differ"
        print(f"R_{i+1}{j+1}{n}  {name:22s} exact = {str(ex):38s} = {sp.N(ex, 25)}   numeric = {mp.nstr(nu, 32)}")
