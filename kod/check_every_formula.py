"""
Second, separately written line-by-line verification of the article.  Written from scratch:
does not import any earlier project code.  Every check prints PASS/FAIL.
"""
import itertools, math
import sympy as sp
import mpmath as mp

RESULTS = []
def check(name, ok, info=""):
    RESULTS.append((name, bool(ok)))
    print(("PASS " if ok else "FAIL ") + name + (("   [" + str(info) + "]") if info != "" else ""), flush=True)

s3 = sp.sqrt(3)
# ---------------- Section 2: geometry ----------------
s = (1 + s3) / 2
al, be = (1 + s3) / 4, (s3 - 1) / 4
A1, A2 = sp.Matrix([s, s]), sp.Matrix([s, -s])
r = [sp.Matrix([s + be, -al]), sp.Matrix([s + al, be]), sp.Matrix([s - be, al]), sp.Matrix([s - al, -be])]
pos = lambda j, n: r[j] + n[0] * A1 + n[1] * A2
EDGES = [(0, 1, (-1, 0)), (0, 1, (0, 0)), (0, 3, (0, 0)), (0, 3, (0, 1)), (1, 2, (0, 0)),
         (1, 2, (0, 1)), (2, 3, (0, 0)), (2, 3, (1, 0)), (0, 2, (0, 1)), (1, 3, (1, 0))]
TT = [(0, 2, (0, 1)), (1, 3, (1, 0))]
for (i, j, n) in EDGES:
    d2 = sp.simplify((pos(j, n) - pos(i, (0, 0))).dot(pos(j, n) - pos(i, (0, 0))))
    check(f"edge ({i+1},{j+1};{n}) has length exactly 1", d2 == 1, d2)
cen = sp.Matrix([s, 0])
check("r1..r4 are the corners of a unit square centred at (s,0)",
      all(sp.simplify((r[k] - cen).dot(r[k] - cen) - sp.Rational(1, 2)) == 0 for k in range(4)) and
      all(sp.simplify((r[(k+1) % 4] - r[k]).dot(r[(k+1) % 4] - r[k]) - 1) == 0 for k in range(4)))
# completeness: no other unit-distance pairs (numerically, generous range)
rf = [tuple(float(c) for c in v) for v in r]
A1f, A2f = (float(s), float(s)), (float(s), -float(s))
def pf(j, n): return (rf[j][0] + n[0]*A1f[0] + n[1]*A2f[0], rf[j][1] + n[0]*A1f[1] + n[1]*A2f[1])
nbr = {}
for i in range(4):
    nbr[i] = []
    for j in range(4):
        for n in itertools.product(range(-3, 4), repeat=2):
            p, q = pf(i, (0, 0)), pf(j, n)
            if abs(math.dist(p, q) - 1) < 1e-9:
                nbr[i].append((j, n))
listed = {i: [] for i in range(4)}
for (i, j, n) in EDGES:
    listed[i].append((j, n)); listed[j].append((i, (-n[0], -n[1])))
check("unit-distance neighbours coincide with the listed edges; every vertex has degree 5",
      all(sorted(nbr[i]) == sorted(listed[i]) and len(nbr[i]) == 5 for i in range(4)))
# triangle-triangle vs triangle-square: count common neighbours of endpoints
def neigh(v):  # v = (j, n)
    j, n = v
    return {(k, (n[0]+m[0], n[1]+m[1])) for (k, m) in listed[j]}
ok = True
for (i, j, n) in EDGES:
    common = neigh((i, (0, 0))) & neigh((j, n))
    want = 2 if (i, j, n) in TT else 1
    ok &= (len(common) == want)
check("TT edges lie in two triangles, the other eight in exactly one triangle", ok)
# cyclic order of faces around vertex 1: angles of the 5 neighbours, gaps of 60 deg -> triangle, 90 -> square
p0 = pf(0, (0, 0))
angs = sorted(math.degrees(math.atan2(pf(j, n)[1]-p0[1], pf(j, n)[0]-p0[0])) % 360 for (j, n) in listed[0])
gaps = [round((angs[(k+1) % 5] - angs[k]) % 360) for k in range(5)]
cyc = ''.join('t' if g == 60 else 's' if g == 90 else '?' for g in gaps)
check("faces around a vertex are t,t,s,t,s in cyclic order", any((cyc*2)[k:k+5] == 'ttsts' for k in range(5)), cyc)

# ---------------- Bloch matrix ----------------
X, Z = sp.symbols('X Z')
M = 5 * sp.eye(4)
for (i, j, n) in EDGES:
    ph = X**n[0] * Z**n[1]
    M[i, j] -= ph; M[j, i] -= 1/ph
Mpaper = sp.Matrix([[5, -1-1/X, -Z, -1-Z],
                    [-1-X, 5, -1-Z, -X],
                    [-1/Z, -1-1/Z, 5, -1-X],
                    [-1-1/Z, -1/X, -1-1/X, 5]])
check("displayed matrix M(theta) equals the matrix built from the edge list", sp.simplify(M - Mpaper) == sp.zeros(4))
check("M is Hermitian on the torus (M^T(1/X,1/Z) = M)", sp.simplify(M.T.subs({X: 1/X, Z: 1/Z}, simultaneous=True) - M) == sp.zeros(4))

# ---------------- Lemma 2 ----------------
a, b = sp.symbols('a b', real=True)
D = a**2 + b**2 - 14*a*b - 36*(a + b) + 84
detM = sp.expand(M.det())
A_ = (X + 1/X)/2; B_ = (Z + 1/Z)/2
check("Lemma 2: det M = 4 D(cos t1, cos t2)  (as Laurent polynomials)", sp.expand(detM - 4*D.subs({a: A_, b: B_})) == 0)
check("Lemma 2: D = 4(a-b)^2 + 3(2-a-b)(14+a+b)", sp.expand(D - (4*(a-b)**2 + 3*(2-a-b)*(14+a+b))) == 0)
C2 = lambda W: (W**2 + 1/W**2)/2
owaidat = 340 - 144*(A_ + B_) - 56*A_*B_ + 2*C2(X) + 2*C2(Z)
check("cosine form of det M (= Owaidat-Asad determinant)", sp.expand(detM - owaidat) == 0)
# D>0 except (1,1): 4(a-b)^2>=0, (2-a-b)>=0 with equality only at a=b=1, (14+a+b)>=12
# the two terms vanish simultaneously only when a = b and (a + b = 2 or a + b = -14); in [-1,1]^2 only (1,1) remains
zeros = sp.solve([a - b, (2 - a - b)*(14 + a + b)], [a, b], dict=True)
zeros_in_box = [z for z in zeros if all(-1 <= z[v] <= 1 for v in (a, b))]
check("D > 0 on [-1,1]^2 minus (1,1) (both terms >= 0; their common zeros inside the box computed: only (1,1))",
      D.subs({a: 1, b: 1}) == 0 and zeros_in_box == [{a: 1, b: 1}] and
      sp.minimum(14 + a + b, a, sp.Interval(-1, 1)).subs(b, -1) == 12, zeros_in_box)

# ---------------- Remark 3 ----------------
P = sp.expand(X**2 * Z**2 * detM)
Pp = sp.Poly(P, X, Z)
mons = Pp.monoms()
from sympy import Point, Polygon, convex_hull
hull = convex_hull(*[Point(*q) for q in mons])
interior = sum(1 for i in range(0, 5) for j in range(0, 5) if hull.encloses_point(Point(i, j)))
check("P(X,Z) is a polynomial; Newton polygon has 5 interior points", all(e >= 0 for q in mons for e in q) and interior == 5, interior)
r30 = sp.sqrt(30)
pts = [(1, 1), (11+2*r30, -1), (11-2*r30, -1), (-1, 11+2*r30), (-1, 11-2*r30)]
H = sp.hessian(P, (X, Z))
ok = True
for (x0, z0) in pts:
    sub = {X: x0, Z: z0}
    ok &= all(sp.simplify(e.subs(sub)) == 0 for e in (P, sp.diff(P, X), sp.diff(P, Z)))
    ok &= sp.simplify(H.det().subs(sub)) != 0
check("Remark 3: the five stated points are ordinary double points of P=0", ok)
G = sp.groebner([P, sp.diff(P, X), sp.diff(P, Z)], X, Z, order='lex')
zpoly = [g for g in G.exprs if g.free_symbols == {Z}][0]
# solve the eliminant, back-substitute into the other Groebner elements, keep the solutions with X Z != 0
sing = set()
for z0 in sp.roots(sp.Poly(zpoly, Z), multiple=True):
    rest = [sp.expand(g.subs(Z, z0)) for g in G.exprs]
    rest = [g for g in rest if g != 0]
    xs = [sol[X] for sol in sp.solve(rest, X, dict=True)] if rest else []
    for x0 in xs:
        sing.add((sp.nsimplify(sp.simplify(x0)), sp.simplify(z0)))
sing_in_torus = [(x0, z0) for (x0, z0) in sing if x0 != 0 and z0 != 0]
def _same(p, q):
    return all(sp.simplify(u - v) == 0 for u, v in zip(p, q))
matched = all(any(_same(p, q) for q in sing_in_torus) for p in pts)
extra = [q for q in sing_in_torus if not any(_same(q, p) for p in pts)]
check("Remark 3: no other singular points in (C*)^2 (Groebner eliminant solved, solutions compared with the five points)",
      matched and not extra and len(sing_in_torus) == 5, (sp.factor(zpoly), len(sing_in_torus), extra))
check("Remark 3: D(-1,b) = (b-11)^2, D(a,-1) = (a-11)^2, D(1,1) = 0",
      sp.expand(D.subs(a, -1) - (b-11)**2) == 0 and sp.expand(D.subs(b, -1) - (a-11)**2) == 0 and D.subs({a: 1, b: 1}) == 0)
m, mu, tau, w, rho = sp.symbols('m mu tau w rho', positive=True)
am = (5 - m**2)/(m**2 - 1)
bp = (11*m**2 + 17 + 16*s3*m)/(m**2 - 1); bm = (11*m**2 + 17 - 16*s3*m)/(m**2 - 1)
check("Remark 3: (a(m), b_+(m)) and (a(m), b_-(m)) lie on the conic; m -> -m swaps b_+ and b_-",
      sp.simplify(D.subs({a: am, b: bp})) == 0 and sp.simplify(D.subs({a: am, b: bm})) == 0 and sp.simplify(bp.subs(m, -m) - bm) == 0)
check("Remark 3: a^2-1 = -8(m^2-3)/(m^2-1)^2", sp.simplify(am**2 - 1 + 8*(m**2-3)/(m**2-1)**2) == 0)
check("Remark 3: b_+^2-1 = 8(sqrt3 m+2)^2 (m+sqrt3)(5m+3sqrt3)/(m^2-1)^2",
      sp.simplify(bp**2 - 1 - 8*(s3*m+2)**2*(m+s3)*(5*m+3*s3)/(m**2-1)**2) == 0)

# ---------------- Lemma 4 (symmetrised integrand) via Laurent averaging ----------------
adj = M.adjugate()
def avg4(expr):
    e = sp.expand(expr)
    return sp.expand((e + e.subs(X, 1/X) + e.subs(Z, 1/Z) + e.subs({X: 1/X, Z: 1/Z}, simultaneous=True))/4)
def numer(i, j, n):
    ph = X**n[0] * Z**n[1]
    return sp.expand(adj[i, i] + adj[j, j] - adj[i, j]/ph - adj[j, i]*ph)
Ntt = avg4(numer(0, 2, (0, 1)))
check("Lemma 4: averaged numerator = 16(8-3a-3b-2ab)", sp.expand(Ntt - 16*(8 - 3*A_ - 3*B_ - 2*A_*B_)) == 0)

# ---------------- eq (3) ----------------
mp.mp.dps = 30
ok = True
for beta in [mp.mpf('1.01'), mp.mpf(3), mp.mpf(49)]:
    lhs = mp.quad(lambda ph: 1/(beta - mp.cos(ph)), [0, mp.pi]) / mp.pi
    ok &= abs(lhs - 1/mp.sqrt(beta**2 - 1)) < mp.mpf(10)**-25
check("eq (3) at beta = 1.01, 3, 49", ok)

# ---------------- Step 1 ----------------
bpl = 7*a + 18 + 4*s3*rho; bmi = 7*a + 18 - 4*s3*rho
red = lambda e: sp.expand(sp.expand(e).subs(rho**2, (a+1)*(a+5)))
check("Step 1: D = (b_+ - b)(b_- - b) with rho^2=(a+1)(a+5)", red((bpl - b)*(bmi - b) - D) == 0)
cp = (2*(a+1)*(7*a+23) + 4*s3*(2*a+3)*rho)/(2*s3*rho)
cm = (2*(a+1)*(7*a+23) - 4*s3*(2*a+3)*rho)/(-2*s3*rho)
lhs = sp.together(cp*(bmi - b) + cm*(bpl - b))
check("Step 1: partial fractions c+/(b+-b) + c-/(b--b) = 4(8-3a-3b-2ab)/D",
      sp.simplify(red(sp.numer(lhs)) - red(4*(8-3*a-3*b-2*a*b)*sp.denom(lhs))) == 0)

# ---------------- Step 2 ----------------
rho_m = 4*m/(m**2-1)
check("Step 2: rho(m)^2 = (a+1)(a+5) and rho(m) > 0", sp.simplify(rho_m**2 - ((am+1)*(am+5))) == 0)
check("Step 2: a(sqrt3) = 1, a(oo) = -1, da/dm < 0 for m > sqrt3",
      sp.simplify(am.subs(m, s3)) == 1 and sp.limit(am, m, sp.oo) == -1 and sp.simplify(sp.diff(am, m) + 8*m/(m**2-1)**2) == 0)
bp_m = bpl.subs({rho: rho_m, a: am}); bm_m = bmi.subs({rho: rho_m, a: am})
ok = True
for sg, bb in ((1, bp_m), (-1, bm_m)):
    ok &= sp.simplify(bb + 1 - 4*(s3*m + sg*2)**2/(m**2-1)) == 0
    ok &= sp.simplify(bb - 1 - 2*(m + sg*s3)*(5*m + sg*3*s3)/(m**2-1)) == 0
check("Step 2: b_pm + 1 and b_pm - 1 formulas", ok)
cp_m = cp.subs({rho: rho_m, a: am}); cm_m = cm.subs({rho: rho_m, a: am})
ok = sp.simplify(cp_m - 2*(s3*m+2)*(m+s3)**2/(s3*m*(m**2-1))) == 0 and \
     sp.simplify(cm_m - 2*(s3*m-2)*(m-s3)**2/(s3*m*(m**2-1))) == 0
check("Step 2: c_pm formula", ok)
ok = True
for sg, cc, bb in ((1, cp_m, bp_m), (-1, cm_m, bm_m)):
    target = (m + sg*s3)**3/(6*m**2*(5*m + sg*3*s3))     # square of claimed term
    ok &= sp.simplify(cc**2/(bb**2 - 1) - target) == 0
check("Step 2: c_pm/sqrt(b_pm^2-1) = (m pm sqrt3)^{3/2}/(sqrt6 m sqrt(5m pm 3sqrt3))  (squares agree; both sides > 0)", ok)
ok = True
for mv in [sp.Rational(18, 10), 3, 7]:
    for sg in (1, -1):
        cc = (cp_m if sg == 1 else cm_m).subs(m, mv); bb = (bp_m if sg == 1 else bm_m).subs(m, mv)
        ok &= sp.N(cc) > 0 and sp.N(bb) > 1
check("Step 2: positivity (c_pm > 0, b_pm > 1) at sample m > sqrt3", ok)

# ---------------- Step 3 ----------------
check("Step 3: 1-a and 1+a formulas", sp.simplify(1 - am - 2*(m**2-3)/(m**2-1)) == 0 and sp.simplify(1 + am - 4/(m**2-1)) == 0)
check("Step 3: (da/dm)^2/(1-a^2) = (2 sqrt2 m/((m^2-1) sqrt(m^2-3)))^2",
      sp.simplify(sp.diff(am, m)**2/(1 - am**2) - 8*m**2/((m**2-1)**2*(m**2-3))) == 0)
ok = True
for sg in (1, -1):
    lhs2 = (1/sp.pi)**2 * 8*m**2/((m**2-1)**2*(m**2-3)) * (m + sg*s3)**3/(6*m**2*(5*m + sg*3*s3))
    rhs2 = (2/(s3*sp.pi))**2 * (m + sg*s3)**2/((m - sg*s3)*(5*m + sg*3*s3)) / (m**2-1)**2
    ok &= sp.simplify(lhs2 - rhs2) == 0
check("Step 3: first display line (squared, termwise)", ok)
ok = True
for sg in (1, -1):
    lhs_m = (2/(s3*sp.pi))**2 * (m + sg*s3)**2/((m - sg*s3)*(5*m + sg*3*s3)) / (m**2-1)**2 * 3   # (dm = sqrt3 dmu)^2
    rhs_mu = (2/sp.pi)**2 * (mu + sg)**2/((mu - sg)*(5*mu + sg*3)) / (3*mu**2-1)**2
    ok &= sp.simplify(lhs_m.subs(m, s3*mu) - rhs_mu) == 0
check("Step 3: mu-substitution (squared, termwise, including Jacobian)", ok)

# ---------------- Step 4 ----------------
ok = True
for sg in (1, -1):
    mut = (tau**2 + 3)/(tau**2 - 5) if sg == 1 else (tau**2 + 3)/(5 - tau**2)
    rad = sp.sqrt((mut - sg)*(5*mut + sg*3))
    rad_claim = (mut - 1)*tau if sg == 1 else (mut + 1)*tau
    ok &= sp.simplify(rad_claim**2 - (mut - sg)*(5*mut + sg*3)) == 0
    ok &= sp.simplify(3*mut**2 - 1 - 2*(tau**4 + 14*tau**2 + 1)/(tau**2 - 5)**2) == 0
    integrand = (mut + sg)/(rad_claim*(3*mut**2 - 1)) * sp.diff(mut, tau)
    target = -2*(tau**2-1)/(tau**4+14*tau**2+1) if sg == 1 else 2*(tau**2-1)/(tau**4+14*tau**2+1)
    ok &= sp.simplify(integrand - target) == 0
check("Step 4: Euler substitutions, 3mu^2-1 identity, transformed integrands", ok)
mut1 = (tau**2 + 3)/(tau**2 - 5); mut2 = (tau**2 + 3)/(5 - tau**2)
check("Step 4: limits (tau=oo->mu=1, tau->sqrt5+ -> mu->oo; tau=1->mu=1, tau->sqrt5- -> mu->oo)",
      sp.limit(mut1, tau, sp.oo) == 1 and sp.limit(mut1, tau, sp.sqrt(5), '+') == sp.oo and
      mut2.subs(tau, 1) == 1 and sp.limit(mut2, tau, sp.sqrt(5), '-') == sp.oo)
check("Step 4: radical positivity: (mu-1)tau > 0 on (sqrt5,oo), (mu+1)tau > 0 on (1,sqrt5)",
      all(sp.N(((mut1 - 1)*tau).subs(tau, t)) > 0 for t in [2.3, 5, 50]) and all(sp.N(((mut2 + 1)*tau).subs(tau, t)) > 0 for t in [1.1, 1.8, 2.2]))

# ---------------- Step 5 ----------------
wt = tau + 1/tau
check("Step 5: (tau^2-1) dtau/(tau^4+14tau^2+1) = dw/(w^2+12)",
      sp.simplify((tau**2-1)/(tau**4+14*tau**2+1) - sp.diff(wt, tau)/(wt**2+12)) == 0)
val = 4/sp.pi * sp.integrate(1/(w**2+12), (w, 2, sp.oo))
check("Step 5: (4/pi) int_2^oo dw/(w^2+12) = 2 sqrt3/9", sp.simplify(val - 2*s3/9) == 0, sp.simplify(val))
check("Step 5: intermediate expression", sp.simplify(4/sp.pi/(2*s3)*(sp.pi/2 - sp.atan(1/s3)) - 2*s3/9) == 0)
I45 = sp.integrate(2*(tau**2-1)/(tau**4+14*tau**2+1), (tau, 1, sp.oo))
I45 = I45.subs(sp.atan(7*s3/3), sp.atan(s3/6) + sp.pi/3)      # arctan identity, checked below
check("Step 4->5: (2/pi) int_1^oo 2(t^2-1)/(t^4+14t^2+1) = 2sqrt3/9", sp.simplify(2/sp.pi*I45 - 2*s3/9) == 0)

# ---------------- printed digits ----------------
mp.mp.dps = 40
check("Theorem 1 digits 0.38490017945975050967", mp.nstr(2*mp.sqrt(3)/9, 21).startswith('0.38490017945975050967'), mp.nstr(2*mp.sqrt(3)/9, 25))
check("Corollary 5 digits 0.40377495513506237258", mp.nstr(mp.mpf(1)/2 - mp.sqrt(3)/18, 21).startswith('0.40377495513506237258'), mp.nstr(mp.mpf(1)/2-mp.sqrt(3)/18, 25))
check("Corollary 5: Rtt + 4 Rts = 2", sp.simplify(2*s3/9 + 4*(sp.Rational(1, 2) - s3/18) - 2) == 0)
check("Corollary 6 numbers 0.6151, 0.5962", round(1 - 2*math.sqrt(3)/9, 4) == 0.6151 and round(0.5 + math.sqrt(3)/18, 4) == 0.5962)
check("Intro 4.8.8 numbers 0.6385, 0.7229",
      round(0.5 + math.sqrt(2)/(4*math.pi)*math.atan(2*math.sqrt(2)), 4) == 0.6385 and
      round(1 - math.sqrt(2)/2 + math.sqrt(2)/math.pi*math.atan(math.sqrt(2)), 4) == 0.7229)
check("Discussion: arctan(7sqrt3/3) - arctan(sqrt3/6) = pi/3", abs(mp.atan(7*mp.sqrt(3)/3) - mp.atan(mp.sqrt(3)/6) - mp.pi/3) < mp.mpf(10)**-35)

# ---------------- Cairo lattice (proposed new section) ----------------
# Built here from scratch as the planar dual: vertices = incentres of the faces of the snub square lattice,
# two of them adjacent iff their distance equals the sum of the inradii (faces sharing an edge).
tri_in, sq_in = 1/(2*s3), sp.Rational(1, 2)
faces = []          # (centre, degree)
for i, j, k in itertools.combinations([(v, n) for v in range(4) for n in itertools.product(range(-1, 2), repeat=2)], 3):
    pts = [pos(*u) for u in (i, j, k)]
    if all(sp.simplify((p - q).dot(p - q) - 1) == 0 for p, q in itertools.combinations(pts, 2)):
        faces.append(((pts[0] + pts[1] + pts[2]) / 3, 3))
faces.append((sp.Matrix([s, 0]), 4)); faces.append((sp.Matrix([0, 0]), 4))
Ainv = sp.Matrix.hstack(A1, A2).inv()
cv = []
for c, dg in faces:     # reduce to the unit cell, remove duplicates
    fr = [sp.nsimplify(x) for x in Ainv * c]
    c0 = sp.simplify(c - sp.floor(fr[0]) * A1 - sp.floor(fr[1]) * A2)
    if not any(sp.simplify((c0 - e).dot(c0 - e)) == 0 for e, _ in cv):
        cv.append((c0, dg))
check("Cairo: 6 vertices per cell (4 triangle centres of degree 3, 2 square centres of degree 4)",
      len(cv) == 6 and sorted(d for _, d in cv) == [3, 3, 3, 3, 4, 4], len(cv))
cvf = [(tuple(float(x) for x in c), d) for c, d in cv]
inr = lambda d: float(tri_in) if d == 3 else 0.5
cedges = []
for (i, (pi_, di)), (j, (pj, dj)) in itertools.product(enumerate(cvf), repeat=2):
    for n in itertools.product(range(-2, 3), repeat=2):
        q = (pj[0] + n[0]*A1f[0] + n[1]*A2f[0], pj[1] + n[0]*A1f[1] + n[1]*A2f[1])
        if (i, n) != (j, (0, 0)) and abs(math.dist(pi_, q) - inr(di) - inr(dj)) < 1e-9 and (j, i, (-n[0], -n[1])) not in cedges:
            cedges.append((i, j, n))
cdeg = [d for _, d in cv]
check("Cairo: 10 edges per cell, vertex degrees match the face sizes",
      len(cedges) == 10 and all(sum((e[0] == v) + (e[1] == v) for e in cedges) == cdeg[v] for v in range(6)))
LAMc = 3 - s3
lens = sorted({sp.nsimplify(sp.simplify(LAMc**2 * ((cv[j][0] + n[0]*A1 + n[1]*A2 - cv[i][0]).dot(cv[j][0] + n[0]*A1 + n[1]*A2 - cv[i][0]))))
               for i, j, n in cedges}, key=float)
check("Cairo: with scale 3-sqrt3, deg3-deg4 edges have length 1 and deg3-deg3 edges length sqrt3-1",
      len(lens) == 2 and sp.simplify(lens[1] - 1) == 0 and sp.simplify(lens[0] - (s3 - 1)**2) == 0, lens)
Lc = sp.diag(*cdeg)
for i, j, n in cedges:
    ph = X**n[0] * Z**n[1]
    Lc[i, j] -= ph; Lc[j, i] -= 1/ph
check("Cairo: det of the Bloch Laplacian equals det M of the snub square lattice (as Laurent polynomials)",
      sp.expand(Lc.det() - detM) == 0)
# Table of the proposed section: value by (degree of source, degree of target, d^2); the d^2 are distinct per source type
r3 = sp.sqrt(3)
CT = {(4, 3, 1): sp.Rational(1, 2) + r3/18, (4, 3, 4 - r3): sp.Rational(3, 2) - 7*r3/18, (4, 4, 3): 4*r3/9,
      (4, 3, 4 + r3): -sp.Rational(1, 2) + r3/2 + 2/sp.pi, (4, 4, 6): 4 - 16*r3/9,
      (4, 3, 7): -sp.Rational(3, 2) + 11*r3/6 - 2/sp.pi, (4, 3, 10 - r3): -sp.Rational(9, 2) + 25*r3/18 + 10/sp.pi,
      (3, 3, 4 - 2*r3): 1 - 2*r3/9, (3, 3, 2): r3/3 + 1/sp.pi, (3, 3, 4): 2 - 2*r3/9 - 2/sp.pi,
      (3, 3, 8 - 2*r3): -2 + 11*r3/9 + 3/sp.pi, (3, 3, 6): 7 - 8*r3/3 - 4/sp.pi, (3, 3, 4 + 2*r3): 2*r3/3}
for (p, q, dd) in list(CT):
    if p == 4 and q == 3:
        CT[(3, 4, dd)] = CT[(p, q, dd)]
dgc = r3/18                         # g_3 - g_4
check("Cairo: Foster 2 R_33 + 8 R_34 = 6 and duality R_33 = 1 - Rtt, R_34 = 1 - Rts",
      sp.simplify(2*CT[(3, 3, 4 - 2*r3)] + 8*CT[(3, 4, 1)] - 6) == 0 and
      sp.simplify(CT[(3, 3, 4 - 2*r3)] - (1 - 2*s3/9)) == 0 and sp.simplify(CT[(3, 4, 1)] - (sp.Rational(1, 2) + s3/18)) == 0)
def cnb(v):
    j, n = v
    out = []
    for p, q, m in cedges:
        if p == j: out.append((q, (n[0] + m[0], n[1] + m[1])))
        if q == j: out.append((p, (n[0] - m[0], n[1] - m[1])))
    return out
def cd2(x, y):
    d = cv[y[0]][0] + y[1][0]*A1 + y[1][1]*A2 - cv[x[0]][0] - x[1][0]*A1 - x[1][1]*A2
    return sp.nsimplify(sp.expand(LAMc**2 * d.dot(d)))
for srcv in (next(v for v in range(6) if cdeg[v] == 4), next(v for v in range(6) if cdeg[v] == 3)):
    x0 = (srcv, (0, 0))
    Rv = {x0: sp.Integer(0)}
    for j in range(6):
        for n in itertools.product(range(-5, 6), repeat=2):
            y = (j, n)
            if y != x0:
                dd = cd2(x0, y)
                if float(dd) < 9 - 1e-12:
                    Rv[y] = CT[(cdeg[srcv], cdeg[j], dd)]
    cnt = bad = 0
    for x in Rv:
        nb = cnb(x)
        if all(w in Rv for w in nb):
            cnt += 1
            kappa = -2 if cdeg[x[0]] == 3 else 4
            bad += sp.simplify(sum(Rv[w] - Rv[x] for w in nb) - (2 if x == x0 else 0) - kappa*dgc) != 0
    check(f"Cairo: Kirchhoff's law (with the term kappa_x (g3-g4), g3-g4 = sqrt3/18) holds exactly from a degree-{cdeg[srcv]} vertex",
          bad == 0 and cnt > 10, f"{cnt} vertices")

print("\nSUMMARY:", sum(ok for _, ok in RESULTS), "passed,", sum(not ok for _, ok in RESULTS), "failed")
