"""
The Cairo pentagonal lattice as the planar dual (Laves tiling V3^2.4.3.4) of the snub square lattice.
Vertices of the Cairo lattice = incentres of the faces of the snub square lattice of Section 2
(4 triangles + 2 squares per unit cell, so 6 vertices per cell; same lattice vectors a1, a2).
Edges = the 10 dual edges per cell (one across each snub square edge).
Requires sympy.   Run:  python cairo_lattice.py   (prints the cell, the edge list and the determinant check)
"""
import itertools, math
import sympy as sp

s3 = sp.sqrt(3)
s = (1 + s3) / 2
al, be = (1 + s3) / 4, (s3 - 1) / 4
A1, A2 = sp.Matrix([s, s]), sp.Matrix([s, -s])
R = [sp.Matrix([s + be, -al]), sp.Matrix([s + al, be]), sp.Matrix([s - be, al]), sp.Matrix([s - al, -be])]
SNUB_EDGES = [(0, 1, (-1, 0)), (0, 1, (0, 0)), (0, 3, (0, 0)), (0, 3, (0, 1)), (1, 2, (0, 0)),
              (1, 2, (0, 1)), (2, 3, (0, 0)), (2, 3, (1, 0)), (0, 2, (0, 1)), (1, 3, (1, 0))]


def spos(j, n):
    return R[j] + n[0] * A1 + n[1] * A2


def _faces():
    """Faces of the snub square lattice meeting cell 0, as lists of (vertex, cell) by tracing the planar embedding."""
    nb = {i: [] for i in range(4)}
    for i, j, n in SNUB_EDGES:
        nb[i].append((j, n)); nb[j].append((i, (-n[0], -n[1])))
    f = lambda v: tuple(float(c) for c in spos(*v))
    def ang(u, v):
        p, q = f(u), f(v)
        return math.atan2(q[1] - p[1], q[0] - p[0])
    def nbrs(v):
        j, n = v
        return [(k, (n[0] + m[0], n[1] + m[1])) for k, m in nb[j]]
    faces = []
    for i in range(4):
        for w in nbrs((i, (0, 0))):
            face, u, v = [(i, (0, 0))], (i, (0, 0)), w
            while v != face[0]:
                face.append(v)
                # next vertex: the neighbour of v making the smallest clockwise turn from the reverse edge (v->u)
                back = ang(v, u)
                cand = sorted(nbrs(v), key=lambda x: (back - ang(v, x)) % (2 * math.pi) or 7)
                u, v = v, cand[0]
            faces.append(face)
    return faces


def build():
    """Returns (pos, deg, edges): 6 Cairo vertices (exact incentres), their degrees, and the 10 edges per cell
    as (i, j, (n1, n2)) joining (i, cell 0) to (j, cell n).  Vertex order: 4 triangle centres, then 2 square centres."""
    Ai = sp.Matrix.hstack(A1, A2).inv()
    reps, keyed = [], {}
    for face in _faces():
        c = sum((spos(*v) for v in face), sp.zeros(2, 1)) / len(face)
        fr = [sp.nsimplify(sp.simplify(x)) for x in Ai * c]
        shift = (sp.floor(fr[0] + sp.Rational(1, 10**6)), sp.floor(fr[1] + sp.Rational(1, 10**6)))
        base = sp.simplify(c - shift[0] * A1 - shift[1] * A2)
        key = tuple(round(float(x), 9) for x in base)
        if key not in keyed:
            keyed[key] = (len(face), base)
    reps = sorted(keyed.values(), key=lambda t: (t[0], float(t[1][0]), float(t[1][1])))
    pos = [p for _, p in reps]
    deg = [k for k, _ in reps]
    # dual edges: two Cairo vertices are adjacent iff their faces share a snub square edge,
    # i.e. iff they are at distance r_i + r_j (inradii) -- checked combinatorially below
    inr = {3: 1 / (2 * s3), 4: sp.Rational(1, 2)}
    edges = []
    pf = [tuple(float(x) for x in p) for p in pos]
    for i, j in itertools.product(range(6), repeat=2):
        for n in itertools.product(range(-2, 3), repeat=2):
            if i == j and n == (0, 0):
                continue
            q = (pf[j][0] + n[0] * float(A1[0]) + n[1] * float(A2[0]), pf[j][1] + n[0] * float(A1[1]) + n[1] * float(A2[1]))
            dd = math.dist(pf[i], q)
            if abs(dd - float(inr[deg[i]] + inr[deg[j]])) < 1e-9:
                if (j, i, (-n[0], -n[1])) not in edges:
                    edges.append((i, j, n))
    return pos, deg, edges


def cpos(pos, j, n):
    return pos[j] + n[0] * A1 + n[1] * A2


X, Z = sp.symbols('X Z')


def bloch_laplacian(deg, edges):
    L = sp.diag(*deg)
    for i, j, (n1, n2) in edges:
        ph = X**n1 * Z**n2
        L[i, j] -= ph
        L[j, i] -= 1 / ph
    return L


# ---------------------------------------------------------------- distances and symmetry classes
# Distances are measured with the degree-3--degree-4 edges (8 of the 10 edges per cell) of length 1.
# In the coordinates above those edges have length 1/(2 sqrt3) + 1/2, so we scale by LAM = 3 - sqrt3;
# the degree-3--degree-3 edges then have length sqrt3 - 1.
LAM = 3 - s3


def dist2(pos, x, y):
    """exact squared distance (scaled units) between Cairo vertices x=(i,n) and y=(j,m)"""
    d = cpos(pos, *y) - cpos(pos, *x)
    return sp.nsimplify(sp.expand(LAM**2 * d.dot(d)))


def symmetries(pos):
    """isometries p -> O p + t (O in the dihedral group of order 8) mapping the Cairo vertex set onto itself"""
    import numpy as np
    P = [np.array([float(c) for c in p]) for p in pos]
    B = np.array([[float(A1[0]), float(A2[0])], [float(A1[1]), float(A2[1])]])
    Bi = np.linalg.inv(B)
    def ident(q):
        for j, p in enumerate(P):
            fr = Bi @ (q - p)
            if np.allclose(fr, np.round(fr), atol=1e-8):
                return j, (int(round(fr[0])), int(round(fr[1])))
        return None
    ops = []
    for k in range(4):
        c, sn = round(math.cos(k * math.pi / 2)), round(math.sin(k * math.pi / 2))
        for refl in (1, -1):
            O = np.array([[c, -sn], [sn, c]]) @ np.array([[1, 0], [0, refl]])
            for j in range(6):
                t = P[j] - O @ P[0]
                if all(ident(O @ p + t) is not None for p in P):
                    ops.append((O, t))
    def apply(op, v):
        O, t = op
        j, n = v
        return ident(O @ (P[j] + n[0] * B[:, 0] + n[1] * B[:, 1]) + t)
    return ops, apply, ident


def targets(pos, src, rmax=3):
    """all Cairo vertices y != src with scaled distance < rmax, grouped in classes of unordered pairs {src, y}
    under the symmetry group.  Returns list of (representative, d^2, members)."""
    ops, apply, _ = symmetries(pos)
    srcv = (src, (0, 0))
    ys = []
    for j in range(6):
        for n in itertools.product(range(-5, 6), repeat=2):
            if (j, n) == srcv:
                continue
            if float(dist2(pos, srcv, (j, n))) < rmax**2 - 1e-12:
                ys.append((j, n))
    def shift(v, w):          # translate v by the cell of w
        return (v[0], (v[1][0] - w[1][0], v[1][1] - w[1][1]))
    def images(y):
        out = set()
        for op in ops:
            gx, gy = apply(op, srcv), apply(op, y)
            if gx[0] == src:                       # g fixes the source up to translation
                out.add(shift(gy, gx))
            if gy[0] == src:                       # g swaps the endpoints
                out.add(shift(gx, gy))
        return out
    yset, classes, seen = set(ys), [], set()
    for y in sorted(ys, key=lambda v: (float(dist2(pos, srcv, v)), v)):
        if y in seen:
            continue
        orb = sorted(images(y) & yset)
        assert y in orb
        seen |= set(orb)
        classes.append((y, dist2(pos, srcv, y), orb))
    return classes


def laurent_to_ab(e, a, b):
    """average a Laurent polynomial in X, Z over X -> 1/X, Z -> 1/Z and write it in a = cos t1, b = cos t2"""
    e = sp.expand(e)
    e = sp.expand((e + e.subs(X, 1/X) + e.subs(Z, 1/Z) + e.subs({X: 1/X, Z: 1/Z}, simultaneous=True)) / 4)
    sh = 60
    out = 0
    for (i, j), c in sp.Poly(sp.expand(e * X**sh * Z**sh), X, Z).terms():
        i -= sh; j -= sh
        if i >= 0 and j >= 0:
            out += c * (2 * sp.chebyshevt(i, a) if i else 1) * (2 * sp.chebyshevt(j, b) if j else 1)
    return sp.expand(out)


_ADJ = None


def numerator(i, j, n):
    """Laurent numerator of the integrand of Eq. (1) for the Cairo pair (i,0)-(j,n): adj entries, det = det M_snub"""
    global _ADJ
    if _ADJ is None:
        pos, deg, edges = build()
        _ADJ = bloch_laplacian(deg, edges).adjugate().applyfunc(sp.expand)
    ph = X**n[0] * Z**n[1]
    return sp.expand(_ADJ[i, i] + _ADJ[j, j] - _ADJ[i, j] / ph - _ADJ[j, i] * ph)


def diag_difference(i, j):
    """Laurent numerator of G_ii - G_jj (times det): its integral is the finite difference g_i - g_j of the
    (divergent) on-site Green's functions, which enters Kirchhoff's law on a lattice with two vertex types"""
    numerator(0, 0, (0, 0))
    return sp.expand(_ADJ[i, i] - _ADJ[j, j])


if __name__ == '__main__':
    pos, deg, edges = build()
    print('Cairo vertices (exact), degrees:')
    for p, d in zip(pos, deg):
        print('   ', [sp.nsimplify(x) for x in p], d)
    print('edges:', edges)
    assert len(edges) == 10 and sorted(deg) == [3, 3, 3, 3, 4, 4]
    L = bloch_laplacian(deg, edges)
    detL = sp.factor(sp.expand(L.det()))
    from exact_snub import DET as DET_SNUB
    print('det L_Cairo =', detL)
    print('det L_Cairo - det M_snub =', sp.simplify(sp.expand(detL - DET_SNUB)))
