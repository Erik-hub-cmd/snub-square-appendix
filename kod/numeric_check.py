"""
Numerical checks for "Exact two-point resistances on the snub square lattice" (Section 6).
  (a) Eq. (1) integrated directly on a midpoint grid (batched 4x4 inversions), no symbolic input.
  (b) Kirchhoff's equations on N x N tori built from the geometric coordinates, Richardson extrapolation.
  (c) The symmetrised integral, Eq. (2), and the one-dimensional integral of Step 4.
  (d) Cairo lattice (dual, 6 vertices per cell): Eq. (1) on a grid and Kirchhoff on tori, for four entries of
      the proposed Cairo table.
Requires numpy, scipy.   Run:  python3 numeric_check.py
"""
import math, itertools
import numpy as np, scipy.sparse as sps, scipy.sparse.linalg as spl
from scipy import integrate

s3 = math.sqrt(3)
EXACT = {'tt': 2*s3/9, 'ts': 0.5 - s3/18, 'square diagonal': s3/9 + 1/math.pi, 'rhombus diagonal': 8*s3/9 - 1}
PAIRS = {'tt': (0, 2, (0, 1)), 'ts': (0, 1, (0, 0)), 'square diagonal': (0, 2, (0, 0)), 'rhombus diagonal': (0, 2, (-1, 0))}
EDGES = [(0, 1, (-1, 0)), (0, 1, (0, 0)), (0, 3, (0, 0)), (0, 3, (0, 1)), (1, 2, (0, 0)),
         (1, 2, (0, 1)), (2, 3, (0, 0)), (2, 3, (1, 0)), (0, 2, (0, 1)), (1, 3, (1, 0))]

# ---------- (a) direct grid integration of Eq. (1) ----------
def grid_resistance(i, j, n, K):
    th = (np.arange(K) + 0.5) * 2*np.pi/K - np.pi
    T1, T2 = np.meshgrid(th, th, indexing='ij')
    M = np.zeros((K, K, 4, 4), complex)
    for k in range(4): M[..., k, k] = 5
    for p, q, (n1, n2) in EDGES:
        ph = np.exp(1j*(n1*T1 + n2*T2))
        M[..., p, q] -= ph; M[..., q, p] -= np.conj(ph)
    G = np.linalg.inv(M)
    ph = np.exp(1j*(n[0]*T1 + n[1]*T2))
    f = G[..., i, i] + G[..., j, j] - G[..., i, j]/ph - G[..., j, i]*ph
    return f.real.mean()

# ---------- (d) Cairo lattice: degrees and edge list of cairo_lattice.py (vertices 0-3 triangle centres, 4-5 square centres)
CDEG = [3, 3, 3, 3, 4, 4]
CEDGES = [(0, 3, (0, -1)), (0, 4, (0, 0)), (0, 5, (0, 0)), (1, 2, (-1, 0)), (1, 4, (0, 1)), (1, 5, (0, 0)),
          (2, 4, (1, 0)), (2, 5, (0, 0)), (3, 4, (1, 1)), (3, 5, (0, 0))]
CEXACT = {'cairo 3-3 edge': (0, 3, (0, -1), 1 - 2*s3/9), 'cairo 3-4 edge': (4, 0, (0, 0), 0.5 + s3/18),
          'cairo 4-4, d^2=3': (4, 5, (-1, -1), 4*s3/9), 'cairo 3-3, d^2=4+2sqrt3': (0, 3, (-1, 0), 2*s3/3)}

def cairo_grid_resistance(i, j, n, K):
    th = (np.arange(K) + 0.5) * 2*np.pi/K - np.pi
    T1, T2 = np.meshgrid(th, th, indexing='ij')
    L = np.zeros((K, K, 6, 6), complex)
    for k in range(6): L[..., k, k] = CDEG[k]
    for p, q, (n1, n2) in CEDGES:
        ph = np.exp(1j*(n1*T1 + n2*T2))
        L[..., p, q] -= ph; L[..., q, p] -= np.conj(ph)
    G = np.linalg.inv(L)
    ph = np.exp(1j*(n[0]*T1 + n[1]*T2))
    return (G[..., i, i] + G[..., j, j] - G[..., i, j]/ph - G[..., j, i]*ph).real.mean()

# ---------- (b) finite tori from geometry ----------
def lattice_from_geometry():
    s = (1 + s3)/2
    A1, A2 = np.array([s, s]), np.array([s, -s])
    pts = []
    for mm, nn in itertools.product(range(-2, 3), repeat=2):
        ang = math.radians(45 + (30 if (mm + nn) % 2 == 0 else -30))
        pts += [np.array([mm*s + math.cos(ang + k*math.pi/2)/math.sqrt(2), nn*s + math.sin(ang + k*math.pi/2)/math.sqrt(2)]) for k in range(4)]
    A = np.array([A1, A2]).T; Ai = np.linalg.inv(A)
    reps = []
    for q in pts:
        fr = Ai @ q; fr -= np.floor(fr + 1e-9); r = A @ fr
        if not any(np.linalg.norm(r - t - c1*A1 - c2*A2) < 1e-6 for t in reps for c1 in (-1, 0, 1) for c2 in (-1, 0, 1)):
            reps.append(r)
    edges = []
    for p, q in itertools.product(range(len(reps)), repeat=2):
        for n1, n2 in itertools.product(range(-2, 3), repeat=2):
            if abs(np.linalg.norm(reps[q] + n1*A1 + n2*A2 - reps[p]) - 1) < 1e-6 and (q, p, (-n1, -n2)) not in edges:
                edges.append((p, q, (n1, n2)))
    return reps, edges

def torus_resistances(N, reps, edges, pairs):
    P = len(reps)
    idx = lambda p, i, j: p + P*((i % N) + N*(j % N))
    rows, cols, vals = [], [], []
    for i, j in itertools.product(range(N), repeat=2):
        for p, q, (n1, n2) in edges:
            u, v = idx(p, i, j), idx(q, i + n1, j + n2)
            rows += [u, v, u, v]; cols += [u, v, v, u]; vals += [1, 1, -1, -1]
    n = P*N*N
    solve = spl.factorized(sps.csr_matrix((vals, (rows, cols)), shape=(n, n))[:-1, :-1].tocsc())
    out = {}
    for name, (p, q, (n1, n2)) in pairs.items():
        u, v = idx(p, N//2, N//2), idx(q, N//2 + n1, N//2 + n2)
        rhs = np.zeros(n - 1); rhs[u] += 1; rhs[v] -= 1
        V = np.append(solve(rhs), 0.0)
        out[name] = V[u] - V[v]
    return out

if __name__ == '__main__':
    print("(a) Eq. (1) on a K x K midpoint grid (error ~ K^-2), Richardson K=400,800:")
    for name, (i, j, n) in PAIRS.items():
        r1, r2 = grid_resistance(i, j, n, 400), grid_resistance(i, j, n, 800)
        print(f"    {name:17s} {(4*r2 - r1)/3:.10f}   exact {EXACT[name]:.10f}")
    reps, edges = lattice_from_geometry()
    # the geometric construction must reproduce the edge list (up to labelling): 4 vertices, 10 edges, degree 5
    assert len(reps) == 4 and len(edges) == 10
    print("(b) tori of N x N cells, Richardson N=160,320:")
    geo_pairs = {}
    # identify the pairs geometrically: tt = unique edge between its two endpoint classes
    tt = [e for e in edges if sum(1 for f in edges if {f[0], f[1]} == {e[0], e[1]}) == 1][0]
    ts = [e for e in edges if e != tt and e[0] != e[1]][0]
    geo_pairs['tt'], geo_pairs['ts'] = tt, ts
    # square diagonal (d^2 = 2) and rhombus diagonal (d^2 = 3): every vertex pair at that squared distance, found
    # from the coordinates; Table 1 has a single class for each, so all such pairs must give the same value
    A1g, A2g = np.array([(1 + s3)/2, (1 + s3)/2]), np.array([(1 + s3)/2, -(1 + s3)/2])
    for name, d2 in (('square diagonal', 2), ('rhombus diagonal', 3)):
        found = [(p, q, (n1, n2)) for p in range(4) for q in range(4) for n1, n2 in itertools.product(range(-3, 4), repeat=2)
                 if abs(np.sum((reps[q] + n1*A1g + n2*A2g - reps[p])**2) - d2) < 1e-6]
        assert found, name
        geo_pairs[name] = found[0]
        for k, f in enumerate(found[1:6]):      # a few more pairs of the same squared distance, same value expected
            geo_pairs[f'{name} (pair {k + 2})'] = f
    r160 = torus_resistances(160, reps, edges, geo_pairs)
    r320 = torus_resistances(320, reps, edges, geo_pairs)
    for name in geo_pairs:
        ex = EXACT[name.split(' (')[0]]
        val = (4*r320[name] - r160[name])/3
        print(f"    {name:26s} {val:.11f}   exact {ex:.11f}   difference {val - ex:+.1e}")
        assert abs(val - ex) < 1e-9, name
    print("(d) Cairo lattice: grid (Richardson K=400,800) and tori (Richardson N=160,320):")
    t160 = torus_resistances(160, [None]*6, CEDGES, {k: (v[0], v[1], v[2]) for k, v in CEXACT.items()})
    t320 = torus_resistances(320, [None]*6, CEDGES, {k: (v[0], v[1], v[2]) for k, v in CEXACT.items()})
    for name, (i, j, n, ex) in CEXACT.items():
        g1, g2 = cairo_grid_resistance(i, j, n, 400), cairo_grid_resistance(i, j, n, 800)
        print(f"    {name:25s} grid {(4*g2 - g1)/3:.10f}  torus {(4*t320[name] - t160[name])/3:.10f}   exact {ex:.10f}")
    print("(c) reduced integrals:")
    f = lambda y, x: 4*(8 - 3*math.cos(x) - 3*math.cos(y) - 2*math.cos(x)*math.cos(y)) / (
        math.cos(x)**2 + math.cos(y)**2 - 14*math.cos(x)*math.cos(y) - 36*(math.cos(x) + math.cos(y)) + 84)
    I2, _ = integrate.dblquad(f, 0, math.pi, 0, math.pi, epsabs=1e-13, epsrel=1e-13)
    I1, _ = integrate.quad(lambda t: 2*(t*t - 1)/(t**4 + 14*t*t + 1), 1, np.inf, epsabs=1e-14)
    print(f"    Eq. (2)           {I2/math.pi**2:.15f}")
    print(f"    Step 4            {2/math.pi*I1:.15f}")
    print(f"    2*sqrt(3)/9       {2*s3/9:.15f}")
