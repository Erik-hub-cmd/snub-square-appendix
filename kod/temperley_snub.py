# Temperleys dubbelgraf till snub square (period 2Z^2). Vit = kantmittpunkt, svart = nod eller ytcentrum.
# Vikt: vit(e)-nod = c_e, vit(e)-yta = 1. Räknar matchningspolynomet (osignerat) och sidkoefficienterna.
import sympy as sp, math, itertools, sys
from fractions import Fraction as Fr
lam = Fr(1)
nodes = [(0,0),(1,0),(0,1),(1,1)]
edges = []  # (p, q) i lyftet, c
for (x,y) in nodes:
    edges.append(((x,y),(x+1,y),Fr(1))); edges.append(((x,y),(x,y+1),Fr(1)))
edges.append(((0,0),(1,1),lam)); edges.append(((2,1),(1,2),lam))
def face_of(p):  # p = punkt i lyftet (Fraction) -> (ytid, translation)
    lx, ly = math.floor(p[0]/2), math.floor(p[1]/2); q = (p[0]-2*lx, p[1]-2*ly)
    i, j = math.floor(q[0]), math.floor(q[1]); fx, fy = q[0]-i, q[1]-j
    if (i,j) == (0,0): fid = 'T00a' if fx > fy else 'T00b'
    elif (i,j) == (1,1): fid = 'T11a' if fx+fy < 1 else 'T11b'
    else: fid = f'S{i}{j}'
    return fid, (lx, ly)
def node_of(p):
    lx, ly = p[0]//2, p[1]//2; return ('v', p[0]-2*lx, p[1]-2*ly), (lx, ly)
eps = Fr(1,1000)
W = []  # vit -> lista (svart, translation, vikt)
for (p,q,c) in edges:
    m = (Fr(p[0]+q[0],2), Fr(p[1]+q[1],2)); d = (q[0]-p[0], q[1]-p[1]); nrm = (-d[1], d[0])
    nb = []
    for v in (p,q): nb.append(node_of(v)+(c,))
    for s in (1,-1): nb.append(face_of((m[0]+s*eps*nrm[0], m[1]+s*eps*nrm[1]))+(Fr(1),))
    W.append(nb)
blacks = sorted({b for nb in W for (b,t,c) in nb}, key=str)
print(len(W), 'vita', len(blacks), 'svarta')
# matchningspolynom via rekursion över vita
from collections import defaultdict
def match(i, used):
    if i == len(W): return {(0,0): Fr(1)}
    out = defaultdict(Fr)
    for (b,t,c) in W[i]:
        if b in used: continue
        for k,v in match(i+1, used|{b}).items(): out[(k[0]+t[0],k[1]+t[1])] += v*c
    return out
P = match(0, frozenset())
for k in sorted(P): print(k, P[k])
