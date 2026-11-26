# Separat skriven omräkning: det Δ (Bareiss över Z[z,w]), Newtonpolygon, sidpolynom, Galois-test.
import math, itertools, sys
from fractions import Fraction
import flint, sympy, mpmath
ctx = flint.fmpz_mpoly_ctx.get(('z','w'), 'lex')
Z, W = ctx.gens()
mpmath.mp.dps = 60

def bloch(V, edges):
    # edges: (i,j,dx,dy,c) betyder kant i (cell 0) -- j (cell (dx,dy))
    rows = [[{} for _ in range(V)] for _ in range(V)]
    def add(i,j,dx,dy,c):
        rows[i][j][(dx,dy)] = rows[i][j].get((dx,dy),0)+c
    for (i,j,dx,dy,c) in edges:
        add(i,i,0,0,c); add(j,j,0,0,c); add(i,j,dx,dy,-c); add(j,i,-dx,-dy,-c)
    M=[]
    for i in range(V):
        mx = max([-d[0] for e in rows[i] for d in e]+[0]); my = max([-d[1] for e in rows[i] for d in e]+[0])
        R=[]
        for j in range(V):
            p = ctx.from_dict({}) if False else 0*Z
            for (dx,dy),c in rows[i][j].items():
                p += c * Z**(dx+mx) * W**(dy+my)
            R.append(p)
        M.append(R)
    return M

def bareiss(M):
    n=len(M); M=[r[:] for r in M]; sign=1; prev=1+0*Z
    for k in range(n-1):
        if M[k][k]==0:
            for r in range(k+1,n):
                if M[r][k]!=0:
                    M[k],M[r]=M[r],M[k]; sign=-sign; break
            else: return 0*Z
        for i in range(k+1,n):
            for j in range(k+1,n):
                num = M[i][j]*M[k][k]-M[i][k]*M[k][j]
                q,r = divmod(num, prev)
                assert r==0
                M[i][j]=q
        prev=M[k][k]
    return sign*M[n-1][n-1]

def hull(pts):
    pts=sorted(set(pts))
    def cr(o,a,b): return (a[0]-o[0])*(b[1]-o[1])-(a[1]-o[1])*(b[0]-o[0])
    lo=[];up=[]
    for p in pts:
        while len(lo)>=2 and cr(lo[-2],lo[-1],p)<=0: lo.pop()
        lo.append(p)
    for p in reversed(pts):
        while len(up)>=2 and cr(up[-2],up[-1],p)<=0: up.pop()
        up.append(p)
    return lo[:-1]+up[:-1]

def analyse(V, edges):
    P = bareiss(bloch(V,edges))
    d = {tuple(int(x) for x in k):int(v) for k,v in zip(P.monoms(), P.coeffs())}
    H = hull(list(d))
    area2 = abs(sum(H[i][0]*H[(i+1)%len(H)][1]-H[(i+1)%len(H)][0]*H[i][1] for i in range(len(H))))
    u = sympy.Symbol('u'); sides=[]
    for i in range(len(H)):
        a,b = H[i],H[(i+1)%len(H)]
        g = math.gcd(b[0]-a[0], b[1]-a[1]); s=((b[0]-a[0])//g,(b[1]-a[1])//g)
        Q = sum(d.get((a[0]+k*s[0],a[1]+k*s[1]),0)*u**k for k in range(g+1))
        sides.append((s, sympy.expand(Q)))
    return d, H, Fraction(area2,2), sides

def galois(Q):
    u=sympy.Symbol('u'); res=[]
    for f,m in sympy.factor_list(Q)[1]:
        deg = sympy.degree(f,u)
        if deg>=2:
            rts = mpmath.polyroots([int(c) for c in sympy.Poly(f,u).all_coeffs()], maxsteps=200, extraprec=200)
            mods = sorted(abs(r) for r in rts)
            distinct = (mods[-1]-mods[0]) > mpmath.mpf(10)**-30
            allreal = all(abs(mpmath.im(r))<1e-30 for r in rts)
            res.append((str(f), deg, distinct, allreal))
    return res
