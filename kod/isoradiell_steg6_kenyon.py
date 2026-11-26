"""
Isoradiell förklaring, steg 6–7: resistanserna på snub square (λ = 1) ur Kenyons formel för det isoradiella Z².

Steg 6. Kenyons Greensfunktion på det isoradiella Z² (tågspårsvinklar α = {0, −Φ}, β = {π/2 − Φ, π/2}, Φ = π/3):
  G(v₁, v₂) = −(1/(8π² i)) ∮ log z · exp_{v₁v₂}(z) dz/z = −(1/(4π)) Σ residyer,
räknad exakt i Q(ζ), ζ = e^{iπ/6}. Resultatet har formen a + b√3 + (c + d√3)/π. Med M_snub⁻¹ = A·M_Z²⁻¹·B⁻¹ (steg 5)
blir varje resistans på snub square en ändlig kombination av dessa värden; den jämförs exakt med runs/tabell.
Steg 7 (--sigma): samma sak med σ(A), σ(B) (√3 → −√3) och Z² med vinkeln −Φ ger samma resistanser.
Kör: python kod/isoradiell_steg6_kenyon.py [antal närmaste, standard 100] [--sigma]
"""
import sys, os, json
sys.stdout.reconfigure(encoding="utf-8")
import sympy as sp
from fractions import Fraction as F
import math
from functools import lru_cache
# --- talkropp Q(ζ), ζ^4 = ζ^2 - 1 ---
class K:
    __slots__=('c',)
    def __init__(s,c): s.c=tuple(F(x) for x in c)
    @staticmethod
    def z(k):
        k%=12; r=K((1,0,0,0))
        for _ in range(k): r=r*K((0,1,0,0))
        return r
    def __add__(s,o): return K([a+b for a,b in zip(s.c,o.c)])
    def __sub__(s,o): return K([a-b for a,b in zip(s.c,o.c)])
    def __neg__(s): return K([-a for a in s.c])
    def __mul__(s,o):
        if not isinstance(o,K): return K([a*o for a in s.c])
        p=[F(0)]*7
        for i,a in enumerate(s.c):
            if a:
                for j,b in enumerate(o.c):
                    if b: p[i+j]+=a*b
        for d in (6,5,4):          # x^d = x^{d-2} - x^{d-4}
            if p[d]: p[d-2]+=p[d]; p[d-4]-=p[d]; p[d]=F(0)
        return K(p[:4])
    __rmul__=__mul__
    def inv(s):
        # solve s*x = 1 (4x4 linear system over Q)
        M=[]
        for j in range(4):
            col=(s*K([1 if t==j else 0 for t in range(4)])).c; M.append(col)
        A=[[M[j][i] for j in range(4)]+[F(1 if i==0 else 0)] for i in range(4)]
        for c in range(4):
            piv=next(r for r in range(c,4) if A[r][c]!=0); A[c],A[piv]=A[piv],A[c]
            pv=A[c][c]; A[c]=[x/pv for x in A[c]]
            for r in range(4):
                if r!=c and A[r][c]!=0:
                    f=A[r][c]; A[r]=[x-f*y for x,y in zip(A[r],A[c])]
        return K([A[i][4] for i in range(4)])
    def re(s):   # (a, b) med Re = a + b√3
        a0,a1,a2,a3=s.c; return (a0+a2/2, a1/2)
    def im(s):
        a0,a1,a2,a3=s.c; return (a1/2+a3, a2/2)
ONE=K((1,0,0,0)); ZERO=K((0,0,0,0)); I=K.z(3)
# --- potensserier i t, koefficienter i K, trunkerade till ordning n ---
def smul(a,b,n): 
    out=[ZERO]*n
    for i in range(min(n,len(a))):
        for j in range(min(n-i,len(b))): out[i+j]=out[i+j]+a[i]*b[j]
    return out
def sinv(a,n):     # 1/a, a[0] != 0
    i0=a[0].inv(); out=[i0]
    for k in range(1,n):
        s=ZERO
        for j in range(1,min(k,len(a)-1)+1): s=s+a[j]*out[k-j]
        out.append(-(s*i0))
    return out
# --- isoradiellt Z²: tågspårsvinklar i enheter av π/6 (Φ = π/3 = 2 enheter) ---
ALPHA={0:0, 1:-2}          # NE-spår u = x - y, klass (u-1/2) mod 2: α = 0, -Φ
BETA={0:1, 1:3}            # NW-spår w = x + y: β = π/2 - Φ, π/2
def steps(p1,p2):
    S=[]
    for (a,b,tab) in ((p1[0]-p1[1],p2[0]-p2[1],ALPHA),(p1[0]+p1[1],p2[0]+p2[1],BETA)):
        if b>a: S+=[(tab[k%2],1) for k in range(a,b)]
        else:   S+=[(tab[k%2],-1) for k in range(b,a)]
    return S          # (vinkel k i enheter π/6, tecken); steg = tecken·ζ^k
@lru_cache(maxsize=None)
def G(p1,p2):
    """Kenyon: G = -(1/(8π² i)) ∮ log z · exp(z) dz/z = -(1/(4π)) Σ Res.  Returnerar (a, b, c, d): G = a + b√3 + (c + d√3)/π"""
    if p1==p2: return (F(0),)*4
    S=steps(p1,p2)
    # poler: steg s·ζ^k ger faktor ((z+ζ^k)/(z-ζ^k))^s, pol i ζ^k om s=+1, i -ζ^k=ζ^{k+6} om s=-1
    pol={}
    fac=[]                       # (e, s): faktor ((z+e)/(z-e))^s
    for k,s in S:
        fac.append((k%12,s))
    # förenkla: (z+e)/(z-e) med e=ζ^k; inversen har pol i ζ^{k+6}
    poles={}
    for k,s in fac:
        pk = k if s==1 else (k+6)%12
        poles[pk]=poles.get(pk,0)+1
    disp_re=sum(s*math.cos(k*math.pi/6) for k,s in S); disp_im=sum(s*math.sin(k*math.pi/6) for k,s in S)
    th0=math.atan2(disp_im,disp_re)
    resP=ZERO; res1=ZERO          # Σ Res = π·resP + res1
    for pk,m in poles.items():
        e=K.z(pk)
        # γ i (θ0-π, θ0+π]
        g=pk*math.pi/6
        while g>th0+math.pi+1e-9: g-=2*math.pi
        while g<=th0-math.pi+1e-9: g+=2*math.pi
        gq=F(round(g/(math.pi/6)),6)          # γ/π
        n=m
        # h(t) = t^m R(e+t)/(e+t)
        h=[ONE]+[ZERO]*(n-1)
        for k,s in fac:
            ek=K.z(k)
            num=[e+ek,ONE]; den=[e-ek,ONE]     # z+ek, z-ek at z=e+t
            if s==-1: num,den=den,num
            # remove the pole factor t
            if den[0].c==ZERO.c: den=den[1:]+[ZERO]   # den = t·1  -> 1 after taking t^m out
            if num[0].c==ZERO.c: num=num[1:]+[ZERO]; raise ValueError('nollställe och pol möts')
            h=smul(h,smul(num,sinv(den,n),n),n)
        h=smul(h,sinv([e,ONE],n),n)
        # log z = iγ + Σ_{k≥1} (-1)^{k+1} t^k/(k e^k)
        ei=e.inv(); lg1=[ZERO]+[(ei_pow:=None) or ZERO for _ in range(n-1)]
        pw=ONE
        for k in range(1,n):
            pw=pw*ei; lg1[k]=pw*F((-1)**(k+1),k)
        # Res = Σ_{k=0}^{m-1} h[m-1-k]·log_k
        resP=resP + (I*gq)*h[m-1]
        for k in range(1,m): res1=res1+h[m-1-k]*lg1[k]
    # G = -(1/(4π))(π resP + res1) = -resP/4 - res1/(4π)
    a,b=resP.re(); c,d=res1.re()
    return (-a/4,-b/4,-c/4,-d/4)

# ---------------- steg 6: kombination via A, B ----------------
X,Z=sp.symbols('X Z'); s3=sp.sqrt(3); q=2+s3
HERE=os.path.dirname(os.path.abspath(__file__))
src=open(os.path.join(HERE,'isoradiell_steg5_laplacian.py'),encoding='utf-8').read()
ns={'sp':sp,'sqrt':sp.sqrt,'X':X,'Z':Z}
exec(src[src.index('A = sp.Matrix'):src.index("\ncheck('")],ns)
A,B=ns['A'],ns['B']
Binv=(B.adjugate()/(576*q**2)).applyfunc(lambda e: sp.expand(sp.radsimp(sp.cancel(e))))
def laurent(e):
    e=sp.expand(e); out={}
    for term in sp.Add.make_args(e):
        c,pw=term.as_coeff_mul()
        px=pz=0; coef=c
        for f in pw:
            b,ex=f.as_base_exp()
            if b==X: px+=int(ex)
            elif b==Z: pz+=int(ex)
            else: coef*=f
        out[(px,pz)]=out.get((px,pz),0)+coef
    return {k:sp.radsimp(v) for k,v in out.items() if sp.simplify(v)!=0}
AL=[[laurent(A[i,j]) for j in range(4)] for i in range(4)]
BL=[[laurent(Binv[i,j]) for j in range(4)] for i in range(4)]
PT={0:(0,0),1:(0,1),2:(1,0),3:(1,1)}   # min cellnod v = 2i+j <-> Z²-punkt (i,j)
def Gsym(i,m,j):   # G_Z2((i, cell m),(j, cell 0)) exakt som sympy-uttryck
    p1=(PT[i][0]+2*m[0],PT[i][1]+2*m[1]); p2=PT[j]
    a,b,c,d=G(p1,p2) if p1<=p2 else G(p2,p1)
    return sp.Rational(a.numerator,a.denominator)+sp.Rational(b.numerator,b.denominator)*s3+(sp.Rational(c.numerator,c.denominator)+sp.Rational(d.numerator,d.denominator)*s3)/sp.pi
def R_snub(a,b,n):
    # w = e_a - X^{-n1} Z^{-n2} e_b ;  p = A^H w ,  r = B^{-1} w
    w=[{} for _ in range(4)]
    w[a][(0,0)]=w[a].get((0,0),0)+1; w[b][(-n[0],-n[1])]=w[b].get((-n[0],-n[1]),0)-1
    p=[{} for _ in range(4)]; r=[{} for _ in range(4)]
    for i in range(4):
        for k in range(4):
            for (ex,ez),c in AL[k][i].items():          # (A^H)_{ik} = conj(A_{ki}) : X^e -> X^{-e}
                for (wx,wz),cw in w[k].items():
                    key=(-ex+wx,-ez+wz); p[i][key]=p[i].get(key,0)+c*cw
            for (ex,ez),c in BL[i][k].items():
                for (wx,wz),cw in w[k].items():
                    key=(ex+wx,ez+wz); r[i][key]=r[i].get(key,0)+c*cw
    tot=0
    for i in range(4):
        for kk,cp in p[i].items():
            for j in range(4):
                for ll,cr in r[j].items():
                    tot+=cp*cr*Gsym(i,(ll[0]-kk[0],ll[1]-kk[1]),j)    # ∫ conj(p_ik) e^{-ik} M^{-1}_ij r_jl e^{il}
    return sp.collect(sp.expand(sp.radsimp(sp.expand(tot))),[sp.pi])   # exakt; nsimplify kunde gissa fel
SIG=(0,2,3,1); OFF=[(0,0),(0,0),(0,-1),(0,-1)]
if __name__=='__main__':
    args=[a for a in sys.argv[1:] if not a.startswith('--')]
    N=int(args[0]) if args else 100
    if '--sigma' in sys.argv:
        AL=[[{k:sp.expand(v.subs(s3,-s3)) for k,v in d.items()} for d in row] for row in AL]
        BL=[[{k:sp.expand(v.subs(s3,-s3)) for k,v in d.items()} for d in row] for row in BL]
        ALPHA.update({0:0,1:2}); BETA.update({0:5,1:3}); G.cache_clear()
        print('Steg 7: σ(A), σ(B) och Z² med vinkeln −Φ')
    d=json.load(open(os.path.join(HERE,'..','runs','tabell','checkpoint.json'), encoding='utf-8'))['results']
    E=sorted((v['dist'],int(k.split(':')[0]),(int(k.split(':')[1]),int(k.split(':')[2])),v['expr'])
             for t in d.values() for k,v in t.items())
    ok=0; utan=0
    for dist,j,n,expr in E[:N]:
        a,b=SIG[0],SIG[j]; nn=(n[0]+OFF[j][0]-OFF[0][0],n[1]+OFF[j][1]-OFF[0][1])
        R=R_snub(a,b,nn); T=sp.sympify(expr)
        same=sp.simplify(R-T)==0; ok+=same
        cpi=sp.expand(R*sp.pi).coeff(sp.pi,0)            # (c + d√3): koefficienten för 1/π
        utan+= sp.expand(cpi).coeff(s3)==0
        print(f'{"OK " if same else "FEL"} R(0 -> {j}, {n}) avstånd {dist:.3f}: Kenyon {R}   tabell {expr}',flush=True)
    print(f'\n{ok} av {N} exakt lika med tabellen; {utan} av {N} utan √3/π-term')
