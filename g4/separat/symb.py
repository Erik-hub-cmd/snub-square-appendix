import sys, math; sys.path.insert(0,'.')
import flint, sympy
import rev
ctx=flint.fmpz_mpoly_ctx.get(('z','w','a','b'),'lex'); Z,W,A,B=ctx.gens()
rev.Z,rev.W,rev.ctx=Z,W,ctx
src=open('arch.py').read().split("for name,(V,E,deg)")[0].replace("from rev import *","")
g={}; exec(src,g)
def run(V,E,wfun,label):
    EE=[(i,j,n1,n2,wfun(i,j,n1,n2)) for (i,j,n1,n2,_) in E]
    P=rev.bareiss(rev.bloch(V,EE))
    a,b,u=sympy.symbols('a b u'); d={}
    for mm,c in zip(P.monoms(),P.coeffs()):
        k=(int(mm[0]),int(mm[1])); d[k]=d.get(k,0)+int(c)*a**int(mm[2])*b**int(mm[3])
    H=rev.hull(list(d))
    for t in range(len(H)):
        s0,s1=H[t],H[(t+1)%len(H)]
        gg=math.gcd(s1[0]-s0[0],s1[1]-s0[1]); s=((s1[0]-s0[0])//gg,(s1[1]-s0[1])//gg)
        Q=sympy.expand(sum(d.get((s0[0]+k*s[0],s0[1]+k*s[1]),0)*u**k for k in range(gg+1)))
        c=sympy.Poly(Q,u).all_coeffs()
        print(label,'sida',s,sympy.factor(Q),' | disk:',sympy.factor(c[1]**2-4*c[0]*c[2]))
for name in ['4.8^2','4.6.12','3^2.4.3.4']:
    V,E,_=g['geo'](*{'4.8^2':g['T']['4.8^2'],}.get(name,(None,)) ) if False else g['T'][name]
    # återskapa positioner
    s2,s3=math.sqrt(2),math.sqrt(3)
    if name=='4.8^2':
        a_=1+s2; g['geo']((a_,0),(0,a_),[(s2/2,0),(-s2/2,0),(0,s2/2),(0,-s2/2)])
    elif name=='4.6.12':
        R=1/(2*math.sin(math.pi/12)); c=3+s3
        g['geo']((c,0),(c/2,c*s3/2),[(R*math.cos(math.radians(15+30*k)),R*math.sin(math.radians(15+30*k))) for k in range(12)])
    else:
        A_=(1+s3)/s2; pts=[]
        for (cx,cy,sg) in [(0,0,1),(A_/2,A_/2,-1)]:
            for k in range(4):
                t=math.radians(45+sg*15+90*k); pts.append((cx+math.cos(t)/s2, cy+math.sin(t)/s2))
        g['geo']((A_,0),(0,A_),pts)
    P,L1,L2=g['LASTPOS']
    def mid(i,j,n1,n2):
        q=(P[j][0]+n1*L1[0]+n2*L2[0], P[j][1]+n1*L1[1]+n2*L2[1]); return ((P[i][0]+q[0])/2,(P[i][1]+q[1])/2)
    def klass(i,j,n1,n2):
        mx,my=mid(i,j,n1,n2)
        best=None
        for k1 in range(-3,4):
            for k2 in range(-3,4):
                cx,cy=k1*L1[0]+k2*L2[0],k1*L1[1]+k2*L2[1]; dd=math.hypot(mx-cx,my-cy)
                if best is None or dd<best[0]: best=(dd,math.degrees(math.atan2(my-cy,mx-cx)))
        return best
    if name=='4.8^2':
        # kvadratkanter: mittpunkt nära kvadratcentrum (avstånd 1/2) -> 1, oktagon-oktagon -> a
        run(V,E,lambda *e: 1+0*Z if abs(klass(*e)[0]-0.5)<1e-6 else A,'4.8^2 [8|8=a]')
    elif name=='4.6.12':
        r12=R*math.cos(math.pi/12)
        def w(*e):
            dd,an=klass(*e)
            if abs(dd-r12)>1e-6: return 1+0*Z      # 4|6
            return B if abs((an/60)-round(an/60))<1e-6 else A   # 4|12=b, 6|12=a
        run(V,E,w,'4.6.12 [4|6=1, 6|12=a, 4|12=b]')
    else:
        # kvadratkanter: mittpunkt på avstånd 1/2 från ett kvadratcentrum; övriga 3|3 = a
        cent=[(0,0),(A_/2,A_/2)]
        def w(i,j,n1,n2):
            mx,my=mid(i,j,n1,n2)
            for (cx,cy) in cent:
                for k1 in range(-3,4):
                    for k2 in range(-3,4):
                        if abs(math.hypot(mx-cx-k1*A_,my-cy-k2*A_)-0.5)<1e-6: return 1+0*Z
            return A
        run(V,E,w,'snub [3|4=1, 3|3=a]')
