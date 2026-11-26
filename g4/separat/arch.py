import math, sys
sys.path.insert(0,'.')
from rev import *
def geo(L1,L2,pts,tol=1e-9):
    # dedupe mod gitter
    det=L1[0]*L2[1]-L1[1]*L2[0]
    def red(p):
        a=( p[0]*L2[1]-p[1]*L2[0])/det; b=(L1[0]*p[1]-L1[1]*p[0])/det
        fa=math.floor(a+1e-9); fb=math.floor(b+1e-9)
        return (round(a-fa,6)%1.0, round(b-fb,6)%1.0), (fa,fb)
    reps={}
    for p in pts:
        k,(fa,fb)=red(p)
        if k not in reps: reps[k]=(p[0]-fa*L1[0]-fb*L2[0], p[1]-fa*L1[1]-fb*L2[1])
    P=list(reps.values()); V=len(P); E=[]
    for i in range(V):
        for j in range(V):
            for n1 in range(-2,3):
                for n2 in range(-2,3):
                    q=(P[j][0]+n1*L1[0]+n2*L2[0], P[j][1]+n1*L1[1]+n2*L2[1])
                    if abs(math.dist(P[i],q)-1)<tol:
                        if (i,j,n1,n2) < (j,i,-n1,-n2): E.append((i,j,n1,n2,1))
    deg=[0]*V
    for e in E: deg[e[0]]+=1; deg[e[1]]+=1
    global LASTPOS; LASTPOS=(P,L1,L2)
    return V,E,deg
s2=math.sqrt(2); s3=math.sqrt(3)
T={}
a=1+s2; T['4.8^2']=geo((a,0),(0,a),[(s2/2,0),(-s2/2,0),(0,s2/2),(0,-s2/2)])
tri=lambda a,b:(a+b/2, b*s3/2)
T['3^4.6']=geo(tri(2,1),tri(-1,3),[tri(a,b) for a in range(-7,8) for b in range(-7,8) if (a-2*b)%7])
R=1/(2*math.sin(math.pi/12)); c=3+s3
T['4.6.12']=geo((c,0),(c/2,c*s3/2),[(R*math.cos(math.radians(15+30*k)),R*math.sin(math.radians(15+30*k))) for k in range(12)])
# snub square: två kvadrater, prova vridningar
A=(1+s3)/s2
for phi in [15,-15,30,-30]:
    pts=[]
    for (cx,cy,sg) in [(0,0,1),(A/2,A/2,-1)]:
        for k in range(4):
            t=math.radians(45+sg*phi+90*k); pts.append((cx+math.cos(t)/s2, cy+math.sin(t)/s2))
    V,E,deg=geo((A,0),(0,A),pts)
    if V==4 and deg==[5]*4: T['3^2.4.3.4']=(V,E,deg); break
for name,(V,E,deg) in T.items():
    d,H,area,sides=analyse(V,E)
    print(name,'V',V,'E',len(E),'grad',sorted(set(deg)),'Area(N)',area,'hörn',H)
    for s,Q in sides: print('   sida',s,sympy.factor(Q), galois(Q))
