import sys,math; sys.path.insert(0,'.')
src=open('arch.py').read().split("T={}")[0]; g={}; exec(src,g)
from rev import analyse
import sympy
s2=math.sqrt(2); a=1+s2; base=[(s2/2,0),(-s2/2,0),(0,s2/2),(0,-s2/2)]
for (m1,m2) in [(2,1),(2,2),(3,1)]:
    pts=[(p[0]+i*a,p[1]+j*a) for p in base for i in range(-3,4) for j in range(-3,4)]
    V,E,deg=g['geo']((m1*a,0),(0,m2*a),pts); d,H,area,sides=analyse(V,E)
    print((m1,m2),'V',V,'E',len(E),[(s,str(sympy.factor(Q))) for s,Q in sides[:2]])
