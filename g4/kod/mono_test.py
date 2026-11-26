import json, math, sys, random
sys.path.insert(0,'.')
from natgalois import *
from fractions import Fraction
d = json.load(open('uniform12.json'))
items=[(t['name'],)+from_json(t) for t in d if t['summary']['k']==1]
s=1/math.sqrt(2); L=1+math.sqrt(2); pos=[s,1j*s,-s,-1j*s]
items.append(('4.8^2',pos,L,1j*L,geom_edges(pos,L,1j*L)))
random.seed(5)
for name,pos,v1,v2,E in items:
    W=[random.randint(1,30) for _ in E]
    T=tracks(pos,v1,v2,E); mu=monodromy(T,W)
    P=det_poly(len(pos),E,W); H,S=sides(P)
    print('==',name)
    for t,m in zip(T,mu): print('   spår h=',t['h'],'μ=',m)
    for sd in S[:len(S)//2]:
        c=sd['poly'].all_coeffs()  # högsta först
        lead=c[0]; const=c[-1]
        print('   sida',sd['dir'],'längd',sd['len'],'produkt av rötter =',sp.Rational(const,lead)*(-1)**sd['len'], ' faktorer', sp.factor(sd['poly'].as_expr()))
