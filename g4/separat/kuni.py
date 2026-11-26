import json, sys, cmath, math, random, time
sys.path.insert(0,'.')
from rev import *
z12=cmath.exp(1j*math.pi/6)
pos=lambda c: sum(c[k]*z12**k for k in range(4))
out=[]
for fn in sys.argv[1:]:
    data=json.load(open(fn))
    for idx,t in enumerate(data):
        V=len(t['verts']); P=[pos(v[0]) for v in t['verts']]; L1=pos(t['v1']); L2=pos(t['v2'])
        bad=0; E=[]
        for (i,j,(dx,dy),cl) in t['edges']:
            if abs(abs(P[j]+dx*L1+dy*L2-P[i])-1)>1e-9: bad+=1
            E.append((i,j,dx,dy,1))
        # gradkontroll mot vertextyp
        deg=[0]*V
        for e in E: deg[e[0]]+=1; deg[e[1]]+=1
        degok = all(deg[k]==len(t['verts'][k][1]) for k in range(V))
        res={}
        random.seed(idx*1000+V)
        for lab,EE in [('enhet',E),('slump',[(i,j,dx,dy,random.randint(1,97)) for (i,j,dx,dy,_) in E])]:
            d,H,area,sides=analyse(V,EE)
            obs=[ (s,g) for s,Q in sides for g in galois(Q) if g[2]]
            irr=[ (s,g) for s,Q in sides for g in galois(Q)]
            res[lab]=dict(area=float(area), obstr=bool(obs), irr=bool(irr), allreal=all(g[3] for _,g in irr),
                          polys=sorted(set(g[0] for _,g in irr)) if lab=='enhet' else None)
        r=dict(fil=fn.split('/')[-1],idx=idx,name=t.get('name'),V=V,E=len(E),badlen=bad,degok=degok,**{k+'_'+a:v for k,dd in res.items() for a,v in dd.items()})
        print(json.dumps(r),flush=True)
