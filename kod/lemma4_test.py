import random, sys
from lemma4_natverk import *
rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 1)
DIRS = [(1,0),(0,1),(1,1),(1,-1),(2,1),(1,2),(2,-1),(1,-2),(3,1),(1,3)]
ok = True; antal = 0
while antal < int(sys.argv[2] if len(sys.argv) > 2 else 15):
    m = rng.randint(2, 5)
    H = [rng.choice(DIRS) for _ in range(m)]
    if sum(a for a,b in H) % 2 or sum(b for a,b in H) % 2: continue
    if len({(a*1.0/math.hypot(a,b), b/math.hypot(a,b)) if (a>0 or (a==0 and b>0)) else (-a/math.hypot(a,b), -b/math.hypot(a,b)) for a,b in H}) < 2: continue
    E_pred = sum(abs(H[i][0]*H[j][1]-H[i][1]*H[j][0]) for i in range(m) for j in range(i+1, m))
    if E_pred > 14: continue
    O = [rng.random() for _ in H]
    n, K = bygg(H, O)
    if n > 8: continue
    cond = [sp.Rational(rng.randint(1, 9), rng.randint(1, 9)) for _ in K]
    T = blochdet(n, K, cond)
    N = hull(list(T)); Zg = zonogon(H)
    lika = sorted(N) == sorted(Zg)
    sp_ok = True; info = []
    for st, g, f in sidpolynom(T):
        lin = all(sp.Poly(fa, w).degree() <= 1 for fa, mm in sp.factor_list(f, w)[1])
        sp_ok &= lin; info.append((st, g, lin))
    ok &= lika and sp_ok and len(K) == E_pred
    antal += 1
    print(f'H={H} noder={n} kanter={len(K)} (pred {E_pred})  N(P)=Z: {lika}  alla sidpolynom linjära över Q: {sp_ok}', '' if lika else f'\n   N={N}\n   Z={Zg}')
print('ALLT OK' if ok else 'AVVIKELSE')
