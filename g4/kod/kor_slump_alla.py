import json, sys, time, pickle
sys.path.insert(0, '.')
from natgalois import from_json
from reducera import Net, random_walks
files = sys.argv[1].split(','); walks = int(sys.argv[2]); ut = sys.argv[3]
klass = {x[1]: x[2] for x in pickle.load(open('../loggar/sammanstallning.pkl', 'rb'))}
items = []
for fn in files: items += json.load(open(fn))
res = []
for i, t in enumerate(items):
    nm = t.get('name') or ' ; '.join(t['types'])
    full = f"[{i}] {nm} (k={t['summary']['k']}, V={t['summary']['V']})"
    cl = klass.get(full, '?')
    if cl == 'minimalt': continue
    G = Net.from_geom(*from_json(t)); t0 = time.time()
    H, msg = random_walks(G, walks=walks if cl != 'Galois-hindrat' else walks // 3, length=60, seed=i)
    ok = msg.startswith('MINIMALT') or msg.startswith('minimalt')
    flag = '   <== MOTSÄGELSE MOT SATSEN' if ok and cl == 'Galois-hindrat' else ''
    print(f"{full:60s} [{cl}]  E={G.E()} -> {H.E()}  {msg}  ({time.time()-t0:.0f}s){flag}", flush=True)
    res.append((full, cl, G.E(), H.E(), ok))
    pickle.dump(res, open(ut, 'wb'))
