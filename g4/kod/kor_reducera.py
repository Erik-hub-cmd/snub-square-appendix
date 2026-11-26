import json, sys, math, time, pickle
sys.path.insert(0, '.')
from natgalois import from_json, geom_edges
from reducera import Net, search
files = sys.argv[1].split(','); which = sys.argv[2]; depth = int(sys.argv[3]); maxst = int(sys.argv[4])
klass = pickle.load(open('../loggar/sammanstallning.pkl', 'rb'))
kl = {}
items = []
for fn in files:
    for t in json.load(open(fn)): items.append(t)
names = {}
res = []
for idx, t in enumerate(items):
    nm = t.get('name') or ' ; '.join(t['types'])
    full = f"[{idx}] {nm} (k={t['summary']['k']}, V={t['summary']['V']})"
    c = [x for x in klass if x[1] == full]
    cl = c[0][2] if c else '?'
    if which != 'alla' and cl != which: continue
    G = Net.from_geom(*from_json(t))
    t0 = time.time()
    H, msg, path = search(G, maxdepth=depth, maxstates=maxst, log=lambda s: None)
    V = len(H.rot); F = len(H.faces())
    print(f'{full:55s} [{cl}]  E={G.E()} -> {H.E()}  {msg}  (Euler {V - H.E() + F}, {len(path)} Y-Δ-drag, {time.time() - t0:.0f}s)', flush=True)
    res.append((full, cl, G.E(), H.E(), msg, len(path)))
pickle.dump(res, open(f'../loggar/reducera_{which.replace(" ", "_")}_{files[0].split(".")[0]}.pkl', 'wb'))
