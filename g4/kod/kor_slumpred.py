import json, sys, time
sys.path.insert(0, '.')
from natgalois import from_json
from reducera import Net, random_walks
files = sys.argv[1].split(','); idxs = [int(x) for x in sys.argv[2].split(',')]; walks = int(sys.argv[3])
items = []
for fn in files: items += json.load(open(fn))
for i in idxs:
    t = items[i]; nm = t.get('name') or ' ; '.join(t['types'])
    G = Net.from_geom(*from_json(t)); t0 = time.time()
    H, msg = random_walks(G, walks=walks, length=60, seed=i)
    print(f"[{i}] {nm}: E={G.E()} -> {H.E()}  {msg}  ({time.time()-t0:.0f}s)", flush=True)
