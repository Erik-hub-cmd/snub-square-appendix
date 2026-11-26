import json, math, cmath, sys, time
sys.path.insert(0, '.')
from natgalois import *
d = json.load(open('uniform12.json'))
arch = [t for t in d if t['summary']['k'] == 1]
items = []
for t in arch:
    pos, v1, v2, E = from_json(t); items.append((t['name'], pos, v1, v2, E))
# 4.8^2 för hand: kvadrat vriden 45°, gitter (1+√2)
s = 1 / math.sqrt(2); L = 1 + math.sqrt(2)
pos = [s, 1j * s, -s, -1j * s]
E = geom_edges(pos, L, 1j * L)
items.append(('4.8^2', pos, L, 1j * L, E))
mode = sys.argv[1] if len(sys.argv) > 1 else 'enhet'
random.seed(1)
for name, pos, v1, v2, E in items:
    W = [1] * len(E) if mode == 'enhet' else [random.randint(1, 97) for _ in E]
    t0 = time.time(); r = analyse(name, pos, v1, v2, E, W)
    print(report(r), f'  ({time.time()-t0:.1f}s)', flush=True)
