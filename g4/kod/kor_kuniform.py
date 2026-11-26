import json, sys, time, random, pickle
sys.path.insert(0, '.')
from natgalois import *
fil, mode, ut = sys.argv[1], sys.argv[2], sys.argv[3]
d = []
for f in fil.split(','):
    d += json.load(open(f))
res = []
random.seed(7)
for idx, t in enumerate(d):
    pos, v1, v2, E = from_json(t)
    W = [1] * len(E) if mode == 'enhet' else [random.randint(1, 97) for _ in E]
    nm = t.get('name') or ' ; '.join(t['types'])
    nm = f"[{idx}] {nm} (k={t['summary']['k']}, V={len(pos)})"
    t0 = time.time()
    try:
        r = analyse(nm, pos, v1, v2, E, W)
        print(report(r), f'  ({time.time()-t0:.1f}s)', flush=True)
        r['sides'] = [{k: (str(v) if k == 'poly' else v) for k, v in s.items() if k != 'factors'} | {'factors': [{kk: vv for kk, vv in f.items() if kk != 'rotter'} for f in s['factors']]} for s in r['sides']]
        res.append(r)
    except Exception as e:
        print('==', nm, 'FEL', repr(e), flush=True)
    pickle.dump(res, open(ut, 'wb'))
