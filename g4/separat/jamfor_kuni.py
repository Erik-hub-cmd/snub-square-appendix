"""Compare the second, separately written computation (kuni.log, from kuni.py) with the first program
(../loggar/u2_enhet.pkl, u3_enhet.pkl, u4_enhet.pkl, from ../kod/kor_kuniform.py), tiling by tiling.

For every tiling: same number of vertices and edges, and the same verdict "Galois-obstructed" with unit conductances
and with random conductances. Prints the counts per k as in Appendix G.4 and a summary line.
Run from this folder after kuni.py and kor_kuniform.py (both are in expected.json)."""
import json, pickle, collections, sys

rows = [json.loads(l) for l in open('kuni.log', encoding='utf-8') if l.strip()]
first = []
for f in ('u2_enhet', 'u3_enhet', 'u4_enhet'):
    first += pickle.load(open(f'../loggar/{f}.pkl', 'rb'))
slump = []
for f in ('u2_slump', 'u3_slump', 'u4_slump'):
    slump += pickle.load(open(f'../loggar/{f}.pkl', 'rb'))
ok = fel = 0
cnt = collections.Counter()
if not (len(rows) == len(first) == len(slump)):
    print(f'FEL olika antal: kuni.log {len(rows)}, första programmet {len(first)}/{len(slump)}'); sys.exit(1)
for k, (r, a, s) in enumerate(zip(rows, first, slump)):
    kk = int(a['name'].split('(k=')[1].split(',')[0])
    same = (r['V'] == a['V'] and r['E'] == a['E'] and r['badlen'] == 0 and r['degok']
            and bool(r['enhet_obstr']) == bool(a['Hinf']) and bool(r['slump_obstr']) == bool(s['Hinf']))
    cls = 'minimalt' if a['min']['minimal'] else ('Galois-hindrat' if a['Hinf'] else 'ej minimalt, inget hinder')
    cnt[(kk, cls)] += 1
    if same: ok += 1
    else:
        fel += 1; print('FEL', a['name'], r)
for key in sorted(cnt): print(f'k={key[0]} {key[1]}: {cnt[key]}')
print(f'{ok} godkända, {fel} fel  (samma V, E och hinderklass med lika och med slumpade vikter i båda programmen)')
