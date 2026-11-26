"""Exporterar datafilerna som artikel 1 listar som bilagor (Appendix H.3) till artikel/ancillary/:
  resistances_lambda1.csv  - de 1763 exakta resistanserna från nod 1 (runs/tabell/checkpoint.json)
  closed_forms_lambda.csv  - a(λ), b(λ), c(λ) för de 53 klasserna (runs/isoradiell_lambda/SLUTNA_FORMER.md)
Nodnumreringen j är 1-baserad som i artikeln (j = 1..4), n = (n1, n2) cellförskjutning.
Kör: .venv\\Scripts\\python kod\\export_artikel1_data.py
"""
import csv, json, os, re

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, '..')
OUT = os.path.join(ROOT, 'artikel', 'ancillary')
os.makedirs(OUT, exist_ok=True)

d = json.load(open(os.path.join(ROOT, 'runs', 'tabell', 'checkpoint.json'), encoding='utf-8'))['results']
rows = []
for part in d.values():
    for k, v in part.items():
        j, n1, n2 = (int(x) for x in k.split(':'))
        rows.append((v['dist'], j + 1, n1, n2, v['expr']))
rows.sort()
with open(os.path.join(OUT, 'resistances_lambda1.csv'), 'w', newline='', encoding='utf-8') as f:
    w = csv.writer(f)
    w.writerow(['j', 'n1', 'n2', 'euclidean_distance', 'R_1j(n) exact (sympy syntax)'])
    for dist, j, n1, n2, e in rows:
        w.writerow([j, n1, n2, f'{dist:.12f}', e])

md = open(os.path.join(ROOT, 'runs', 'isoradiell_lambda', 'SLUTNA_FORMER.md'), encoding='utf-8').read()
out = []
for line in md.splitlines():
    m = re.match(r'\| (\d+):(-?\d+):(-?\d+) \| (.*?) \| (.*?) \| (.*?) \|', line)
    if m:
        j, n1, n2 = int(m.group(1)) + 1, int(m.group(2)), int(m.group(3))
        clean = lambda s: s.strip().strip('$')
        out.append((j, n1, n2, clean(m.group(4)), clean(m.group(5)), clean(m.group(6))))
with open(os.path.join(OUT, 'closed_forms_lambda.csv'), 'w', newline='', encoding='utf-8') as f:
    w = csv.writer(f)
    w.writerow(['j', 'n1', 'n2', 'a(lambda) LaTeX', 'b(lambda) LaTeX', 'c(lambda) LaTeX'])
    w.writerows(out)
print(len(rows), 'resistanser,', len(out), 'slutna former')
