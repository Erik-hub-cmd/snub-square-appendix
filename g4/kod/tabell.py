import pickle, sympy as sp
U = sp.Symbol('u')
rows = ['| # | tessellation | V | E | minimalt | sidpolynom (lika vikter, irreducibla faktorer av grad ≥ 2) | kropp | klass |', '|---|---|---|---|---|---|---|---|']
for f in ['u2_enhet', 'u3_enhet', 'u4_enhet']:
    R = pickle.load(open(f'../loggar/{f}.pkl', 'rb'))
    for r in R:
        irr = []; fields = set()
        for s in r['sides'][:len(r['sides']) // 2]:
            for fa in s['factors']:
                if fa['deg'] >= 2:
                    irr.append(fa['f'])
                    if fa['deg'] == 2:
                        p = sp.Poly(sp.sympify(fa['f']), U); a, b, c = p.all_coeffs()
                        D = sp.factorint(b * b - 4 * a * c); sq = 1
                        for pr, e in D.items():
                            if e % 2: sq *= pr
                        fields.add(f'Q(√{sq})')
                    else: fields.add(f'grad {fa["deg"]}')
        irr = sorted(set(irr))
        cls = 'minimalt' if r['min']['minimal'] else ('**Galois-hindrat**' if r['Hinf'] else ('hinder egen torus' if r['H0'] else 'icke-minimalt, inget hinder'))
        name = r['name'].split(' (k=')[0]
        k = r['name'].split('(k=')[1].split(',')[0]
        rows.append(f"| {k}-{name.split(']')[0][1:]} | {name.split('] ')[1]} | {r['V']} | {r['E']} | {'ja' if r['min']['minimal'] else 'nej'} | {'; '.join(irr) or '–'} | {', '.join(sorted(fields)) or '–'} | {cls} |")
open('../TABELL_ALLA.md', 'w').write('# Galois-testet på alla 1-, 2-, 3- och 4-uniforma tessellationer i projektet (lika vikter)\n\nGenererad av `kod/tabell.py` ur `loggar/u*_enhet.pkl`. Kolumnen # är k-index i indatafilen. 4.8² (inte i Z[ζ₁₂]-datan) finns i `loggar/arkimedes_enhet.log`.\n\n' + '\n'.join(rows) + '\n')
print(len(rows) - 2)
