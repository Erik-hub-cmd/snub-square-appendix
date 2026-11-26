"""
Isoradiell förklaring, steg 9.1: kontroll av ALLA 159 värden i runs/isoradiell_lambda (λ = 1/2, 2, 3; 53 noder var)
mot den snävare formen
    R = a + b/π + c·√d·Φ/π,   a, b, c ∈ Q,   d = λ(λ+2),  Φ = arccos(1/(1+λ)) = arctan √d.
För varje värde:
  (1) PSLQ mot hela basen (1, √d, 1/π, √d/π, Φ/π, √dΦ/π) (som i nattkörningen) – termerna √d, √d/π, Φ/π ska vara 0;
  (2) PSLQ mot den snäva basen (1, 1/π, √dΦ/π) – samma a, b, c ska komma ut;
  (3) a + b/π + c√dΦ/π räknas om på 100 siffror och jämförs med båda värdena (75 och 60 siffror).
De två värdena kommer från samma kod (kod/jobb_isoradiell_lambda.py) i två precisioner: de skyddar mot avrundningsfel,
inte mot systematiska fel (granskning 4, fynd 5). Separata kontroller: isoradiell_steg9_kenyon_t.py och granskning4_k3_stickprov.py.
Skriver hela tabellen (159 rader) till runs/isoradiell_lambda/KONTROLL_159.md.
Kör: python kod/isoradiell_steg9_kontroll159.py
"""
import sys, os, json
sys.stdout.reconfigure(encoding="utf-8")
import mpmath as mp
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
RUN = os.path.join(HERE, '..', 'runs', 'isoradiell_lambda')
LAMBDAS = ['1/2', '2', '3']


def basis(lv):
    s = mp.sqrt(lv * (lv + 2)); Phi = mp.acos(1 / (1 + lv))
    return s, Phi


def pslq(v, B, dps):
    mp.mp.dps = dps - 3
    rel = mp.pslq([v] + B, tol=mp.mpf(10) ** (-(dps - 8)), maxcoeff=10 ** 9, maxsteps=10 ** 6)
    if rel is None or rel[0] == 0:
        return None
    return [F(-r, rel[0]) for r in rel[1:]]


def main():
    res = json.load(open(os.path.join(RUN, 'checkpoint.json'), encoding='utf-8'))['results']
    data = {}
    for cid, r in res.items():
        ls, dps, tg = cid.split('|')
        data.setdefault((ls, tg), {})[int(dps)] = r
    out = {}
    L = ['# Kontroll av alla 159 värden: snäv form R = a + b/π + c·√d·Φ/π (steg 9.1)', '',
         'd = λ(λ+2), Φ = arccos(1/(1+λ)) = arctan √d. Kolumnen "fel 75/60" är |R − (a + b/π + c√dΦ/π)| för de två',
         'beräkningarna (75 resp. 60 siffror), uttryckt som antal korrekta decimaler. De två beräkningarna är samma kod',
         '(kod/jobb_isoradiell_lambda.py) i två precisioner: de skyddar mot avrundningsfel, inte mot systematiska fel.', '']
    tot_ok = 0
    for ls in LAMBDAS:
        lv = mp.mpf(F(ls).numerator) / F(ls).denominator
        rows = sorted(((k[1], d) for k, d in data.items() if k[0] == ls), key=lambda kd: (kd[1][75]['dist'], kd[0]))
        L += [f'## λ = {ls}  (d = {F(ls) * (F(ls) + 2)})', '',
              '| nod j:n1:n2 | avstånd | a | b | c | 6-bas: √d, √d/π, Φ/π | fel 75/60 (siffror) | OK |',
              '|---|---|---|---|---|---|---|---|']
        nok = 0
        for tg, d in rows:
            mp.mp.dps = 90
            v75, v60 = mp.mpf(d[75]['value']), mp.mpf(d[60]['value'])
            s, Phi = basis(lv)
            full = pslq(v75, [mp.mpf(1), s, 1 / mp.pi, s / mp.pi, Phi / mp.pi, s * Phi / mp.pi], 75)
            mp.mp.dps = 90
            s, Phi = basis(lv)
            snav = pslq(v75, [mp.mpf(1), 1 / mp.pi, s * Phi / mp.pi], 75)
            ok = full is not None and snav is not None
            if ok:
                ok = full[1] == 0 and full[3] == 0 and full[4] == 0 and [full[0], full[2], full[5]] == snav
            e75 = e60 = 0
            if snav is not None:
                mp.mp.dps = 110
                s, Phi = basis(lv)
                a, b, c = (mp.mpf(x.numerator) / x.denominator for x in snav)
                val = a + b / mp.pi + c * s * Phi / mp.pi
                e75 = int(-mp.log10(abs(val - v75) + mp.mpf(10) ** -105))
                e60 = int(-mp.log10(abs(val - v60) + mp.mpf(10) ** -105))
                ok = ok and e75 >= 70 and e60 >= 58
            nok += ok
            out[f'{ls}|{tg}'] = [str(x) for x in snav] if snav else None
            ft = 'ingen' if full is None else f'{full[1]}, {full[3]}, {full[4]}'
            A_, B_, C_ = (str(x) for x in snav) if snav else ('?', '?', '?')
            L.append(f'| {tg} | {d[75]["dist"]:.3f} | {A_} | {B_} | {C_} | {ft} | {e75}/{e60} | {"ja" if ok else "NEJ"} |')
        L += ['', f'**λ = {ls}: {nok} av {len(rows)} i den snäva formen.**', '']
        print(f'λ = {ls}: {nok} av {len(rows)} i den snäva formen', flush=True)
        tot_ok += nok
    L.insert(5, f'**Totalt: {tot_ok} av {len(data)} värden i den snäva formen (alla tre villkor uppfyllda).**\n')
    open(os.path.join(RUN, 'KONTROLL_159.md'), 'w', encoding='utf-8').write('\n'.join(L) + '\n')
    json.dump(out, open(os.path.join(RUN, 'koefficienter_159.json'), 'w', encoding='utf-8'), indent=1, ensure_ascii=False)
    print(f'\n{tot_ok} av {len(data)} godkända, {len(data) - tot_ok} fel')


if __name__ == '__main__':
    main()
