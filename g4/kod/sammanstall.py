import pickle, sys, collections
L = '../loggar/'
def load(p): return pickle.load(open(L + p, 'rb'))
out = []
for grp, fe, fs in [('1-uniforma (arkimediska, utom 4.8²)', 'u2_enhet.pkl', 'u2_slump.pkl'), ('3-uniforma', 'u3_enhet.pkl', 'u3_slump.pkl'), ('4-uniforma', 'u4_enhet.pkl', 'u4_slump.pkl')]:
    E, S = load(fe), load(fs)
    for re_, rs in zip(E, S):
        assert re_['name'] == rs['name']
        k = int(re_['name'].split('(k=')[1].split(',')[0])
        g = {1: '1-uniforma', 2: '2-uniforma', 3: '3-uniforma', 4: '4-uniforma'}[k]
        cls = 'minimalt' if re_['min']['minimal'] else ('Galois-hindrat' if re_['Hinf'] else ('hinder bara på egen torus' if re_['H0'] else 'ej minimalt, inget hinder'))
        cls_s = 'minimalt' if rs['min']['minimal'] else ('Galois-hindrat' if rs['Hinf'] else ('hinder bara på egen torus' if rs['H0'] else 'ej minimalt, inget hinder'))
        out.append((g, re_['name'], cls, cls_s, re_['N_eq_Z']))
cnt = collections.Counter((g, c, cs) for g, n, c, cs, nz in out)
for k in sorted(cnt): print(k, cnt[k])
print()
for g, n, c, cs, nz in out:
    if c != cs: print('OLIKA enhet/slump:', n, c, '|', cs)
print()
for g, n, c, cs, nz in out:
    if g in ('1-uniforma', '2-uniforma'): print(f'{n:60s} enhet: {c:28s} slump: {cs:28s} N=Z:{nz}')
pickle.dump(out, open(L + 'sammanstallning.pkl', 'wb'))
