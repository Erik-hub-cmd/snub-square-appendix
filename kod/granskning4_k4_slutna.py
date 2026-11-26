r"""Separat kontroll av Sats B (Sats 3 i arbetsanteckningarna): slutna former a(λ), b(λ), c(λ) för tre noder UTANFÖR de 53, räknade med projektets
steg 9-funktioner (import, ingen fil skrivs), utvärderade i λ = 1/2, 3 och 1.
Se GRANSKNING4.md, avsnitt 5. Kör: .venv\Scripts\python kod\granskning4_k4_slutna.py"""
import sys, os
sys.stdout.reconfigure(encoding="utf-8")
sys.argv = ['x', '0']
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import isoradiell_steg9_kenyon_t as s9
import sympy as sp, json
# tillagt: jämförelse med den exakta tabellen i λ = 1 (runs/tabell/checkpoint.json), där R = a + c·√3/3 + b/π
TAB = {}
for chunk in json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'runs', 'tabell', 'checkpoint.json'),
                            encoding='utf-8'))['results'].values():
    TAB.update(chunk)
nok = 0; nfel = 0
for nd in ['2:-2:2', '1:-2:-2', '0:-2:2']:
    j, n1, n2 = map(int, nd.split(':'))
    nn = (n1 + s9.OFF[j][0], n2 + s9.OFF[j][1])
    r0, r1, r2 = s9.R_snub(s9.SIG[0], s9.SIG[j], nn)
    sym = s9.tau(r0) == r0 and s9.tau(r1) == r1 and s9.tau(r2) == -r2
    fa, fb, fc = (s9.to_lambda(f) for f in (r0, r1, r2 / s9.KAPPA))
    print(f'{nd}: τ-symmetri {sym}; a = {fa}; b = {fb}; c = {fc}')
    for lv in (sp.Rational(1, 2), sp.Integer(3), sp.Integer(1)):
        print(f'   λ = {lv}: a, b, c = {[sp.cancel(f.subs(s9.lam_s, lv)) for f in (fa, fb, fc)]}', flush=True)
    a1, b1, c1 = (sp.cancel(f.subs(s9.lam_s, 1)) for f in (fa, fb, fc))
    tab = [sp.Rational(x) for x in TAB[nd]['coeffs']]   # (1, √3, 1/π)
    ok = sym and tab == [a1, c1 / 3, b1]
    nok += ok; nfel += not ok
    print(f"{'OK ' if ok else 'FEL'} {nd}: τ-symmetri och λ = 1 lika med tabellen ({TAB[nd]['expr']})", flush=True)
print(f'{nok} godkända, {nfel} fel')
