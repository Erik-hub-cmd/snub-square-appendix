r"""Separat kontroll 3 av Sats B: Kenyons formel (egen numerisk konturintegral) på Z²_t.
(a) rombinbäddningens konduktanser = Z²_t:s vikter; (b) ΔG = δ; (c) symmetri; (d) t = 1 ger vanliga Z²;
(e) −G(v0,v1) − G(v1,v0) = resistansen ur grundintegralen för Z²_t; för t > 1 och t < 1.
Se GRANSKNING4.md, avsnitt 3. Kör: .venv\Scripts\python kod\granskning4_k2_kenyon.py"""
import sys
sys.stdout.reconfigure(encoding="utf-8")
import mpmath as mp, sympy as sp
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from granskning4_egen import *

ok = []
def check(name, cond):
    ok.append(bool(cond)); print(('OK  ' if cond else 'FEL ') + name, flush=True)

DPS = 25
for tv in [mp.mpf(4) + mp.sqrt(15), 1 / (mp.mpf(4) + mp.sqrt(15)), mp.mpf(2)]:
    mp.mp.dps = DPS; Ph = Phi_of(tv)
    print(f'--- t = {mp.nstr(tv, 12)}, Φ = {mp.nstr(Ph, 12)}')
    # (a)
    bad = 0
    for p in [(0, 0), (1, 0), (0, 1), (1, 1), (-3, 2)]:
        for r, w in neighbours(p, None, hz(tv)):
            if abs(cond_geom(p, r, Ph) - w) > mp.mpf(10) ** (-DPS + 3): bad += 1
    check('(a) tan(halv rombvinkel i primalnoden) = vikterna 1/t, t (horisontellt, x+y jämnt/udda) och 1 (vertikalt)', bad == 0)
    # (b)
    worst = 0
    for y in [(0, 0), (1, 0), (2, 1), (-1, 2), (3, -1)]:
        L = mp.fsum(w * (kenyon((0, 0), y, tv, DPS) - kenyon((0, 0), r, tv, DPS)) for r, w in neighbours(y, None, hz(tv)))
        worst = max(worst, abs(L - (1 if y == (0, 0) else 0)))
    check(f'(b) ΔG = δ i 5 noder (största avvikelse {mp.nstr(worst, 3)})', worst < mp.mpf(10) ** (-DPS + 5))
    # (c)
    asym = max(abs(kenyon(p, q, tv, DPS) - kenyon(q, p, tv, DPS)) for p, q in [((0, 0), (2, 1)), ((1, 0), (-2, 3)), ((0, 1), (3, 3))])
    check(f'(c) G(v0, v1) = G(v1, v0) (avvikelse {mp.nstr(asym, 3)}), och G är reellt', asym < mp.mpf(10) ** (-DPS + 5)
          and abs(mp.im(kenyon((0, 0), (3, 2), tv, DPS))) < mp.mpf(10) ** (-DPS + 5))
    # (e)
    M = bloch(None, hz(sp.Symbol('t')))
    gr = Grund(M, tval=tv)
    for (p, q) in [((0, 0), (1, 0)), ((0, 0), (2, 1)), ((1, 1), (-2, 3)), ((0, 1), (4, 2))]:
        a, ca = (TYPES.index((p[0] % 2, p[1] % 2)), (p[0] // 2, p[1] // 2))
        b, cb = (TYPES.index((q[0] % 2, q[1] % 2)), (q[0] // 2, q[1] // 2))
        Rf = gr.R(a, b, (cb[0] - ca[0], cb[1] - ca[1]), dps=DPS)
        Rk = -kenyon(p, q, tv, DPS) - kenyon(q, p, tv, DPS)
        check(f'(e) R{p}–{q}: grundintegral {mp.nstr(Rf, 18)}  Kenyon {mp.nstr(mp.re(Rk), 18)}', abs(Rf - Rk) < mp.mpf(10) ** (-DPS + 6))
# (d) vanliga Z²
mp.mp.dps = DPS
g1, g2, g3 = kenyon((0, 0), (1, 0), 1, DPS), kenyon((0, 0), (1, 1), 1, DPS), kenyon((0, 0), (2, 0), 1, DPS)
check(f'(d) t = 1: G(e1) = {mp.nstr(mp.re(g1), 15)} (−1/4), G(1,1) = {mp.nstr(mp.re(g2), 15)} (−1/π), G(2,0) = {mp.nstr(mp.re(g3), 15)} (−(1 − 2/π))',
      abs(g1 + mp.mpf(1) / 4) < 1e-20 and abs(g2 + 1 / mp.pi) < 1e-20 and abs(g3 + 1 - 2 / mp.pi) < 1e-20)
print(f'\n{sum(ok)} godkända, {len(ok) - sum(ok)} fel')
