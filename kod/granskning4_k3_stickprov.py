r"""Separat kontroll 5 av Sats B: snub square med λ = 1/2 och λ = 3.
(i) egen grundintegral (egen Bloch-matris);  (ii) Sats 3-värdet a + b/π + c√dΦ/π med a, b, c ur RESULTAT.md resp.
anteckningens slutna former;  (iii) bevisets kedja numeriskt: Σ conj(p)·r·G_Kenyon med A(t), B(t) (steg 8) och MIN
Kenyon-konturintegral, för t > 1 och (några noder) för 1/t;  (iv) för noder utanför de 53: PSLQ mot (1, 1/π, √dΦ/π).
Se GRANSKNING4.md, avsnitt 5. Kör: .venv\Scripts\python kod\granskning4_k3_stickprov.py (ca 6 min)"""
import sys, os, time
sys.stdout.reconfigure(encoding="utf-8")
import sympy as sp, mpmath as mp
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from granskning4_egen import *

ROOT = os.path.dirname(os.path.abspath(__file__))
SIG = (0, 2, 3, 1); OFF = [(0, 0), (0, 0), (0, -1), (0, -1)]
lam_t = (t - 1)**2 / (2 * t)

src = open(os.path.join(ROOT, 'isoradiell_steg8_allmant_lambda.py'), encoding='utf-8').read()
ns = {'sp': sp, 'X': X, 'Z': Z, 't': t}
exec(src[src.index('A = sp.Matrix'):src.index("\ncheck('")], ns)
A, B = ns['A'], ns['B']
D = 4 * (t + 1)**4 * (t**2 + 1)**2 / (t - 1)**4
AH = A.subs({X: 1 / X, Z: 1 / Z}, simultaneous=True).T          # A*: koefficienterna är reella för reellt t
Binv = B.adjugate() / D


def zmap(j, n):
    """artikelnod (j, n1, n2) -> (Z²-modellens typ, cell) enligt anteckningens steg 6 (SIG, OFF)."""
    return SIG[j], (n[0] + OFF[j][0], n[1] + OFF[j][1])


def laurent2(e):
    e = sp.cancel(sp.together(e)); num, den = sp.fraction(e)
    pd = sp.Poly(den, X, Z); assert len(pd.terms()) == 1
    (dx, dz), dc = pd.terms()[0]
    return {(ex - dx, ez - dz): sp.cancel(c / dc) for (ex, ez), c in sp.Poly(sp.expand(num), X, Z).terms()}


def canon(p0, p1):
    """Flytta p0 till (0,0) (x+y jämnt) eller (1,0) (x+y udda); translationer (a, b) med a + b jämnt bevarar
    Z²_t och rombinbäddningen."""
    sh = (p0[0], p0[1]) if (p0[0] + p0[1]) % 2 == 0 else (p0[0] - 1, p0[1])
    return (p0[0] - sh[0], p0[1] - sh[1]), (p1[0] - sh[0], p1[1] - sh[1])


CACHE = {}
def Gk(p0, p1, tv, dps):
    key = (canon(p0, p1), str(tv))
    if key not in CACHE:
        c0, c1 = canon(p0, p1)
        CACHE[key] = mp.re(kenyon(c0, c1, tv, dps))
    return CACHE[key]


def chain(a, b, n, tv, dps=22):
    w = sp.zeros(4, 1); w[a] += 1; w[b] -= X**(-n[0]) * Z**(-n[1])
    p = AH * w; r = Binv * w
    mp.mp.dps = dps + 10
    tvs = sp.Float(mp.nstr(tv, dps + 15), dps + 15)
    P = [{k: mp.mpf(str(sp.N(v.subs(t, tvs), dps + 12))) for k, v in laurent2(p[i]).items()} for i in range(4)]
    Rr = [{k: mp.mpf(str(sp.N(v.subs(t, tvs), dps + 12))) for k, v in laurent2(r[j]).items()} for j in range(4)]
    tot = mp.mpf(0)
    for i in range(4):
        for kk, cp in P[i].items():
            for j in range(4):
                for ll, cr in Rr[j].items():
                    m = (ll[0] - kk[0], ll[1] - kk[1])
                    pi_ = (TYPES[i][0] + 2 * m[0], TYPES[i][1] + 2 * m[1]); pj = TYPES[j]
                    tot += cp * cr * Gk(pj, pi_, tv, dps)       # G̃(i@m, j@0) = Kenyons G(v0 = j@0, v1 = i@m) (+ konstant)
    return tot


# förväntade koefficienter (1 ; √d ; 1/π ; √d/π ; Φ/π ; √dΦ/π) ur runs/isoradiell_lambda/RESULTAT.md
EXP = {
 '1/2': {'1:-1:0': '1/2;0;0;0;0;-1/5', '2:0:1': '0;0;0;0;0;8/5', '2:0:0': '0;0;1;0;0;4/5', '2:-1:0': '-1/2;0;0;0;0;18/5',
         '0:0:1': '5/2;0;-1;0;0;-26/5', '1:0:-1': '-7/4;0;5;0;0;14/5', '0:1:1': '1/2;0;-20;0;0;22', '0:-1:1': '0;0;-6;0;0;44/5',
         '2:0:2': '-5;0;21;0;0;-16/5'},
 '3': {'1:-1:0': '1/2;0;0;0;0;-1/10', '2:0:1': '0;0;0;0;0;2/15', '2:0:0': '0;0;1;0;0;1/15', '2:-1:0': '-3;0;0;0;0;32/15',
       '0:0:1': '5;0;-1;0;0;-13/5', '1:0:-1': '-31/2;0;10;0;0;79/10', '0:1:1': '3;0;-20/3;0;0;-2/9', '0:-1:1': '0;0;-8/3;0;0;38/45',
       '2:0:2': '-10;0;23/3;0;0;223/45'}}
UTANFOR = ['2:-2:2', '1:-2:-2', '0:-2:2']           # avstånd 4,11; 5,04; 5,46 (utanför de 53 noderna med avstånd < 4)
MED_INV = {'1:-1:0', '2:0:1', '0:1:1', '2:-2:2'}     # här räknas kedjan också med 1/t

ok = []
def check(name, cond):
    ok.append(bool(cond)); print(('OK  ' if cond else 'FEL ') + name, flush=True)

for ls in ['1/2', '3']:
    lv = sp.Rational(ls)
    mp.mp.dps = 50
    lm = mp.mpf(lv.p) / lv.q
    sd = mp.sqrt(lm * (lm + 2)); Ph = mp.acos(1 / (1 + lm)); tv = 1 + lm + sd
    basis = [mp.mpf(1), sd, 1 / mp.pi, sd / mp.pi, Ph / mp.pi, sd * Ph / mp.pi]
    Ms = bloch(lv, lambda x, y: 1)
    gr = Grund(Ms)
    print(f'=== λ = {ls}: t = {mp.nstr(tv, 15)}, Φ = arccos(1/(1+λ)) = {mp.nstr(Ph, 15)}, 2·atan(t) − π/2 = {mp.nstr(2*mp.atan(tv) - mp.pi/2, 15)}')
    for nd in list(EXP[ls]) + UTANFOR:
        j, n1, n2 = map(int, nd.split(':'))
        b, cell = zmap(j, (n1, n2)); a = SIG[0]
        t0 = time.time()
        Rg = gr.R(a, b, cell, dps=32)
        mp.mp.dps = 32
        line = f'{nd}: grundintegral {mp.nstr(Rg, 25)}'
        if nd in EXP[ls]:
            c = [Fr(x) for x in EXP[ls][nd].split(';')]
            Rs = mp.fsum(mp.mpf(x.numerator) / x.denominator * bb for x, bb in zip(c, basis))
            check(line + f' | Sats 3-värde {mp.nstr(Rs, 25)}', abs(Rg - Rs) < mp.mpf(10)**-27)
        else:
            mp.mp.dps = 32
            rel = mp.pslq([Rg, 1, 1 / mp.pi, sd * Ph / mp.pi], tol=mp.mpf(10)**-27, maxcoeff=10**7, maxsteps=10**6)
            if rel and rel[0]:
                cc = [Fr(-x, rel[0]) for x in rel[1:]]
                print(f'   PSLQ (1, 1/π, √dΦ/π): {cc}')
            rel6 = mp.pslq([Rg] + basis, tol=mp.mpf(10)**-27, maxcoeff=10**6, maxsteps=10**6)
            check(line + f' | PSLQ snäv bas: {rel}, 6-bas: {rel6}', rel is not None and rel[0] != 0)
            UT_VAL = globals().setdefault('UT_VAL', {}); UT_VAL[(ls, nd)] = (Rg, rel)
        Rc = chain(a, b, cell, tv)
        mp.mp.dps = 22
        check(f'   bevisets kedja (A, B, Kenyon) med t: {mp.nstr(Rc, 18)}', abs(Rc - Rg) < mp.mpf(10)**-17)
        if nd in MED_INV:
            Rc2 = chain(a, b, cell, 1 / tv)
            check(f'   … och med 1/t (annan A, B, annat Z²): {mp.nstr(Rc2, 18)}', abs(Rc2 - Rg) < mp.mpf(10)**-17)
        print(f'   ({time.time() - t0:.0f} s, {len(CACHE)} G-värden i cache)', flush=True)
print(f'\n{sum(ok)} godkända, {len(ok) - sum(ok)} fel')
