r"""Separat kontroll 1-2 av Sats A: egen konstruktion av M_snub(λ) och M_Z²(t) (utan projektets lap()/z2e()),
sedan M_snub·A = B·M_Z², det A, det B, Laurent-egenskap, nämnare i t, nollmoden.
A, B läses (som data) ur kod/isoradiell_steg8_allmant_lambda.py och kod/isoradiell_steg5_laplacian.py.
Se GRANSKNING4.md, avsnitt 1–2. Kör: .venv\Scripts\python kod\granskning4_k1_identitet.py"""
import sys, os, random
sys.stdout.reconfigure(encoding="utf-8")
import sympy as sp

ROOT = os.path.dirname(os.path.abspath(__file__))
X, Z, t = sp.symbols('X Z t')
TYPES = [(0, 0), (0, 1), (1, 0), (1, 1)]          # projektets ordning: index 2*(x mod 2) + (y mod 2)


def typ(p):
    return TYPES.index((p[0] % 2, p[1] % 2)), (p[0] // 2, p[1] // 2)   # (typ, cell)


def has_diag(p, r):
    """Diagonal i rutan med nedre vänstra hörnet c: huvuddiagonal om c ≡ (0,0), bidiagonal om c ≡ (1,1) (mod 2)."""
    for c in [(min(p[0], r[0]), min(p[1], r[1]))]:
        if abs(p[0] - r[0]) != 1 or abs(p[1] - r[1]) != 1: return False
        main = (p[0] - r[0]) == (p[1] - r[1])
        cm = (c[0] % 2, c[1] % 2)
        return (cm == (0, 0) and main) or (cm == (1, 1) and not main)


def neighbours(p, lam, h):
    out = []
    for d in [(1, 0), (-1, 0), (0, 1), (0, -1)]:
        r = (p[0] + d[0], p[1] + d[1])
        if d[1] == 0: x = min(p[0], r[0]); w = h(x, p[1])
        else: w = 1
        out.append((r, w))
    if lam is not None:
        for d in [(1, 1), (-1, -1), (1, -1), (-1, 1)]:
            r = (p[0] + d[0], p[1] + d[1])
            if has_diag(p, r): out.append((r, lam))
    return out


def bloch(lam, h):
    """(Mf)_a = Σ_nb w (f_a − X^{cx} Z^{cy} f_b), cell (cx, cy) för grannen; samma konvention som projektets lap()."""
    M = sp.zeros(4, 4)
    for a, p in enumerate(TYPES):
        for r, w in neighbours(p, lam, h):
            b, c = typ(r)
            M[a, a] += w; M[a, b] -= w * X**c[0] * Z**c[1]
    return M


ok = []
def check(name, cond):
    ok.append(bool(cond)); print(('OK  ' if cond else 'FEL ') + name, flush=True)

# --- 0. gittret är snub square: hörnfigur 3.3.4.3.4 i alla fyra nodtyper, en diagonal per nod
def vertex_figure(p):
    # rutorna runt p i cyklisk ordning (NE, NW, SW, SE); en ruta med diagonal genom p ger 2 trianglar, annars 1 triangel/1 kvadrat
    seq = []
    for c in [(p[0], p[1]), (p[0] - 1, p[1]), (p[0] - 1, p[1] - 1), (p[0], p[1] - 1)]:
        cm = (c[0] % 2, c[1] % 2)
        if cm == (0, 0) or cm == (1, 1):
            corners = [(c[0], c[1]), (c[0] + 1, c[1] + 1)] if cm == (0, 0) else [(c[0] + 1, c[1]), (c[0], c[1] + 1)]
            seq += ['3', '3'] if p in corners else ['3']
        else:
            seq.append('4')
    return seq
figs = [vertex_figure(p) for p in TYPES]
def cyc_eq(s, target):
    n = len(s); return any((s[i:] + s[:i]) == target or (s[i:] + s[:i])[::-1] == target for i in range(n))
check('hörnfigur 3.3.4.3.4 i alla fyra typer: ' + str(figs), all(cyc_eq(f, list('33434')) for f in figs))
lam_t = (t - 1)**2 / (2 * t)
Ms = bloch(lam_t, lambda x, y: 1)
Mz = bloch(None, lambda x, y: 1 / t if (x + y) % 2 == 0 else t)
check('grad 5 i snub (diagonalsumma = 4 + λ)', all(sp.simplify(Ms.subs({X: 1, Z: 1})[a, a] - 4 - lam_t) == 0 for a in range(4)))
check('radsummor 0 i X = Z = 1 (Laplacian)', Ms.subs({X: 1, Z: 1}) * sp.ones(4, 1) == sp.zeros(4, 1) and
      sp.simplify(Mz.subs({X: 1, Z: 1}) * sp.ones(4, 1)) == sp.zeros(4, 1))
check('hermitesk: M(1/X,1/Z)^T = M', sp.simplify(Ms.subs({X: 1/X, Z: 1/Z}, simultaneous=True).T - Ms) == sp.zeros(4, 4))

# --- 1. läs A, B (steg 8) som data
src = open(os.path.join(ROOT, 'isoradiell_steg8_allmant_lambda.py'), encoding='utf-8').read()
ns = {'sp': sp, 'X': X, 'Z': Z, 't': t}
exec(src[src.index('A = sp.Matrix'):src.index("\ncheck('")], ns)
A, B = ns['A'], ns['B']
E = (Ms * A - B * Mz).applyfunc(lambda e: sp.cancel(sp.together(e)))
check('M_snub·A = B·M_Z² identiskt i X, Z, t (egen konstruktion av båda Laplacianerna)', E == sp.zeros(4, 4))
# fasta rationella punkter som extra kontroll, på båda sidor om t = 1 (förr slumpvisa utan frö)
for sub in ({X: sp.Rational(3, 7), Z: sp.Rational(11, 5), t: sp.Rational(2, 5)},
            {X: sp.Rational(13, 4), Z: sp.Rational(2, 9), t: sp.Rational(9, 4)},
            {X: sp.Rational(5, 17), Z: sp.Rational(19, 6), t: sp.Rational(7, 3)}):
    check(f'  … och numeriskt exakt i {sub}', (Ms * A - B * Mz).subs(sub) == sp.zeros(4, 4))
D = 4 * (t + 1)**4 * (t**2 + 1)**2 / (t - 1)**4
dA, dB = sp.factor(sp.cancel(A.det())), sp.factor(sp.cancel(B.det()))
check(f'det A = {dA} (konstant i X, Z)', sp.cancel(dA - D) == 0)
check(f'det B = {dB} (konstant i X, Z)', sp.cancel(dB - D) == 0)
check('D ≠ 0 för alla t > 0, t ≠ 1 (nollställen bara i t = −1, ±i)', set(sp.solve(sp.numer(sp.together(D)), t)) <= {-1, sp.I, -sp.I})

def is_laurent_S(M):
    good = True
    for e in M:
        n, d = sp.fraction(sp.factor(sp.cancel(sp.together(e))))
        for f, _ in sp.factor_list(d)[1]:
            fx = f.free_symbols
            if fx & {X, Z}:
                if not (f in (X, Z)): good = False
            elif t in fx:
                if sp.expand(f - t) != 0 and sp.expand(f - t + 1) != 0 and sp.expand(f - t - 1) != 0 and sp.expand(f - t**2 - 1) != 0:
                    good = False; print('   otillåten faktor', f)
    return good
Binv = (B.adjugate() / D)
Ainv = (A.adjugate() / D)
check('A, B, A⁻¹, B⁻¹ har Laurent-element i X, Z med nämnare bara t, t−1, t+1, t²+1 (dvs. i S)',
      all(is_laurent_S(M) for M in (A, B, Ainv, Binv)))
check('M_snub⁻¹ = A·M_Z²⁻¹·B⁻¹ (kontroll i en rationell punkt)',
      sp.simplify((A * Mz.inv() * Binv - Ms.inv()).subs({X: sp.Rational(2, 3), Z: sp.Rational(5, 7), t: 3})) == sp.zeros(4, 4))

# --- 2. nollmoden: 1ᵀ B(1,1) ∝ 1ᵀ och A(1,1)·1 ∝ 1, så att Σ p_i(1,1) = Σ r_j(1,1) = 0
one = sp.ones(4, 1)
B1 = (one.T * B.subs({X: 1, Z: 1})).applyfunc(sp.cancel)
A1 = (A.subs({X: 1, Z: 1}) * one).applyfunc(sp.cancel)
print('   1ᵀB(1,1) =', list(B1), '\n   A(1,1)·1 =', list(A1))
check('1ᵀ B(1,1) = c·1ᵀ och A(1,1)·1 = c′·1 med c, c′ ≠ 0', len(set(B1)) == 1 and len(set(A1)) == 1 and B1[0] != 0 and A1[0] != 0)

# --- 3. λ = 1 (steg 5)
src5 = open(os.path.join(ROOT, 'isoradiell_steg5_laplacian.py'), encoding='utf-8').read()
ns5 = {'sp': sp, 'sqrt': sp.sqrt, 'X': X, 'Z': Z}
exec(src5[src5.index('A = sp.Matrix'):src5.index("\ncheck('")], ns5)
A5, B5 = ns5['A'], ns5['B']; q = 2 + sp.sqrt(3)
Ms1, Mz1 = Ms.subs(t, q), Mz.subs(t, q)
E5 = (Ms1 * A5 - B5 * Mz1).applyfunc(lambda e: sp.radsimp(sp.cancel(sp.expand(sp.together(e) * X**3 * Z**3))))
check('λ = 1 (steg 5): M_snub·A = B·M_Z² med egen konstruktion', all(sp.simplify(e) == 0 for e in E5))
for nm, M in (('A', A5), ('B', B5)):
    d = sp.simplify(sp.radsimp(sp.cancel(sp.together(M.det()))))
    check(f'λ = 1: det {nm} = {d} = 576 q²', sp.simplify(d - 576 * q**2) == 0)
d8 = sp.radsimp(D.subs(t, q))
print('   steg 8:s determinant i t = 2+√3:', sp.simplify(d8), ' (steg 5 använder en annan normering av A, B)')
print(f'\n{sum(ok)} godkända, {len(ok) - sum(ok)} fel')
