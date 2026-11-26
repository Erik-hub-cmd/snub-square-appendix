# What each script checks

For every script: the statement of the article it supports, what the script computes, the output that counts as
a pass (the same patterns are in `expected.json` and are checked by `run_all.py` and by CI), and points that a
reader checking the claim should look at (marked **check**).
Numbering follows the current version of the article (version 13, November 2026). The list of ancillary files in
Appendix H of the article names the same scripts.

The scripts are aids to verification: every proof in the article is written out by hand, and a script checks
identities or values that the proof uses or that the text states. Script comments and output are in Swedish
(`godkända` = passed, `fel` = failed, `stämmer` = agrees). A failed check is printed on a line starting with `FEL`
or `FAIL`; `run_all.py` fails the run if any such line appears, if a script exits with a non-zero code, or if an
expected line is missing.

Layout: `kod/` holds the scripts and is the working directory for every run except those of Appendix G.4, which are
in `g4/kod/` (first program) and `g4/separat/` (second, separately written program); `runs/` holds stored numerical
results from earlier high-precision runs (read as data only, see the end of this file); `artikel/` holds the
matrices printed in Appendix A and the two data files.

---

## Theorem A and Appendix A (the snub square network and the isoradial square lattice $\mathbb Z^2_t$)

### `isoradiell_steg8_allmant_lambda.py`
- Claim: for every $\lambda>0$, with $t=1+\lambda+\sqrt{\lambda(\lambda+2)}$, there are $A,B$ with entries Laurent
  polynomials over $\mathbb Q(t)$ such that $M^\square_{\lambda(t)}A=BM_t$ and $\det A=\det B=4(t+1)^4(t^2+1)^2/(t-1)^4$;
  $A,B$ are unique with the normalisation $A_{11}=1$; the constants $\alpha,\beta$ of Eq. (alphabeta).
- Computes: solves the linear system for $A,B$ symbolically in $t$, expands the identity, the determinants, the rank
  at $t=2$ and $t=5/7$ (uniqueness), and $\alpha,\beta$.
- Pass: `9 godkända, 0 fel`.
- **check**: the edge list of the snub square network (function `lap`, `z2e`) is that of Section 2 / Figure 1.

### `isoradiell_steg5_laplacian.py`
- Claim: the case $\lambda=1$ ($t=q=2+\sqrt3$) of the identity over $\mathbb Q(\sqrt3)$, $\det A=\det B=576q^2$.
- Pass: `3 godkända, 0 fel`.

### `artikel1_kontroll.py`
- Claims: (T1) the matrices printed in Appendix A (`artikel/ab_matrices.tex`) are exactly those of the previous
  script; (T2) Lemma 2.2 / Proposition 2.3 (the cell numbering and the graph isomorphism to $\mathbb Z^2$ with
  diagonals); (T3) $\mathbb Z^2_t$ is isoradial with the half-angles and conductances stated; (T4) the partial
  fractions, factorisations and Euler substitutions of Appendix E for general $\lambda$; (T5) Table 2 at $\lambda=1$
  gives Table 1, and at $\lambda=1/2,2,3$ the stored high-precision values; (T6–T7) Theorem F from Kenyon's
  log-det formula, and $R_{TT}(\lambda)$ from Kenyon's $cR=2\theta/\pi$.
- Pass: `23 godkända, 0 fel`.
- **check**: T1 compares generated text with the file verbatim, so it fails if Appendix A is edited by hand.

### `granskning4_k1_identitet.py`
- Claim: as for `isoradiell_steg8_allmant_lambda.py`, by a separate construction of both Bloch matrices that does
  not use that script's `lap()`/`z2e()`; only $A,B$ are read from it as data. Besides the identity in $t$, it checks
  the identity at three fixed rational points, $t=2/5$ (below 1), $t=9/4$ and $t=7/3$ (above 1), as Appendix A says
  ("at rational points on both sides of $t=1$").
- Pass: `17 godkända, 0 fel`.
- Changed in this copy: the three points used to be drawn at random without a seed, so the number of checks was 16
  or 17 and both points could lie on the same side of $t=1$.

## Remark 3.5 (no electrical moves on the own torus and on small covers)

### `isoradiell_steg1_elektrisk.py`
- Claim: on its own torus, sequences of star–triangle moves with series/parallel reductions and removal of dead
  branches never reach a network of the form of $\mathbb Z^2$; the Bloch determinant changes by the expected factor
  in every move.
- Computes: all reachable graphs up to depth 7.
- Pass: `... in varje drag: True` and `... med Z²-form: 0`.

### `torus_sok.py`
- Claim: the exhaustive search of Remark 3.5 over all networks reachable by the moves on nine tori, including the
  $2\times1$ cover; on none of them is the minimal number of edges reached.
- Computes: breadth-first search over all reachable combinatorial maps up to isomorphism and gauge.
- Pass: `stämmer med satsen: JA` (no `NEJ`, no `AVBRUTEN`), for the seven tori of the default run, all searched
  exhaustively. "Satsen" is an internal formula: the minimum number of edges reachable is $[L:\Lambda](8+1/N_a+1/N_b)$,
  where $N_a,N_b$ are read off the cover. `JA` means that the search finds exactly this minimum. Since the formula
  is larger than $8[L:\Lambda]$, the number of edges of a minimal network, `JA` implies the statement of the remark.
  The output agrees with Remark 3.5: 10 edges (against 8) on the own torus, 18 (against 16) on the $2\times1$ cover.
  The nine tori, and the results for $2\times2$ and $4\times1$, go further than the remark, which names only the own
  torus and the $2\times1$ cover.
- Long runs (CI job `long-tori`, started by hand or monthly): `--fall 2x2 --kompakt` must give
  `326379 olika kartor (uttömmande)` and minimum 36 edges; `--fall 4x1 --kompakt` must give
  `2074204 olika kartor (uttömmande)` and minimum 34 edges. On the author's computer they took 15 min and 94 min.
- **check**: the moves implemented are star–triangle moves in both directions, series and parallel reductions and
  removal of dead branches; the search uses no splittings and no added pendant edges or loops (the inverse moves in
  (M) other than $\Delta\to Y$). This matches the remark, which speaks of reductions only.

## Corollary 3.2 and the second proof of Corollary C

### `uppfoljning_isoradiell.py`
- Claim: the snub square network with weight $\lambda$ and $\mathbb Z^2_t$ have the same Bloch determinant,
  identically in $X,Z$; the resulting formula for $R_{TT}(\lambda)$.
- Pass: `11 godkända, 0 fel`.

## Theorem 3.6 and Appendix G (no sequence of moves reaches a minimal network)

### `galois_sidor.py`
- Claim: Lemma G.1, the side polynomials of $\det M^\square_\lambda$ and of $\det M_t$.
- Computes: the four side polynomials of each determinant and their discriminants.
- Pass: snub square discriminant `16*lambda*(lambda + 1)**2*(lambda + 2)` (not a square in $\mathbb Q(\lambda)$),
  and the side polynomial of $\mathbb Z^2_t$ is `-(t**2 - w)*(t**2*w - 1)/t**2`.
- Reads the first part of `isoradiell_steg8_allmant_lambda.py` and executes it (to share the matrices).

### `galois_tackning.py` (uses `galois_lemma.py` as a library)
- Claim: Lemma G.2 on five tori (the own torus and four covers: $2\times1$, the diagonal cover $\langle a,b\rangle$,
  $3\times1$ and $2\times2$): the side polynomials of the snub square network at $\lambda=1$ have
  no rational roots ($w^2-14w+1$, $w^2-194w+1$, $w^2-2702w+1$ and their powers), while the square lattice on the
  same tori has only rational linear factors.
- Pass: the three quadratics appear, no line `rationella rötter: n` with $n>0$, no `False`.

### `galois_lemma.py`
- Library for `galois_tackning.py`, also runnable: (1) on $\mathbb Z^2$ with eight symbolic conductances every side
  polynomial splits into linear factors; (3) $w^2-14w+1$ is irreducible over $\mathbb Q$. Pass: `ALLT OK`.
- Changed in this copy: part (2), random walks from $\mathbb Z^2$, was removed. $\mathbb Z^2$ has no vertices of
  degree 3 and no triangles, so the walks made 0 moves and tested only the starting lattice. Lemma G.5 is checked by
  the two `lemma4_*` scripts below. The script is not cited in the article.

### `lemma4_test.py` (with `lemma4_natverk.py`)
- Claim: Lemma G.5, geodesic case: for minimal networks the Newton polygon is the zonogon of the train tracks and
  every side polynomial splits into linear factors over the field of the conductances.
- Computes: random geodesic minimal networks with rational conductances; the exact Bloch determinant; Newton
  polygon; factorisation of the side polynomials over $\mathbb Q$.
- Pass: `ALLT OK`; the default run tests 15 networks.
- **check**: the generated networks are minimal (tracks simple, crossing numbers $|h_i\wedge h_j|$).

### `lemma4_pl.py`
- Claim: Lemma G.5 for non-geodesic (piecewise linear) minimal networks; the roots equal the train-track
  monodromies (alternating conductance products) computed in the network itself.
- Pass: `ALLT OK (20 nät` for seed 31 and `ALLT OK (6 nät` for seed 41 with `--bara-ickegeo` (only
  combinatorially non-geodesic networks), and no line with `rötter=monodromi(i nätet självt):False`. The script's own
  verdict `ALLT OK` covers the Newton polygon, the splitting and the number of edges but not the comparison with the
  monodromies, so that comparison is checked by this pattern.
- Control: with `--icke-minimal` (seed 5) every network must be reported as `AVVIKELSE (4 nät`, which shows that
  the test can fail.
- **check**: the minimality test in the code matches the definition used in the article.

### `lemma5_yd.py`
- Claim: Remark G.6 (monodromies are invariant under star–triangle moves, and each side polynomial is a constant
  times $\prod_T(u+\varepsilon\mu_T^\sigma)$).
- Computes: hexagonal lattice on $k\times k$ tori ($k=2,3,4$) with random rational conductances, random moves.
- Pass: `fel 0` for each $k$; the default runs make 12, 21 and 40 moves.
- **check**: the script makes $Y\to\Delta$ moves only, while Remark G.6 is about star–triangle moves in both
  directions. Invariance under $\Delta\to Y$ follows because it is the inverse move, but it is not run.
- The numbers of networks and moves in Appendix H.1 (57, 448 with 54 non-geodesic, 410 moves) come from larger runs
  made while the lemma was being checked; the default runs here are smaller (15, 26 with 6 non-geodesic, 73 moves).

### `temperley_snub.py`
- Claim: the Temperley graph in Remark 3.8: the corner coefficients of the dimer determinant give side
  polynomials $w^2\pm14w+1$.
- Pass: lines `(2, 0) 1`, `(1, 1) 14`, `(0, 2) 1`.

## Appendix G.4 (the criterion for general networks; further lattices)

The first program is in `g4/kod/` (`natgalois.py`: train tracks, minimality test, exact Bloch determinant by
evaluation at roots of unity modulo large primes with Chinese remaindering under a Hadamard bound, Newton polygon,
side polynomials, the Galois test). The tiling data are `g4/kod/uniform12.json` (Archimedean and 2-uniform),
`uniform3.json` with `uniform3_saknad.json`, `uniform4.json` with `uniform4_saknad.json`, in coordinates over
$\mathbb Z[\zeta_{12}]$. The second program is in `g4/separat/` (`rev.py`: fraction-free elimination over
$\mathbb Z[z^{\pm1},w^{\pm1}]$; `arch.py`: its own geometric construction of the Archimedean lattices).

- `kor_arkimedes.py enhet` / `slump` (Proposition G.10): with unit and with random conductances, $3^4.6$,
  $3^2.4.3.4$ ($\Gamma$), $4.6.12$ and $4.8^2$ are obstructed (`GALOIS-HINDER: egen torus JA, alla torusar JA`), and
  $4.8^2$ has the side polynomial $u^2+6u+1$; the kagome lattice $3.6.3.6$ and $3.12^2$ are not minimal and not
  obstructed; the other five are minimal.
- `kor_kuniform.py` (six runs: 1-/2-, 3- and 4-uniform, unit and random conductances), then `sammanstall.py`: the
  counts of Appendix G.4 per $k$ (minimal / obstructed / neither): $k=1$ in the data 5/3/2 (plus $4.8^2$, obstructed),
  $k=2$ 6/13/1, $k=3$ 15/42/5, $k=4$ 30/106/14, and unit and random conductances give the same class everywhere.
  **check**: the 3-uniform data have 62 entries with one duplicate (entries 21 and 60, same $V$, $E$ and side
  polynomials), so the article counts $61=15+42+4$; one of the 151 four-uniform tilings is missing from the data, as
  the article says.
- `tabell.py`: writes `g4/TABELL_ALLA.md` (all 242 rows with their irreducible factors); CI checks that it is
  unchanged.
- `oberoende_det.py` (a check of the first program's determinants with SymPy, symbolically): agreement for eight
  tilings.
- `symb_klass.py`: side polynomials with one weight per edge class; for $4.8^2$ $-a(au^2+(2a+4)u+a)$ with
  discriminant $16a^2(a+1)$, as in the proof of Proposition G.10.
- `kor_slump_alla.py` (numerical search for reductions, Appendix G.4, last paragraph): every non-minimal tiling
  without obstruction reaches a minimal network, and no obstructed one does (`MOTSÄGELSE` would mark a
  contradiction with Theorem G.9). All three runs (1-/2-, 3- and 4-uniform) are part of the quick CI job
  (about 1, 4 and 12 minutes). **check**: a failed search proves nothing; the sequences of moves are not saved.
- Second program: `arch.py`, `overtorus.py` ($4.8^2$ on three covers), `symb.py` (weighted side polynomials),
  `kuni.py` (all 242 tilings from the same data files, with a check that every edge has length 1 and that the
  degrees match the vertex types). `jamfor_kuni.py` compares `kuni.py` with the first program tiling by tiling (same
  $V$, $E$ and verdict with both kinds of conductances). Pass: `242 godkända, 0 fel`.
- **check**: both programs read the same data files for the $k$-uniform tilings, so an error in the data would not
  be caught; for the Archimedean lattices the second program builds the geometry itself.
- **check**: the second program tests the side polynomials on the own torus only (`overtorus.py` adds three covers
  for $4.8^2$), while the first program's verdict is "obstructed on every torus". They agree here because no tiling is
  obstructed on the own torus only. The statement for all tori rests on the first program and on the argument of
  Theorem G.9. The second program computes no train tracks, so the class "minimal" comes from the first program
  (`natgalois.py`) alone. `jamfor_kuni.py` pairs the tilings by their order in the data files and checks that $V$ and
  $E$ agree.
- `kor_arkimedes.py enhet` also must report `minimalt=JA` for $3^6$, $4^4$, $6^3$, $3^3.4^2$, $3.4.6.4$ and
  `minimalt=nej` for the other six, and each `kor_kuniform.py` run must reach the last tiling of its data
  (index 29, 61 and 149).

### Appendix G.5: the Lean formalisation
- In the separate repository https://github.com/Erik-hub-cmd/lean-galois-snub, not here.

## Theorem B (closed forms for general $\lambda$)

### `isoradiell_steg9_kontroll159.py` (run first)
- Claim: the 159 stored values ($\lambda=1/2,2,3$, 53 node classes each) have the form
  $a+b/\pi+c\sqrt d\,\Phi/\pi$ with rational $a,b,c$.
- Computes: PSLQ against the full and the narrow basis; recomputes the value at 100 digits.
- Pass: `159 av 159 godkända, 0 fel`. Writes `runs/isoradiell_lambda/KONTROLL_159.md`.
- **check**: the stored values come from one program in two precisions; they protect against rounding errors, not
  systematic ones (the next script is the exact check).

### `isoradiell_steg9_kenyon_t.py`
- Claims: Lemma 4.4 (Kenyon's Green function on $\mathbb Z^2_t$ has the form $g_0+g_1/\pi+g_2\Phi/\pi$ with
  $g_i\in\mathbb Q(t)$; the branch of the logarithm; $\Delta G=\delta$ as an identity in $t$) and, for the 53 classes
  at distance less than 4, the symmetry under $t\mapsto1/t$ that gives Theorem B, with exact agreement with the
  stored values at $\lambda=1/2,2,3$ and with the table at $\lambda=1$.
- Pass: `9 godkända, 0 fel`. Writes `runs/isoradiell_lambda/SLUTNA_FORMER.md`.

### `isoradiell_steg6_kenyon.py`
- Claim: at $\lambda=1$ the resistances from Kenyon's formula and $A,B$ equal the exact table (100 nearest nodes);
  with `--sigma` the Galois-conjugate computation gives the same values.
- Pass: `100 av 100 exakt lika med tabellen; 100 av 100 utan √3/π-term`, with and without `--sigma`.

### `granskning4_egen.py` (library), `granskning4_k2_kenyon.py`, `granskning4_k3_stickprov.py`, `granskning4_k4_slutna.py`
- Written separately from the scripts above: their own Bloch matrices and a numerical contour integral for
  Kenyon's Green function, without residues.
- `granskning4_k2_kenyon.py`: Kenyon's formula on $\mathbb Z^2_t$ (conductances of the rhombic embedding,
  $\Delta G=\delta$, symmetry, $t=1$, resistances from the basic integral, $t>1$ and $t<1$). Pass: `22 godkända, 0 fel`.
- `granskning4_k3_stickprov.py`: at $\lambda=1/2$ and $3$, direct integration against the closed forms and the
  chain of the proof numerically; PSLQ for nodes outside the 53. Pass: `56 godkända, 0 fel`.
- `granskning4_k4_slutna.py`: closed forms for three nodes outside the 53, the symmetry under $t\mapsto1/t$, and at
  $\lambda=1$ exact agreement $a+c\sqrt3/3+b/\pi$ = the stored exact value. Pass: `3 godkända, 0 fel`, no
  `τ-symmetri False`. (The `OK (K1)`–`(K4)` lines in its output come from `isoradiell_steg9_kenyon_t.py`, which it
  imports; they are not about these three nodes.)
- Changed in this copy: the comparison with the table at $\lambda=1$ and the summary line were added; before, the
  script printed the values without comparing them.

## Theorem F (spanning-tree constant for general $\lambda$)

### `kenyon_logdet_kontroll.py`
- Claim: Kenyon's log-det formula for $\mathbb Z^2_t$ gives $m(\lambda)$ of Theorem F; numerical comparison.
- Pass: `10 av 10 godkända`.

### `bevis_viktad.py`
- Claim: the derivation of $R_{TT}(\lambda)$, $R_{TS}(\lambda)$, $m(\lambda)$ in Theorem F and Appendix E, every step
  exactly in $\lambda$ or numerically to at least 50 digits.
- Pass: `46 godkända, 0 fel`.

### `bevis_spannkonstant.py`
- Claim: $\lambda=1$: $m=\tfrac{40}{3}G/\pi+\tfrac43\log(2+\sqrt3)$.
- Pass: `19 godkända, 0 fel`.

## Appendix D and Theorem G (all resistances in $\mathbb Q+\mathbb Q\sqrt3+\mathbb Q/\pi$; Cairo lattice)

### `bevis_formodan8_direkt.py` (needs python-flint)
- Claim: the residue and Galois argument of Appendix D for every pair: the parametrisation, orders of $Z$ in the
  adjugate, poles only at the primitive 24th roots of unity and 0, $R=A-B/(2\pi)$ with $A\in\mathbb Q(\sqrt3)$,
  $B\in\mathbb Q$; exact comparison with all 1763 stored values.
- Pass: `SUMMARY: 10 passed, 0 failed`.
- **check**: the identities of rational functions are verified at (degree bound + 1) points; the degree bounds are
  derived in the script.

### `bevis_cairo.py` (needs python-flint)
- Claim: the same for the Cairo lattice (Theorem G, Lemma 6.1, Table 4), including Kirchhoff's law with the term
  $g_3-g_4=\sqrt3/18$.
- Pass: `7 godkända, 0 fel`.

### Libraries `k24.py`, `residue_tools.py` (need python-flint), `cairo_lattice.py`, `exact_snub.py`
- `k24.py`: exact arithmetic in $\mathbb Q(e^{i\pi/12})$. `residue_tools.py`: the residue machinery for any lattice
  with the same spectral curve. CI checks that they load.
- `cairo_lattice.py`: the Cairo lattice as the dual; pass: `det L_Cairo - det M_snub = 0` (Lemma 6.1).

## Appendix E (the first reduction; Table 4)

### `exact_snub.py`
- Claim: the reduction of Appendix E to integrals of rational functions, with exact integration and an
  mpmath cross-check of Eq. (1); the values of Table 1.
- Pass: $R_{13}(0,0)=\sqrt3/9+1/\pi$ and $R_{12}(-2,0)=-157\sqrt3/18-18/\pi+43/2$ appear, and in all 11 rows
  the exact value (printed to 25 digits) and the numerical value agree within $10^{-24}$ (checked by `run_all.py`).
- With argument `cairo`: Table 4; pass: `g3 - g4 = sqrt(3)/18`, $R_{14}(-1,0)=2\sqrt3/3$, and in all 18 rows
  $|\text{exact}-\text{numeric}|\le10^{-30}$.

### `snub_viktad.py`, `bevis_viktad.py`
- `snub_viktad.py`: the weighted determinant; pass: `λ = 1 ger 4 D(a, b): True`.

## Section 6.2 (sandpile and degree distribution)

### `sandhog_y.py` (needs python-flint)
- Claim: Proposition 6.4 and Table 5: the height-one probability and the degree distribution of the uniform
  spanning tree, exactly, for the Cairo lattice and the snub square lattice, with exact agreement with Goodchild's
  values.
- Pass: `11 godkända, 0 fel` (run without the optional argument `sim`).

## Theorem H (asymptotics)

### `asymptotik.py` (needs python-flint)
- Claim: Theorem H, the residues, the constant $\sqrt3/36$, harmonic positions, the $1/r$ and $1/r^2$ terms, and
  the error on all 1763 stored values.
- Pass: `SUMMARY: 18 passed, 0 failed`.

## Numerical checks (Appendix H)

### `high_precision_rtt.py`
- Direct evaluation of Eq. (1) for $R_{TT}$ at 80 digits and the PSLQ search that found $2\sqrt3/9$.
- Pass: `PSLQ ... [-9, 0, 2, 0, 0]`, `relation 9 R_tt - 2 sqrt(3) = 0 confirmed`, and agreement of at least 53
  significant digits.

### `second_implementation.py`
- A second implementation that uses neither the substitutions of Appendix E nor residues; compared with Table 1.
- Pass: `worst agreement:` at least 66 digits, as the article states.

### `numeric_check.py`
- Eq. (1) on a grid, Kirchhoff's equations on tori with Richardson extrapolation, Eq. (2) and the 1D integral,
  four Cairo entries.
- Pass: Eq. (2) and Step 4 both print `0.384900179459751`; every printed torus difference is at most $10^{-10}$
  (the article: "better than $10^{-10}$"; the largest is $6.4\cdot10^{-11}$); the grid values agree with the exact
  ones to the 10 printed digits ($2\cdot10^{-10}$), and the Cairo torus values within $3\cdot10^{-10}$
  (checked by `run_all.py`).

### `check_every_formula.py`
- 61 symbolic checks, written without importing other project code, of Lemma 2.1, Remark 3.3, the case $\lambda=1$
  of Appendix E, Table 1 and the Cairo lattice. Pass: `SUMMARY: 61 passed, 0 failed`.

### `granskning2_kontroll.py`
- Builds both lattices from the geometry without project code and computes resistances directly from Eq. (1);
  14 checks including 40 random stored values. Pass: `SUMMARY: 14 godkända, 0 fel` (about 7 minutes).

## Data

### `export_artikel1_data.py`
- Writes `artikel/ancillary/resistances_lambda1.csv` (the 1763 exact resistances from vertex 1) and
  `artikel/ancillary/closed_forms_lambda.csv` (the coefficients $a,b,c$ for the 53 classes).
- Pass: `1763 resistanser, 53 slutna former`, and CI checks with `git diff` that the regenerated files are
  identical to the committed ones (likewise the two files under `runs/isoradiell_lambda/` written by the Theorem B
  scripts).

## Stored data in `runs/`

These files are test data, not input that a proof relies on: each is compared exactly with a separate
derivation in the scripts named. The programs that wrote them (for example `jobb_isoradiell_lambda.py`, named in the
docstring of `isoradiell_steg9_kontroll159.py`) are not part of this repository.

- `runs/tabell/checkpoint.json`: the 1763 exact resistances at $\lambda=1$ with their high-precision numerical
  values, from an earlier computation (not rerun here; Appendix D's proof script compares all of them exactly with
  the residue formula, and `granskning2_kontroll.py` recomputes 40 of them from Eq. (1)).
- `runs/isoradiell_lambda/checkpoint.json`, `koefficienter_159.json`: values at $\lambda=1/2,2,3$ and their PSLQ
  coefficients (checked exactly by `isoradiell_steg9_kenyon_t.py`).
- `runs/cairo/tabell.json`: the Cairo table from an earlier computation (compared exactly by `bevis_cairo.py`, C4;
  `exact_snub.py cairo` computes the same values but neither writes nor reads this file).
