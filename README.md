# snub-square-appendix

Ancillary scripts and data for E. Ersson, *Exact two-point resistances on the snub square lattice via a hidden
isoradial square lattice* (2026). Appendix H of the article lists which script checks which statement;
`STATEMENTS.md` gives the same list with the output that counts as a pass for every script.

    pip install -r requirements.txt
    python run_all.py --skip-long     # every script except two long exhaustive searches (about 1 hour)
    python run_all.py --only-long     # the torus searches on the 2x2 and 4x1 covers (about 2 hours)
    python run_all.py lemma4_test.py  # a single script (by file name or by the tag in expected.json)

The entries for Appendix G.4 (tags beginning with `g4_`) must run in the order of `expected.json`: `sammanstall.py`,
`tabell.py` and `jamfor_kuni.py` read files written by the entries before them. Run them with `--skip-long`, or name
all of them in one call.

`run_all.py` runs each script from `kod/`, writes its output to `ci_logs/`, and compares it with the patterns in
`expected.json`; it exits with code 1 if a script fails, prints a failed check, or misses an expected line.
GitHub Actions runs `run_all.py --skip-long` on every push, checks that the regenerated data files are identical
to the committed ones, and runs the long searches when started by hand and once a month.

- `kod/` – the scripts (comments and output in Swedish: `godkända` = passed, `fel` = failed).
- `artikel/ab_matrices.tex` – the matrices $A, B$ printed in Appendix A (checked verbatim by `artikel1_kontroll.py`).
- `artikel/ancillary/` – `resistances_lambda1.csv` (the 1763 exact resistances from vertex 1) and
  `closed_forms_lambda.csv` (the coefficients $a, b, c$ for the 53 classes at distance less than 4).
- `runs/` – stored results of earlier high-precision computations, read as data (see `STATEMENTS.md`).
- `g4/` – the computations of Appendix G.4: `g4/kod/` (first program and the tiling data), `g4/separat/` (second
  program, written separately, and the comparison `jamfor_kuni.py`), `g4/TABELL_ALLA.md` (the full table, regenerated
  by `tabell.py`), `g4/loggar/` (intermediate files written during a run).

## Notes on the scripts

Files whose names begin with `granskning` ('review') and the scripts described as separately written were written
as internal checks in separate sessions with AI assistance, without reading the code they check. They are not an
independent review. Where the article says "independent" about code, it means code written separately that does not
import the other scripts.

Comments in the scripts sometimes name working notes (`forslag/*.md`, `GRANSKNING*.md`, `ekvivalens_teori.md`) or
scripts (for example `jobb_isoradiell_lambda.py`) that are not published, and some give run commands for the
author's Windows setup (`.venv\Scripts\python kod\...`). Here every script is run as `python <script>` from its
folder, as `run_all.py` does.

These copies differ from the author's originals as follows (each change is described in `STATEMENTS.md`):
`galois_sidor.py` and `isoradiell_steg6_kenyon.py` open files with `encoding='utf-8'`; `granskning4_k1_identitet.py`
uses three fixed test points instead of unseeded random ones; `granskning4_k4_slutna.py` compares its values with the
stored table and prints a summary; part (2) of `galois_lemma.py` is removed; and some docstrings were reworded
("separately written" instead of "independent", article numbering).

Tested with Python 3.13 and the versions in `requirements.txt` (on Windows 11 and Linux).
The Lean formalisation of Theorem 3.6 is a separate repository: https://github.com/Erik-hub-cmd/lean-galois-snub

Licence: Apache-2.0.
