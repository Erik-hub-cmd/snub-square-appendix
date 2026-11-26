"""
Section 6, item 1 of "Exact two-point resistances on the snub square lattice":
direct evaluation of Eq. (1) for R_tt with 80-digit working precision (theta2 integral by residues, theta1 integral
by tanh-sinh quadrature, both in exact_snub.numeric_resistance), comparison with 2 sqrt(3)/9, and the PSLQ search that
found the closed form.
Requires: sympy, mpmath, exact_snub.py.   Run:  python3 high_precision_rtt.py
"""
import mpmath as mp
import sympy as sp
import exact_snub as es

DPS = 80
rtt = es.numeric_resistance(0, 2, (0, 1), dps=DPS)           # R_tt = R_13(0,1)
mp.mp.dps = DPS
exact = 2*mp.sqrt(3)/9
digits = int(-mp.log10(abs(rtt - exact)))
print("R_tt numeric :", mp.nstr(rtt, 60))
print("2*sqrt(3)/9  :", mp.nstr(exact, 60))
print("agreement    :", digits, "significant digits")
assert digits >= 44, "agreement below the 44 digits claimed in Section 6"

mp.mp.dps = digits - 4                                       # stay inside the digits that are correct
rel = mp.pslq([rtt, 1, mp.sqrt(3), 1/mp.pi, mp.sqrt(3)/mp.pi], maxcoeff=10**6, maxsteps=10**6)
print("PSLQ on (R_tt, 1, sqrt3, 1/pi, sqrt3/pi):", rel)
assert rel == [9, 0, -2, 0, 0] or rel == [-9, 0, 2, 0, 0], "PSLQ relation is not 9 R_tt - 2 sqrt3 = 0"
print("relation 9 R_tt - 2 sqrt(3) = 0 confirmed")
