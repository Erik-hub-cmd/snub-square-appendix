"""Kontroll (2026-10-19, disposition artikel 1): spännträdskonstanten m(λ) ur Kenyons log det-formel.

Det isoradiella Z²_t (cell 2×2: 4 vertikala kanter med halvvinkel π/4, 2 horisontella med π/4 + Φ/2 och 2 med π/4 − Φ/2)
har samma Bloch-determinant som snub square med vikt λ. Kenyons formel för isoradiella grafer,
log det₁ Δ = Σ_e (2/π)(θ log tan θ + L(θ) + L(π/2 − θ)), L(x) = ½ Cl₂(2x) (Lobachevskij),
ska därför ge m(λ) = ∫∫ log det M_λ dθ/(2π)², dvs. korollariet i artikelns avsnitt 7.
Jämförs här med korollariets formel och med direkt numerisk integration. Kräver numpy, scipy, mpmath.
Kör: python kenyon_logdet_kontroll.py
"""
import numpy as np
from scipy import integrate
from mpmath import clsin, mp, catalan, log, acos, sqrt, pi, tan

mp.dps = 30
ok = []


def L(x):
    return clsin(2, 2*x)/2


def f(th):
    return (2/pi)*(th*log(tan(th)) + L(th) + L(pi/2 - th))


for lam in [0.25, 0.5, 1.0, 3.0, 10.0]:
    lam = mp.mpf(lam); kap = sqrt(lam*(lam + 2)); Phi = acos(1/(1 + lam)); t = 1 + lam + kap
    kenyon = 4*f(pi/4) + 2*f(pi/4 + Phi/2) + 2*f(pi/4 - Phi/2)
    korollarium = 8*catalan/pi + 4/pi*(Phi*log(t) + clsin(2, pi/2 + Phi) + clsin(2, pi/2 - Phi))
    lf = float(lam)
    D = lambda a, b: a*a + b*b - (4*lf**2 + 8*lf + 2)*a*b - 4*(lf + 2)**2*(a + b) + 4*(3*lf + 4)*(lf + 2)
    val, _ = integrate.dblquad(lambda x, y: np.log(4*D(np.cos(x), np.cos(y))), 0, np.pi, 0, np.pi,
                               epsabs=1e-12, epsrel=1e-12)
    direkt = val/np.pi**2
    c1 = abs(kenyon - korollarium) < mp.mpf(10)**-25
    c2 = abs(float(kenyon) - direkt) < 1e-10
    ok += [c1, c2]
    print(f"λ = {lf}: Kenyon {float(kenyon):.15f}  korollarium {float(korollarium):.15f}  direkt {direkt:.15f}  "
          + ("OK" if c1 and c2 else "FEL"))
print(f"{sum(ok)} av {len(ok)} godkända")
