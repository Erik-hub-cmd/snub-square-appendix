# Slutna former för allmänt λ (steg 9.2)

R(nod 1 → nod) = a(λ) + b(λ)/π + c(λ)·κ·Φ/π,  κ = √(λ(λ+2)),  Φ = arccos(1/(1+λ)) = arctan κ.
Räknat exakt ur Kenyons formel för det isoradiella Z² och A(t), B(t) från steg 8 (kod/isoradiell_steg9_kenyon_t.py).
Kolumnerna "λ = 1/2, 2, 3" och "λ = 1" säger om värdena stämmer exakt med PSLQ (KONTROLL_159.md) och runs/tabell.

| nod j:n1:n2 | a(λ) | b(λ) | c(λ) | λ = 1/2, 2, 3 | λ = 1 |
|---|---|---|---|---|---|
| 0:-2:0 | $2 \left(\lambda + 2\right) \left(5 \lambda^{2} + 14 \lambda + 10\right)$ | $- \frac{2 \left(10 \lambda^{3} + 38 \lambda^{2} + 39 \lambda + 5\right)}{\lambda}$ | $- \frac{2 \left(10 \lambda^{5} + 58 \lambda^{4} + 119 \lambda^{3} + 97 \lambda^{2} + 19 \lambda - 5\right)}{\lambda^{2} \left(\lambda + 2\right)}$ | ja | ja |
| 0:-1:-1 | $\lambda$ | $- \frac{4 \left(\lambda + 2\right)}{\lambda}$ | $- \frac{2 \left(\lambda^{2} - 2 \lambda - 2\right)}{\lambda^{2}}$ | ja | ja |
| 0:-1:0 | $\lambda + 2$ | $-1$ | $- \frac{2 \lambda^{2} + 6 \lambda + 3}{\lambda \left(\lambda + 2\right)}$ | ja | ja |
| 0:-1:1 | $0$ | $- \frac{2 \left(\lambda + 1\right)}{\lambda}$ | $\frac{2 \left(\lambda^{2} + 3 \lambda + 1\right)}{\lambda^{2} \left(\lambda + 2\right)}$ | ja | ja |
| 0:0:-2 | $2 \left(\lambda + 2\right) \left(5 \lambda^{2} + 14 \lambda + 10\right)$ | $- \frac{2 \left(10 \lambda^{3} + 38 \lambda^{2} + 39 \lambda + 5\right)}{\lambda}$ | $- \frac{2 \left(10 \lambda^{5} + 58 \lambda^{4} + 119 \lambda^{3} + 97 \lambda^{2} + 19 \lambda - 5\right)}{\lambda^{2} \left(\lambda + 2\right)}$ | ja | ja |
| 0:0:-1 | $\lambda + 2$ | $-1$ | $- \frac{2 \lambda^{2} + 6 \lambda + 3}{\lambda \left(\lambda + 2\right)}$ | ja | ja |
| 0:0:1 | $\lambda + 2$ | $-1$ | $- \frac{2 \lambda^{2} + 6 \lambda + 3}{\lambda \left(\lambda + 2\right)}$ | ja | ja |
| 0:0:2 | $2 \left(\lambda + 2\right) \left(5 \lambda^{2} + 14 \lambda + 10\right)$ | $- \frac{2 \left(10 \lambda^{3} + 38 \lambda^{2} + 39 \lambda + 5\right)}{\lambda}$ | $- \frac{2 \left(10 \lambda^{5} + 58 \lambda^{4} + 119 \lambda^{3} + 97 \lambda^{2} + 19 \lambda - 5\right)}{\lambda^{2} \left(\lambda + 2\right)}$ | ja | ja |
| 0:1:-1 | $0$ | $- \frac{2 \left(\lambda + 1\right)}{\lambda}$ | $\frac{2 \left(\lambda^{2} + 3 \lambda + 1\right)}{\lambda^{2} \left(\lambda + 2\right)}$ | ja | ja |
| 0:1:0 | $\lambda + 2$ | $-1$ | $- \frac{2 \lambda^{2} + 6 \lambda + 3}{\lambda \left(\lambda + 2\right)}$ | ja | ja |
| 0:1:1 | $\lambda$ | $- \frac{4 \left(\lambda + 2\right)}{\lambda}$ | $- \frac{2 \left(\lambda^{2} - 2 \lambda - 2\right)}{\lambda^{2}}$ | ja | ja |
| 0:2:0 | $2 \left(\lambda + 2\right) \left(5 \lambda^{2} + 14 \lambda + 10\right)$ | $- \frac{2 \left(10 \lambda^{3} + 38 \lambda^{2} + 39 \lambda + 5\right)}{\lambda}$ | $- \frac{2 \left(10 \lambda^{5} + 58 \lambda^{4} + 119 \lambda^{3} + 97 \lambda^{2} + 19 \lambda - 5\right)}{\lambda^{2} \left(\lambda + 2\right)}$ | ja | ja |
| 1:-2:-1 | $\frac{4 \lambda^{2} + 8 \lambda + 1}{2}$ | $- \frac{4 \left(\lambda + 2\right)^{2}}{\lambda}$ | $- \frac{8 \lambda^{4} + 24 \lambda^{3} - 3 \lambda^{2} - 56 \lambda - 32}{2 \lambda^{2} \left(\lambda + 2\right)}$ | ja | ja |
| 1:-2:0 | $\frac{6 \lambda^{2} + 20 \lambda + 17}{2}$ | $- 6 \left(\lambda + 2\right)$ | $- \frac{12 \lambda^{3} + 52 \lambda^{2} + 69 \lambda + 24}{2 \lambda \left(\lambda + 2\right)}$ | ja | ja |
| 1:-2:1 | $\frac{2 \lambda^{2} + 4 \lambda + 1}{2}$ | $- \frac{2 \left(\lambda - 2\right) \left(\lambda + 2\right)}{\lambda}$ | $- \frac{4 \lambda^{4} + 12 \lambda^{3} + 9 \lambda^{2} + 8 \lambda + 16}{2 \lambda^{2} \left(\lambda + 2\right)}$ | ja | ja |
| 1:-1:-1 | $- \frac{2 \lambda^{2} + 4 \lambda + 1}{2}$ | $2 \left(\lambda + 2\right)$ | $\frac{4 \lambda^{2} + 12 \lambda + 7}{2 \left(\lambda + 2\right)}$ | ja | ja |
| 1:-1:0 | $\frac{1}{2}$ | $0$ | $- \frac{1}{2 \left(\lambda + 2\right)}$ | ja | ja |
| 1:-1:1 | $- \frac{1}{2}$ | $0$ | $\frac{3 \lambda + 8}{2 \lambda \left(\lambda + 2\right)}$ | ja | ja |
| 1:-1:2 | $- \frac{14 \lambda^{2} + 52 \lambda + 49}{2}$ | $\frac{2 \left(\lambda + 2\right) \left(7 \lambda + 6\right)}{\lambda}$ | $\frac{28 \lambda^{4} + 132 \lambda^{3} + 187 \lambda^{2} + 48 \lambda - 48}{2 \lambda^{2} \left(\lambda + 2\right)}$ | ja | ja |
| 1:0:-1 | $- \frac{2 \lambda^{2} + 4 \lambda + 1}{2}$ | $2 \left(\lambda + 2\right)$ | $\frac{4 \lambda^{2} + 12 \lambda + 7}{2 \left(\lambda + 2\right)}$ | ja | ja |
| 1:0:0 | $\frac{1}{2}$ | $0$ | $- \frac{1}{2 \left(\lambda + 2\right)}$ | ja | ja |
| 1:0:1 | $- \frac{1}{2}$ | $0$ | $\frac{3 \lambda + 8}{2 \lambda \left(\lambda + 2\right)}$ | ja | ja |
| 1:0:2 | $- \frac{14 \lambda^{2} + 52 \lambda + 49}{2}$ | $\frac{2 \left(\lambda + 2\right) \left(7 \lambda + 6\right)}{\lambda}$ | $\frac{28 \lambda^{4} + 132 \lambda^{3} + 187 \lambda^{2} + 48 \lambda - 48}{2 \lambda^{2} \left(\lambda + 2\right)}$ | ja | ja |
| 1:1:-1 | $\frac{4 \lambda^{2} + 8 \lambda + 1}{2}$ | $- \frac{4 \left(\lambda + 2\right)^{2}}{\lambda}$ | $- \frac{8 \lambda^{4} + 24 \lambda^{3} - 3 \lambda^{2} - 56 \lambda - 32}{2 \lambda^{2} \left(\lambda + 2\right)}$ | ja | ja |
| 1:1:0 | $\frac{6 \lambda^{2} + 20 \lambda + 17}{2}$ | $- 6 \left(\lambda + 2\right)$ | $- \frac{12 \lambda^{3} + 52 \lambda^{2} + 69 \lambda + 24}{2 \lambda \left(\lambda + 2\right)}$ | ja | ja |
| 1:1:1 | $\frac{2 \lambda^{2} + 4 \lambda + 1}{2}$ | $- \frac{2 \left(\lambda - 2\right) \left(\lambda + 2\right)}{\lambda}$ | $- \frac{4 \lambda^{4} + 12 \lambda^{3} + 9 \lambda^{2} + 8 \lambda + 16}{2 \lambda^{2} \left(\lambda + 2\right)}$ | ja | ja |
| 2:-2:0 | $- 2 \left(\lambda + 2\right) \left(2 \lambda^{2} + 4 \lambda + 1\right)$ | $\frac{8 \lambda^{3} + 24 \lambda^{2} + 21 \lambda + 8}{\lambda}$ | $\frac{8 \lambda^{5} + 40 \lambda^{4} + 64 \lambda^{3} + 28 \lambda^{2} - 11 \lambda - 8}{\lambda^{2} \left(\lambda + 2\right)}$ | ja | ja |
| 2:-2:1 | $- 2 \left(\lambda + 1\right)^{2} \left(\lambda + 2\right)$ | $\frac{2 \left(\lambda + 1\right) \left(2 \lambda^{2} + 4 \lambda + 1\right)}{\lambda}$ | $\frac{2 \left(2 \lambda^{5} + 10 \lambda^{4} + 17 \lambda^{3} + 11 \lambda^{2} + 2 \lambda - 1\right)}{\lambda^{2} \left(\lambda + 2\right)}$ | ja | ja |
| 2:-1:-1 | $- 2 \left(\lambda + 2\right) \left(2 \lambda^{2} + 4 \lambda + 1\right)$ | $\frac{8 \lambda^{3} + 24 \lambda^{2} + 21 \lambda + 8}{\lambda}$ | $\frac{8 \lambda^{5} + 40 \lambda^{4} + 64 \lambda^{3} + 28 \lambda^{2} - 11 \lambda - 8}{\lambda^{2} \left(\lambda + 2\right)}$ | ja | ja |
| 2:-1:0 | $- \lambda$ | $0$ | $\frac{2 \left(\lambda + 1\right)^{2}}{\lambda \left(\lambda + 2\right)}$ | ja | ja |
| 2:-1:1 | $0$ | $1$ | $\frac{1}{\lambda \left(\lambda + 2\right)}$ | ja | ja |
| 2:-1:2 | $- 2 \left(\lambda + 2\right)$ | $\frac{2 \left(\lambda + 1\right)}{\lambda}$ | $\frac{2 \left(2 \lambda^{3} + 7 \lambda^{2} + 6 \lambda - 1\right)}{\lambda^{2} \left(\lambda + 2\right)}$ | ja | ja |
| 2:0:-1 | $- 2 \left(\lambda + 1\right)^{2} \left(\lambda + 2\right)$ | $\frac{2 \left(\lambda + 1\right) \left(2 \lambda^{2} + 4 \lambda + 1\right)}{\lambda}$ | $\frac{2 \left(2 \lambda^{5} + 10 \lambda^{4} + 17 \lambda^{3} + 11 \lambda^{2} + 2 \lambda - 1\right)}{\lambda^{2} \left(\lambda + 2\right)}$ | ja | ja |
| 2:0:0 | $0$ | $1$ | $\frac{1}{\lambda \left(\lambda + 2\right)}$ | ja | ja |
| 2:0:1 | $0$ | $0$ | $\frac{2}{\lambda \left(\lambda + 2\right)}$ | ja | ja |
| 2:0:2 | $- 2 \left(\lambda + 2\right)$ | $\frac{5 \lambda + 8}{\lambda}$ | $\frac{4 \lambda^{3} + 12 \lambda^{2} + 5 \lambda - 8}{\lambda^{2} \left(\lambda + 2\right)}$ | ja | ja |
| 2:1:0 | $- 2 \left(\lambda + 2\right)$ | $\frac{2 \left(\lambda + 1\right)}{\lambda}$ | $\frac{2 \left(2 \lambda^{3} + 7 \lambda^{2} + 6 \lambda - 1\right)}{\lambda^{2} \left(\lambda + 2\right)}$ | ja | ja |
| 2:1:1 | $- 2 \left(\lambda + 2\right)$ | $\frac{5 \lambda + 8}{\lambda}$ | $\frac{4 \lambda^{3} + 12 \lambda^{2} + 5 \lambda - 8}{\lambda^{2} \left(\lambda + 2\right)}$ | ja | ja |
| 2:1:2 | $- \lambda$ | $- \frac{24 \left(\lambda + 2\right)}{\lambda^{2}}$ | $\frac{2 \left(\lambda^{4} + 2 \lambda^{3} + 9 \lambda^{2} + 28 \lambda + 24\right)}{\lambda^{3} \left(\lambda + 2\right)}$ | ja | ja |
| 3:-1:-1 | $\frac{4 \lambda^{2} + 8 \lambda + 1}{2}$ | $- \frac{4 \left(\lambda + 2\right)^{2}}{\lambda}$ | $- \frac{8 \lambda^{4} + 24 \lambda^{3} - 3 \lambda^{2} - 56 \lambda - 32}{2 \lambda^{2} \left(\lambda + 2\right)}$ | ja | ja |
| 3:-1:0 | $- \frac{2 \lambda^{2} + 4 \lambda + 1}{2}$ | $2 \left(\lambda + 2\right)$ | $\frac{4 \lambda^{2} + 12 \lambda + 7}{2 \left(\lambda + 2\right)}$ | ja | ja |
| 3:-1:1 | $- \frac{2 \lambda^{2} + 4 \lambda + 1}{2}$ | $2 \left(\lambda + 2\right)$ | $\frac{4 \lambda^{2} + 12 \lambda + 7}{2 \left(\lambda + 2\right)}$ | ja | ja |
| 3:-1:2 | $\frac{4 \lambda^{2} + 8 \lambda + 1}{2}$ | $- \frac{4 \left(\lambda + 2\right)^{2}}{\lambda}$ | $- \frac{8 \lambda^{4} + 24 \lambda^{3} - 3 \lambda^{2} - 56 \lambda - 32}{2 \lambda^{2} \left(\lambda + 2\right)}$ | ja | ja |
| 3:0:-1 | $\frac{6 \lambda^{2} + 20 \lambda + 17}{2}$ | $- 6 \left(\lambda + 2\right)$ | $- \frac{12 \lambda^{3} + 52 \lambda^{2} + 69 \lambda + 24}{2 \lambda \left(\lambda + 2\right)}$ | ja | ja |
| 3:0:0 | $\frac{1}{2}$ | $0$ | $- \frac{1}{2 \left(\lambda + 2\right)}$ | ja | ja |
| 3:0:1 | $\frac{1}{2}$ | $0$ | $- \frac{1}{2 \left(\lambda + 2\right)}$ | ja | ja |
| 3:0:2 | $\frac{6 \lambda^{2} + 20 \lambda + 17}{2}$ | $- 6 \left(\lambda + 2\right)$ | $- \frac{12 \lambda^{3} + 52 \lambda^{2} + 69 \lambda + 24}{2 \lambda \left(\lambda + 2\right)}$ | ja | ja |
| 3:1:-1 | $\frac{2 \lambda^{2} + 4 \lambda + 1}{2}$ | $- \frac{2 \left(\lambda - 2\right) \left(\lambda + 2\right)}{\lambda}$ | $- \frac{4 \lambda^{4} + 12 \lambda^{3} + 9 \lambda^{2} + 8 \lambda + 16}{2 \lambda^{2} \left(\lambda + 2\right)}$ | ja | ja |
| 3:1:0 | $- \frac{1}{2}$ | $0$ | $\frac{3 \lambda + 8}{2 \lambda \left(\lambda + 2\right)}$ | ja | ja |
| 3:1:1 | $- \frac{1}{2}$ | $0$ | $\frac{3 \lambda + 8}{2 \lambda \left(\lambda + 2\right)}$ | ja | ja |
| 3:1:2 | $\frac{2 \lambda^{2} + 4 \lambda + 1}{2}$ | $- \frac{2 \left(\lambda - 2\right) \left(\lambda + 2\right)}{\lambda}$ | $- \frac{4 \lambda^{4} + 12 \lambda^{3} + 9 \lambda^{2} + 8 \lambda + 16}{2 \lambda^{2} \left(\lambda + 2\right)}$ | ja | ja |
| 3:2:0 | $- \frac{14 \lambda^{2} + 52 \lambda + 49}{2}$ | $\frac{2 \left(\lambda + 2\right) \left(7 \lambda + 6\right)}{\lambda}$ | $\frac{28 \lambda^{4} + 132 \lambda^{3} + 187 \lambda^{2} + 48 \lambda - 48}{2 \lambda^{2} \left(\lambda + 2\right)}$ | ja | ja |
| 3:2:1 | $- \frac{14 \lambda^{2} + 52 \lambda + 49}{2}$ | $\frac{2 \left(\lambda + 2\right) \left(7 \lambda + 6\right)}{\lambda}$ | $\frac{28 \lambda^{4} + 132 \lambda^{3} + 187 \lambda^{2} + 48 \lambda - 48}{2 \lambda^{2} \left(\lambda + 2\right)}$ | ja | ja |
