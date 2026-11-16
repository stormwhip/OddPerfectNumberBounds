# Deterministic algorithm research directions

Date: 2026-11-16. These are hypotheses and experimental plans, not claims that a more powerful algorithm has already been implemented. The benchmark is [QBoundVerifyAPRefined.py](QBoundVerifyAPRefined.py), which certifies q > 2.20 × 10^20 using the existing finite-sieve data.

## 1. Highest-priority immediate direction: modulus-specific prime-counting errors

The previous majorant has 97 terms. The original verifier used the same coefficient 1/840 for every term. Bennett–Martin–O'Bryant–Rechnitzer, *Explicit bounds for primes in arithmetic progressions* (Illinois J. Math. 62, 2018), Theorem 1.3 and Table 1, provide a sharper constant for q=7; their equation (1.15) provides another for q=1.

Changing only these two terms and preserving all other bounds produces:

- The old error constant: C = 433/840.
- The improved one: C = 21575999/42000000.
- At Z = 57,623,000,000,000,000,000,000, q > 220,000,000,000,000,000,000 is confirmed by exact rational arithmetic.

The next experiment is to retrieve the **published certified c_pi(m) constants and their validity thresholds for every modulus m actually present in the 97 terms**, compute the new total bound, and jointly re-optimize the polynomial and Z. The constants must be rounded *outward* and checked for all x >= B.

Source: https://personal.math.ubc.ca/~bennett/BeMaObRe-IJM-2018.pdf

## 2. Novel sieve-weight family: Selberg upper-bound (Lambda-squared) sieve

Let f(r)=r²+r+1 and P be a set of excluded primes. For squarefree d composed of elements of P, choose rational lambda_d with lambda_1=1. Define

G_lambda(r) = (sum_{d | f(r)} lambda_d)^2.

Then G_lambda(r) >= 1 whenever f(r) avoids every prime in P, and G_lambda(r) >= 0 otherwise, so it is a valid pointwise majorant. Expand the square using lcm(d,e):

delta(lambda) = sum_{d,e} lambda_d lambda_e * prod_{l | lcm(d,e)} 2/(l-1).

The explicit error bound can be obtained termwise with modulus m=lcm(d,e), provided each modulus and x threshold is permitted by the cited AP theorem. The search objective is an upper bound of the form A*delta(lambda)+E*C(lambda). Require exact rational coefficients and a separate verifier.

**Important baseline**: For the same legal monomial pool and the same error coefficients, a Lambda-squared polynomial is already representable by the full polynomial LP. Therefore it cannot improve upon the corresponding *globally optimized* LP without changing the allowed moduli, weight structure, or analytic upper bound. Its main potential is searching a much larger admissible domain efficiently.

General background: Selberg's upper-bound sieve (Lambda-squared construction), a classical sieve method; inspect exact explicit remainder estimates before applying it to this new application.

## 3. Improve the tail bound rather than only the LP

- Optimize the allowed cutoff Z with exact interval arithmetic.
- Use explicit error constants specific to each AP modulus and validity range.
- Explore an interval-dependent majorant G_j for B_j < r <= B_{j+1}; account for **all** endpoint error terms, which may offset the main-term gain.
- Investigate direct estimates for unions of residue classes, replacing individual triangle inequalities where a stronger published theorem is available.

Searches that merely increase the polynomial degree within a fixed monomial pool are not expected to improve bounds when the corresponding polynomial optimization has already reached its optimum; future experiments should target stronger admissible weights or analytic estimates.

## 4. Expand the finite congruence sieve efficiently

Increase the finite cutoff B using CRT wheel pre-sieving, segmented Eratosthenes, cached forbidden residues, and independent interval checks. Use reciprocal sums with integer upward rounding. Benchmark **new certified q per CPU-hour** rather than raw sieve speed alone.

An exploratory scan of the existing 4,943 row prefixes, optimizing Z numerically for each, preferred the largest recorded endpoint B = 5.029e13; this suggests that additional correctly computed finite-sieve intervals may be productive. It does not predict a numeric improvement beyond the saved data.

## 5. Independent verification & publication criteria

Before announcing a new number-theory theorem: (1) check the analytic theorem's valid x and modulus range for every term; (2) generate an integer/rational certificate for coefficients, cutoffs, finite sieve totals, and prime-counting errors; (3) rerun the finite sieve independently on a different code path; (4) archive exact commands, software versions, and raw outputs.

All numerical improvements should be backed by independently checkable integer or rational certificates.
