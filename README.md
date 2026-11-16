# Odd Perfect Number Bounds

Code and numerical data supporting an explicit lower bound for the exceptional exponent in restrained odd perfect numbers.

## Files

- **[prefix.cpp](prefix.cpp)** — Multithreaded segmented congruence sieve. Counts surviving prime candidates in a specified interval and computes an integer upper bound for their reciprocal sum.
- **[sieve.tsv](sieve.tsv)** — Recorded results from 4,943 sieve intervals, including candidate counts and scaled reciprocal-sum bounds.
- **[verify.py](verify.py)** — Standalone exact-arithmetic verifier. Checks the 97-term sieve majorant, validates the recorded sieve totals, and verifies the numerical inequality yielding \(q>2.20\times10^{20}\).

## Run

```bash
python verify.py
g++ -std=c++17 -O2 -pthread prefix.cpp -o prefix
./prefix 50290000000000 16 50280000000000 3000000
```

The Python verifier uses only the standard library. It checks the saved sieve records; a complete independent rerun of the finite sieve is a separate verification step.

The structural reduction and excluded-prime conditions are taken from Pascal Ochem and Joshua Zelinsky, [*On odd perfect numbers with exactly one even exponent greater than 2*](https://arxiv.org/abs/2607.19746).
