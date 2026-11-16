# Odd Perfect Number Bounds

Code and numerical data supporting an explicit lower bound for the exceptional exponent in restrained odd perfect numbers.

## Files

- **[prefix.cpp](prefix.cpp)** — Multithreaded segmented congruence sieve. Counts surviving prime candidates in a specified interval and computes an integer upper bound for their reciprocal sum.
- **[sieve.tsv](sieve.tsv)** — Recorded results from 4,943 sieve intervals, including candidate counts and scaled reciprocal-sum bounds.
- **[chain.cpp](chain.cpp)** — Two-step factor-chain sieve. Tests the candidates left by the full congruence sieve and groups candidates that require a specified special prime.
- **[chain.tsv](chain.tsv)** — Counts, scaled reciprocal sums and witnesses for the two-step sieve on \((3\,000\,000,10^9]\).
- **[verify.py](verify.py)** — Exact-arithmetic verifier. Checks the 97-term majorant and certifies \(q>2.20\times10^{20}\) from `sieve.tsv`, or \(q>1.14\times10^{22}\) with `chain.tsv`.

## Run

```bash
python verify.py
g++ -std=c++17 -O2 -pthread prefix.cpp -o prefix
./prefix 50290000000000 16 50280000000000 3000000
g++ -std=c++17 -O2 -pthread chain.cpp -o chain
./chain 1000000000 20 3000000 10000 > chain.tsv
python verify.py chain
```

The Python verifier uses only the standard library. It checks the saved sieve records; a complete independent rerun of the finite sieve is a separate verification step.

The structural reduction and excluded-prime conditions are taken from Pascal Ochem and Joshua Zelinsky, [*On odd perfect numbers with exactly one even exponent greater than 2*](https://arxiv.org/abs/2607.19746).
