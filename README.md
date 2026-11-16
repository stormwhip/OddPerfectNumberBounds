# Odd Perfect Number Bounds / 受限奇完全数指数下界研究

Repository with **two separated research stages**. The verified mathematical bound currently recorded is **q > 2.20 × 10^20**. The original article (Ochem–Zelinsky, 2026) proved **q ≥ 2.3 × 10^10**.

## Stage 1 — Pure mathematical and computational result (no ML)

[Stage 1](Stage%201/) ([阶段说明](Stage%201/README.md)) contains all material needed for the non-ML result:

- [Mathematical report](Stage%201/ResearchReport.md): proof, finite sieve, optimized 97-term inclusion–exclusion majorant.
- [Latest exact verifier](Stage%201/QBoundVerifyAPRefined.py): q > 2.20 × 10^20 using modulus-specific prime-counting errors.
- [Deterministic research plan](Stage%201/DeterministicResearchPlan.md): error constants, Selberg sieve and finite-sieve improvements.
- [Earlier exact verifier](Stage%201/QBoundVerifyOptimized.py): verifies 4096 Boolean assignments, 97 terms, rational tail bound and 4943 recorded sieve rows.
- [Initial verifier](Stage%201/QBoundVerify.py), [C++ finite sieve](Stage%201/QBoundPrefix.cpp) and [recorded output](Stage%201/QBoundSieveRuns.tsv).
- [Original paper and companion data](Stage%201/OriginalWork/): Ochem and Zelinsky's archived source, factor-tree data and attributed PDF, kept intact.

Run from the repository root:

```bash
cd "Stage 1"
python QBoundVerifyAPRefined.py
python QBoundVerifyOptimized.py
g++ -std=c++17 -O3 -pthread QBoundPrefix.cpp -o QBoundPrefix
./QBoundPrefix 50290000000000 16 50280000000000 3000000
```

Latest verifier: `CERTIFIED: q > 220000000000000000000`. Earlier verifier: `CERTIFIED: q > 219300000000000000000`. The displayed final finite-sieve interval produces 36,314,737 candidates and scaled reciprocal sum 722,196,535,949. The TSV stores 4,943 contiguous intervals; the whole finite sieve awaits a fully independent rerun before paper submission.

## Stage 2 — Learning-guided search: experiments and planned extensions

[Stage 2](Stage%202/) ([阶段说明](Stage%202/README.md)) contains:

- [Research direction and plan](Stage%202/LearningGuidedSievePlan.md).
- [Initial pilot report](Stage%202/experiments/PilotReport.md).
- [Local ML-vs-rule experiment](Stage%202/experiments/learning_guided_pilot.py) and [full LP baseline](Stage%202/experiments/lp_full_search.py).

Current pilot: **LightGBM 72.70%**, deterministic rule **75.06%**, random selection **54.49%** of full LP objective improvement across four seeds and eighty hold-out instances. The model beat random choice, but *did not* beat the deterministic rule. No new q theorem has been obtained from ML.

Run from the root (requires numpy, scipy, scikit-learn, lightgbm):

```bash
python "Stage 2/experiments/lp_full_search.py" 12 100000
python "Stage 2/experiments/learning_guided_pilot.py" --train 60 --test 20 --budget 5 --seed 42
```

Experimental floating-point LP objectives are distinct from Stage 1's rigorous integer/rational proof. No paid model or third-party API is used.

## Original work and attribution

Pascal Ochem and Joshua Zelinsky, *On odd perfect numbers with exactly one even exponent greater than 2*, [arXiv:2607.19746v1](https://arxiv.org/abs/2607.19746v1). The [archived PDF](Stage%201/OriginalWork/RestrainedOddPerfectNumbers.pdf) is distributed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/); the authors' archived source and data are unchanged.

## Latest Stage 2 exploration: larger exact-constraint LPs

[Extended numerical LP exploration](Stage%202/experiments/ExtendedSearchReport.md): with degree ≤4 and modulus ≤100,000, full numerical LP over 12/14/15 primes (214/279/310 columns, 4,095/16,383/32,767 Boolean constraints) returns the same 97-term polynomial. The current small-modulus LP family appears saturated; this does not prove global optimality for more general analytic sieves.

Code: [extended_lp_search.py](Stage%202/experiments/extended_lp_search.py); [results JSON](Stage%202/experiments/extended_search_summary.json). Research will next investigate alternative mathematically admissible sieve weights and genuinely large-scale column/constraint generation.

### Follow-up: deterministic constraint generation

[Constraint-generation report](Stage%202/experiments/ConstraintGenerationReport.md), [implementation](Stage%202/experiments/constraint_generation.py), [recorded summary](Stage%202/experiments/constraint_generation_summary.json). In a 16-prime problem with 65,535 Boolean states and 340 candidate monomials, just 339 initial active constraints found the same 97-term LP result; enumeration checked all states and found no violations. This is a **numerical deterministic search result**, not an improved number-theory certificate. The new ML baseline should exceed this approach on genuinely difficult instances.
