# Odd Perfect Number Bounds

Certified analytic and computational bounds for the exceptional even exponent of a restrained odd perfect number, extending Ochem and Zelinsky (2026).

## Current result

**q > 219,300,000,000,000,000,000 (2.193 × 10^20)**, conditional on the original paper's reduction and exclusions and the recorded finite sieve data.

This improves the prior repository bound **q > 7.705 × 10^19** by approximately **2.85×**. The original article proved q ≥ 2.3 × 10^10. The new 97-term majorant was discovered through **deterministic optimization**; no machine-learning model has been trained.

The finite sieve covers 4,943 contiguous intervals up to B=50,290,000,000,000, with 189,326,366,793 candidates and scaled reciprocal sum 86,154,121,237,231,668. The improved certificate uses 12 excluded primes and checks every one of 4,096 Boolean assignments, with delta=64959667/169646400 and C=433/840.

## Quick reproduction (Python standard library only)

Run the improved verifier from the repository root:

    python QBoundVerifyOptimized.py

Expected terminal line:

    CERTIFIED: q > 219300000000000000000

The program checks the upper-majorant inequality, 97 permitted terms, exact coefficients, entire saved TSV continuity and totals, and the final exact rational inequality. It imports the original QBoundVerify.py.

Reproduce the earlier bound:

    python QBoundVerify.py 50290000000000 18909040000000000000000 77050000000000000000 189326366793 86154121237231668

## Finite sieve

Build using GCC/MinGW-w64 and C++17:

    g++ -std=c++17 -O3 -pthread QBoundPrefix.cpp -o QBoundPrefix
    ./QBoundPrefix 50290000000000 16 50280000000000 3000000

The last interval reproduces 36,314,737 candidates and scaled reciprocal upper sum 722,196,535,949. The recorded TSV has 4,943 interval rows. **A fully independent rerun of every interval remains outstanding** before submission.

## Project files and research direction

- [Research report](ResearchReport.md): mathematical proof, original and optimized bounds, and Question 9.
- [Optimized verifier](QBoundVerifyOptimized.py): fully exact new certificate.
- [Original verifier](QBoundVerify.py), [C++ segmented sieve](QBoundPrefix.cpp), and [recorded runs](QBoundSieveRuns.tsv).
- [Learning-guided certified sieve search: detailed research plan](research/LearningGuidedSievePlan.md): mathematical LP formulation, learning-guided column selection, non-ML baselines, exact verification, experimental design, metrics, and six-week milestones.

The **learning-guided method is proposed research, not an implemented ML result**. The plan uses local/open-source tools with no paid model/API calls.

## Original work and attribution

Pascal Ochem and Joshua Zelinsky, *On odd perfect numbers with exactly one even exponent greater than 2*, [arXiv:2607.19746v1](https://arxiv.org/abs/2607.19746v1). Their unchanged [PDF](OriginalWork/RestrainedOddPerfectNumbers.pdf) is provided under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). Their archived companion code/data are in [OriginalWork](OriginalWork), unchanged by this update.

The report also gives Ruize Zhao's Pell-equation prime example addressing the paper's Question 9.
