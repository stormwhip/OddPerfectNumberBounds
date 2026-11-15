# Odd Perfect Number Bounds

This repository contains a computational extension of the exponent bound for restrained odd perfect numbers. Under the reduction, small-prime exclusions, and Lemma 8 of Ochem and Zelinsky, the recorded sieve and rational verification give **q > 7.705 x 10^19** at **B = 5.029 x 10^13**. The [research report](ResearchReport.md) contains the proof, data summary, and full reproduction instructions.

## Verify the recorded bound

Python 3.11 is sufficient. The verifier uses the standard library and exact rational arithmetic.

    python QBoundVerify.py 50290000000000 18909040000000000000000 77050000000000000000 189326366793 86154121237231668

The last two arguments are the sum of the candidate counts and the sum of the integer reciprocal bounds in [QBoundSieveRuns.tsv](QBoundSieveRuns.tsv). Its intervals are lower-exclusive and upper-inclusive. Reciprocal bounds use a scale of 10^18.

## Reproduce or extend the sieve

Build with GCC or MinGW-w64 and C++17. The sieve accepts upper endpoints through 10^14 and accumulates reciprocal bounds in unsigned 128-bit integers.

    g++ -std=c++17 -O3 -pthread QBoundPrefix.cpp -o QBoundPrefix.exe

Each row can be reproduced with the arguments upper endpoint, worker count, lower endpoint, and exclusion-prime limit. For example, the last complete interval was computed with:

    .\QBoundPrefix.exe 50290000000000 16 50280000000000 3000000

The program prints the candidate count and scaled reciprocal bound on separate lines. For this interval the outputs are 36314737 and 722196535949.

## Original work and attribution

Pascal Ochem and Joshua Zelinsky, *On odd perfect numbers with exactly one even exponent greater than 2*, [arXiv:2607.19746v1](https://arxiv.org/abs/2607.19746v1). The [bundled paper](OriginalWork/RestrainedOddPerfectNumbers.pdf) is provided under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). Its PDF contents are unchanged. The original programs and data are from the [authors' companion page](https://www.lirmm.fr/~ochem/opn/restrained/). They are reproduced with their contents intact and their filenames adapted to the repository convention.

| Companion filename | Repository file |
| --- | --- |
| all.cpp | [All.cpp](OriginalWork/Source/All.cpp) |
| pairs.py | [Pairs.py](OriginalWork/Source/Pairs.py) |
| loopless.cpp | [Loopless.cpp](OriginalWork/Source/Loopless.cpp) |
| tree.cpp | [Tree.cpp](OriginalWork/Source/Tree.cpp) |
| tree.txt | [Tree.txt](OriginalWork/Source/Tree.txt) |
| checkfacts.txt | [CheckFacts.txt](OriginalWork/Source/CheckFacts.txt) |
| 3p2.txt | [ThreeSquared.txt](OriginalWork/Source/ThreeSquared.txt) |

When running the original tree program, pass -db CheckFacts.txt from the OriginalWork/Source directory to select the renamed data file.

The appendix of the [research report](ResearchReport.md) includes Ruize Zhao's Pell-equation example for Question 9.
