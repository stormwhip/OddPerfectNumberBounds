"""Rigorous refinement of the OPN exponent bound using modulus-specific AP errors.

Bennett, Martin, O'Bryant, Rechnitzer: Illinois J. Math. 62 (2018),
Theorem 1.3 / Table 1 (q=7), and equation (1.15) (q=1).
The cited q=7 constant is < 0.0004857 at x >= 1,765,715,753;
we conservatively use 0.000486. For q=1, use 0.0008375
at x >= 1,474,279,333. Our B is larger than both thresholds.

The other terms use the original 1/840 bound. This verifier relies on
the recorded finite sieve data, and uses exact Fraction arithmetic.
"""
from fractions import Fraction
from math import ceil, floor
from pathlib import Path
import QBoundVerify as base
import QBoundVerifyOptimized as prior

B = prior.B
Z = 57_623_000_000_000_000_000_000
Q = 220_000_000_000_000_000_000


def main():
    delta, old_C = prior.certificate()
    assert old_C == Fraction(433, 840)
    # Polynomial coefficients: 1 for the empty set and -1 for {7}.
    # The singleton {7} corresponds to two reduced residue classes.
    c_pi_1 = Fraction(8375, 10_000_000)
    c_pi_7 = Fraction(486, 1_000_000)
    C = old_C - (Fraction(1,840) - c_pi_1) - 2*(Fraction(1,840) - c_pi_7)
    assert C == Fraction(21_575_999, 42_000_000)
    count, scaled_sum, rows = prior.finite_sieve_totals(
        Path(__file__).with_name('QBoundSieveRuns.tsv')
    )
    assert rows == 4943
    assert count == 189_326_366_793
    assert scaled_sum == 86_154_121_237_231_668

    base.tail_count_constants = lambda: (delta, C)
    upper, threshold = base.verify(B, Z, Q, count, scaled_sum)
    scale = 10**15
    up, low = ceil(upper * scale), floor(threshold * scale)
    print(f'delta={delta}, C={C}, Z={Z}')
    print(f'upper <= {up//scale}.{up%scale:015d}')
    print(f'log(4/3) >= {low//scale}.{low%scale:015d}')
    if not upper < threshold:
        raise ArithmeticError('Uncertified input; numerical inequality failed')
    print(f'CERTIFIED: q > {Q}')


if __name__ == '__main__':
    main()
