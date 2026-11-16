"""Verify a q bound using the output of QBoundPrefix and rational arithmetic.

Use B = upper endpoint of the combined sieve, Z > B, and a proposed q bound.
Count the survivors in (3,000,000, B] and sum their integer reciprocal bounds.
Arguments: B Z q survivors scaled_sum
"""

import argparse
from fractions import Fraction
from itertools import combinations
from math import ceil, floor, prod


SCALE = 10**18


def log_unit_interval(value):
    """Bound log(value) for 1 <= value <= 2 by a series with a geometric tail."""
    ratio = (value - 1) / (value + 1)
    square = ratio * ratio
    power = ratio
    lower = Fraction(0)
    terms = 24
    for index in range(terms):
        lower += 2 * power / (2 * index + 1)
        power *= square
    upper = lower + 2 * power / ((2 * terms + 1) * (1 - square))
    return lower, upper


def log_bounds(value):
    value = Fraction(value)
    exponent = 0
    while value > 2:
        value /= 2
        exponent += 1
    lower, upper = log_unit_interval(value)
    log2_lower, log2_upper = log_unit_interval(Fraction(2))
    lower += exponent * log2_lower
    upper += exponent * log2_upper
    return Fraction(floor(lower * SCALE), SCALE), Fraction(ceil(upper * SCALE), SCALE)


def tail_count_constants():
    """Majorize the survivor indicator using exact exclusion at 7 and Bonferroni."""
    other_primes = (13, 19, 31, 37, 43, 61, 67)
    density = Fraction(0)
    error = Fraction(0)
    for size in range(3):
        for subset in combinations(other_primes, size):
            for factors in (subset, (7, *subset)):
                modulus = prod(factors)
                # The empty term counts all primes, using Bennett et al. (1.15).
                # Theorem 1.3 covers every other modulus, up to 7*61*67.
                divisor = 840 if modulus <= 10000 else 160
                density += (-1)**len(factors) * prod(
                    Fraction(2, p - 1) for p in factors
                )
                error += Fraction(2**len(factors), divisor)
    return density, error


def verify(B, Z, q, survivors, scaled_sum):
    log_B_lower, _ = log_bounds(B)
    _, log_Z_upper = log_bounds(Z)
    _, log_ratio_upper = log_bounds(log_Z_upper / log_B_lower)

    density, error_coefficient = tail_count_constants()
    b_ratio = Fraction(B, B - 1)

    # Lower bound on the integral of 1/log(t) over [B,Z] by one
    # integration by parts, bounding the integral of 1/log(t)^2 below.
    li_lower = (
        Z / log_Z_upper - B / log_B_lower + (Z - B) / log_Z_upper**2
    )
    tail_main = density * (b_ratio * log_ratio_upper - li_lower / (Z - 1))
    tail_error = error_coefficient * (
        b_ratio / log_B_lower**2
        + b_ratio**2 * (1 / log_B_lower - 1 / log_Z_upper)
    )

    # Lemma 8 allows at most 2q+2 ordinary primes. Each gets a baseline
    # 1/(Z-1); the remaining weight is counted only for candidates below Z.
    ordinary = Fraction(scaled_sum, SCALE) + tail_main + tail_error
    ordinary += Fraction(2 * q + 2 - survivors, Z - 1)
    total = ordinary + Fraction(1, 2999999)
    target, _ = log_bounds(Fraction(4, 3))
    return total, target


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("B", "Z", "q", "survivors", "scaled_sum"):
        parser.add_argument(name, type=int)
    args = parser.parse_args()
    if not (8 * 10**9 <= args.B < args.Z and args.q > 0
            and args.survivors >= 0 and args.scaled_sum >= 0):
        parser.error("require 8e9 <= B < Z, positive q and nonnegative sieve results")

    total, target = verify(**vars(args))
    # Print outward-rounded decimals while comparing the exact fractions.
    digits = 10**12
    upper_total = ceil(total * digits)
    lower_target = floor(target * digits)
    print(f"upper total: {upper_total // digits}.{upper_total % digits:012d}")
    print(f"lower target: {lower_target // digits}.{lower_target % digits:012d}")
    if total >= target:
        print("This input does not certify the proposed bound.")
        return 1
    print(f"Certified: q > {args.q}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
