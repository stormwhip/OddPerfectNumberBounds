"""Verify the explicit exponent bound for restrained odd perfect numbers.

Run: python verify.py
Requires: sieve.tsv in the same directory.
Uses Python's standard library and exact rational arithmetic throughout
the final certificate. The saved sieve data require an independent full
rerun before submission.
"""
from fractions import Fraction
from itertools import combinations
from math import ceil, floor, prod
from pathlib import Path
import csv

SCALE = 10**18
P = (7, 13, 19, 31, 37, 43, 61, 67, 73, 79, 97, 103)
TRIPLES = [(7, a, b) for a, b in combinations(P[1:], 2) if 7 * a * b <= 10_000]
B = 50_290_000_000_000
Z = 57_623_000_000_000_000_000_000
Q = 220_000_000_000_000_000_000


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



def certificate():
    terms = [(subset, (-1) ** size)
             for size in range(3) for subset in combinations(P, size)]
    terms += [(triple, -1) for triple in TRIPLES]
    assert len(terms) == 97 and len(TRIPLES) == 18
    assert all(prod(subset) <= 10_000 for subset, _ in terms)

    # Check every Boolean assignment, independently of floating-point values.
    for mask in range(1 << len(P)):
        bits = {p: bool(mask & (1 << j)) for j, p in enumerate(P)}
        k = mask.bit_count()
        G = 1 - k + k * (k - 1) // 2
        G -= sum(all(bits[p] for p in triple) for triple in TRIPLES)
        assert G >= (1 if mask == 0 else 0)

    delta = sum((sign * prod(Fraction(2, p - 1) for p in subset)
                 for subset, sign in terms), Fraction(0))
    error = sum((Fraction(2 ** len(subset), 840)
                 for subset, sign in terms), Fraction(0))
    assert delta == Fraction(64959667, 169646400)
    assert error == Fraction(433, 840)
    return delta, error


def finite_sieve_totals(filename):
    with open(filename, newline='', encoding='utf-8') as fh:
        rows = list(csv.DictReader(fh, delimiter='\t'))
    assert rows and int(rows[0]['lower_exclusive']) == 3_000_000
    count = total = 0
    previous = 3_000_000
    for row in rows:
        assert int(row['lower_exclusive']) == previous
        previous = int(row['upper_inclusive'])
        count += int(row['survivors'])
        total += int(row['scaled_reciprocal_upper'])
    assert previous == B
    return count, total, len(rows)



def verify(B, Z, q, survivors, scaled_sum, density, error_coefficient):
    log_B_lower, _ = log_bounds(B)
    _, log_Z_upper = log_bounds(Z)
    _, log_ratio_upper = log_bounds(log_Z_upper / log_B_lower)

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
    density, old_error = certificate()
    assert old_error == Fraction(433, 840)

    # Bennett et al.: per-residue-class AP errors for moduli 1 and 7.
    # The remaining 95 monomials retain the uniform 1/840 coefficient.
    c_pi_1 = Fraction(8375, 10_000_000)
    c_pi_7 = Fraction(486, 1_000_000)
    error = (old_error
             - (Fraction(1, 840) - c_pi_1)
             - 2 * (Fraction(1, 840) - c_pi_7))
    assert error == Fraction(21_575_999, 42_000_000)

    survivors, scaled_sum, rows = finite_sieve_totals(
        Path(__file__).with_name("sieve.tsv")
    )
    assert rows == 4943
    assert survivors == 189_326_366_793
    assert scaled_sum == 86_154_121_237_231_668

    upper, target = verify(B, Z, Q, survivors, scaled_sum, density, error)
    scale = 10**15
    upper_rounded = ceil(upper * scale)
    lower_rounded = floor(target * scale)

    print(f"Validated {1 << len(P)} Boolean assignments and {rows} sieve rows.")
    print(f"97 majorant terms; {len(TRIPLES)} selected triple intersections.")
    print(f"delta = {density}; C = {error}")
    print(f"upper <= {upper_rounded // scale}.{upper_rounded % scale:015d}")
    print(f"log(4/3) >= {lower_rounded // scale}.{lower_rounded % scale:015d}")
    if not upper < target:
        raise ArithmeticError("Final inequality does not certify the bound.")
    print(f"CERTIFIED: q > {Q}")


if __name__ == "__main__":
    main()
