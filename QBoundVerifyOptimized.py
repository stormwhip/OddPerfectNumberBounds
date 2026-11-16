"""Reproduce a certified improvement of the restrained OPN exponent lower bound.

This script uses only Python's standard library and the accompanying
QBoundVerify.py and QBoundSieveRuns.tsv from the user's original report.
It does not need paid APIs or new large-scale sieving.
"""
from fractions import Fraction
from itertools import combinations
from math import prod
from pathlib import Path
import csv
import QBoundVerify as original

P = (7, 13, 19, 31, 37, 43, 61, 67, 73, 79, 97, 103)
TRIPLES = [(7, a, b) for a, b in combinations(P[1:], 2) if 7 * a * b <= 10_000]
B = 50_290_000_000_000
Z = 57_433_500_000_000_000_000_000
Q = 219_300_000_000_000_000_000


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


if __name__ == '__main__':
    delta, error = certificate()
    source = Path(__file__).with_name('QBoundSieveRuns.tsv')
    survivors, scaled_sum, nrows = finite_sieve_totals(source)
    assert survivors == 189_326_366_793
    assert scaled_sum == 86_154_121_237_231_668
    original.tail_count_constants = lambda: (delta, error)
    upper, threshold = original.verify(B, Z, Q, survivors, scaled_sum)
    print(f'Validated {1 << len(P)} Boolean assignments and {nrows} sieve rows.')
    print(f'G has 97 terms, with {len(TRIPLES)} selected triple intersections.')
    print(f'delta = {delta}; C = {error}')
    print(f'B = {B}; Z = {Z}; proposed q = {Q}')
    from math import ceil, floor
    scale = 10 ** 15
    upper_rounded = ceil(upper * scale)
    lower_rounded = floor(threshold * scale)
    print(f'upper   <= {upper_rounded // scale}.{upper_rounded % scale:015d}')
    print(f'target  >= {lower_rounded // scale}.{lower_rounded % scale:015d}')
    print(f'margin  ~= {float(threshold - upper):.12g}')
    assert upper < threshold
    print(f'CERTIFIED: q > {Q}')
