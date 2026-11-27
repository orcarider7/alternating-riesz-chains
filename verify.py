"""Reproduce inverse-square examples with Python 3 and no dependencies."""

from fractions import Fraction
from math import fsum, isclose, pi, sin, sqrt


def periodic_energy(gaps, mass):
    """Per-site signed energy, including every periodic image, for p = 2."""
    length = fsum(gaps)
    count = len(gaps)
    positions = [0.0]
    for gap in gaps[:-1]:
        positions.append(positions[-1] + gap)
    charges = [mass if i % 2 == 0 else -1 for i in range(count)]

    # The pair kernel is the exact sum over integer translates.
    pair_sum = fsum(
        charges[i] * charges[j] / sin(pi * (positions[j] - positions[i]) / length) ** 2
        for i in range(count)
        for j in range(i + 1, count)
    )
    self_images = fsum(charge**2 for charge in charges) / 3
    return -(pi / length) ** 2 * (2 * pair_sum + self_images) / count


def uniform_coefficient(mass):
    return pi**2 / 24 * (6 * mass - 1 - mass**2)


def open_block_coefficient(count, mass):
    """Exact coefficient of the inverse-square energy on an open unit-gap block."""
    return sum(
        Fraction((count - offset) * (2 * mass if offset % 2 else -(1 + mass**2)),
                 count * offset**2)
        for offset in range(1, count)
    )


def main():
    boundary = 3 + 2 * sqrt(2)
    gaps = (0.4, 1.1, 0.8, 1.2, 0.9, 1.6)
    spacing = fsum(gaps) / len(gaps)
    for mass in (1.0, 5.0, boundary, 1 / boundary):
        uniform = periodic_energy((spacing,) * len(gaps), mass)
        nonuniform = periodic_energy(gaps, mass)
        coefficient = uniform_coefficient(mass)
        gap_bound = coefficient * fsum(d**-2 for d in gaps) / len(gaps)
        assert isclose(uniform, coefficient * spacing**-2, abs_tol=1e-12)
        assert nonuniform > gap_bound - 1e-12
        assert nonuniform > uniform
        print(f"m={mass:.9f} | equal={uniform:.9f} | gap bound={gap_bound:.9f}"
              f" | unequal={nonuniform:.9f}")

    # Rational enclosures on both sides of the finite-block sign change.
    c160 = open_block_coefficient(160, 6)
    c162 = open_block_coefficient(162, 6)
    assert Fraction(2347, 10**6) < c160 < Fraction(2349, 10**6)
    assert Fraction(-1801, 10**6) < c162 < Fraction(-1798, 10**6)
    print(f"m=6 | C_160={float(c160):.9f} | C_162={float(c162):.9f}")

    count, mass = 256, 6
    coefficient = open_block_coefficient(count, mass)
    assert coefficient < 0  # This sign check is rational arithmetic.
    print(f"m={mass}, N={count} | open-block coefficient={float(coefficient):.9f}")
    energies = []
    for short_gap in (1 / 4, 1 / 16, 1 / 64):
        gaps = (short_gap,) * (count - 1) + (count - (count - 1) * short_gap,)
        energy = periodic_energy(gaps, mass)
        energies.append(energy)
        print(f"gap={short_gap:.6f} | E={energy:.9f}"
              f" | gap^2 E={short_gap**2 * energy:.9f}")
    assert energies[2] < energies[1] < energies[0] < 0
    assert abs(short_gap**2 * energies[-1] - float(coefficient)) < 1e-3


if __name__ == "__main__":
    main()
