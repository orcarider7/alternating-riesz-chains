"""Check four enclosures of width 10^-7 with exact integer arithmetic."""

import gzip
import json
from fractions import Fraction
from pathlib import Path

from check_collapse_onset import require, verify_flow


# The integer ticks are exact charge ratios with denominator ten million.
LOWER_TICKS = {73: 60049142, 75: 60004364, 77: 59961805, 79: 59921303}
TICK_SCALE = 10**7
ENERGY_SCALE = 10**15


def energy_upper(reference, unit, charge):
    """Bound the signed open energy at x = reference / unit from above."""
    p, d = charge.numerator, charge.denominator
    charges = [p if i % 2 == 0 else -d for i in range(len(reference))]
    terms = (
        (ENERGY_SCALE * unit**2 * charges[i] * charges[j])
        // (d**2 * (reference[j] - reference[i])**2)
        for j in range(len(reference)) for i in range(j)
    )
    return -sum(terms)


def check(n, ticks):
    path = Path(__file__).with_name("charge-windows") / f"open-{n}.json.gz"
    with gzip.open(path, "rt", encoding="utf-8") as stream:
        data = json.load(stream)
    lower, upper = Fraction(ticks, TICK_SCALE), Fraction(ticks + 1, TICK_SCALE)
    require(Fraction(data["lower"]) == lower and Fraction(data["upper"]) == upper, "wrong interval")
    _, _, records, utilization = verify_flow(data, n, lower)
    bound = energy_upper(data["reference"], data["unit"], upper)
    require(bound < 0, "upper witness must have negative energy")
    return records, utilization, bound


def print_tick(tick):
    whole, remainder = divmod(tick, TICK_SCALE)
    return f"{whole}.{remainder:07d}"


if __name__ == "__main__":
    for n, tick in LOWER_TICKS.items():
        count, fraction, bound = check(n, tick)
        print(f"Period {n + 1}: ({print_tick(tick)}, {print_tick(tick + 1)}), {count} flow records")
        print(f"  Maximum capacity fraction {fraction} < 1.  At upper endpoint, 10^15 S <= {bound}.")
