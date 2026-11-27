"""Exact certificates for the first collapsing even period. Python standard library."""

import json
from collections import defaultdict
from fractions import Fraction
from pathlib import Path


def require(condition, message):
    if not condition:
        raise ValueError(message)


def verify_lower(data):
    """Check that aggregate source flows serve every sink without exceeding capacity."""
    reference = data["reference"]
    unit, scale = data["unit"], data["scale"]
    n = len(reference)
    require(n == 75 and all(type(v) is int for v in reference), "expected 75 integer coordinates")
    require(all(a < b for a, b in zip(reference, reference[1:])), "reference order")
    require(all(type(v) is int and v > 0 for v in (unit, scale)), "invalid scale")
    normalization = scale * unit**3
    balance = defaultdict(int)
    load = defaultdict(int)
    for source, u, v, amount in data["flow"]:
        require(all(type(k) is int for k in (source, u, v, amount)), "noninteger flow")
        require(0 <= source <= u < v < n and (v-u) % 2 == 1 and amount > 0, "invalid arc")
        balance[source, u] -= amount
        balance[source, v] += amount
        load[u, v] += amount

    demands = 0
    for source in range(n):
        for sink in range(source + 1, n):
            absorption = balance[source, sink]
            same = (sink-source) % 2 == 0
            weight = (36 if source % 2 == 0 else 1) if same else 0
            require(absorption * (reference[sink]-reference[source])**3 >= weight * normalization,
                    f"unmet demand at source {source}, sink {sink}")
            demands += int(same)

    mixed = [(u, v) for u in range(n) for v in range(u + 1, n, 2)]
    used = max(load[u, v] * (reference[v]-reference[u])**3 for u, v in mixed)
    capacity = 6 * normalization
    require(used < capacity, "insufficient mixed-edge capacity")
    return demands, len(mixed), len(data["flow"]), Fraction(used, capacity)


def verify_upper():
    """Verify a strictly negative signed pair sum for the 77-site open block."""
    radius = 38
    numerators = [2*k*radius**4 + k**5 for k in range(-radius, radius + 1)]
    z = [2*h if j % 2 == 0 else numerators[j-1] + numerators[j+1]
         for j, h in enumerate(numerators)]
    denominator = 2 * radius**5
    require(all(a < b for a, b in zip(z, z[1:])), "witness order")
    charges = [6 if j % 2 == 0 else -1 for j in range(len(z))]
    million = 10**6
    upper = -sum((million * denominator**2 * charges[i] * charges[j]) // (z[j]-z[i])**2
                 for j in range(len(z)) for i in range(j))
    require(upper < 0, "witness energy must be negative")
    return upper


if __name__ == "__main__":
    cert = Path(__file__).with_name("flow-75.json")
    result = verify_lower(json.loads(cert.read_text(encoding="utf-8")))
    print("Lower certificate: {} demands, {} mixed edges, {} integer records.".format(*result[:3]))
    print(f"Maximum capacity fraction = {result[3]} < 1.")
    print(f"Upper certificate: 1000000 S_77 <= {verify_upper()} < 0.")
