"""Generate an alternative lower certificate from the published reference using SciPy.

Run python generate_flow.py NEW_FILE.json. The routing can depend on the
SciPy/HiGHS version. Validation uses only exact integer arithmetic.
"""

import argparse
import json
from collections import defaultdict
from pathlib import Path

from scipy.optimize import linprog
from scipy.sparse import coo_matrix

from check_collapse_onset import verify_lower


def generate(reference, unit, scale):
    n = len(reference)
    mixed = [(u, v) for u in range(n) for v in range(u + 1, n, 2)]
    arcs = []
    equality_rows, equality_cols, equality_values = [], [], []
    capacity_rows, capacity_cols, capacity_values = [], [], []
    demands = []
    # Source flow can terminate at any later site of the same parity.
    for source in range(n - 2):
        last = n - 1 - (n - 1 - source) % 2
        weight = 36 if source % 2 == 0 else 1
        first_row = len(demands)
        for sink in range(source + 1, last + 1):
            length = (reference[sink] - reference[source]) / unit
            demands.append(weight / length**3 if (sink - source) % 2 == 0 else 0)
        for edge, (u, v) in enumerate(mixed):
            if u < source or v > last:
                continue
            col = len(arcs)
            arcs.append((source, u, v))
            capacity_rows.append(edge)
            capacity_cols.append(col)
            capacity_values.append(1)
            equality_rows.append(first_row + v - source - 1)
            equality_cols.append(col)
            equality_values.append(1)
            if u > source:
                equality_rows.append(first_row + u - source - 1)
                equality_cols.append(col)
                equality_values.append(-1)

    # Minimize the largest fraction of capacity used on any mixed edge.
    max_load_col = len(arcs)
    cols = max_load_col + 1
    for edge, (u, v) in enumerate(mixed):
        capacity_rows.append(edge)
        capacity_cols.append(max_load_col)
        capacity_values.append(-6 * (unit / (reference[v] - reference[u]))**3)
    capacity = coo_matrix((capacity_values, (capacity_rows, capacity_cols)),
                          shape=(len(mixed), cols)).tocsc()
    balance = coo_matrix((equality_values, (equality_rows, equality_cols)),
                         shape=(len(demands), cols)).tocsc()
    result = linprog([0] * len(arcs) + [1], A_ub=capacity, b_ub=[0] * len(mixed),
                     A_eq=balance, b_eq=demands, bounds=(0, None), method="highs-ds",
                     options={"primal_feasibility_tolerance": 1e-9,
                              "dual_feasibility_tolerance": 1e-9})
    if not result.success:
        raise RuntimeError(result.message)

    flow = defaultdict(int)
    absorbed = defaultdict(int)
    for (source, u, v), value in zip(arcs, result.x):
        amount = round(scale * float(value))
        if amount > 0:
            flow[source, u, v] = amount
            absorbed[source, u] -= amount
            absorbed[source, v] += amount

    # An incoming unit on the nearest edge repairs a deficit at sink and
    # moves it one site left. A right-to-left pass preserves repaired sinks.
    normalization = scale * unit**3
    for source in range(n - 2):
        weight = 36 if source % 2 == 0 else 1
        for sink in range(n - 1, source, -1):
            demand = 0
            if (sink - source) % 2 == 0:
                span = (reference[sink] - reference[source])**3
                demand = (weight * normalization + span - 1) // span
            extra = demand - absorbed[source, sink]
            if extra > 0:
                flow[source, sink - 1, sink] += extra
                absorbed[source, sink] += extra
                absorbed[source, sink - 1] -= extra
    return {"unit": unit, "scale": scale, "reference": reference,
            "flow": [[*arc, amount] for arc, amount in sorted(flow.items())]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path, help="new certificate file to create")
    parser.add_argument("--reference", type=Path, default=Path(__file__).with_name("flow-75.json"),
                        help="JSON file supplying the reference coordinates and integer scales")
    args = parser.parse_args()
    if args.output.exists():
        parser.error("output already exists")
    source = json.loads(args.reference.read_text(encoding="utf-8"))
    data = generate(source["reference"], source["unit"], source["scale"])
    demands, edges, records, fraction = verify_lower(data)
    with args.output.open("x", encoding="utf-8") as stream:
        json.dump(data, stream, separators=(",", ":"))
        stream.write("\n")
    print(f"Verified {demands} demands, {edges} mixed edges and {records} flow records.")
    print(f"Maximum capacity fraction = {fraction}")


if __name__ == "__main__":
    main()
