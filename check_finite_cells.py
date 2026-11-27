"""Reproduce the numerical checks in four-sites.md and inverse-square.md.

Run with Python 3.8 or later. No third-party libraries or network required.
"""

import cmath
from fractions import Fraction
from math import fsum, isclose, pi, sin, tan


def check_close(actual, expected, tolerance=2e-9):
    if not isclose(actual, expected, rel_tol=tolerance, abs_tol=tolerance):
        raise AssertionError(f"{actual!r} differs from {expected!r}")


def inverse_square(s, length):
    return (pi / length / sin(pi * s / length)) ** 2


def shifted_image_sum(s, length, p):
    """Compute S(s)-S(L/2), pairing opposite images to cancel slow tails."""
    center = length / 2
    terms = [s ** (-p) - center ** (-p)]
    for k in range(1, 2501):
        base = k * length
        terms.extend(((base + s) ** (-p) - (base + center) ** (-p),
                      (base - s) ** (-p) - (base - center) ** (-p)))
    return fsum(terms)


def pair_energy(positions, charges, kernel):
    return fsum(-charges[i] * charges[j] * kernel(positions[j] - positions[i])
                for j in range(len(positions)) for i in range(j))


def check_four_sites():
    length = 4
    uniform = (0, 1, 2, 3)
    for p in (1.5, 2, 3.25):
        def kernel(s):
            return shifted_image_sum(s, length, p)

        # Comparing two configurations cancels both self-images and the
        # constant subtracted from every distinct pair in shifted_image_sum.
        for ratio in (1, 2 ** (p + 1), 2 ** (p + 1) + 1):
            charges = (ratio, -1, ratio, -1)
            reference = pair_energy(uniform, charges, kernel)
            for x in (0.7, 1.4, 2.6):
                centered = (0, x / 2, x, (length + x) / 2)
                actual = pair_energy(centered, charges, kernel)
                expected = ratio * (2 ** (p + 1) - ratio) * kernel(x)
                check_close(actual - reference, expected, 2e-8)
                shifted = (0, x * 0.4, x, centered[3] + 0.17)
                assert pair_energy(shifted, charges, kernel) > actual

    def kernel(s):
        return inverse_square(s, length)

    for sites in (uniform, (0, 0.5, 1, 2.5)):
        pair = pair_energy(sites, (8, -1, 8, -1), kernel)
        total = pair / 2 - (8 ** 2 + 1) * pi ** 2 / (6 * length ** 2)
        check_close(total, -17 * pi ** 2 / 24)
    delta = 1e-4
    for ratio in (6, 9):
        trimer = (0, delta, 2 * delta, length / 2)
        scaled = delta ** 2 * pair_energy(trimer, (ratio, -1, ratio, -1), kernel)
        check_close(scaled, 2 * ratio - ratio ** 2 / 4, 1e-7)
    print("Four sites: midpoint identity checked at p=1.5, 2 and 3.25.")
    print(f"At p=2, M=8, both exact energy checks give {-17 * pi ** 2 / 24:.10f}.")


def hessian_block(n, m, mode):
    """Construct the two-site Bloch block from periodic pair potentials."""
    theta = 2 * pi * mode / n
    block = [[0j, 0j], [0j, 0j]]
    for i, qi in enumerate((m, -1)):
        for j in range(n):
            if i == j:
                continue  # same-site images give a constant under displacement
            qj = m if j % 2 == 0 else -1
            kernel = inverse_square(j - i, n)
            # Second derivative of (pi/N)^2 csc^2(pi*s/N).
            curvature = -qi * qj * 2 * kernel * (3 * kernel - 2 * (pi / n) ** 2)
            block[i][i] += curvature
            block[i][j % 2] -= curvature * cmath.exp(1j * theta * (j - i))
    normalization = 12 * pi ** 4 / 96
    for row in block:
        for v in row:
            check_close(v.imag, 0)
    return [[v.real / normalization for v in row] for row in block]


def check_hessian():
    m = 6
    margins = {}
    for n in (4, 6, 94, 96):
        for mode in range(n // 4 + 1):
            t = Fraction(2 * mode, n)
            u = t * (1 - t)
            e = 2 * u ** 2
            d = 2 * t ** 2 * (3 - 2 * t)
            block = hessian_block(n, m, mode)
            check_close(block[0][0], m - m ** 2 * e)
            check_close(block[1][1], m - e)
            check_close(block[0][1], -m * (1 - d))
            check_close(block[0][1], block[1][0])
            determinant = block[0][0] * block[1][1] - block[0][1] ** 2
            margin = 6 + 8 * u + 2 * u ** 2 - (m + Fraction(1, m))
            check_close(determinant, m ** 2 * e * margin)
            if mode == 1:
                margins[n] = margin
                assert (determinant > 0) == (margin > 0)
    assert margins[94] == Fraction(23183, 29278086)
    assert margins[96] == Fraction(-7007, 2654208)
    print(f"Hessian: H_94 - 37/6 = {margins[94]}, H_96 - 37/6 = {margins[96]}.")


def check_cotangent_clusters():
    coefficients = {}
    for n in (6, 8, 98, 100):
        for m in (Fraction(1, 2), Fraction(6)):
            a, b = (m - 1) / 2, (m + 1) / 2
            exact = -Fraction(n, 16) * (a ** 2 * (n - 1) ** 2
                                        - b ** 2 * (Fraction(n ** 2, 2) - 1))
            charges = [float(m) if i % 2 == 0 else -1 for i in range(n)]
            for shift in (0.23, 0.5):
                sites = [-1 / tan(pi * (j + shift) / n) for j in range(n)]
                direct = pair_energy(sites, charges, lambda s: 1 / s ** 2)
                check_close(direct, exact)
            if m == 6:
                coefficients[n] = exact
    assert coefficients[98] == Fraction(147, 4)
    assert coefficients[100] == Fraction(-925, 8)
    assert 2 * coefficients[100] / 100 == Fraction(-37, 16)
    print(f"Cotangent clusters: S_98 = {coefficients[98]}, S_100 = {coefficients[100]}.")


def check_rational_clusters():
    """Certify collapse using ordered rational points and integer arithmetic."""
    # Numerical minimization suggested these integers. They specify the
    # certificate exactly and no optimizer is needed to verify it.
    base = (
        -2055656, -1665239, -1386664, -1295986, -1203764, -1136963,
        -1069087, -1011194, -952559, -903460, -853865, -811668,
        -769108, -730925, -692292, -656624, -620935, -587911,
        -554682, -523002, -491190, -461106, -430929, -402262,
        -373515, -345934, -318298, -291720, -265102, -239131,
        -213125, -187580, -162011, -136821, -111619, -86702,
        -61775, -36992, -12211, 12469, 37130, 61849,
        86557, 111509, 136446, 161752, 187034, 212882,
        238691, 265210, 291689, 319018, 346281, 374372,
        402375, 431468, 460492, 491169, 521732, 553984,
        586119, 620892, 655446, 693461, 731426, 772249,
        812946, 857931, 902636, 954215, 1005218, 1065333,
        1124763, 1211505, 1297176, 1446015, 1579906, 2057713,
    )
    for n, expected in ((78, -6134970), (96, -6135499), (98, -6135404)):
        z = base + tuple(10**9 * k for k in range(1, n - 78 + 1))
        assert len(z) == n and all(a < b for a, b in zip(z, z[1:]))
        q = [6 if j % 2 == 0 else -1 for j in range(n)]
        # y_i = z_i/10^6. Each term is the ceiling of 10^6 times its
        # signed pair contribution, so their sum bounds 10^6 S from above.
        upper = sum(-((q[i] * q[j] * 10**18) // (z[j] - z[i])**2)
                    for j in range(n) for i in range(j))
        assert upper == expected and upper < 0
        print(f"Rational cluster: N={n}, 10^6 S <= {upper}.")


if __name__ == "__main__":
    check_four_sites()
    check_hessian()
    check_cotangent_clusters()
    check_rational_clusters()
