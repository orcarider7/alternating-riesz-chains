"""Enclose the period-six threshold with rational arithmetic."""

from fractions import Fraction
from math import ceil, floor


def polynomial(u):
    """Evaluate the stationary zero-energy equation after clearing denominators."""
    return (
        256 * u**2 * (3 + u) * (1 + 3 * u)
        - 128 * (1 + 3 * u) * (1 - u)**3
        + 9 * (1 - u)**6
    )


def charge_ratio(u):
    return 16 * (1 + 3 * u) / (9 * (1 - u)**3)


def enclose():
    left, right = Fraction(1, 4), Fraction(13, 50)
    if not polynomial(left) < 0 < polynomial(right):
        raise ValueError("The polynomial does not change sign on the initial bracket")
    while charge_ratio(right) - charge_ratio(left) > Fraction(1, 10**8):
        middle = (left + right) / 2
        value = polynomial(middle)
        if value == 0:
            return charge_ratio(middle), charge_ratio(middle)
        if value < 0:
            left = middle
        else:
            right = middle
    return charge_ratio(left), charge_ratio(right)


def decimal_endpoint(scaled_integer, places):
    quotient, remainder = divmod(scaled_integer, 10**places)
    return f"{quotient}.{remainder:0{places}d}"


if __name__ == "__main__":
    lower, upper = enclose()
    places = 11
    scale = 10**places
    left_integer = floor(lower * scale)
    right_integer = ceil(upper * scale)
    if Fraction(right_integer - left_integer, scale) > Fraction(1, 10**8):
        raise ValueError("The printed enclosure is wider than 10^-8")
    print("Exact rational enclosure:", lower, upper)
    print("Outward-rounded decimal enclosure:",
          decimal_endpoint(left_integer, places), decimal_endpoint(right_integer, places))
