#!/usr/bin/env python3
"""Independent checks for the fixed Student-t density-ratio e-value."""

import csv
import math
import sys
from pathlib import Path


SHIFTS = (2.0, 4.0, 6.0)
SCALES = (2.0, 4.0, 8.0)


def logsumexp(values):
    maximum = max(values)
    return maximum + math.log(sum(math.exp(value - maximum) for value in values))


def log_t_density(value, df):
    return (
        math.lgamma(0.5 * (df + 1.0))
        - math.lgamma(0.5 * df)
        - 0.5 * math.log(df * math.pi)
        - 0.5 * (df + 1.0) * math.log1p(value * value / df)
    )


def log_evalue(value, df):
    terms = []
    for shift in SHIFTS:
        terms.append(
            logsumexp(
                (log_t_density(value - shift, df) + math.log(0.5),
                 log_t_density(value + shift, df) + math.log(0.5))
            )
        )
    for scale in SCALES:
        terms.append(log_t_density(value / scale, df) - math.log(scale))
    return logsumexp(terms) - math.log(len(terms)) - log_t_density(value, df)


def beta_fraction(a, b, x):
    qab = a + b
    qap = a + 1.0
    qam = a - 1.0
    c = 1.0
    d = 1.0 - qab * x / qap
    d = 1e-300 if abs(d) < 1e-300 else d
    d = 1.0 / d
    result = d
    for iteration in range(1, 201):
        twice = 2 * iteration
        aa = iteration * (b - iteration) * x / ((qam + twice) * (a + twice))
        d = 1.0 + aa * d
        d = 1e-300 if abs(d) < 1e-300 else d
        c = 1.0 + aa / c
        c = 1e-300 if abs(c) < 1e-300 else c
        d = 1.0 / d
        result *= d * c
        aa = -(a + iteration) * (qab + iteration) * x / (
            (a + twice) * (qap + twice)
        )
        d = 1.0 + aa * d
        d = 1e-300 if abs(d) < 1e-300 else d
        c = 1.0 + aa / c
        c = 1e-300 if abs(c) < 1e-300 else c
        d = 1.0 / d
        delta = d * c
        result *= delta
        if abs(delta - 1.0) < 3e-14:
            break
    return result


def regularized_beta(a, b, x):
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    front = math.exp(
        math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
        + a * math.log(x) + b * math.log1p(-x)
    )
    if x < (a + 1.0) / (a + b + 2.0):
        return front * beta_fraction(a, b, x) / a
    return 1.0 - front * beta_fraction(b, a, 1.0 - x) / b


def two_sided_p(value, df):
    x = df / (df + value * value)
    return regularized_beta(0.5 * df, 0.5, x)


def statistic_from_p(p_value, df):
    lower = 0.0
    upper = 2.0
    while two_sided_p(upper, df) > p_value:
        upper *= 2.0
    for _ in range(100):
        middle = 0.5 * (lower + upper)
        if two_sided_p(middle, df) > p_value:
            lower = middle
        else:
            upper = middle
    return 0.5 * (lower + upper)


def check_null_expectation(df, intervals=40000):
    width = math.pi / intervals
    total = 0.0
    for index in range(intervals):
        angle = -0.5 * math.pi + (index + 0.5) * width
        value = math.tan(angle)
        log_integrand = log_t_density(value, df) + log_evalue(value, df)
        total += math.exp(log_integrand) * (1.0 + value * value)
    return total * width


def check_output(path):
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    row = next(item for item in rows if item["Protein"] == "P_FACTOR")
    for leg, count_name in (("XM", "factor_nA"), ("MY", "factor_nB")):
        df = float(row[count_name]) - 2.0
        p_value = float(row[f"factor_p_{leg}"])
        observed = float(row[f"factor_log_e_{leg}"])
        statistic = statistic_from_p(p_value, df)
        expected = log_evalue(statistic, df)
        if abs(observed - expected) > 2e-5:
            raise SystemExit(
                f"{leg} log e-value mismatch: observed={observed}, expected={expected}"
            )


def main():
    if len(sys.argv) != 2:
        raise SystemExit("usage: check_student_t_evalue.py RESULT.mediation")
    check_output(Path(sys.argv[1]))
    for df in (1.0, 2.0, 6.0, 30.0):
        expectation = check_null_expectation(df)
        if abs(expectation - 1.0) > 5e-4:
            raise SystemExit(f"null expectation failed for df={df}: {expectation}")


if __name__ == "__main__":
    main()
