from __future__ import annotations

import math

import numpy as np
from scipy.stats import norm

from coverage_coefficients import (
    METHOD_BOOTSTRAP_T,
    METHOD_PERCENTILE,
    METHOD_T,
    coefficient_lookup,
    direct_coefficients,
    finite_b_design,
    predict_coverage,
)


def test_finite_b_baseline_and_symmetry() -> None:
    for B in (399, 999):
        _, checks = direct_coefficients(0.05, B)
        assert checks["baseline_integral_abs_error"] < 1e-10
        assert checks["acceptance_symmetry_max_error"] < 1e-12
        assert abs(checks["percentile_R_minus_half_coefficient"]) < 1e-10
        assert checks["max_second_order_quad_error"] < 1e-9


def test_symbolic_component_values_at_large_B() -> None:
    rows, _ = direct_coefficients(0.05, 20001)
    values = coefficient_lookup(rows)
    z = norm.ppf(0.975)
    phi = norm.pdf(z)

    expected = {
        (METHOD_PERCENTILE, "A0"): -z * (z**2 + 3) * phi / 2,
        (METHOD_PERCENTILE, "A4"): z * (7 * z**2 - 13) * phi / 12,
        (METHOD_PERCENTILE, "A33"): -z
        * (3 * z**4 + 6 * z**2 - 11)
        * phi
        / 12,
        (METHOD_BOOTSTRAP_T, "A0"): 0.0,
        (METHOD_BOOTSTRAP_T, "A4"): -z * (2 * z**2 + 1) * phi / 3,
        (METHOD_BOOTSTRAP_T, "A33"): z * (2 * z**2 + 1) * phi / 2,
    }
    # B is large but finite, so the binomial acceptance function is still a
    # smoothed version of the limiting indicator.  The remaining discrepancy
    # is a finite-B effect, not quadrature error.
    for key, target in expected.items():
        assert math.isclose(values[key], target, rel_tol=7e-4, abs_tol=3e-5)


def test_t_normal_benchmark_is_exact_to_this_order() -> None:
    prediction = predict_coverage(
        METHOD_T,
        R=16,
        alpha=0.05,
        B=399,
        lambda3_sq=0.0,
        lambda4=0.0,
    )
    assert np.isclose(prediction, 0.95)


if __name__ == "__main__":
    test_finite_b_baseline_and_symmetry()
    test_symbolic_component_values_at_large_B()
    test_t_normal_benchmark_is_exact_to_this_order()
    print("All direct-coefficient checks passed.")
