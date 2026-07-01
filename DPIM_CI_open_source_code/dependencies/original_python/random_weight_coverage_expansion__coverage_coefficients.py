"""Direct finite-B coverage coefficients for three two-sided intervals.

This module is deliberately standalone.  It does not import or call any legacy
coverage-calibration code in the parent repository.

The percentile-bootstrap and bootstrap-t coefficients are the one-dimensional
finite-B integrals obtained from the Hall-type Edgeworth polynomials for a
studentized or nonstudentized sample mean.  The t-interval coefficients use the
Student critical value, so the distribution-free constant term cancels and the
normal outer-law benchmark has exact coverage.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np
import sympy as sp
from scipy.integrate import quad
from scipy.special import ndtr
from scipy.stats import binom, norm


METHOD_T = "t distribution"
METHOD_PERCENTILE = "percentile bootstrap"
METHOD_BOOTSTRAP_T = "bootstrap-t"


@dataclass(frozen=True)
class FiniteBDesign:
    alpha: float
    B: int
    k_minus: int
    k_plus: int
    C0B: float


def finite_b_design(alpha: float, B: int) -> FiniteBDesign:
    if not 0.0 < alpha < 1.0:
        raise ValueError("alpha must lie in (0, 1).")
    if B < 2:
        raise ValueError("B must be at least 2.")
    k_minus = max(1, math.floor(alpha * (B + 1) / 2))
    k_plus = min(B, math.ceil((1 - alpha / 2) * (B + 1)))
    if k_plus <= k_minus:
        raise ValueError("Invalid finite-B endpoint ranks.")
    return FiniteBDesign(
        alpha=alpha,
        B=B,
        k_minus=k_minus,
        k_plus=k_plus,
        C0B=(k_plus - k_minus) / (B + 1),
    )


def edgeworth_polynomials() -> dict[str, sp.Expr]:
    """Return the independently specified Hall-polynomial convention.

    The CDF convention is
        Phi(z) + R^(-1/2) p1(z) phi(z) + R^(-1) p2(z) phi(z)
    for the standardized mean and the same expression with q1, q2 for the
    studentized mean.  gamma is standardized skewness and kappa is excess
    kurtosis of one outer curve ordinate.
    """

    z, gamma, kappa = sp.symbols("z gamma kappa", real=True)
    p1 = -gamma * (z**2 - 1) / 6
    p2 = -z * (
        kappa * (z**2 - 3) / 24
        + gamma**2 * (z**4 - 10 * z**2 + 15) / 72
    )
    q1 = gamma * (2 * z**2 + 1) / 6
    q2 = z * (
        kappa * (z**2 - 3) / 12
        - gamma**2 * (z**4 + 2 * z**2 - 3) / 18
        - (z**2 + 3) / 4
    )

    # The random plug-in cumulant contribution in the bootstrap expansion.
    a_percentile = (p1 / gamma) * (kappa - sp.Rational(3, 2) * gamma**2)
    a_bootstrap_t = (q1 / gamma) * (kappa - sp.Rational(3, 2) * gamma**2)

    inverse_p2 = p1 * sp.diff(p1, z) - z * p1**2 / 2 - p2
    r1_percentile = sp.expand(p1 - q1)
    r2_percentile = sp.expand(
        q2
        + inverse_p2
        - z * p1**2 / 2
        + p1 * (z * q1 - sp.diff(q1, z))
        - a_percentile * z
    )

    inverse_q2 = q1 * sp.diff(q1, z) - z * q1**2 / 2 - q2
    r2_bootstrap_t = sp.expand(
        q2
        + inverse_q2
        - z * q1**2 / 2
        + q1 * (z * q1 - sp.diff(q1, z))
        - a_bootstrap_t * z
    )

    return {
        "z": z,
        "gamma": gamma,
        "kappa": kappa,
        "p1": sp.factor(p1),
        "p2": sp.factor(p2),
        "q1": sp.factor(q1),
        "q2": sp.factor(q2),
        "r1_percentile": sp.factor(r1_percentile),
        "r2_percentile": sp.factor(r2_percentile),
        "r2_bootstrap_t": sp.factor(r2_bootstrap_t),
    }


def coefficient_polynomials() -> dict[tuple[str, str], sp.Expr]:
    """Return r_{m,j}(z) for A0, A4 and A33."""

    poly = edgeworth_polynomials()
    z = poly["z"]
    gamma = poly["gamma"]
    kappa = poly["kappa"]

    def split(expr: sp.Expr) -> tuple[sp.Expr, sp.Expr, sp.Expr]:
        expanded = sp.expand(expr)
        constant = expanded.subs({gamma: 0, kappa: 0})
        kappa_part = sp.diff(expanded, kappa)
        gamma_sq_part = sp.diff(expanded, gamma, 2) / 2
        return tuple(sp.factor(v) for v in (constant, kappa_part, gamma_sq_part))

    p0, p4, p33 = split(poly["r2_percentile"])
    bt0, bt4, bt33 = split(poly["r2_bootstrap_t"])
    return {
        (METHOD_PERCENTILE, "A0"): p0,
        (METHOD_PERCENTILE, "A4"): p4,
        (METHOD_PERCENTILE, "A33"): p33,
        (METHOD_BOOTSTRAP_T, "A0"): bt0,
        (METHOD_BOOTSTRAP_T, "A4"): bt4,
        (METHOD_BOOTSTRAP_T, "A33"): bt33,
    }


def finite_b_acceptance(z_value: float, design: FiniteBDesign) -> float:
    """P{k_minus <= Bin(B, Phi(z)) <= k_plus - 1}."""

    u = ndtr(z_value)
    return float(
        binom.cdf(design.k_plus - 1, design.B, u)
        - binom.cdf(design.k_minus - 1, design.B, u)
    )


def _finite_b_integral(
    polynomial: sp.Expr, z: sp.Symbol, design: FiniteBDesign
) -> tuple[float, float]:
    p = sp.lambdify(z, polynomial, "numpy")
    dp = sp.lambdify(z, sp.diff(polynomial, z), "numpy")

    def integrand(x: float) -> float:
        density_correction = dp(x) - x * p(x)
        return finite_b_acceptance(x, design) * density_correction * norm.pdf(x)

    value, error = quad(
        integrand,
        -9.0,
        9.0,
        epsabs=1e-13,
        epsrel=1e-11,
        limit=500,
        points=[-2.5, -1.96, 0.0, 1.96, 2.5],
    )
    return float(value), float(error)


def direct_coefficients(alpha: float, B: int) -> tuple[list[dict], dict]:
    """Calculate all coefficients without fitting coverage observations."""

    design = finite_b_design(alpha, B)
    poly = edgeworth_polynomials()
    z_symbol = poly["z"]
    component_poly = coefficient_polynomials()
    rows: list[dict] = []

    max_quad_error = 0.0
    for method in (METHOD_PERCENTILE, METHOD_BOOTSTRAP_T):
        for name in ("A0", "A4", "A33"):
            polynomial = component_poly[(method, name)]
            value, error = _finite_b_integral(polynomial, z_symbol, design)
            max_quad_error = max(max_quad_error, error)
            rows.append(
                {
                    "method": method,
                    "coefficient": name,
                    "value": value,
                    "integration_error": error,
                    "polynomial": str(polynomial),
                }
            )

    z_critical = norm.ppf(1 - alpha / 2)
    phi = norm.pdf(z_critical)
    t_values = {
        "A0": 0.0,
        "A4": 2 * z_critical * phi * (z_critical**2 - 3) / 12,
        "A33": -2
        * z_critical
        * phi
        * (z_critical**4 + 2 * z_critical**2 - 3)
        / 18,
    }
    for name, value in t_values.items():
        rows.append(
            {
                "method": METHOD_T,
                "coefficient": name,
                "value": float(value),
                "integration_error": 0.0,
                "polynomial": "closed_form_student_critical",
            }
        )

    baseline_integral, baseline_error = quad(
        lambda x: finite_b_acceptance(x, design) * norm.pdf(x),
        -9.0,
        9.0,
        epsabs=1e-13,
        epsrel=1e-11,
        limit=500,
    )
    r1_gamma = sp.factor(poly["r1_percentile"] / poly["gamma"])
    first_order, first_order_error = _finite_b_integral(
        r1_gamma, z_symbol, design
    )
    symmetry_grid = np.linspace(-7, 7, 281)
    symmetry_error = max(
        abs(
            finite_b_acceptance(x, design)
            - finite_b_acceptance(-x, design)
        )
        for x in symmetry_grid
    )
    checks = {
        "alpha": alpha,
        "B": B,
        "k_minus": design.k_minus,
        "k_plus": design.k_plus,
        "C0B": design.C0B,
        "baseline_integral": baseline_integral,
        "baseline_integral_abs_error": abs(baseline_integral - design.C0B),
        "baseline_quad_error": baseline_error,
        "acceptance_symmetry_max_error": symmetry_error,
        "percentile_R_minus_half_coefficient": first_order,
        "percentile_R_minus_half_quad_error": first_order_error,
        "max_second_order_quad_error": max_quad_error,
        "symbolic": {
            name: str(poly[name])
            for name in (
                "p1",
                "p2",
                "q1",
                "q2",
                "r1_percentile",
                "r2_percentile",
                "r2_bootstrap_t",
            )
        },
    }
    return rows, checks


def coefficient_lookup(rows: list[dict]) -> dict[tuple[str, str], float]:
    return {
        (row["method"], row["coefficient"]): float(row["value"])
        for row in rows
    }


def predict_coverage(
    method: str,
    R: int,
    alpha: float,
    B: int,
    lambda3_sq: np.ndarray | float,
    lambda4: np.ndarray | float,
) -> np.ndarray:
    rows, _ = direct_coefficients(alpha, B)
    coefficients = coefficient_lookup(rows)
    baseline = 1 - alpha if method == METHOD_T else finite_b_design(alpha, B).C0B
    A = (
        coefficients[(method, "A0")]
        + coefficients[(method, "A4")] * np.asarray(lambda4)
        + coefficients[(method, "A33")] * np.asarray(lambda3_sq)
    )
    return baseline + A / R

