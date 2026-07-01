from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import sympy as sp
from matplotlib.lines import Line2D
from scipy.integrate import quad
from scipy.special import ndtr
from scipy.stats import beta, binom, chi2, norm


METHOD_T = "t distribution"
METHOD_P = "percentile bootstrap"
METHOD_BT = "bootstrap-t"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Direct, non-fitted calculation of finite-B coverage coefficients."
    )
    parser.add_argument("results_root", type=Path)
    parser.add_argument("output_root", type=Path)
    parser.add_argument("--alpha", type=float, default=0.05)
    return parser.parse_args()


def finite_b_ranks(alpha: float, bootstrap_replicates: int) -> tuple[int, int]:
    k_minus = max(1, math.floor(alpha / 2 * (bootstrap_replicates + 1)))
    k_plus = min(
        bootstrap_replicates,
        math.ceil((1 - alpha / 2) * (bootstrap_replicates + 1)),
    )
    if k_plus <= k_minus:
        raise ValueError("Invalid finite-B endpoint ranks.")
    return k_minus, k_plus


def symbolic_polynomials() -> dict[str, sp.Expr]:
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

    # Hall's random plug-in cumulant term E(T_n Delta_n).
    a_percentile = (p1 / gamma) * (kappa - sp.Rational(3, 2) * gamma**2)
    a_bootstrap_t = (q1 / gamma) * (kappa - sp.Rational(3, 2) * gamma**2)

    s2_percentile = p1 * sp.diff(p1, z) - z * p1**2 / 2 - p2
    r1_percentile = sp.expand(p1 - q1)
    r2_percentile = sp.expand(
        q2
        + s2_percentile
        - z * p1**2 / 2
        + p1 * (z * q1 - sp.diff(q1, z))
        - a_percentile * z
    )

    s2_bootstrap_t = q1 * sp.diff(q1, z) - z * q1**2 / 2 - q2
    r2_bootstrap_t = sp.expand(
        q2
        + s2_bootstrap_t
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


def polynomial_components(expr: sp.Expr, z: sp.Symbol, gamma: sp.Symbol, kappa: sp.Symbol):
    expanded = sp.expand(expr)
    constant = expanded.subs({gamma: 0, kappa: 0})
    kappa_part = sp.diff(expanded, kappa)
    gamma_sq_part = sp.diff(expanded, gamma, 2) / 2
    return tuple(sp.factor(v) for v in (constant, kappa_part, gamma_sq_part))


def acceptance_probability(
    z_value: float, bootstrap_replicates: int, k_minus: int, k_plus: int
) -> float:
    u = ndtr(z_value)
    return float(
        binom.cdf(k_plus - 1, bootstrap_replicates, u)
        - binom.cdf(k_minus - 1, bootstrap_replicates, u)
    )


def finite_b_integral(
    polynomial: sp.Expr,
    z: sp.Symbol,
    bootstrap_replicates: int,
    k_minus: int,
    k_plus: int,
) -> tuple[float, float]:
    p = sp.lambdify(z, polynomial, "numpy")
    dp = sp.lambdify(z, sp.diff(polynomial, z), "numpy")

    def integrand(x: float) -> float:
        pi_b = acceptance_probability(x, bootstrap_replicates, k_minus, k_plus)
        return float(pi_b * (dp(x) - x * p(x)) * norm.pdf(x))

    value, error = quad(
        integrand,
        -9.0,
        9.0,
        epsabs=1e-13,
        epsrel=1e-11,
        limit=500,
        points=[-2.5, -1.96, 0.0, 1.96, 2.5],
    )
    return value, error


def calculate_coefficients(alpha: float, bootstrap_replicates: int):
    polynomials = symbolic_polynomials()
    z_symbol = polynomials["z"]
    gamma = polynomials["gamma"]
    kappa = polynomials["kappa"]
    k_minus, k_plus = finite_b_ranks(alpha, bootstrap_replicates)
    c0b = (k_plus - k_minus) / (bootstrap_replicates + 1)

    p_components = polynomial_components(
        polynomials["r2_percentile"], z_symbol, gamma, kappa
    )
    bt_components = polynomial_components(
        polynomials["r2_bootstrap_t"], z_symbol, gamma, kappa
    )

    rows = []
    integration_errors = {}
    for method, components in (
        (METHOD_P, p_components),
        (METHOD_BT, bt_components),
    ):
        for label, polynomial in zip(("A0", "A4", "A33"), components, strict=True):
            value, error = finite_b_integral(
                polynomial, z_symbol, bootstrap_replicates, k_minus, k_plus
            )
            rows.append(
                {
                    "method": method,
                    "coefficient": label,
                    "value": value,
                    "integration_error": error,
                    "polynomial": str(polynomial),
                }
            )
            integration_errors[f"{method}:{label}"] = error

    z_critical = norm.ppf(1 - alpha / 2)
    phi = norm.pdf(z_critical)
    t_a4 = 2 * z_critical * phi * (z_critical**2 - 3) / 12
    t_a33 = -2 * z_critical * phi * (
        z_critical**4 + 2 * z_critical**2 - 3
    ) / 18
    rows.extend(
        [
            {
                "method": METHOD_T,
                "coefficient": "A0",
                "value": 0.0,
                "integration_error": 0.0,
                "polynomial": "0",
            },
            {
                "method": METHOD_T,
                "coefficient": "A4",
                "value": t_a4,
                "integration_error": 0.0,
                "polynomial": "closed_form",
            },
            {
                "method": METHOD_T,
                "coefficient": "A33",
                "value": t_a33,
                "integration_error": 0.0,
                "polynomial": "closed_form",
            },
        ]
    )

    r1_gamma = sp.factor(polynomials["r1_percentile"] / gamma)
    a1_value, a1_error = finite_b_integral(
        r1_gamma, z_symbol, bootstrap_replicates, k_minus, k_plus
    )

    baseline_value, baseline_error = quad(
        lambda x: acceptance_probability(
            x, bootstrap_replicates, k_minus, k_plus
        )
        * norm.pdf(x),
        -9.0,
        9.0,
        epsabs=1e-13,
        epsrel=1e-11,
        limit=500,
    )
    symmetry_grid = np.linspace(-7, 7, 281)
    symmetry_error = max(
        abs(
            acceptance_probability(x, bootstrap_replicates, k_minus, k_plus)
            - acceptance_probability(-x, bootstrap_replicates, k_minus, k_plus)
        )
        for x in symmetry_grid
    )

    checks = {
        "B": bootstrap_replicates,
        "alpha": alpha,
        "k_minus": k_minus,
        "k_plus": k_plus,
        "C0B": c0b,
        "baseline_integral": baseline_value,
        "baseline_integral_abs_error": abs(baseline_value - c0b),
        "baseline_quad_error": baseline_error,
        "acceptance_symmetry_max_error": symmetry_error,
        "percentile_R_minus_half_coefficient": a1_value,
        "percentile_R_minus_half_quad_error": a1_error,
        "max_second_order_quad_error": max(integration_errors.values()),
        "symbolic": {
            name: str(polynomials[name])
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
    return pd.DataFrame(rows), checks


def coefficient_map(coefficients: pd.DataFrame) -> dict[tuple[str, str], float]:
    return {
        (row.method, row.coefficient): float(row.value)
        for row in coefficients.itertuples(index=False)
    }


def build_predictions(
    coverage: pd.DataFrame,
    moments: pd.DataFrame,
    coefficients: pd.DataFrame,
) -> pd.DataFrame:
    cmap = coefficient_map(coefficients)
    records = []
    moment_by_model = moments.set_index("model")
    for row in coverage.itertuples(index=False):
        moment = moment_by_model.loc[row.model]
        substitutions = {
            "fixed_weight_factorization": (
                float(moment.kappa_h * moment.mean_rho4_w),
                float(moment.gamma_h**2 * moment.mean_rho3_sq_w),
            ),
            "empirical_outer_law_sensitivity": (
                float(moment.empirical_pool_excess_kurtosis),
                float(moment.empirical_pool_skewness**2),
            ),
        }
        for moment_source, (kappa_effective, gamma_sq_effective) in substitutions.items():
            a0 = cmap[(row.method, "A0")]
            a4 = cmap[(row.method, "A4")]
            a33 = cmap[(row.method, "A33")]
            coefficient_a = a0 + a4 * kappa_effective + a33 * gamma_sq_effective
            predicted = float(row.formula_baseline + coefficient_a / row.R)
            residual = float(row.coverage - predicted)
            records.append(
                {
                    "model": row.model,
                    "y0": row.y0,
                    "h": row.h,
                    "method": row.method,
                    "R": int(row.R),
                    "M": int(row.M),
                    "B": int(row.B),
                    "target_mode": row.target_mode,
                    "moment_source": moment_source,
                    "baseline": float(row.formula_baseline),
                    "A0": a0,
                    "A4": a4,
                    "A33": a33,
                    "kappa_effective": kappa_effective,
                    "gamma_sq_effective": gamma_sq_effective,
                    "coefficient_A": coefficient_a,
                    "observed_coverage": float(row.coverage),
                    "coverage_mcse": float(row.coverage_mcse),
                    "predicted_coverage": predicted,
                    "residual": residual,
                    "residual_mcse_units": residual / max(float(row.coverage_mcse), 1e-15),
                    "fallback_rate": float(row.fallback_rate),
                    "interval_inf_rate": float(row.interval_inf_rate),
                }
            )
    return pd.DataFrame.from_records(records)


def summarize_predictions(predictions: pd.DataFrame) -> pd.DataFrame:
    records = []
    grouping = predictions.groupby(["model", "method", "moment_source"], sort=False)
    for keys, group in grouping:
        model, method, moment_source = keys
        residual = group["residual"].to_numpy()
        z_residual = group["residual_mcse_units"].to_numpy()
        records.append(
            {
                "model": model,
                "method": method,
                "moment_source": moment_source,
                "coefficient_A": group["coefficient_A"].iloc[0],
                "mean_abs_error": np.mean(np.abs(residual)),
                "root_mean_square_error": np.sqrt(np.mean(residual**2)),
                "max_abs_error": np.max(np.abs(residual)),
                "max_abs_residual_mcse_units": np.max(np.abs(z_residual)),
                "points_within_1_96_mcse": np.mean(np.abs(z_residual) <= 1.96),
                "max_fallback_rate": group["fallback_rate"].max(),
                "max_interval_inf_rate": group["interval_inf_rate"].max(),
            }
        )
    return pd.DataFrame.from_records(records)


def holm_adjusted_p_values(p_values: np.ndarray) -> np.ndarray:
    order = np.argsort(p_values)
    adjusted = np.empty_like(p_values, dtype=float)
    running_max = 0.0
    count = len(p_values)
    for rank, index in enumerate(order):
        candidate = min(1.0, (count - rank) * float(p_values[index]))
        running_max = max(running_max, candidate)
        adjusted[index] = running_max
    return adjusted


def build_direct_validation(predictions: pd.DataFrame, alpha: float = 0.05) -> pd.DataFrame:
    direct = predictions[
        predictions["moment_source"] == "fixed_weight_factorization"
    ].copy()
    direct["hits"] = np.rint(
        direct["observed_coverage"] * direct["M"]
    ).astype(int)
    direct["null_mcse"] = np.sqrt(
        direct["predicted_coverage"]
        * (1 - direct["predicted_coverage"])
        / direct["M"]
    )
    direct["z_null"] = direct["residual"] / direct["null_mcse"]
    direct["pointwise_p_value"] = 2 * norm.sf(np.abs(direct["z_null"]))

    lower = []
    upper = []
    for hits, trials in zip(direct["hits"], direct["M"], strict=True):
        lower.append(
            0.0
            if hits == 0
            else float(beta.ppf(alpha / 2, hits, trials - hits + 1))
        )
        upper.append(
            1.0
            if hits == trials
            else float(beta.ppf(1 - alpha / 2, hits + 1, trials - hits))
        )
    direct["exact_95_lower"] = lower
    direct["exact_95_upper"] = upper
    direct["direct_inside_exact_95"] = (
        (direct["predicted_coverage"] >= direct["exact_95_lower"])
        & (direct["predicted_coverage"] <= direct["exact_95_upper"])
    )
    direct["holm_adjusted_p_all_48"] = holm_adjusted_p_values(
        direct["pointwise_p_value"].to_numpy()
    )
    return direct.sort_values(["model", "method", "R"]).reset_index(drop=True)


def summarize_direct_validation(direct: pd.DataFrame) -> pd.DataFrame:
    records = []
    for (model, method), group in direct.groupby(["model", "method"], sort=False):
        pearson_q = float(np.sum(group["z_null"] ** 2))
        degrees_of_freedom = int(len(group))
        records.append(
            {
                "model": model,
                "method": method,
                "coefficient_A": group["coefficient_A"].iloc[0],
                "point_count": degrees_of_freedom,
                "mean_abs_error": float(np.mean(np.abs(group["residual"]))),
                "root_mean_square_error": float(
                    np.sqrt(np.mean(group["residual"] ** 2))
                ),
                "max_abs_error": float(np.max(np.abs(group["residual"]))),
                "max_abs_z_null": float(np.max(np.abs(group["z_null"]))),
                "all_direct_inside_exact_95": bool(
                    group["direct_inside_exact_95"].all()
                ),
                "minimum_pointwise_p_value": float(
                    group["pointwise_p_value"].min()
                ),
                "minimum_holm_adjusted_p_all_48": float(
                    group["holm_adjusted_p_all_48"].min()
                ),
                "pearson_Q_zero_fit": pearson_q,
                "pearson_df": degrees_of_freedom,
                "pearson_p_value_zero_fit": float(
                    chi2.sf(pearson_q, degrees_of_freedom)
                ),
            }
        )
    return pd.DataFrame.from_records(records)


def summarize_uncentered(uncentered: pd.DataFrame) -> pd.DataFrame:
    return (
        uncentered.groupby(["model", "method"], sort=False)
        .agg(
            min_coverage=("coverage", "min"),
            max_coverage=("coverage", "max"),
            mean_coverage=("coverage", "mean"),
            max_abs_error_to_baseline=("coverage_error", lambda x: np.max(np.abs(x))),
        )
        .reset_index()
    )


def use_validation_chart_theme() -> None:
    sns.set_theme(
        style="whitegrid",
        rc={
            "figure.facecolor": "#FCFCFD",
            "axes.facecolor": "#FFFFFF",
            "axes.edgecolor": "#D7DBE7",
            "axes.labelcolor": "#1F2430",
            "axes.titlecolor": "#1F2430",
            "grid.color": "#E6E8F0",
            "grid.linewidth": 0.8,
            "font.family": "sans-serif",
            "font.sans-serif": [
                "Aptos",
                "Segoe UI",
                "DejaVu Sans",
                "Arial",
                "sans-serif",
            ],
        },
    )


def plot_observed_vs_direct(direct: pd.DataFrame, output_root: Path) -> None:
    use_validation_chart_theme()
    models = list(direct["model"].drop_duplicates())
    methods = [METHOD_T, METHOD_P, METHOD_BT]
    titles = {
        METHOD_T: "t distribution",
        METHOD_P: "percentile bootstrap",
        METHOD_BT: "bootstrap-t",
    }
    fig, axes = plt.subplots(
        len(models), len(methods), figsize=(13.2, 7.2), sharex=True, sharey=True
    )
    if len(models) == 1:
        axes = np.asarray([axes])

    for i, model in enumerate(models):
        for j, method in enumerate(methods):
            ax = axes[i, j]
            group = direct[
                (direct["model"] == model) & (direct["method"] == method)
            ].sort_values("R")
            observed = group["observed_coverage"].to_numpy()
            lower_error = observed - group["exact_95_lower"].to_numpy()
            upper_error = group["exact_95_upper"].to_numpy() - observed
            ax.errorbar(
                group["R"],
                observed,
                yerr=np.vstack([lower_error, upper_error]),
                fmt="o",
                ms=4.7,
                color="#1F2430",
                ecolor="#7A828F",
                elinewidth=1.0,
                capsize=2.5,
                label="experiment: exact 95% binomial interval",
                zorder=3,
            )
            sns.lineplot(
                data=group,
                x="R",
                y="predicted_coverage",
                ax=ax,
                color="#5477C4",
                marker="s",
                markersize=4.0,
                linewidth=1.3,
                label="direct formula (no fit)",
                zorder=2,
            )
            ax.axhline(
                group["baseline"].iloc[0],
                color="#464C55",
                linestyle=":",
                linewidth=1.0,
                label="formula baseline",
            )
            if ax.legend_ is not None:
                ax.legend_.remove()
            ax.set_xscale("log", base=2)
            ax.set_xticks(group["R"], [str(value) for value in group["R"]])
            ax.set_ylim(0.90, 0.968)
            ax.set_title(titles[method], fontsize=10.5)
            if j == 0:
                ax.set_ylabel(f"{model}\ncoverage probability")
            else:
                ax.set_ylabel("")
            if i == len(models) - 1:
                ax.set_xlabel("outer replicate count R")
            else:
                ax.set_xlabel("")
            sns.despine(ax=ax)

    legend_handles = [
        Line2D(
            [0],
            [0],
            color="#5477C4",
            marker="s",
            markersize=5,
            linewidth=1.3,
            label="direct formula (no fit)",
        ),
        Line2D(
            [0],
            [0],
            color="#464C55",
            linestyle=":",
            linewidth=1.0,
            label="formula baseline",
        ),
        Line2D(
            [0],
            [0],
            color="#7A828F",
            marker="o",
            markerfacecolor="#1F2430",
            markeredgecolor="#1F2430",
            markersize=5,
            linewidth=1.0,
            label="experiment: exact 95% binomial interval",
        ),
    ]
    fig.legend(
        handles=legend_handles,
        loc="upper center",
        ncol=3,
        frameon=False,
        bbox_to_anchor=(0.5, 0.955),
    )
    fig.suptitle(
        "Experimental coverage versus directly calculated coverage",
        fontsize=14,
        fontweight="semibold",
        color="#1F2430",
        y=0.995,
    )
    fig.text(
        0.5,
        0.965,
        "Centered probability-weighted DPIM; B=999, M=1200; "
        "intervals are exact Clopper-Pearson intervals for experiment counts",
        ha="center",
        va="top",
        fontsize=9,
        color="#6F768A",
    )
    fig.tight_layout(rect=(0, 0, 1, 0.89))
    for extension in ("png", "pdf", "svg"):
        fig.savefig(
            output_root / f"observed_vs_direct_coverage.{extension}",
            dpi=320,
            bbox_inches="tight",
        )
    plt.close(fig)


def plot_standardized_residuals(direct: pd.DataFrame, output_root: Path) -> None:
    use_validation_chart_theme()
    models = list(direct["model"].drop_duplicates())
    methods = [METHOD_T, METHOD_P, METHOD_BT]
    titles = {
        METHOD_T: "t distribution",
        METHOD_P: "percentile bootstrap",
        METHOD_BT: "bootstrap-t",
    }
    fig, axes = plt.subplots(
        len(models), len(methods), figsize=(13.2, 6.5), sharex=True, sharey=True
    )
    if len(models) == 1:
        axes = np.asarray([axes])
    for i, model in enumerate(models):
        for j, method in enumerate(methods):
            ax = axes[i, j]
            group = direct[
                (direct["model"] == model) & (direct["method"] == method)
            ].sort_values("R")
            sns.lineplot(
                data=group,
                x="R",
                y="z_null",
                ax=ax,
                color="#5477C4",
                marker="o",
                markersize=4.5,
                linewidth=1.0,
            )
            ax.axhline(0.0, color="#464C55", linewidth=0.9)
            ax.axhline(1.96, color="#CC6F47", linestyle="--", linewidth=1.0)
            ax.axhline(-1.96, color="#CC6F47", linestyle="--", linewidth=1.0)
            ax.set_xscale("log", base=2)
            ax.set_xticks(group["R"], [str(value) for value in group["R"]])
            ax.set_ylim(-2.35, 2.35)
            ax.set_title(titles[method], fontsize=10.5)
            if j == 0:
                ax.set_ylabel(f"{model}\nstandardized residual")
            else:
                ax.set_ylabel("")
            if i == len(models) - 1:
                ax.set_xlabel("outer replicate count R")
            else:
                ax.set_xlabel("")
            sns.despine(ax=ax)
    fig.suptitle(
        "No-fit standardized residuals remain within pointwise 95% limits",
        fontsize=14,
        fontweight="semibold",
        color="#1F2430",
        y=0.99,
    )
    fig.text(
        0.5,
        0.955,
        r"Residual = (experimental coverage - direct coverage) / "
        r"sqrt(C_direct(1-C_direct)/M); orange lines are +/-1.96",
        ha="center",
        va="top",
        fontsize=9,
        color="#6F768A",
    )
    fig.tight_layout(rect=(0, 0, 1, 0.89))
    for extension in ("png", "pdf", "svg"):
        fig.savefig(
            output_root / f"standardized_direct_residuals.{extension}",
            dpi=320,
            bbox_inches="tight",
        )
    plt.close(fig)


def plot_direct_comparison(predictions: pd.DataFrame, output_root: Path) -> None:
    models = list(predictions["model"].drop_duplicates())
    methods = [METHOD_T, METHOD_P, METHOD_BT]
    titles = {
        METHOD_T: "t distribution",
        METHOD_P: "percentile bootstrap",
        METHOD_BT: "bootstrap-t",
    }
    colors = {
        "fixed_weight_factorization": "#2667A8",
        "empirical_outer_law_sensitivity": "#D97706",
    }
    fig, axes = plt.subplots(
        len(models), len(methods), figsize=(13.2, 7.0), sharex=True, sharey=True
    )
    if len(models) == 1:
        axes = np.asarray([axes])

    for i, model in enumerate(models):
        for j, method in enumerate(methods):
            ax = axes[i, j]
            group = predictions[
                (predictions["model"] == model) & (predictions["method"] == method)
            ]
            observed = group[
                group["moment_source"] == "fixed_weight_factorization"
            ].sort_values("R")
            x = 1.0 / observed["R"].to_numpy()
            y = observed["observed_coverage"].to_numpy() - observed["baseline"].to_numpy()
            yerr = 1.96 * observed["coverage_mcse"].to_numpy()
            ax.errorbar(
                x,
                y,
                yerr=yerr,
                fmt="o",
                ms=4.5,
                color="#202124",
                ecolor="#9AA0A6",
                capsize=2,
                label="observed (95% MC interval)",
                zorder=3,
            )
            for source, label, linestyle in (
                ("fixed_weight_factorization", "direct: weighted formula", "-"),
                (
                    "empirical_outer_law_sensitivity",
                    "direct: empirical-cumulant sensitivity",
                    "--",
                ),
            ):
                source_group = group[group["moment_source"] == source].sort_values("R")
                ax.plot(
                    1.0 / source_group["R"],
                    source_group["predicted_coverage"] - source_group["baseline"],
                    linestyle,
                    lw=1.8,
                    color=colors[source],
                    label=label,
                )
            ax.axhline(0.0, color="#6B7280", lw=0.9)
            ax.grid(axis="y", color="#E5E7EB", lw=0.7)
            ax.set_title(titles[method], fontsize=10.5)
            if j == 0:
                ax.set_ylabel(f"{model}\ncoverage error")
            if i == len(models) - 1:
                ax.set_xlabel(r"$1/R$")
            ax.tick_params(labelsize=8.5)

    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper center",
        ncol=3,
        frameon=False,
        bbox_to_anchor=(0.5, 1.01),
    )
    fig.suptitle(
        "Direct second-order coverage prediction versus centered experiment",
        fontsize=13,
        y=1.045,
    )
    fig.tight_layout()
    for extension in ("png", "pdf"):
        fig.savefig(
            output_root / f"direct_coverage_comparison.{extension}",
            dpi=320,
            bbox_inches="tight",
        )
    plt.close(fig)


def plot_coefficient_decomposition(
    predictions: pd.DataFrame, output_root: Path
) -> None:
    fixed = predictions[
        predictions["moment_source"] == "fixed_weight_factorization"
    ].drop_duplicates(["model", "method"])
    fixed = fixed.copy()
    fixed["constant_term"] = fixed["A0"]
    fixed["kurtosis_term"] = fixed["A4"] * fixed["kappa_effective"]
    fixed["skewness_sq_term"] = fixed["A33"] * fixed["gamma_sq_effective"]

    models = list(fixed["model"].drop_duplicates())
    fig, axes = plt.subplots(1, len(models), figsize=(11.2, 4.4), sharey=True)
    if len(models) == 1:
        axes = [axes]
    palette = ["#6B7280", "#2667A8", "#D97706"]
    for ax, model in zip(axes, models, strict=True):
        group = fixed[fixed["model"] == model].set_index("method").loc[
            [METHOD_T, METHOD_P, METHOD_BT]
        ]
        x = np.arange(3)
        width = 0.21
        for offset, column, label, color in zip(
            (-width, 0.0, width),
            ("constant_term", "kurtosis_term", "skewness_sq_term"),
            (r"$A^{(0)}$", r"$\kappa_{\rm eff} A^{(4)}$", r"$\gamma_{\rm eff}^2 A^{(33)}$"),
            palette,
            strict=True,
        ):
            values = group[column].to_numpy()
            ax.bar(
                x + offset,
                values,
                width=width,
                color=color,
                edgecolor="white",
                linewidth=0.6,
                label=label,
            )
        ax.scatter(
            x,
            group["coefficient_A"],
            marker="D",
            s=34,
            color="#202124",
            label=r"total $A_m$",
            zorder=4,
        )
        ax.axhline(0.0, color="#202124", lw=0.9)
        ax.set_xticks(x, ["t", "percentile", "bootstrap-t"])
        ax.set_title(model)
        ax.grid(axis="y", color="#E5E7EB", lw=0.7)
        ax.set_ylabel(r"direct coefficient $A_m$")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper center",
        ncol=4,
        frameon=False,
        bbox_to_anchor=(0.5, 0.94),
    )
    fig.suptitle("Probability-weighted direct coefficient decomposition", y=1.01)
    fig.tight_layout(rect=(0, 0, 1, 0.86))
    for extension in ("png", "pdf"):
        fig.savefig(
            output_root / f"direct_coefficient_decomposition.{extension}",
            dpi=320,
            bbox_inches="tight",
        )
    plt.close(fig)


def main() -> None:
    args = parse_args()
    results_root = args.results_root.resolve()
    output_root = args.output_root.resolve()
    output_root.mkdir(parents=True, exist_ok=True)

    coverage = pd.read_csv(results_root / "coverage_results.csv")
    moments = pd.read_csv(results_root / "moment_factorization.csv")
    uncentered = pd.read_csv(results_root / "coverage_results_uncentered.csv")
    bootstrap_replicates = int(coverage["B"].iloc[0])
    if coverage["B"].nunique() != 1:
        raise ValueError("This document expects one common B value.")

    coefficients, checks = calculate_coefficients(args.alpha, bootstrap_replicates)
    if checks["baseline_integral_abs_error"] > 1e-10:
        raise RuntimeError("Finite-B baseline integral failed.")
    if checks["acceptance_symmetry_max_error"] > 1e-12:
        raise RuntimeError("Finite-B acceptance polynomial is not symmetric.")
    if abs(checks["percentile_R_minus_half_coefficient"]) > 1e-10:
        raise RuntimeError("Two-sided percentile R^(-1/2) cancellation failed.")
    if checks["max_second_order_quad_error"] > 1e-9:
        raise RuntimeError("Finite-B coefficient quadrature did not converge.")

    predictions = build_predictions(coverage, moments, coefficients)
    summary = summarize_predictions(predictions)
    direct_validation = build_direct_validation(predictions, alpha=args.alpha)
    direct_validation_summary = summarize_direct_validation(direct_validation)
    uncentered_summary = summarize_uncentered(uncentered)

    coefficients.to_csv(output_root / "direct_coefficients.csv", index=False)
    predictions.to_csv(output_root / "direct_predictions.csv", index=False)
    summary.to_csv(output_root / "direct_prediction_summary.csv", index=False)
    direct_validation.to_csv(
        output_root / "observed_vs_direct_validation.csv", index=False
    )
    direct_validation_summary.to_csv(
        output_root / "observed_vs_direct_validation_summary.csv", index=False
    )
    uncentered_summary.to_csv(
        output_root / "uncentered_coverage_summary.csv", index=False
    )
    moments.to_csv(output_root / "moment_inputs.csv", index=False)
    with (output_root / "derivation_self_checks.json").open(
        "w", encoding="utf-8"
    ) as handle:
        json.dump(checks, handle, ensure_ascii=False, indent=2)

    plot_direct_comparison(predictions, output_root)
    plot_observed_vs_direct(direct_validation, output_root)
    plot_standardized_residuals(direct_validation, output_root)
    plot_coefficient_decomposition(predictions, output_root)

    manifest = {
        "source_results_root": str(results_root),
        "output_root": str(output_root),
        "calculation_policy": "direct_non_fitted",
        "coverage_source": "coverage_results.csv (centered_expansion)",
        "primary_validation": (
            "experiment hit proportions versus fixed-weight probability-weighted "
            "direct formula; exact binomial intervals and zero-fit Pearson tests"
        ),
        "sensitivity_source": "empirical outer-law cumulants, not fitted coefficients",
        "files": sorted(path.name for path in output_root.iterdir()),
    }
    with (output_root / "manifest.json").open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, ensure_ascii=False, indent=2)


if __name__ == "__main__":
    main()
