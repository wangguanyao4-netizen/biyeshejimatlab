from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import beta

from coverage_coefficients import (
    METHOD_T,
    coefficient_lookup,
    direct_coefficients,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("source_root", type=Path)
    parser.add_argument("output_root", type=Path)
    return parser.parse_args()


def exact_interval(hits: np.ndarray, trials: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    lower = np.where(hits == 0, 0.0, beta.ppf(0.025, hits, trials - hits + 1))
    upper = np.where(
        hits == trials, 1.0, beta.ppf(0.975, hits + 1, trials - hits)
    )
    return lower, upper


def main() -> None:
    args = parse_args()
    source_root = args.source_root.resolve()
    output_root = args.output_root.resolve()
    output_root.mkdir(parents=True, exist_ok=True)

    coverage = pd.read_csv(source_root / "coverage_results.csv")
    moments = pd.read_csv(source_root / "moment_factorization.csv").set_index("model")
    if coverage["B"].nunique() != 1:
        raise ValueError("Expected one common B in the linear/nonlinear source.")
    B = int(coverage["B"].iloc[0])
    coefficient_rows, checks = direct_coefficients(0.05, B)
    coefficients = coefficient_lookup(coefficient_rows)

    records = []
    for row in coverage.itertuples(index=False):
        moment = moments.loc[row.model]
        sources = {
            "complete outer cumulants": (
                float(moment.empirical_pool_excess_kurtosis),
                float(moment.empirical_pool_skewness) ** 2,
            ),
            "fixed-weight factorization": (
                float(moment.factorized_pool_excess_kurtosis),
                float(moment.gamma_h) ** 2 * float(moment.mean_rho3_sq_w),
            ),
        }
        for source, (lambda4, lambda3_sq) in sources.items():
            A0 = coefficients[(row.method, "A0")]
            A4 = coefficients[(row.method, "A4")]
            A33 = coefficients[(row.method, "A33")]
            coefficient_A = A0 + A4 * lambda4 + A33 * lambda3_sq
            predicted = row.formula_baseline + coefficient_A / row.R
            records.append(
                {
                    **row._asdict(),
                    "moment_source": source,
                    "lambda4_input": lambda4,
                    "lambda3_sq_input": lambda3_sq,
                    "A0": A0,
                    "A4": A4,
                    "A33": A33,
                    "coefficient_A": coefficient_A,
                    "predicted_coverage": predicted,
                    "residual": row.coverage - predicted,
                }
            )

    predictions = pd.DataFrame(records)
    predictions["hit_count"] = np.rint(
        predictions["coverage"] * predictions["M"]
    ).astype(int)
    lower, upper = exact_interval(
        predictions["hit_count"].to_numpy(), predictions["M"].to_numpy()
    )
    predictions["coverage_exact95_lower"] = lower
    predictions["coverage_exact95_upper"] = upper
    predictions["prediction_inside_exact95"] = (
        (predictions["predicted_coverage"] >= lower)
        & (predictions["predicted_coverage"] <= upper)
    )
    predictions["residual_mcse_units"] = predictions["residual"] / np.maximum(
        predictions["coverage_mcse"], 1e-15
    )

    summary_rows = []
    for keys, group in predictions.groupby(
        ["model", "method", "moment_source"], sort=False
    ):
        model, method, source = keys
        residual = group["residual"].to_numpy()
        summary_rows.append(
            {
                "model": model,
                "method": method,
                "moment_source": source,
                "point_count": len(group),
                "coefficient_A": group["coefficient_A"].iloc[0],
                "mean_abs_error": np.mean(np.abs(residual)),
                "root_mean_square_error": np.sqrt(np.mean(residual**2)),
                "max_abs_error": np.max(np.abs(residual)),
                "share_prediction_inside_exact95": group[
                    "prediction_inside_exact95"
                ].mean(),
                "max_abs_residual_mcse_units": np.max(
                    np.abs(group["residual_mcse_units"])
                ),
                "max_fallback_rate": group["fallback_rate"].max(),
                "max_interval_inf_rate": group["interval_inf_rate"].max(),
            }
        )
    summary = pd.DataFrame(summary_rows)

    predictions.to_csv(output_root / "formula_predictions.csv", index=False)
    summary.to_csv(output_root / "formula_validation_summary.csv", index=False)
    pd.DataFrame(coefficient_rows).assign(B=B).to_csv(
        output_root / "direct_coefficients.csv", index=False
    )
    (output_root / "coefficient_self_checks.json").write_text(
        json.dumps(checks, indent=2), encoding="utf-8"
    )
    manifest = {
        "source_root": str(source_root),
        "calculation_policy": "direct_non_fitted",
        "coverage_target": "centered_expansion",
        "note": (
            "The t-interval uses the Student critical value. Percentile and "
            "bootstrap-t use finite-B order-statistic baselines."
        ),
    }
    (output_root / "manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )


if __name__ == "__main__":
    main()
