from __future__ import annotations

import argparse
import json
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import beta

from coverage_coefficients import (
    METHOD_BOOTSTRAP_T,
    METHOD_PERCENTILE,
    METHOD_T,
    coefficient_lookup,
    direct_coefficients,
)


METHOD_ORDER = [METHOD_T, METHOD_PERCENTILE, METHOD_BOOTSTRAP_T]
METHOD_LABEL = {
    METHOD_T: "t distribution",
    METHOD_PERCENTILE: "percentile bootstrap",
    METHOD_BOOTSTRAP_T: "bootstrap-t",
}
COLORS = {
    "observed": "#222831",
    "formula": "#2F5D9E",
    "baseline": "#6F7785",
    "band": "#E8ECF2",
    "grid": "#E3E7EF",
    "ink": "#222831",
    "threshold_light": "#D8E4F2",
    "threshold_mid": "#7FA3C7",
    "threshold_dark": "#2F5D9E",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Recompute Kirchhoff plate coverage predictions using the original "
            "paper fixed-weight factorized inputs only."
        )
    )
    parser.add_argument(
        "--source-root",
        type=Path,
        default=None,
        help="Existing plate_outer_cumulant_formal_* directory.",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=None,
        help="Directory for original-formula CSVs, figures, and change log.",
    )
    return parser.parse_args()


def latest_plate_formal_root(package_root: Path) -> Path:
    candidates = sorted(
        (package_root / "results").glob("plate_outer_cumulant_formal_*"),
        key=lambda path: path.stat().st_mtime,
        reverse=True,
    )
    if not candidates:
        raise FileNotFoundError("No plate_outer_cumulant_formal_* result found.")
    return candidates[0]


def exact_binomial_interval(hits: np.ndarray, trials: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    lower = np.where(
        hits == 0,
        0.0,
        beta.ppf(0.025, hits, trials - hits + 1),
    )
    upper = np.where(
        hits == trials,
        1.0,
        beta.ppf(0.975, hits + 1, trials - hits),
    )
    return lower, upper


def build_original_predictions(
    source_root: Path,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, dict]:
    coverage = pd.read_csv(source_root / "coverage_results.csv")
    moments = pd.read_csv(source_root / "outer_moments_active_grid.csv")

    coefficient_rows: list[dict] = []
    coefficient_checks: dict[str, dict] = {}
    coefficient_maps: dict[int, dict] = {}
    for B in sorted(coverage["B"].unique()):
        rows, checks = direct_coefficients(0.05, int(B))
        coefficient_rows.extend({"B": int(B), **row} for row in rows)
        coefficient_checks[str(int(B))] = checks
        coefficient_maps[int(B)] = coefficient_lookup(rows)

    merged = coverage.merge(
        moments,
        on="response_m",
        how="left",
        validate="many_to_one",
    )

    records: list[dict] = []
    for row in merged.itertuples(index=False):
        coefficients = coefficient_maps[int(row.B)]
        A0 = coefficients[(row.method, "A0")]
        A4 = coefficients[(row.method, "A4")]
        A33 = coefficients[(row.method, "A33")]

        # Conditional fixed-probability-weight factorization:
        # lambda_3^2 = gamma_h^2 * rho_3(w)^2,
        # lambda_4   = kappa_h   * rho_4(w).
        # The plate run uses independent randomized weights per curve, so the
        # reproducible deterministic input is the pilot average of rho moments.
        lambda3_sq = row.factorized_lambda3_sq
        lambda4 = row.factorized_lambda4
        coefficient_A = A0 + A4 * lambda4 + A33 * lambda3_sq
        predicted = row.formula_baseline + coefficient_A / row.R

        record = row._asdict()
        record.update(
            {
                "formula_policy": "original_fixed_weight_factorization",
                "A0": A0,
                "A4": A4,
                "A33": A33,
                "lambda3_sq_input": lambda3_sq,
                "lambda4_input": lambda4,
                "coefficient_A": coefficient_A,
                "predicted_coverage": predicted,
                "residual": row.coverage - predicted,
            }
        )
        records.append(record)

    predictions = pd.DataFrame.from_records(records)
    predictions["relative_density"] = (
        predictions["reference_density_per_m"]
        / predictions["reference_density_per_m"].max()
    )
    lower, upper = exact_binomial_interval(
        predictions["hit_count"].to_numpy(),
        predictions["M"].to_numpy(),
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
    predictions["abs_lambda3_over_sqrt_R"] = (
        np.sqrt(np.maximum(predictions["lambda3_sq_input"], 0.0))
        / np.sqrt(predictions["R"])
    )
    predictions["abs_lambda4_over_R"] = np.abs(predictions["lambda4_input"]) / predictions["R"]

    summary = summarize(predictions, ["method", "R", "B"])
    region_rows = []
    for threshold in (0.01, 0.05, 0.10):
        region = predictions[predictions["relative_density"] >= threshold]
        group_summary = summarize(region, ["method", "R", "B"])
        group_summary.insert(0, "relative_density_threshold", threshold)
        region_rows.append(group_summary)
    region_summary = pd.concat(region_rows, ignore_index=True)
    coefficient_frame = pd.DataFrame(coefficient_rows)

    return predictions, summary, region_summary, coefficient_frame, coefficient_checks


def summarize(predictions: pd.DataFrame, keys: list[str]) -> pd.DataFrame:
    rows = []
    for values, group in predictions.groupby(keys, sort=False):
        if not isinstance(values, tuple):
            values = (values,)
        residual = group["residual"].to_numpy()
        row = dict(zip(keys, values, strict=True))
        row.update(
            {
                "grid_point_count": len(group),
                "mean_observed_coverage": group["coverage"].mean(),
                "mean_predicted_coverage": group["predicted_coverage"].mean(),
                "mean_signed_error": residual.mean(),
                "mean_abs_error": np.mean(np.abs(residual)),
                "root_mean_square_error": np.sqrt(np.mean(residual**2)),
                "max_abs_error": np.max(np.abs(residual)),
                "share_prediction_inside_exact95": group[
                    "prediction_inside_exact95"
                ].mean(),
                "max_abs_residual_mcse_units": np.max(
                    np.abs(group["residual_mcse_units"])
                ),
                "max_abs_lambda3_over_sqrt_R": group[
                    "abs_lambda3_over_sqrt_R"
                ].max(),
                "max_abs_lambda4_over_R": group["abs_lambda4_over_R"].max(),
                "max_bootstrap_t_diagnostic_rate": group[
                    "bootstrap_t_diagnostic_rate"
                ].max(),
                "max_interval_inf_rate": group["interval_inf_rate"].max(),
            }
        )
        rows.append(row)
    return pd.DataFrame(rows)


def set_theme() -> None:
    sns.set_theme(
        style="whitegrid",
        rc={
            "figure.facecolor": "#FFFFFF",
            "savefig.facecolor": "#FFFFFF",
            "axes.facecolor": "#FFFFFF",
            "axes.edgecolor": COLORS["ink"],
            "axes.labelcolor": COLORS["ink"],
            "axes.titlecolor": COLORS["ink"],
            "axes.linewidth": 0.85,
            "text.color": COLORS["ink"],
            "grid.color": COLORS["grid"],
            "grid.linewidth": 0.65,
            "font.family": "serif",
            "font.serif": ["Times New Roman", "STIXGeneral", "DejaVu Serif"],
            "mathtext.fontset": "stix",
            "pdf.fonttype": 42,
        },
    )


def save_figure(fig: plt.Figure, output_root: Path, name: str) -> None:
    for extension in ("pdf", "png"):
        fig.savefig(
            output_root / f"{name}.{extension}",
            dpi=320,
            bbox_inches="tight",
            pad_inches=0.035,
        )
    plt.close(fig)


def plot_formula_inputs(predictions: pd.DataFrame, output_root: Path) -> None:
    data = predictions.drop_duplicates("response_m").sort_values("response_m")
    set_theme()
    fig, axes = plt.subplots(2, 2, figsize=(11.0, 7.0), sharex=True)
    series = [
        ("kernel_gamma", r"kernel skewness $\gamma_h(y)$"),
        ("kernel_excess_kurtosis", r"kernel excess kurtosis $\kappa_h(y)$"),
        ("lambda3_sq_input", r"$\gamma_h^2(y)\rho_3^2(w)$"),
        ("lambda4_input", r"$\kappa_h(y)\rho_4(w)$"),
    ]
    for ax, (column, label) in zip(axes.ravel(), series, strict=True):
        ax.plot(data["response_m"], data[column], color=COLORS["formula"], linewidth=1.2)
        ax.axhline(0, color=COLORS["baseline"], linewidth=0.8)
        ax.set_ylabel(label)
        ax.ticklabel_format(axis="x", style="sci", scilimits=(-3, 3))
        sns.despine(ax=ax)
    axes[1, 0].set_xlabel("plate-center deflection response (m)")
    axes[1, 1].set_xlabel("plate-center deflection response (m)")
    fig.tight_layout()
    save_figure(fig, output_root, "original_formula_inputs")


def plot_coverage_profiles(predictions: pd.DataFrame, output_root: Path) -> None:
    B = int(predictions["B"].max())
    r_values = [int(predictions["R"].min()), int(predictions["R"].max())]
    set_theme()
    fig, axes = plt.subplots(3, 2, figsize=(7.35, 5.25), sharex=True, sharey=True)
    fig.subplots_adjust(
        left=0.09, right=0.99, bottom=0.11, top=0.88, hspace=0.34, wspace=0.15
    )
    for i, method in enumerate(METHOD_ORDER):
        for j, R in enumerate(r_values):
            ax = axes[i, j]
            group = predictions[
                (predictions["method"] == method)
                & (predictions["R"] == R)
                & (predictions["B"] == B)
            ].sort_values("response_m")
            ax.plot(
                group["response_m"],
                group["coverage"],
                color=COLORS["observed"],
                linewidth=0.95,
                label="centered experiment",
                zorder=3,
            )
            ax.fill_between(
                group["response_m"],
                group["coverage_exact95_lower"],
                group["coverage_exact95_upper"],
                color=COLORS["band"],
                alpha=0.65,
                linewidth=0,
                label="exact 95% interval",
                zorder=1,
            )
            ax.plot(
                group["response_m"],
                group["predicted_coverage"],
                color=COLORS["formula"],
                linewidth=1.12,
                label="factorized formula",
                zorder=4,
            )
            ax.axhline(
                group["formula_baseline"].iloc[0],
                color=COLORS["baseline"],
                linestyle=":",
                linewidth=0.9,
            )
            ax.set_title(f"{METHOD_LABEL[method]}, R={R}, B={B}", fontsize=7.8)
            ax.ticklabel_format(axis="x", style="sci", scilimits=(-3, 3))
            if j == 0:
                ax.set_ylabel("pointwise coverage")
            if i == 2:
                ax.set_xlabel("plate-center deflection response (m)")
            sns.despine(ax=ax)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper center",
        ncol=3,
        frameon=False,
        bbox_to_anchor=(0.5, 0.995),
        fontsize=7.6,
        handlelength=2.4,
        columnspacing=1.0,
    )
    save_figure(fig, output_root, "original_formula_coverage_profiles")


def plot_region_rmse(region_summary: pd.DataFrame, output_root: Path) -> None:
    data = (
        region_summary.groupby(["relative_density_threshold", "method"], as_index=False)
        .agg(root_mean_square_error=("root_mean_square_error", "mean"))
    )
    data["method"] = pd.Categorical(data["method"], categories=METHOD_ORDER, ordered=True)
    set_theme()
    fig, ax = plt.subplots(figsize=(6.15, 3.35))
    fig.subplots_adjust(left=0.13, right=0.985, bottom=0.20, top=0.78)
    sns.barplot(
        data=data,
        x="method",
        y="root_mean_square_error",
        hue="relative_density_threshold",
        palette=[
            COLORS["threshold_light"],
            COLORS["threshold_mid"],
            COLORS["threshold_dark"],
        ],
        edgecolor=COLORS["ink"],
        linewidth=0.7,
        ax=ax,
    )
    ax.set_xlabel("")
    ax.set_ylabel("mean pointwise coverage RMSE")
    ax.yaxis.set_major_formatter(mticker.FormatStrFormatter("%.3f"))
    ax.legend(
        title="density threshold",
        frameon=False,
        ncol=3,
        loc="upper center",
        bbox_to_anchor=(0.5, 1.25),
        fontsize=8,
        title_fontsize=8,
    )
    ax.margins(y=0.10)
    sns.despine(ax=ax)
    save_figure(fig, output_root, "original_formula_region_rmse")


def plot_r_scaling(predictions: pd.DataFrame, output_root: Path) -> None:
    B = int(predictions["B"].max())
    data = (
        predictions[predictions["B"] == B]
        .groupby(["method", "R"], as_index=False)
        .agg(
            observed=("coverage", "mean"),
            predicted=("predicted_coverage", "mean"),
            baseline=("formula_baseline", "mean"),
        )
    )
    set_theme()
    fig, axes = plt.subplots(1, 3, figsize=(12.6, 4.2), sharey=True)
    for ax, method in zip(axes, METHOD_ORDER, strict=True):
        group = data[data["method"] == method].sort_values("R")
        ax.plot(
            1 / group["R"],
            group["observed"] - group["baseline"],
            marker="o",
            color=COLORS["observed"],
            linewidth=1.0,
            label="centered experiment",
        )
        ax.plot(
            1 / group["R"],
            group["predicted"] - group["baseline"],
            marker="s",
            color=COLORS["formula"],
            linewidth=1.2,
            label="factorized formula",
        )
        ax.axhline(0, color=COLORS["baseline"], linewidth=0.8)
        ax.set_title(METHOD_LABEL[method])
        ax.set_xlabel(r"$1/R$")
        sns.despine(ax=ax)
    axes[0].set_ylabel("grid-mean coverage minus baseline")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, frameon=False, ncol=2, loc="upper center")
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    save_figure(fig, output_root, "original_formula_r_scaling")


def write_change_log(
    output_root: Path,
    source_root: Path,
    summary: pd.DataFrame,
    region_summary: pd.DataFrame,
) -> None:
    overall = (
        summary.groupby("method", as_index=False)
        .agg(
            mean_mae=("mean_abs_error", "mean"),
            mean_rmse=("root_mean_square_error", "mean"),
            inside_exact95=("share_prediction_inside_exact95", "mean"),
            max_diag=("max_bootstrap_t_diagnostic_rate", "max"),
        )
    )
    core = (
        region_summary[region_summary["relative_density_threshold"] == 0.05]
        .groupby("method", as_index=False)
        .agg(
            mean_mae=("mean_abs_error", "mean"),
            mean_rmse=("root_mean_square_error", "mean"),
            inside_exact95=("share_prediction_inside_exact95", "mean"),
            max_l3=("max_abs_lambda3_over_sqrt_R", "max"),
            max_l4=("max_abs_lambda4_over_R", "max"),
        )
    )
    lines = [
        "# Conditional fixed-probability-weight formula rerun",
        "",
        f"- Source centered plate coverage data: `{source_root}`.",
        "- No coverage experiment was rerun; this rerun recomputes formula predictions, summaries, and figures from the saved centered thin-plate data.",
        "- Main change: use only the conditional fixed-probability-weight factorized inputs, not the complete outer-cumulant inputs.",
        "- Formula input: `lambda3_sq = gamma_h^2 * mean(rho3(w)^2)` and `lambda4 = kappa_h * mean(rho4(w))`.",
        "- Reason for `mean(rho3(w)^2)`: the plate pilot has randomized Voronoi weights per curve; the deterministic fixed-weight analogue for the second-order coverage term is the pilot average of the squared weight factor.",
        "- Centering convention is unchanged: the 600-curve validation pool is translated pointwise to the MC reference curve before coverage counting.",
        "- Figures remove the complete outer-cumulant curve to make the factorized formula the sole displayed prediction.",
        "",
        "## Overall summary",
        "",
        frame_to_markdown(overall),
        "",
        "## Core density region summary",
        "",
        "The core region uses reference density at least 5% of its peak.",
        "",
        frame_to_markdown(core),
        "",
        "## Output files",
        "",
        "- `original_formula_predictions.csv`: pointwise formula predictions and observed centered coverage.",
        "- `original_formula_summary.csv`: method/R/B summaries over the active grid.",
        "- `original_formula_density_region_summary.csv`: summaries after density-threshold filtering.",
        "- `original_formula_coefficients.csv`: direct finite-B coefficients used by the formulas.",
        "- `original_formula_coverage_profiles.pdf`: pointwise coverage profiles.",
        "- `original_formula_inputs.pdf`: kernel cumulants and factorized formula inputs.",
        "- `original_formula_region_rmse.pdf`: density-region RMSE comparison.",
        "- `original_formula_r_scaling.pdf`: grid-mean coverage error versus 1/R.",
    ]
    (output_root / "ORIGINAL_FORMULA_CHANGELOG.md").write_text(
        "\n".join(lines) + "\n",
        encoding="utf-8",
    )


def frame_to_markdown(frame: pd.DataFrame) -> str:
    columns = list(frame.columns)
    rows = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join("---" for _ in columns) + " |",
    ]
    for record in frame.itertuples(index=False):
        values = []
        for value in record:
            if isinstance(value, (float, np.floating)):
                values.append(f"{float(value):.5f}")
            else:
                values.append(str(value))
        rows.append("| " + " | ".join(values) + " |")
    return "\n".join(rows)


def main() -> None:
    package_root = Path(__file__).resolve().parent
    args = parse_args()
    source_root = args.source_root.resolve() if args.source_root else latest_plate_formal_root(package_root)
    if args.output_root:
        output_root = args.output_root.resolve()
    else:
        stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        output_root = package_root / "results" / f"plate_original_formula_validation_{stamp}"
    output_root.mkdir(parents=True, exist_ok=False)

    predictions, summary, region_summary, coefficients, checks = build_original_predictions(
        source_root
    )
    predictions.to_csv(output_root / "original_formula_predictions.csv", index=False)
    summary.to_csv(output_root / "original_formula_summary.csv", index=False)
    region_summary.to_csv(
        output_root / "original_formula_density_region_summary.csv", index=False
    )
    coefficients.to_csv(output_root / "original_formula_coefficients.csv", index=False)
    (output_root / "original_formula_coefficient_self_checks.json").write_text(
        json.dumps(checks, indent=2),
        encoding="utf-8",
    )

    plot_formula_inputs(predictions, output_root)
    plot_coverage_profiles(predictions, output_root)
    plot_region_rmse(region_summary, output_root)
    plot_r_scaling(predictions, output_root)
    write_change_log(output_root, source_root, summary, region_summary)

    manifest = {
        "schema_version": "plate_original_formula_validation_v1",
        "source_root": str(source_root),
        "output_root": str(output_root),
        "calculation_policy": "direct_non_fitted",
        "formula_input": "original_fixed_weight_factorization",
        "lambda3_sq": "gamma_h^2 * mean(rho3(w)^2)",
        "lambda4": "kappa_h * mean(rho4(w))",
        "complete_outer_cumulants_used": False,
        "centered_coverage_data_reused": True,
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "files": sorted(path.name for path in output_root.iterdir()),
    }
    (output_root / "analysis_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    print(output_root)


if __name__ == "__main__":
    main()
