from __future__ import annotations

import argparse
import json
from pathlib import Path

import h5py
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
SOURCE_ORDER = ["complete outer cumulants", "fixed-weight factorization"]
COLORS = {
    "observed": "#1F2430",
    "complete outer cumulants": "#5477C4",
    "fixed-weight factorization": "#CC6F47",
    "baseline": "#7A828F",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("results_root", type=Path)
    return parser.parse_args()


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


def read_bootstrap_moments(mat_path: Path, active_count: int) -> tuple[np.ndarray, np.ndarray]:
    with h5py.File(mat_path, "r") as handle:
        lambda3 = np.asarray(handle["lambda3Boot"])
        lambda4 = np.asarray(handle["lambda4Boot"])
    if lambda3.shape[1] == active_count:
        return lambda3, lambda4
    if lambda3.shape[0] == active_count:
        return lambda3.T, lambda4.T
    raise ValueError(
        f"Unexpected bootstrap moment shape {lambda3.shape}; active count is {active_count}."
    )


def build_predictions(
    results_root: Path,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict]:
    coverage = pd.read_csv(results_root / "coverage_results.csv")
    moments = pd.read_csv(results_root / "outer_moments_active_grid.csv")
    lambda3_boot, lambda4_boot = read_bootstrap_moments(
        results_root / "plate_outer_cumulant_results.mat", len(moments)
    )

    coefficient_rows: list[dict] = []
    coefficient_checks: dict[str, dict] = {}
    coefficient_maps: dict[int, dict] = {}
    for B in sorted(coverage["B"].unique()):
        rows, checks = direct_coefficients(0.05, int(B))
        for row in rows:
            coefficient_rows.append({"B": int(B), **row})
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
        source_values = {
            "complete outer cumulants": (row.outer_lambda4, row.outer_lambda3_sq),
            "fixed-weight factorization": (
                row.factorized_lambda4,
                row.factorized_lambda3_sq,
            ),
        }
        for source, (lambda4, lambda3_sq) in source_values.items():
            coefficient_A = A0 + A4 * lambda4 + A33 * lambda3_sq
            predicted = row.formula_baseline + coefficient_A / row.R
            record = row._asdict()
            record.update(
                {
                    "moment_source": source,
                    "A0": A0,
                    "A4": A4,
                    "A33": A33,
                    "lambda4_input": lambda4,
                    "lambda3_sq_input": lambda3_sq,
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

    full_mask = predictions["moment_source"] == "complete outer cumulants"
    for (method, R, B), index in predictions.loc[full_mask].groupby(
        ["method", "R", "B"], sort=False
    ).groups.items():
        coefficients = coefficient_maps[int(B)]
        A0 = coefficients[(method, "A0")]
        A4 = coefficients[(method, "A4")]
        A33 = coefficients[(method, "A33")]
        baseline = predictions.loc[index, "formula_baseline"].iloc[0]
        boot_prediction = baseline + (
            A0 + A4 * lambda4_boot + A33 * lambda3_boot**2
        ) / R
        predictions.loc[index, "formula_boot95_lower"] = np.quantile(
            boot_prediction, 0.025, axis=0
        )
        predictions.loc[index, "formula_boot95_upper"] = np.quantile(
            boot_prediction, 0.975, axis=0
        )

    summaries = []
    for keys, group in predictions.groupby(
        ["method", "R", "B", "moment_source"], sort=False
    ):
        method, R, B, source = keys
        residual = group["residual"].to_numpy()
        summaries.append(
            {
                "method": method,
                "R": R,
                "B": B,
                "moment_source": source,
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
                "max_bootstrap_t_diagnostic_rate": group[
                    "bootstrap_t_diagnostic_rate"
                ].max(),
                "max_interval_inf_rate": group["interval_inf_rate"].max(),
            }
        )

    region_rows = []
    for threshold in (0.01, 0.05, 0.10):
        region = predictions[predictions["relative_density"] >= threshold]
        for keys, group in region.groupby(
            ["method", "R", "B", "moment_source"], sort=False
        ):
            method, R, B, source = keys
            residual = group["residual"].to_numpy()
            region_rows.append(
                {
                    "relative_density_threshold": threshold,
                    "method": method,
                    "R": R,
                    "B": B,
                    "moment_source": source,
                    "grid_point_count": len(group),
                    "mean_abs_error": np.mean(np.abs(residual)),
                    "root_mean_square_error": np.sqrt(np.mean(residual**2)),
                    "max_abs_error": np.max(np.abs(residual)),
                    "share_prediction_inside_exact95": group[
                        "prediction_inside_exact95"
                    ].mean(),
                    "max_abs_lambda3_over_sqrt_R": np.max(
                        np.abs(group["outer_lambda3"]) / np.sqrt(R)
                    ),
                    "max_abs_lambda4_over_R": np.max(
                        np.abs(group["outer_lambda4"]) / R
                    ),
                }
            )

    return predictions, pd.DataFrame(summaries), pd.DataFrame(region_rows), {
        "coefficient_rows": coefficient_rows,
        "coefficient_checks": coefficient_checks,
    }


def set_theme() -> None:
    sns.set_theme(
        style="whitegrid",
        rc={
            "figure.facecolor": "#FCFCFD",
            "axes.facecolor": "#FFFFFF",
            "axes.edgecolor": "#D7DBE7",
            "axes.labelcolor": "#1F2430",
            "text.color": "#1F2430",
            "grid.color": "#E6E8F0",
            "grid.linewidth": 0.8,
            "font.family": "sans-serif",
            "font.sans-serif": ["Aptos", "Segoe UI", "DejaVu Sans"],
        },
    )


def save_figure(fig: plt.Figure, root: Path, name: str) -> None:
    for extension in ("pdf", "png"):
        fig.savefig(
            root / f"{name}.{extension}",
            dpi=320,
            bbox_inches="tight",
        )
    plt.close(fig)


def plot_cumulants(predictions: pd.DataFrame, root: Path) -> None:
    moments = predictions.drop_duplicates("response_m").sort_values("response_m")
    set_theme()
    fig, axes = plt.subplots(2, 1, figsize=(9.2, 7.0), sharex=True)
    for ax, outer, factorized, label in (
        (
            axes[0],
            "outer_lambda3",
            "factorized_lambda3",
            r"standardized skewness $\lambda_3$",
        ),
        (
            axes[1],
            "outer_lambda4",
            "factorized_lambda4",
            r"excess kurtosis $\lambda_4$",
        ),
    ):
        ax.plot(
            moments["response_m"],
            moments[outer],
            color=COLORS["complete outer cumulants"],
            linewidth=1.3,
            label="complete outer cumulants",
        )
        ax.plot(
            moments["response_m"],
            moments[factorized],
            color=COLORS["fixed-weight factorization"],
            linewidth=1.2,
            linestyle="--",
            label="fixed-weight factorization",
        )
        ax.axhline(0, color=COLORS["baseline"], linewidth=0.8)
        ax.set_ylabel(label)
        ax.ticklabel_format(axis="x", style="sci", scilimits=(-3, 3))
        sns.despine(ax=ax)
    axes[0].legend(frameon=False, ncol=2, loc="upper center")
    axes[1].set_xlabel("plate-center deflection response (m)")
    fig.suptitle("Random-weight cross structure changes the outer cumulants", y=0.995)
    fig.text(
        0.5,
        0.945,
        "Independent 200-curve pilot; probability-weighted RQMC/Voronoi DPIM",
        ha="center",
        color="#6F768A",
        fontsize=9,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    save_figure(fig, root, "outer_vs_factorized_cumulants")


def plot_profiles(predictions: pd.DataFrame, root: Path) -> None:
    B = int(predictions["B"].max())
    r_values = [int(predictions["R"].min()), int(predictions["R"].max())]
    set_theme()
    fig, axes = plt.subplots(3, 2, figsize=(12.0, 10.2), sharex=True, sharey=True)
    for i, method in enumerate(METHOD_ORDER):
        for j, R in enumerate(r_values):
            ax = axes[i, j]
            group = predictions[
                (predictions["method"] == method)
                & (predictions["R"] == R)
                & (predictions["B"] == B)
            ].sort_values("response_m")
            observed = group[group["moment_source"] == SOURCE_ORDER[0]]
            ax.plot(
                observed["response_m"],
                observed["coverage"],
                color=COLORS["observed"],
                linewidth=1.0,
                label="experiment",
            )
            ax.fill_between(
                observed["response_m"],
                observed["coverage_exact95_lower"],
                observed["coverage_exact95_upper"],
                color="#E2E5EA",
                alpha=0.65,
                linewidth=0,
                label="experiment exact 95% interval",
            )
            for source, linestyle in zip(SOURCE_ORDER, ("-", "--"), strict=True):
                part = group[group["moment_source"] == source]
                ax.plot(
                    part["response_m"],
                    part["predicted_coverage"],
                    color=COLORS[source],
                    linestyle=linestyle,
                    linewidth=1.3,
                    label=source,
                )
            ax.axhline(
                observed["formula_baseline"].iloc[0],
                color=COLORS["baseline"],
                linestyle=":",
                linewidth=0.9,
            )
            ax.set_title(f"{method}, R={R}, B={B}", fontsize=10)
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
        ncol=4,
        frameon=False,
        bbox_to_anchor=(0.5, 0.97),
    )
    fig.suptitle("Direct coverage prediction versus centered plate experiment", y=0.998)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    save_figure(fig, root, "plate_coverage_profiles")


def plot_rmse(summary: pd.DataFrame, root: Path) -> None:
    overall = (
        summary.groupby(["method", "moment_source"], as_index=False)
        .agg(root_mean_square_error=("root_mean_square_error", "mean"))
    )
    overall["method"] = pd.Categorical(
        overall["method"], categories=METHOD_ORDER, ordered=True
    )
    overall["moment_source"] = pd.Categorical(
        overall["moment_source"], categories=SOURCE_ORDER, ordered=True
    )
    overall = overall.sort_values(["method", "moment_source"])
    set_theme()
    fig, ax = plt.subplots(figsize=(8.8, 4.8))
    palette = {source: COLORS[source] for source in SOURCE_ORDER}
    sns.barplot(
        data=overall,
        x="method",
        y="root_mean_square_error",
        hue="moment_source",
        palette=palette,
        edgecolor="#464C55",
        linewidth=0.7,
        ax=ax,
    )
    ax.set_xlabel("")
    ax.set_ylabel("mean pointwise coverage RMSE")
    ax.legend(frameon=False, ncol=2, loc="upper center")
    ax.yaxis.set_major_formatter(mticker.FormatStrFormatter("%.3f"))
    ax.set_title("Complete outer cumulants do not automatically dominate at finite pilot size")
    fig.text(
        0.5,
        0.91,
        "Average of configuration-level RMSE across R and B; lower is better",
        ha="center",
        color="#6F768A",
        fontsize=9,
    )
    sns.despine(ax=ax)
    fig.tight_layout(rect=(0, 0, 1, 0.88))
    save_figure(fig, root, "prediction_rmse_by_moment_source")


def plot_r_scaling(predictions: pd.DataFrame, root: Path) -> None:
    B = int(predictions["B"].max())
    data = (
        predictions[predictions["B"] == B]
        .groupby(["method", "R", "moment_source"], as_index=False)
        .agg(
            observed=("coverage", "mean"),
            predicted=("predicted_coverage", "mean"),
            baseline=("formula_baseline", "mean"),
        )
    )
    set_theme()
    fig, axes = plt.subplots(1, 3, figsize=(12.6, 4.2), sharey=True)
    for ax, method in zip(axes, METHOD_ORDER, strict=True):
        group = data[data["method"] == method]
        observed = group[group["moment_source"] == SOURCE_ORDER[0]].sort_values("R")
        ax.plot(
            1 / observed["R"],
            observed["observed"] - observed["baseline"],
            marker="o",
            color=COLORS["observed"],
            linewidth=1.0,
            label="experiment",
        )
        for source, linestyle in zip(SOURCE_ORDER, ("-", "--"), strict=True):
            part = group[group["moment_source"] == source].sort_values("R")
            ax.plot(
                1 / part["R"],
                part["predicted"] - part["baseline"],
                color=COLORS[source],
                linestyle=linestyle,
                linewidth=1.3,
                label=source,
            )
        ax.axhline(0, color=COLORS["baseline"], linewidth=0.8)
        ax.set_title(method)
        ax.set_xlabel(r"$1/R$")
        sns.despine(ax=ax)
    axes[0].set_ylabel("grid-mean coverage minus baseline")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, frameon=False, ncol=3, loc="upper center")
    fig.suptitle(f"Coverage error scaling for B={B}", y=1.01)
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    save_figure(fig, root, "coverage_error_scaling")


def write_html_report(results_root: Path, summary: pd.DataFrame) -> None:
    overall = (
        summary.groupby(["method", "moment_source"], as_index=False)
        .agg(
            mean_rmse=("root_mean_square_error", "mean"),
            mean_mae=("mean_abs_error", "mean"),
            inside_exact95=("share_prediction_inside_exact95", "mean"),
            max_diagnostic=("max_bootstrap_t_diagnostic_rate", "max"),
        )
    )
    table_html = overall.to_html(index=False, float_format=lambda value: f"{value:.5f}")
    html = f"""<!doctype html>
<html lang="zh-CN">
<head>
<meta charset="utf-8">
<title>随机权重外层累积量覆盖展开验证</title>
<style>
body {{ font-family: "Segoe UI", sans-serif; margin: 42px auto; max-width: 1040px; color: #1f2430; line-height: 1.65; }}
h1, h2 {{ font-family: Georgia, "Times New Roman", serif; }}
img {{ width: 100%; margin: 12px 0 28px; }}
table {{ border-collapse: collapse; width: 100%; font-size: 14px; }}
th, td {{ border-bottom: 1px solid #d7dbe7; padding: 8px 10px; text-align: right; }}
th:first-child, td:first-child, th:nth-child(2), td:nth-child(2) {{ text-align: left; }}
.note {{ background: #f4f5f7; border-left: 4px solid #5477c4; padding: 12px 16px; }}
</style>
</head>
<body>
<h1>随机权重外层累积量覆盖展开验证</h1>
<p class="note">本报告以完整外层曲线为统计单位。固定权重因子化结果仅作为特殊情形对照，不作为随机 Voronoi/RQMC 的一般定理。</p>
<h2>累积量差异</h2>
<img src="outer_vs_factorized_cumulants.png" alt="外层累积量比较">
<h2>逐点覆盖率</h2>
<img src="plate_coverage_profiles.png" alt="逐点覆盖率比较">
<h2>误差汇总</h2>
{table_html}
<img src="prediction_rmse_by_moment_source.png" alt="预测误差汇总">
<h2>R 阶机制</h2>
<img src="coverage_error_scaling.png" alt="覆盖误差随1/R变化">
<p>严格理论边界、推导和引用见同目录 LaTeX 报告。数值结论必须与 pilot 累积量不确定性、有限曲线池误差和 Edgeworth 余项同时解读。</p>
</body>
</html>"""
    (results_root / "report.html").write_text(html, encoding="utf-8")


def main() -> None:
    args = parse_args()
    results_root = args.results_root.resolve()
    predictions, summary, region_summary, coefficient_data = build_predictions(
        results_root
    )
    predictions.to_csv(results_root / "formula_predictions.csv", index=False)
    summary.to_csv(results_root / "formula_validation_summary.csv", index=False)
    region_summary.to_csv(
        results_root / "formula_validation_by_density_region.csv", index=False
    )
    pd.DataFrame(coefficient_data["coefficient_rows"]).to_csv(
        results_root / "direct_coefficients.csv", index=False
    )
    (results_root / "coefficient_self_checks.json").write_text(
        json.dumps(coefficient_data["coefficient_checks"], indent=2),
        encoding="utf-8",
    )
    plot_cumulants(predictions, results_root)
    plot_profiles(predictions, results_root)
    plot_rmse(summary, results_root)
    plot_r_scaling(predictions, results_root)
    write_html_report(results_root, summary)
    manifest = {
        "calculation_policy": "direct_non_fitted",
        "primary_formula_input": "independent-pilot complete outer cumulants",
        "control_formula_input": "fixed-weight factorization",
        "coverage_target": "MC reference after pointwise translation of the validation pool",
        "files": sorted(path.name for path in results_root.iterdir()),
    }
    (results_root / "analysis_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
