from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.patches import Circle, Polygon, Rectangle


TOKENS = {
    "surface": "#FFFFFF",
    "ink": "#1A1A1A",
    "muted": "#666666",
    "grid": "#D9D9D9",
    "gold": "#1F4E79",
    "green": "#2E7D32",
    "blue": "#1F4E79",
    "red": "#B22222",
    "band": "#ECECEC",
    "light_blue": "#6BAED6",
}

METHOD_STYLES = {
    "Student-t": {"color": TOKENS["blue"], "marker": "o", "linestyle": "-"},
    "t distribution": {"color": TOKENS["blue"], "marker": "o", "linestyle": "-"},
    "Percentile bootstrap": {
        "color": TOKENS["green"],
        "marker": "s",
        "linestyle": "--",
    },
    "percentile bootstrap": {
        "color": TOKENS["green"],
        "marker": "s",
        "linestyle": "--",
    },
    "Bootstrap-t": {"color": TOKENS["red"], "marker": "^", "linestyle": "-."},
    "bootstrap-t": {"color": TOKENS["red"], "marker": "^", "linestyle": "-."},
    "Percentile max-deviation band": {
        "color": TOKENS["green"],
        "marker": "s",
        "linestyle": "--",
    },
    "Bootstrap-t max-stat band": {
        "color": TOKENS["blue"],
        "marker": "^",
        "linestyle": "-.",
    },
    "percentile max-deviation band": {
        "color": TOKENS["green"],
        "marker": "s",
        "linestyle": "--",
    },
    "bootstrap-t max-stat band": {
        "color": TOKENS["blue"],
        "marker": "^",
        "linestyle": "-.",
    },
}

DISPLAY_NAMES = {
    "Student-t": "Student t",
    "t distribution": "Student t",
    "Percentile bootstrap": "percentile bootstrap",
    "percentile bootstrap": "percentile bootstrap",
    "Bootstrap-t": "bootstrap-t",
    "bootstrap-t": "bootstrap-t",
    "Percentile max-deviation band": "percentile bootstrap band",
    "Bootstrap-t max-stat band": "bootstrap-t band",
    "percentile max-deviation band": "percentile bootstrap band",
    "bootstrap-t max-stat band": "bootstrap-t band",
}


def use_theme() -> None:
    sns.set_theme(
        style="whitegrid",
        rc={
            "figure.facecolor": TOKENS["surface"],
            "savefig.facecolor": TOKENS["surface"],
            "axes.facecolor": TOKENS["surface"],
            "axes.edgecolor": TOKENS["ink"],
            "axes.labelcolor": TOKENS["ink"],
            "axes.titlecolor": TOKENS["ink"],
            "grid.color": TOKENS["grid"],
            "grid.linewidth": 0.45,
            "font.family": "serif",
            "font.serif": ["Times New Roman", "STIXGeneral", "DejaVu Serif"],
            "mathtext.fontset": "stix",
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.linewidth": 0.75,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.major.width": 0.75,
            "ytick.major.width": 0.75,
            "pdf.fonttype": 42,
        },
    )


def save_pdf(fig: plt.Figure, output_dir: Path, stem: str) -> None:
    fig.savefig(output_dir / f"{stem}.pdf", bbox_inches="tight", pad_inches=0.025)
    plt.close(fig)


def legend_above(ax: plt.Axes, ncol: int = 2, y: float = 1.18) -> None:
    ax.legend(
        frameon=False,
        fontsize=8,
        ncol=ncol,
        loc="upper center",
        bbox_to_anchor=(0.5, y),
        handlelength=2.1,
        columnspacing=0.85,
    )


def plot_method_lines(
    ax: plt.Axes,
    frame: pd.DataFrame,
    x_col: str,
    y_col: str,
    methods: list[str],
    baseline_col: str | None = None,
) -> None:
    for method in methods:
        part = frame.loc[frame["method"] == method].sort_values(x_col)
        style = METHOD_STYLES[method]
        ax.plot(
            part[x_col],
            part[y_col],
            label=DISPLAY_NAMES[method],
            color=style["color"],
            marker=style["marker"],
            linestyle=style["linestyle"],
            linewidth=1.15,
            markersize=4.1,
            markerfacecolor="white",
            markeredgewidth=1.0,
            zorder=3,
        )
    if baseline_col is not None:
        baseline = float(frame[baseline_col].iloc[0])
        ax.axhline(
            baseline,
            color=TOKENS["ink"],
            linestyle=":",
            linewidth=1.0,
            label=rf"$C_{{0,B}}={baseline:.2f}$",
        )


def build_finite_b_figure(formal_root: Path, output_dir: Path) -> None:
    frame = pd.read_csv(
        formal_root
        / "_manuscript_assets"
        / "chart_data"
        / "E1_finite_B_clean.csv"
    )
    fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.45))
    fig.subplots_adjust(left=0.09, right=0.985, bottom=0.18, top=0.96, wspace=0.27)

    ax = axes[0]
    ax.errorbar(
        frame["B"],
        frame["simulated_coverage"],
        yerr=frame["simulation_ci_half_width_95"],
        fmt="o",
        color=TOKENS["blue"],
        markerfacecolor="white",
        capsize=3,
        linewidth=1.0,
        label="simulation",
    )
    ax.plot(
        frame["B"],
        frame["C0B"],
        color=TOKENS["ink"],
        marker="s",
        markersize=4,
        linewidth=1.1,
        label=r"$C_{0,B}$",
    )
    ax.axhline(0.95, color=TOKENS["muted"], linestyle=":", linewidth=0.9)
    ax.set_xscale("log")
    ax.set_xlabel(r"Bootstrap replicate count $B$")
    ax.set_ylabel("Coverage probability")
    ax.legend(frameon=False, fontsize=8)

    ax = axes[1]
    ax.axhline(0, color=TOKENS["ink"], linewidth=0.9)
    ax.plot(
        frame["B"],
        frame["grid_error_vs_nominal"],
        color=TOKENS["red"],
        marker="o",
        markerfacecolor="white",
        linewidth=1.1,
        label=r"$C_{0,B}-(1-\alpha)$",
    )
    ax.plot(
        frame["B"],
        frame["simulation_error_vs_C0B_signed"],
        color=TOKENS["blue"],
        marker="s",
        markerfacecolor="white",
        linestyle="--",
        linewidth=1.1,
        label=r"$\widehat C-C_{0,B}$",
    )
    ax.set_xscale("log")
    ax.set_xlabel(r"Bootstrap replicate count $B$")
    ax.set_ylabel("Coverage difference")
    ax.legend(frameon=False, fontsize=8)
    save_pdf(fig, output_dir, "finite_B_order_statistic")


def load_confirmation(path: Path) -> pd.DataFrame:
    methods = ["Student-t", "Percentile bootstrap", "Bootstrap-t"]
    frame = pd.read_csv(path)
    frame = frame.loc[frame["method"].isin(methods)].copy()
    frame = frame.loc[frame["selection_reason"] == "tuning_selected_stable"]
    return (
        frame.groupby(["method", "R"], as_index=False)
        .agg(
            mean_coverage=("coverage", "mean"),
            mean_formula_baseline=("formula_baseline", "mean"),
        )
        .sort_values(["method", "R"])
    )


def build_pointwise_figure(
    source_path: Path, output_dir: Path, stem: str
) -> None:
    methods = ["Student-t", "Percentile bootstrap", "Bootstrap-t"]
    frame = load_confirmation(source_path)
    fig, ax = plt.subplots(figsize=(6.15, 3.35))
    fig.subplots_adjust(left=0.13, right=0.985, bottom=0.19, top=0.82)
    plot_method_lines(
        ax,
        frame,
        "R",
        "mean_coverage",
        methods,
        baseline_col="mean_formula_baseline",
    )
    ax.set_xticks(sorted(frame["R"].unique()))
    ax.set_xlabel(r"Outer replicate count $R$")
    ax.set_ylabel("Mean pointwise coverage")
    ax.margins(y=0.08)
    legend_above(ax, ncol=2)
    save_pdf(fig, output_dir, stem)


def build_weighted_density_figure(formal_root: Path, output_dir: Path) -> None:
    frame = pd.read_csv(
        formal_root
        / "_manuscript_assets"
        / "chart_data"
        / "E5_coverage_by_R_method.csv"
    )
    methods = ["Student-t", "Percentile bootstrap", "Bootstrap-t"]
    frame = frame.loc[frame["method"].isin(methods)].copy()
    fig, ax = plt.subplots(figsize=(6.15, 3.35))
    fig.subplots_adjust(left=0.13, right=0.985, bottom=0.19, top=0.82)
    plot_method_lines(
        ax,
        frame,
        "R",
        "mean_coverage",
        methods,
        baseline_col="mean_formula_baseline",
    )
    ax.set_xticks(sorted(frame["R"].unique()))
    ax.set_xlabel(r"Outer replicate count $R$")
    ax.set_ylabel("Mean density-ordinate coverage")
    ax.margins(y=0.08)
    legend_above(ax, ncol=2)
    save_pdf(fig, output_dir, "linear_weighted_density_coverage")


def build_band_figure(formal_root: Path, output_dir: Path) -> None:
    frame = pd.read_csv(
        formal_root
        / "_manuscript_assets"
        / "chart_data"
        / "E7_band_coverage_by_R_method.csv"
    )
    methods = ["Percentile max-deviation band", "Bootstrap-t max-stat band"]
    frame = frame.loc[frame["method"].isin(methods)].copy()
    fig, ax = plt.subplots(figsize=(6.15, 3.35))
    fig.subplots_adjust(left=0.13, right=0.985, bottom=0.19, top=0.82)
    plot_method_lines(
        ax,
        frame,
        "R",
        "mean_coverage",
        methods,
        baseline_col="mean_formula_baseline",
    )
    ax.set_xticks(sorted(frame["R"].unique()))
    ax.set_xlabel(r"Outer replicate count $R$")
    ax.set_ylabel("Finite-grid simultaneous coverage")
    ax.margins(y=0.08)
    legend_above(ax, ncol=1, y=1.16)
    save_pdf(fig, output_dir, "linear_finite_grid_band_coverage")


def build_direct_coverage_figure(
    direct_assets: Path, output_dir: Path, model: str
) -> None:
    frame = pd.read_csv(direct_assets / "observed_vs_direct_validation.csv")
    frame = frame.loc[frame["model"] == model].copy()
    methods = ["t distribution", "percentile bootstrap", "bootstrap-t"]
    fig, axes = plt.subplots(1, 3, figsize=(8.35, 2.95), sharey=True)
    fig.subplots_adjust(left=0.075, right=0.99, bottom=0.20, top=0.80, wspace=0.16)
    legend_handles = None
    legend_labels = None
    for index, (ax, method) in enumerate(zip(axes, methods)):
        part = frame.loc[frame["method"] == method].sort_values("R")
        style = METHOD_STYLES[method]
        lower = part["observed_coverage"] - part["exact_95_lower"]
        upper = part["exact_95_upper"] - part["observed_coverage"]
        ax.errorbar(
            part["R"],
            part["observed_coverage"],
            yerr=np.vstack([lower, upper]),
            fmt=style["marker"],
            color=TOKENS["ink"],
            markerfacecolor="white",
            capsize=2.4,
            linewidth=0.9,
            markersize=4,
            label="experiment",
        )
        ax.plot(
            part["R"],
            part["predicted_coverage"],
            color=style["color"],
            linestyle=style["linestyle"],
            marker=style["marker"],
            markerfacecolor="white",
            linewidth=1.2,
            markersize=4,
            label="direct calculation",
        )
        ax.axhline(
            float(part["baseline"].iloc[0]),
            color=TOKENS["muted"],
            linestyle=":",
            linewidth=0.9,
        )
        ax.text(
            0.03,
            1.03,
            f"({chr(97 + index)}) {DISPLAY_NAMES[method]}",
            transform=ax.transAxes,
            ha="left",
            va="bottom",
            fontsize=9,
        )
        ax.set_xlabel(r"$R$")
        ax.set_xticks([16, 48, 96, 192])
        if index == 0:
            legend_handles, legend_labels = ax.get_legend_handles_labels()
    axes[0].set_ylabel("Coverage probability")
    if legend_handles is not None and legend_labels is not None:
        fig.legend(
            legend_handles,
            legend_labels,
            frameon=False,
            fontsize=7.6,
            ncol=2,
            loc="upper center",
            bbox_to_anchor=(0.5, 0.995),
        )
    save_pdf(fig, output_dir, f"{model}_direct_formula_coverage")


def build_direct_residual_figure(
    direct_assets: Path, output_dir: Path, model: str
) -> None:
    frame = pd.read_csv(direct_assets / "observed_vs_direct_validation.csv")
    frame = frame.loc[frame["model"] == model].copy()
    methods = ["t distribution", "percentile bootstrap", "bootstrap-t"]
    fig, ax = plt.subplots(figsize=(6.15, 3.35))
    fig.subplots_adjust(left=0.13, right=0.985, bottom=0.20, top=0.82)
    for method in methods:
        part = frame.loc[frame["method"] == method].sort_values("R")
        style = METHOD_STYLES[method]
        ax.plot(
            part["R"],
            part["z_null"],
            label=DISPLAY_NAMES[method],
            color=style["color"],
            marker=style["marker"],
            markerfacecolor="white",
            linestyle=style["linestyle"],
            linewidth=1.15,
            markersize=4.3,
        )
    ax.axhline(0, color=TOKENS["ink"], linewidth=0.8)
    ax.axhline(1.96, color=TOKENS["muted"], linestyle=":", linewidth=0.9)
    ax.axhline(-1.96, color=TOKENS["muted"], linestyle=":", linewidth=0.9)
    ax.set_xlabel(r"Outer replicate count $R$")
    ax.set_ylabel("Standardized direct-calculation residual")
    legend_above(ax, ncol=2)
    save_pdf(fig, output_dir, f"{model}_direct_formula_residuals")


def build_uncentered_figure(
    direct_root: Path, output_dir: Path, model: str
) -> None:
    frame = pd.read_csv(direct_root / "coverage_results_uncentered.csv")
    methods = ["t distribution", "percentile bootstrap", "bootstrap-t"]
    frame = frame.loc[
        (frame["model"] == model) & frame["method"].isin(methods)
    ].copy()
    fig, ax = plt.subplots(figsize=(6.4, 3.65))
    fig.subplots_adjust(left=0.13, right=0.98, bottom=0.19, top=0.97)
    plot_method_lines(
        ax,
        frame,
        "R",
        "coverage",
        methods,
        baseline_col="formula_baseline",
    )
    ax.set_xlabel(r"Outer replicate count $R$")
    ax.set_ylabel("End-to-end coverage probability")
    ax.legend(frameon=False, fontsize=8, ncol=2)
    save_pdf(fig, output_dir, f"{model}_uncentered_coverage")


def build_beam_schematic(output_dir: Path) -> None:
    fig, ax = plt.subplots(figsize=(6.2, 1.9))
    fig.subplots_adjust(left=0.03, right=0.97, bottom=0.12, top=0.92)
    ax.plot([0.7, 6.3], [0, 0], color=TOKENS["ink"], linewidth=4)
    for x in np.linspace(1.0, 6.0, 11):
        ax.annotate(
            "",
            xy=(x, 0.10),
            xytext=(x, 0.92),
            arrowprops={"arrowstyle": "-|>", "color": TOKENS["blue"], "lw": 0.9},
        )
    left_support = Polygon(
        [[0.7, -0.05], [0.43, -0.48], [0.97, -0.48]],
        closed=True,
        facecolor="none",
        edgecolor=TOKENS["ink"],
        linewidth=1.0,
    )
    ax.add_patch(left_support)
    right_support = Polygon(
        [[6.3, -0.05], [6.03, -0.42], [6.57, -0.42]],
        closed=True,
        facecolor="none",
        edgecolor=TOKENS["ink"],
        linewidth=1.0,
    )
    ax.add_patch(right_support)
    ax.add_patch(Circle((6.15, -0.50), 0.08, fill=False, color=TOKENS["ink"], lw=0.9))
    ax.add_patch(Circle((6.45, -0.50), 0.08, fill=False, color=TOKENS["ink"], lw=0.9))
    ax.plot([0.25, 1.15], [-0.58, -0.58], color=TOKENS["ink"], linewidth=0.8)
    ax.plot([5.85, 6.75], [-0.62, -0.62], color=TOKENS["ink"], linewidth=0.8)
    ax.plot(3.5, 0, marker="o", markersize=6, color=TOKENS["red"])
    ax.text(3.5, -0.28, r"$x_0=L/2$", ha="center", va="top", fontsize=10)
    ax.text(3.5, 1.02, r"$q(x,\boldsymbol{\Theta})$", ha="center", fontsize=10)
    ax.annotate(
        "",
        xy=(0.7, -0.84),
        xytext=(6.3, -0.84),
        arrowprops={"arrowstyle": "<->", "color": TOKENS["ink"], "lw": 0.8},
    )
    ax.text(3.5, -0.80, r"$L$", ha="center", va="bottom", fontsize=10)
    ax.set_xlim(0.05, 6.95)
    ax.set_ylim(-1.05, 1.2)
    ax.axis("off")
    save_pdf(fig, output_dir, "linear_euler_bernoulli_beam_model")


def build_plate_schematic(output_dir: Path) -> None:
    fig, ax = plt.subplots(figsize=(4.15, 3.65))
    fig.subplots_adjust(left=0.03, right=0.97, bottom=0.04, top=0.97)
    rect = Rectangle((0, 0), 1, 1, fill=False, edgecolor=TOKENS["ink"], linewidth=1.7)
    ax.add_patch(rect)
    for value in np.linspace(0.1, 0.9, 9):
        ax.plot([value, value], [0, 1], color=TOKENS["grid"], linewidth=0.45)
        ax.plot([0, 1], [value, value], color=TOKENS["grid"], linewidth=0.45)
    for x in np.linspace(0.08, 0.92, 7):
        ax.add_patch(
            Polygon(
                [[x, -0.01], [x - 0.025, -0.07], [x + 0.025, -0.07]],
                closed=True,
                facecolor="none",
                edgecolor=TOKENS["muted"],
                linewidth=0.55,
            )
        )
        ax.add_patch(
            Polygon(
                [[x, 1.01], [x - 0.025, 1.07], [x + 0.025, 1.07]],
                closed=True,
                facecolor="none",
                edgecolor=TOKENS["muted"],
                linewidth=0.55,
            )
        )
    for y in np.linspace(0.08, 0.92, 7):
        ax.add_patch(
            Polygon(
                [[-0.01, y], [-0.07, y - 0.025], [-0.07, y + 0.025]],
                closed=True,
                facecolor="none",
                edgecolor=TOKENS["muted"],
                linewidth=0.55,
            )
        )
        ax.add_patch(
            Polygon(
                [[1.01, y], [1.07, y - 0.025], [1.07, y + 0.025]],
                closed=True,
                facecolor="none",
                edgecolor=TOKENS["muted"],
                linewidth=0.55,
            )
        )
    ax.plot(0.5, 0.5, marker="o", markersize=6, color=TOKENS["red"])
    ax.annotate(
        r"$q_0$",
        xy=(0.5, 0.53),
        xytext=(0.68, 0.78),
        arrowprops={"arrowstyle": "-|>", "color": TOKENS["blue"], "lw": 1.0},
        fontsize=11,
    )
    ax.text(0.5, 0.43, r"$w_c$", ha="center", va="top", fontsize=10)
    ax.annotate(
        "",
        xy=(0, -0.14),
        xytext=(1, -0.14),
        arrowprops={"arrowstyle": "<->", "color": TOKENS["ink"], "lw": 0.8},
    )
    ax.text(0.5, -0.12, r"$a=1\,\mathrm{m}$", ha="center", va="bottom", fontsize=9.5)
    ax.annotate(
        "",
        xy=(-0.14, 0),
        xytext=(-0.14, 1),
        arrowprops={"arrowstyle": "<->", "color": TOKENS["ink"], "lw": 0.8},
    )
    ax.text(-0.12, 0.5, r"$b=1\,\mathrm{m}$", rotation=90, ha="left", va="center", fontsize=9.5)
    ax.set_xlim(-0.22, 1.12)
    ax.set_ylim(-0.22, 1.12)
    ax.set_aspect("equal")
    ax.axis("off")
    save_pdf(fig, output_dir, "plate_sfem_model")


def build_plate_density_reference_figure(plate_density_root: Path, output_dir: Path) -> None:
    frame = pd.read_csv(plate_density_root / "selected_density_curve.csv")
    fig, ax = plt.subplots(figsize=(6.15, 3.35))
    fig.subplots_adjust(left=0.14, right=0.985, bottom=0.20, top=0.82)
    ax.plot(
        frame["response_m"],
        frame["reference_density_per_m"],
        color=TOKENS["ink"],
        linewidth=1.35,
        label="MC reference",
        zorder=3,
    )
    ax.plot(
        frame["response_m"],
        frame["pool_mean_density_per_m"],
        color=TOKENS["blue"],
        linewidth=1.25,
        linestyle="--",
        label="RQMC pool mean",
        zorder=3,
    )
    ax.set_xlabel(r"Plate-center deflection response (m)")
    ax.set_ylabel(r"Probability density (m$^{-1}$)")
    ax.ticklabel_format(axis="x", style="sci", scilimits=(-3, 3))
    legend_above(ax, ncol=2)
    save_pdf(fig, output_dir, "plate_density_reference_vs_pool")


def build_plate_diagnostics(plate_root: Path, output_dir: Path) -> None:
    active = pd.read_csv(plate_root / "_summary" / "active_grid_summary.csv")
    active = active.loc[
        (active["experiment"] == "E7_plate_SFEM")
        & (active["active_definition"] == "active_ref_abs_ge_1e-12")
    ].sort_values("h")
    fig, ax = plt.subplots(figsize=(6.4, 3.65))
    fig.subplots_adjust(left=0.13, right=0.86, bottom=0.19, top=0.97)
    ax.plot(
        active["h"],
        active["bootstrap_t_coverage"],
        color=TOKENS["blue"],
        marker="o",
        markerfacecolor="white",
        linewidth=1.15,
        markersize=4.2,
        label="bootstrap-t coverage",
    )
    ax.axhline(0.95, color=TOKENS["ink"], linestyle=":", linewidth=0.9)
    ax.set_xscale("log")
    ax.set_xlabel(r"Bandwidth $h$")
    ax.set_ylabel("Active-grid pointwise coverage")
    count_ax = ax.twinx()
    count_ax.plot(
        active["h"],
        active["active_grid_count"],
        color=TOKENS["red"],
        marker="s",
        markerfacecolor="white",
        linestyle="--",
        linewidth=1.0,
        markersize=3.8,
        label="active grid points",
    )
    count_ax.set_ylabel("Active grid-point count")
    count_ax.spines["top"].set_visible(False)
    lines = ax.get_lines()[:1] + count_ax.get_lines()
    labels = [line.get_label() for line in lines]
    ax.legend(lines, labels, frameon=False, fontsize=8, loc="best")
    save_pdf(fig, output_dir, "plate_e7_active_grid_coverage")

    band = pd.read_csv(plate_root / "_summary" / "active_grid_summary.csv")
    band = band.loc[
        (band["experiment"] == "E8_simultaneous_band")
        & (band["active_definition"] == "active_ref_abs_ge_1e-12")
        & band["simultaneous_band_coverage"].notna()
    ].sort_values("h")
    fig, ax = plt.subplots(figsize=(6.4, 3.65))
    fig.subplots_adjust(left=0.13, right=0.86, bottom=0.19, top=0.97)
    ax.plot(
        band["h"],
        band["simultaneous_band_coverage"],
        color=TOKENS["blue"],
        marker="o",
        markerfacecolor="white",
        linewidth=1.15,
        markersize=4.2,
        label="simultaneous coverage",
    )
    ax.axhline(0.95, color=TOKENS["ink"], linestyle=":", linewidth=0.9)
    ax.set_xscale("log")
    ax.set_xlabel(r"Bandwidth $h$")
    ax.set_ylabel("Active-grid simultaneous coverage")
    count_ax = ax.twinx()
    count_ax.plot(
        band["h"],
        band["active_grid_count"],
        color=TOKENS["red"],
        marker="s",
        markerfacecolor="white",
        linestyle="--",
        linewidth=1.0,
        markersize=3.8,
        label="active grid points",
    )
    count_ax.set_ylabel("Active grid-point count")
    count_ax.spines["top"].set_visible(False)
    lines = ax.get_lines()[:1] + count_ax.get_lines()
    labels = [line.get_label() for line in lines]
    ax.legend(lines, labels, frameon=False, fontsize=8, loc="best")
    save_pdf(fig, output_dir, "plate_e8_simultaneous_band")


def build_nr_coverage_figure(nr_root: Path, output_dir: Path) -> None:
    frame = pd.read_csv(nr_root / "nr_coverage_results.csv")
    frame = frame.loc[frame["target_mode"] == "end_to_end"].copy()
    models = [("linear", "Linear response"), ("nonlinear", "Nonlinear response")]
    methods = ["t distribution", "percentile bootstrap", "bootstrap-t"]
    n_values = sorted(frame["n"].unique())
    colors = [TOKENS["blue"], TOKENS["red"], TOKENS["green"], TOKENS["muted"]]
    fig, axes = plt.subplots(2, 3, figsize=(10.6, 6.2), sharex=True, sharey=True)
    fig.subplots_adjust(
        left=0.075, right=0.985, bottom=0.12, top=0.93, hspace=0.28, wspace=0.18
    )
    for row, (model, model_label) in enumerate(models):
        for col, method in enumerate(methods):
            ax = axes[row, col]
            subset = frame.loc[
                (frame["model"] == model) & (frame["method"] == method)
            ]
            for color, n_value in zip(colors, n_values):
                part = subset.loc[subset["n"] == n_value].sort_values("R")
                ax.errorbar(
                    part["R"],
                    part["coverage"],
                    yerr=1.96 * part["coverage_mcse"],
                    color=color,
                    marker="o",
                    markerfacecolor="white",
                    linewidth=1.05,
                    markersize=3.8,
                    capsize=2.2,
                    label=rf"$n={n_value}$",
                )
            ax.axhline(0.95, color=TOKENS["ink"], linestyle=":", linewidth=0.9)
            ax.set_ylim(0.60, 0.99)
            ax.set_xticks([16, 32, 64, 128])
            if row == 0:
                ax.set_title(DISPLAY_NAMES[method], fontsize=10.5)
            if row == 1:
                ax.set_xlabel(r"Outer curve count $R$")
            if col == 0:
                ax.set_ylabel(f"{model_label}\ncoverage")
            if row == 0 and col == 0:
                ax.legend(frameon=False, fontsize=7.5, loc="lower left")
    save_pdf(fig, output_dir, "nr_end_to_end_coverage")


def build_grid_size_band_figure(nr_root: Path, output_dir: Path) -> None:
    frame = pd.read_csv(nr_root / "finite_grid_band_results.csv")
    fig, ax = plt.subplots(figsize=(6.5, 3.75))
    fig.subplots_adjust(left=0.13, right=0.98, bottom=0.19, top=0.97)
    g_styles = {
        5: {"linestyle": "-", "markerfacecolor": "white"},
        9: {"linestyle": "--", "markerfacecolor": TOKENS["surface"]},
    }
    for method in ["percentile max-deviation band", "bootstrap-t max-stat band"]:
        for g_value in [5, 9]:
            part = frame.loc[
                (frame["method"] == method) & (frame["G"] == g_value)
            ].sort_values("R")
            style = METHOD_STYLES[method]
            ax.errorbar(
                part["R"],
                part["coverage"],
                yerr=1.96 * part["coverage_mcse"],
                color=style["color"],
                marker=style["marker"],
                markerfacecolor=g_styles[g_value]["markerfacecolor"],
                linestyle=g_styles[g_value]["linestyle"],
                linewidth=1.15,
                markersize=4.2,
                capsize=2.5,
                label=f"{DISPLAY_NAMES[method]}, $G={g_value}$",
            )
    ax.axhline(0.95, color=TOKENS["ink"], linestyle=":", linewidth=0.9)
    ax.set_xticks([16, 32, 64, 128])
    ax.set_ylim(0.90, 0.98)
    ax.set_xlabel(r"Outer curve count $R$")
    ax.set_ylabel("Finite-grid simultaneous coverage")
    ax.legend(frameon=False, fontsize=7.6, loc="upper center", bbox_to_anchor=(0.5, 1.18), ncol=2)
    save_pdf(fig, output_dir, "linear_finite_grid_G5_G9_coverage")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--formal-root", type=Path, required=True)
    parser.add_argument("--direct-root", type=Path, required=True)
    parser.add_argument("--direct-assets", type=Path, required=True)
    parser.add_argument("--plate-root", type=Path, required=True)
    parser.add_argument("--plate-density-root", type=Path, required=False)
    parser.add_argument("--nr-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    use_theme()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for old_file in args.output_dir.glob("*"):
        if old_file.is_file():
            old_file.unlink()

    build_pointwise_figure(
        args.formal_root
        / "evidence_suite"
        / "E2_scalar_normal_R_order_weighted"
        / "confirmation_summary.csv",
        args.output_dir,
        "linear_pointwise_coverage",
    )
    build_pointwise_figure(
        args.formal_root
        / "evidence_suite"
        / "E4_standard_normal_nonlinear_tail_weighted"
        / "confirmation_summary.csv",
        args.output_dir,
        "nonlinear_pointwise_coverage",
    )
    build_weighted_density_figure(args.formal_root, args.output_dir)
    build_grid_size_band_figure(args.nr_root, args.output_dir)
    for model in ["linear", "nonlinear"]:
        build_direct_coverage_figure(args.direct_assets, args.output_dir, model)
        build_direct_residual_figure(args.direct_assets, args.output_dir, model)
    build_beam_schematic(args.output_dir)
    build_plate_schematic(args.output_dir)
    if args.plate_density_root is not None:
        build_plate_density_reference_figure(args.plate_density_root, args.output_dir)


if __name__ == "__main__":
    main()
