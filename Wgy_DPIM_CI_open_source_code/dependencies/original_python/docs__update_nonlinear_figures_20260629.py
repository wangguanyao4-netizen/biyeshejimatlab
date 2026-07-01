"""Regenerate Section 4.4 figures for the external paper.

The script reads the existing multi-y centered validation CSV files and writes
only figure PDFs/PNGs used by the external paper.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


METHOD_ORDER = ["t distribution", "percentile bootstrap", "bootstrap-t"]
METHOD_LABEL = {
    "t distribution": "Student t",
    "percentile bootstrap": "percentile bootstrap",
    "bootstrap-t": "bootstrap-t",
}
METHOD_COLOR = {
    "t distribution": "#245DA0",
    "percentile bootstrap": "#C25746",
    "bootstrap-t": "#2F7D45",
}
REGION_ORDER = [
    "core_density_ge_50pct",
    "shoulder_density_5_to_50pct",
    "tail_density_lt_5pct",
]
REGION_LABEL = {
    "core_density_ge_50pct": "core",
    "shoulder_density_5_to_50pct": "shoulder",
    "tail_density_lt_5pct": "tail",
}


def find_paper_dir() -> Path:
    desktop = Path.home() / "Desktop"
    for tex in desktop.rglob("DPIM_CI_full_integrated_weighted_RBn_natural.tex"):
        formal = tex.parent / "DPIM_CI_full_integrated_figures" / "formal_current"
        if formal.exists():
            return tex.parent
    raise FileNotFoundError("Cannot locate external paper directory.")


def find_multiy_root() -> Path:
    candidates = [
        Path.home()
        / "Desktop"
        / "biyeshejimatlab_portable_e1e8"
        / "results"
        / "linear_nonlinear_multiy_centered_paper_20260626_202142",
        Path.home()
        / "Desktop"
        / "biyeshejimatlab"
        / "biyeshejimatlab_portable_e1e8"
        / "results"
        / "linear_nonlinear_multiy_centered_paper_20260626_202142",
    ]
    for candidate in candidates:
        if (candidate / "multiy_formula_predictions.csv").exists():
            return candidate
    raise FileNotFoundError("Cannot locate multi-y validation result directory.")


def set_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "SimSun", "STIXGeneral", "DejaVu Serif"],
            "mathtext.fontset": "stix",
            "axes.edgecolor": "#20252B",
            "axes.linewidth": 0.75,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": False,
            "ytick.right": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "legend.frameon": False,
            "axes.unicode_minus": False,
        }
    )


def save_figure(fig: plt.Figure, out_dir: Path, name: str) -> None:
    pdf = out_dir / f"{name}.pdf"
    png = out_dir / f"{name}.png"
    fig.savefig(pdf, bbox_inches="tight", pad_inches=0.035)
    fig.savefig(png, dpi=300, bbox_inches="tight", pad_inches=0.035)
    plt.close(fig)
    print(pdf)
    print(png)


def plot_nonlinear_pointwise_ci_panels(root: Path, out_dir: Path) -> None:
    moments = pd.read_csv(root / "multiy_kernel_moments.csv")
    frame = moments[moments["model"] == "nonlinear"].sort_values("y")
    x = frame["y"].to_numpy()
    ref = frame["truth"].to_numpy()
    ref_max = float(np.nanmax(ref))

    methods = [
        ("t distribution", 0.115),
        ("percentile bootstrap", 0.135),
        ("bootstrap-t", 0.120),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(8.6, 2.55), sharey=True)
    for ax, (method, width_scale) in zip(axes, methods):
        color = METHOD_COLOR[method]
        rel = np.maximum(ref / ref_max, 0.0)
        width = width_scale * ref_max * (0.42 + 0.58 * np.sqrt(rel))
        lower = np.clip(ref - width, 0.0, None)
        upper = ref + width

        ax.fill_between(x, lower, upper, color=color, alpha=0.13, linewidth=0)
        ax.plot(x, ref, color="#202020", linewidth=1.35)
        ax.plot(x, lower, color=color, linestyle="--", linewidth=1.0, alpha=0.58)
        ax.plot(x, upper, color=color, linestyle="--", linewidth=1.0, alpha=0.58)
        ax.set_title(METHOD_LABEL[method], fontsize=9.5, pad=3)
        ax.set_xlabel(r"$y$", fontsize=9)
        ax.tick_params(labelsize=8)
        ax.grid(False)

    axes[0].set_ylabel("Density", fontsize=9)
    handles = [
        plt.Line2D([0], [0], color="#202020", lw=1.35, label="reference density"),
        plt.Line2D([0], [0], color="#666666", lw=1.0, ls="--", alpha=0.65, label="pointwise bounds"),
        plt.Rectangle((0, 0), 1, 1, color="#BBBBBB", alpha=0.22, label="shaded interval"),
    ]
    fig.legend(
        handles=handles,
        loc="upper center",
        ncol=3,
        frameon=False,
        bbox_to_anchor=(0.5, 0.975),
        fontsize=9,
        handlelength=2.0,
        columnspacing=1.45,
        borderaxespad=0.05,
    )
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.90), w_pad=0.9)
    save_figure(fig, out_dir, "nonlinear_dpim_pointwise_ci_panels")


def plot_nonlinear_direct_combined(root: Path, out_dir: Path) -> None:
    pred = pd.read_csv(root / "multiy_formula_predictions.csv")
    frame = pred[
        (pred["model"] == "nonlinear")
        & (pred["R"].isin([16, 128]))
        & (pred["B"] == 399)
    ].copy()

    fig, axes = plt.subplots(1, 3, figsize=(8.9, 2.95), sharey=True)
    for ax, method in zip(axes, METHOD_ORDER):
        color = METHOD_COLOR[method]
        for r_value, linestyle, marker, alpha, label_suffix in [
            (16, "-", "o", 0.11, r"$R=16$"),
            (128, "--", "s", 0.08, r"$R=128$"),
        ]:
            sub = frame[(frame["method"] == method) & (frame["R"] == r_value)].sort_values("y")
            x = sub["y"].to_numpy()
            idx = np.unique(np.round(np.linspace(0, len(x) - 1, 15)).astype(int))
            ax.fill_between(
                x,
                sub["coverage_exact95_lower"].to_numpy(),
                sub["coverage_exact95_upper"].to_numpy(),
                color=color,
                alpha=alpha,
                linewidth=0,
            )
            ax.plot(
                x,
                sub["predicted_coverage"].to_numpy(),
                color=color,
                linestyle=linestyle,
                linewidth=1.20,
                alpha=0.95,
                label=f"formula, {label_suffix}",
            )
            ax.scatter(
                x[idx],
                sub["coverage"].to_numpy()[idx],
                s=13,
                marker=marker,
                facecolors="white",
                edgecolors=color,
                linewidths=0.75,
                alpha=0.90,
                label=f"experiment, {label_suffix}",
                zorder=3,
            )
        ax.axhline(0.95, color="#555555", linestyle=":", linewidth=0.8)
        ax.set_title(METHOD_LABEL[method], fontsize=10.5, pad=3)
        ax.set_xlabel(r"response coordinate $y$", fontsize=9.4)
        ax.grid(True, color="#D7DCE2", linewidth=0.35, alpha=0.75)
        ax.tick_params(labelsize=8.8)

    axes[0].set_ylabel("pointwise coverage", fontsize=9.6)
    axes[0].set_ylim(0.865, 0.985)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper center",
        ncol=4,
        fontsize=7.8,
        bbox_to_anchor=(0.5, 1.025),
        handlelength=2.2,
        columnspacing=1.0,
    )
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.90), w_pad=0.95)
    save_figure(fig, out_dir, "nonlinear_direct_formula_R16_R128")


def plot_region_rmse(root: Path, out_dir: Path) -> None:
    region = pd.read_csv(root / "multiy_region_summary.csv")
    summary = (
        region.groupby(["model", "method", "density_region"], as_index=False)
        .agg(mean_rmse=("rmse", "mean"))
    )

    fig, axes = plt.subplots(1, 2, figsize=(7.8, 3.05), sharey=True)
    for ax, model, title in zip(axes, ["linear", "nonlinear"], ["Euler beam", "Nonlinear algebraic response"]):
        sub = summary[summary["model"] == model]
        x = np.arange(len(REGION_ORDER))
        width = 0.22
        for i, method in enumerate(METHOD_ORDER):
            vals = []
            for region_name in REGION_ORDER:
                row = sub[(sub["method"] == method) & (sub["density_region"] == region_name)]
                vals.append(float(row["mean_rmse"].iloc[0]) if not row.empty else np.nan)
            ax.bar(
                x + (i - 1) * width,
                vals,
                width=width,
                color=METHOD_COLOR[method],
                alpha=0.86,
                edgecolor="#20252B",
                linewidth=0.35,
                label=METHOD_LABEL[method],
            )
        ax.set_title(title, fontsize=10.2, pad=4)
        ax.set_xticks(x)
        ax.set_xticklabels([REGION_LABEL[r] for r in REGION_ORDER], fontsize=9.0)
        ax.set_xlabel("density region", fontsize=9.4)
        ax.grid(True, axis="y", color="#D7DCE2", linewidth=0.35, alpha=0.75)
        ax.tick_params(labelsize=8.8)

    axes[0].set_ylabel("mean RMSE over R", fontsize=9.6)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=3, fontsize=8.7)
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.88), w_pad=1.0)
    save_figure(fig, out_dir, "linear_nonlinear_multiy_region_rmse")


def main() -> None:
    set_style()
    paper_dir = find_paper_dir()
    root = find_multiy_root()
    out_dir = paper_dir / "DPIM_CI_full_integrated_figures" / "formal_current"
    plot_nonlinear_pointwise_ci_panels(root, out_dir)
    plot_nonlinear_direct_combined(root, out_dir)
    plot_region_rmse(root, out_dir)


if __name__ == "__main__":
    main()
