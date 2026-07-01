from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


MULTIY_ROOT = Path(
    r"C:\Users\Wangg\Desktop\biyeshejimatlab_portable_e1e8\results"
    r"\linear_nonlinear_multiy_centered_paper_20260626_202142"
)
PLATE_ROOT = Path(
    r"C:\Users\Wangg\Desktop\biyeshejimatlab_portable_e1e8\results"
    r"\plate_original_formula_validation_20260624_figurestyle"
)
FIG_ROOT = Path(
    r"C:\Users\Wangg\Desktop\期刊论文\DPIM_CI_full_integrated_figures"
    r"\formal_current"
)

METHODS = ["t distribution", "percentile bootstrap", "bootstrap-t"]
METHOD_LABELS = {
    "t distribution": "Student t",
    "percentile bootstrap": "percentile bootstrap",
    "bootstrap-t": "bootstrap-t",
}
COLORS = {16: "#1F4E79", 128: "#B22222"}


def set_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "STIXGeneral", "DejaVu Serif"],
            "mathtext.fontset": "stix",
            "axes.edgecolor": "#1A1A1A",
            "axes.linewidth": 0.75,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def plot_profile_panel(ax, data: pd.DataFrame, x_col: str, y_label: str, method: str) -> None:
    sub = data[data["method"] == method].copy()
    for r in [16, 128]:
        sr = sub[sub["R"] == r].sort_values(x_col)
        if sr.empty:
            continue
        x = sr[x_col].to_numpy()
        ax.fill_between(
            x,
            sr["coverage_exact95_lower"].to_numpy(),
            sr["coverage_exact95_upper"].to_numpy(),
            color=COLORS[r],
            alpha=0.10,
            linewidth=0,
        )
        ax.plot(
            x,
            sr["coverage"].to_numpy(),
            color=COLORS[r],
            linewidth=1.15,
            marker="o",
            markersize=2.2,
            label=f"{METHOD_LABELS[method]}, R={r}",
        )
        ax.plot(
            x,
            sr["predicted_coverage"].to_numpy(),
            color=COLORS[r],
            linewidth=1.05,
            linestyle="--",
        )
    ax.axhline(0.95, color="#666666", linewidth=0.70, linestyle=":")
    ax.set_title(METHOD_LABELS[method], fontsize=9)
    ax.set_xlabel(y_label)
    ax.set_ylabel("pointwise coverage")
    ax.grid(True, color="#D9D9D9", linewidth=0.35, alpha=0.70)


def rebuild_multiy_profiles() -> None:
    data = pd.read_csv(MULTIY_ROOT / "multiy_formula_predictions.csv")
    for model, xlim, outfile in [
        ("linear", (-3.05, 3.05), "linear_multiy_centered_formula_profiles.pdf"),
        ("nonlinear", (1.45, 3.05), "nonlinear_multiy_centered_formula_profiles.pdf"),
    ]:
        sub = data[(data["model"] == model) & (data["B"] == 399)].copy()
        fig, axes = plt.subplots(1, 3, figsize=(7.0, 2.55), sharey=True)
        for ax, method in zip(axes, METHODS):
            plot_profile_panel(ax, sub, "y", "response coordinate $y$", method)
            ax.set_xlim(*xlim)
            ax.set_ylim(0.84, 1.005)
        handles, labels = axes[0].get_legend_handles_labels()
        fig.legend(handles, labels, loc="upper center", ncol=2, frameon=False, fontsize=8)
        fig.text(
            0.52,
            0.035,
            "solid line: centered experiment; dashed line: factorized formula; shaded band: exact 95% interval",
            ha="center",
            fontsize=8,
        )
        fig.tight_layout(rect=[0.02, 0.12, 0.995, 0.84], w_pad=1.1)
        fig.savefig(FIG_ROOT / outfile, bbox_inches="tight")
        plt.close(fig)


def rebuild_plate_profiles() -> None:
    data = pd.read_csv(PLATE_ROOT / "original_formula_predictions.csv")
    sub = data[(data["B"] == 999) & (data["active_grid"] == 1)].copy()
    fig, axes = plt.subplots(1, 3, figsize=(7.0, 2.55), sharey=True)
    for ax, method in zip(axes, METHODS):
        plot_profile_panel(
            ax,
            sub[sub["method"] == method],
            "response_m",
            "plate-center deflection response (m)",
            method,
        )
        ax.set_ylim(0.58, 1.02)
    handles, labels = axes[0].get_legend_handles_labels()
    labels = [label + ", B=999" for label in labels]
    fig.legend(handles, labels, loc="upper center", ncol=2, frameon=False, fontsize=8)
    fig.text(
        0.52,
        0.035,
        "solid line: centered experiment; dashed line: factorized formula; shaded band: exact 95% interval",
        ha="center",
        fontsize=8,
    )
    fig.tight_layout(rect=[0.02, 0.12, 0.995, 0.84], w_pad=1.1)
    fig.savefig(FIG_ROOT / "plate_original_formula_coverage_profiles.pdf", bbox_inches="tight")
    plt.close(fig)


def rebuild_plate_region_rmse() -> None:
    data = pd.read_csv(PLATE_ROOT / "original_formula_density_region_summary.csv")
    pivot = data.pivot_table(
        index="relative_density_threshold",
        columns="method",
        values="root_mean_square_error",
        aggfunc="mean",
    ).sort_index()
    x = np.arange(len(pivot.index))
    width = 0.22
    fig, ax = plt.subplots(figsize=(4.8, 2.55))
    for i, method in enumerate(METHODS):
        ax.bar(
            x + (i - 1) * width,
            pivot[method].to_numpy(),
            width=width,
            color=["#1F4E79", "#B22222", "#2E7D32"][i],
            label=METHOD_LABELS[method],
            edgecolor="#1A1A1A",
            linewidth=0.35,
        )
    ax.set_xticks(x)
    ax.set_xticklabels([f"{v:g}" for v in pivot.index])
    ax.set_xlabel("density threshold")
    ax.set_ylabel("mean pointwise coverage RMSE")
    ax.grid(True, axis="y", color="#D9D9D9", linewidth=0.35, alpha=0.70)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, 1.16), ncol=3, frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG_ROOT / "plate_original_formula_region_rmse.pdf", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    set_style()
    FIG_ROOT.mkdir(parents=True, exist_ok=True)
    rebuild_multiy_profiles()
    rebuild_plate_profiles()
    rebuild_plate_region_rmse()


if __name__ == "__main__":
    main()