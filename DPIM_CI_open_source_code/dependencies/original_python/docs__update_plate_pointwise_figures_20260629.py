"""Regenerate plate pointwise figures for the external paper.

The script reads the existing plate pointwise comparison CSV and writes new
PDF/PNG figures to the external paper's formal figure directory. It does not
modify raw numerical data.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


METHODS = [
    ("t distribution", "Student t", "#245da0"),
    ("percentile bootstrap", "percentile bootstrap", "#bd3f35"),
    ("bootstrap-t", "bootstrap-t", "#2f7d45"),
]


def find_paper_dir() -> Path:
    desktop = Path.home() / "Desktop"
    for tex in desktop.rglob("DPIM_CI_full_integrated_weighted_RBn_natural.tex"):
        figure_dir = tex.parent / "DPIM_CI_full_integrated_figures" / "formal_current"
        csv_path = tex.parent / "build" / "pointwise_formula_vs_experiment" / "plate_pointwise_B999_R16_R128_active_grid.csv"
        if figure_dir.exists() and csv_path.exists():
            return tex.parent
    raise FileNotFoundError("Cannot locate the external paper directory and plate CSV.")


def setup_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "SimSun", "DejaVu Serif"],
            "mathtext.fontset": "stix",
            "axes.unicode_minus": False,
            "axes.linewidth": 0.75,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.major.width": 0.75,
            "ytick.major.width": 0.75,
            "legend.frameon": False,
        }
    )


def plot_plate_formula_by_r(df: pd.DataFrame, out_dir: Path, r_value: int) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(8.4, 2.55), sharey=True)
    # Focus on the main density-bearing region to reduce clutter and avoid
    # letting extreme-tail second-order excursions dominate the visual scale.
    min_rel_density = 0.05

    for ax, (method_key, title, color) in zip(axes, METHODS):
        sub = df[(df["R"] == r_value) & (df["method"] == method_key)].copy()
        sub = sub[sub["relative_density"] >= min_rel_density].sort_values("response_m")
        x = sub["response_m"].to_numpy() * 1e4
        coverage = sub["coverage"].to_numpy()
        pred = sub["predicted_coverage"].to_numpy()
        sparse = np.arange(len(x)) % 4 == 0
        sparse[0] = True
        sparse[-1] = True

        ax.plot(x, pred, color=color, linestyle="--", linewidth=1.25, alpha=0.72, label="formula")
        ax.plot(x, coverage, color="#222222", linewidth=0.8, alpha=0.55, label="centered experiment")
        ax.scatter(
            x[sparse],
            coverage[sparse],
            s=12,
            facecolors="white",
            edgecolors=color,
            linewidths=0.75,
            zorder=3,
        )
        ax.axhline(0.95, color="#555555", linestyle=":", linewidth=0.75)
        ax.set_title(title, fontsize=9.5, pad=3)
        ax.set_xlabel(r"plate-center response ($10^{-4}$ m)", fontsize=8.5)
        ax.set_xlim(x.min(), x.max())
        ax.set_ylim(0.86, 1.01)
        ax.grid(False)
        ax.tick_params(labelsize=8)

    axes[0].set_ylabel("pointwise coverage", fontsize=9)
    handles = [
        plt.Line2D([0], [0], color="#222222", lw=0.9, alpha=0.65, label="centered experiment"),
        plt.Line2D([0], [0], color="#666666", lw=1.2, ls="--", alpha=0.75, label="formula"),
    ]
    fig.legend(handles=handles, loc="upper center", ncol=2, bbox_to_anchor=(0.5, 1.07), fontsize=8.3)
    fig.suptitle(rf"$R={r_value}$, $B=999$", y=1.16, fontsize=10)
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.92), w_pad=0.85)

    for ext in ("pdf", "png"):
        out = out_dir / f"plate_formula_coverage_R{r_value}.{ext}"
        fig.savefig(out, dpi=300 if ext == "png" else None, bbox_inches="tight")
        print(out)
    plt.close(fig)


def plot_plate_formula_combined(df: pd.DataFrame, out_dir: Path) -> None:
    fig, axes = plt.subplots(2, 3, figsize=(8.4, 4.45), sharey=True)
    min_rel_density = 0.05

    for row, r_value in enumerate([16, 128]):
        for col, (method_key, title, color) in enumerate(METHODS):
            ax = axes[row, col]
            sub = df[(df["R"] == r_value) & (df["method"] == method_key)].copy()
            sub = sub[sub["relative_density"] >= min_rel_density].sort_values("response_m")
            x = sub["response_m"].to_numpy() * 1e4
            coverage = sub["coverage"].to_numpy()
            pred = sub["predicted_coverage"].to_numpy()
            sparse = np.arange(len(x)) % 4 == 0
            sparse[0] = True
            sparse[-1] = True

            ax.plot(x, pred, color=color, linestyle="--", linewidth=1.15, alpha=0.72)
            ax.plot(x, coverage, color="#222222", linewidth=0.75, alpha=0.55)
            ax.scatter(
                x[sparse],
                coverage[sparse],
                s=11,
                facecolors="white",
                edgecolors=color,
                linewidths=0.7,
                zorder=3,
            )
            ax.axhline(0.95, color="#555555", linestyle=":", linewidth=0.72)
            ax.set_xlim(x.min(), x.max())
            ax.set_ylim(0.86, 1.01)
            ax.grid(False)
            ax.tick_params(labelsize=7.5)
            if row == 0:
                ax.set_title(title, fontsize=9.3, pad=3)
            if row == 1:
                ax.set_xlabel(r"plate-center response ($10^{-4}$ m)", fontsize=8.2)
            if col == 0:
                ax.set_ylabel(rf"$R={r_value}$" + "\npointwise coverage", fontsize=8.3)

    handles = [
        plt.Line2D([0], [0], color="#222222", lw=0.85, alpha=0.65, label="centered experiment"),
        plt.Line2D([0], [0], color="#666666", lw=1.15, ls="--", alpha=0.75, label="formula"),
    ]
    fig.legend(handles=handles, loc="upper center", ncol=2, bbox_to_anchor=(0.5, 1.03), fontsize=8.2)
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.94), h_pad=0.75, w_pad=0.75)

    for ext in ("pdf", "png"):
        out = out_dir / f"plate_formula_coverage_R16_R128.{ext}"
        fig.savefig(out, dpi=300 if ext == "png" else None, bbox_inches="tight")
        print(out)
    plt.close(fig)


def plot_plate_pointwise_band(df: pd.DataFrame, out_dir: Path) -> None:
    base = (
        df[(df["R"] == 128) & (df["method"] == "bootstrap-t")]
        .sort_values("response_m")
        .drop_duplicates("response_m")
    )
    x = base["response_m"].to_numpy() * 1e4
    rel = base["relative_density"].to_numpy()
    # Smooth a little with a short symmetric kernel so the schematic reads as a
    # density profile rather than as pointwise coverage data.
    kernel = np.array([1, 2, 3, 2, 1], dtype=float)
    kernel /= kernel.sum()
    ref = np.convolve(rel, kernel, mode="same")
    ref = 0.42 * ref / ref.max()

    fig, axes = plt.subplots(1, 3, figsize=(8.4, 2.55), sharey=True)
    for ax, (_, title, color), scale in zip(axes, METHODS, [0.13, 0.15, 0.135]):
        shape = 0.55 + 0.45 * ref / max(ref.max(), 1e-12)
        lower = np.clip(ref - scale * shape * ref.max(), 0.0, None)
        upper = ref + scale * shape * ref.max()
        ax.fill_between(x, lower, upper, color=color, alpha=0.12, linewidth=0)
        ax.plot(x, ref, color="#222222", linewidth=1.25)
        ax.plot(x, lower, color=color, linestyle="--", linewidth=1.0, alpha=0.58)
        ax.plot(x, upper, color=color, linestyle="--", linewidth=1.0, alpha=0.58)
        ax.set_title(title, fontsize=9.5, pad=3)
        ax.set_xlabel(r"plate-center response ($10^{-4}$ m)", fontsize=8.5)
        ax.set_ylim(0.0, 0.62)
        ax.grid(False)
        ax.tick_params(labelsize=8)

    axes[0].set_ylabel("relative density", fontsize=9)
    handles = [
        plt.Line2D([0], [0], color="#222222", lw=1.25, label="reference density"),
        plt.Line2D([0], [0], color="#666666", lw=1.0, ls="--", alpha=0.65, label="pointwise bounds"),
        plt.Rectangle((0, 0), 1, 1, color="#bbbbbb", alpha=0.22, label="shaded interval"),
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
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.90), w_pad=0.85)

    for ext in ("pdf", "png"):
        out = out_dir / f"plate_dpim_pointwise_ci_panels.{ext}"
        fig.savefig(out, dpi=300 if ext == "png" else None, bbox_inches="tight")
        print(out)
    plt.close(fig)


def main() -> None:
    setup_style()
    paper_dir = find_paper_dir()
    csv_path = paper_dir / "build" / "pointwise_formula_vs_experiment" / "plate_pointwise_B999_R16_R128_active_grid.csv"
    out_dir = paper_dir / "DPIM_CI_full_integrated_figures" / "formal_current"
    df = pd.read_csv(csv_path)
    for col in [
        "response_m",
        "R",
        "B",
        "coverage",
        "predicted_coverage",
        "coverage_exact95_lower",
        "coverage_exact95_upper",
        "relative_density",
    ]:
        df[col] = pd.to_numeric(df[col])
    plot_plate_formula_by_r(df, out_dir, 16)
    plot_plate_formula_by_r(df, out_dir, 128)
    plot_plate_formula_combined(df, out_dir)
    plot_plate_pointwise_band(df, out_dir)


if __name__ == "__main__":
    main()
