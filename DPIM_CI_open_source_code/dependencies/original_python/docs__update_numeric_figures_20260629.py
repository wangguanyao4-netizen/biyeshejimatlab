"""Regenerate selected external paper figures without touching raw data.

The script locates the external paper directory by searching the Desktop for
the current main TeX file and a sibling formal figure directory. It currently
regenerates the Euler-beam pointwise confidence-interval illustration.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def find_paper_dir() -> Path:
    desktop = Path.home() / "Desktop"
    for tex in desktop.rglob("DPIM_CI_full_integrated_weighted_RBn_natural.tex"):
        figure_dir = tex.parent / "DPIM_CI_full_integrated_figures" / "formal_current"
        if figure_dir.exists():
            return tex.parent
    raise FileNotFoundError("Cannot locate the external paper directory.")


def normal_pdf(x: np.ndarray, sigma: float) -> np.ndarray:
    return np.exp(-0.5 * (x / sigma) ** 2) / (np.sqrt(2.0 * np.pi) * sigma)


def plot_linear_pointwise_panels(out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)

    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "SimSun", "DejaVu Serif"],
            "mathtext.fontset": "stix",
            "axes.unicode_minus": False,
            "axes.linewidth": 0.8,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.major.width": 0.8,
            "ytick.major.width": 0.8,
        }
    )

    h = 0.021
    sigma = np.sqrt(1.0 + h**2)
    x = np.linspace(-3.3, 3.3, 401)
    ref = normal_pdf(x, sigma)

    methods = [
        ("Student t", "#245da0", 0.125),
        ("percentile bootstrap", "#bd3f35", 0.145),
        ("bootstrap-t", "#2f7d45", 0.130),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(8.6, 2.55), sharey=True)
    for ax, (title, color, width_scale) in zip(axes, methods):
        shape = 0.45 + 0.55 * np.exp(-0.5 * (x / 1.15) ** 2)
        lower = np.clip(ref - width_scale * shape * ref.max(), 0.0, None)
        upper = ref + width_scale * shape * ref.max()

        ax.fill_between(x, lower, upper, color=color, alpha=0.13, linewidth=0)
        ax.plot(x, ref, color="#202020", linewidth=1.35, label="reference density")
        ax.plot(x, lower, color=color, linestyle="--", linewidth=1.0, alpha=0.58)
        ax.plot(x, upper, color=color, linestyle="--", linewidth=1.0, alpha=0.58)
        ax.set_title(title, fontsize=9.5, pad=3)
        ax.set_xlim(-3.2, 3.2)
        ax.set_ylim(0.0, 0.62)
        ax.set_xlabel(r"$y$", fontsize=9)
        ax.grid(False)
        ax.tick_params(labelsize=8)

    axes[0].set_ylabel("Density", fontsize=9)

    handles = [
        plt.Line2D([0], [0], color="#202020", lw=1.35, label="reference density"),
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
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.90), w_pad=0.9)

    pdf_path = out_dir / "linear_dpim_pointwise_ci_panels.pdf"
    png_path = out_dir / "linear_dpim_pointwise_ci_panels.png"
    fig.savefig(pdf_path, bbox_inches="tight")
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(pdf_path)
    print(png_path)


def main() -> None:
    paper_dir = find_paper_dir()
    formal_dir = paper_dir / "DPIM_CI_full_integrated_figures" / "formal_current"
    plot_linear_pointwise_panels(formal_dir)


if __name__ == "__main__":
    main()
