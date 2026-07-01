from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


SOURCE = Path(
    r"C:\Users\Wangg\Desktop\biyeshejimatlab_portable_e1e8\results"
    r"\linear_nonlinear_multiy_centered_paper_20260626_202142"
    r"\multiy_region_summary.csv"
)
OUTPUT = Path(
    r"C:\Users\Wangg\Desktop\期刊论文\DPIM_CI_full_integrated_figures"
    r"\formal_current\linear_nonlinear_multiy_region_rmse.pdf"
)


METHOD_LABELS = {
    "t distribution": "Student t",
    "percentile bootstrap": "percentile bootstrap",
    "bootstrap-t": "bootstrap-t",
}

REGION_ORDER = [
    "core_density_ge_50pct",
    "shoulder_density_5_to_50pct",
    "tail_density_lt_5pct",
]
REGION_LABELS = {
    "core_density_ge_50pct": "core",
    "shoulder_density_5_to_50pct": "shoulder",
    "tail_density_lt_5pct": "tail",
}

COLORS = {
    ("linear", "t distribution"): "#1F4E79",
    ("linear", "percentile bootstrap"): "#B22222",
    ("linear", "bootstrap-t"): "#2E7D32",
    ("nonlinear", "t distribution"): "#8FB7D9",
    ("nonlinear", "percentile bootstrap"): "#E69595",
    ("nonlinear", "bootstrap-t"): "#A8D5A2",
}
HATCHES = {
    "linear": "",
    "nonlinear": "///",
}


def main() -> None:
    frame = pd.read_csv(SOURCE)
    summary = (
        frame.groupby(["model", "method", "density_region"], as_index=False)
        .agg(mean_rmse=("rmse", "mean"))
    )

    series = [
        ("linear", "t distribution"),
        ("linear", "percentile bootstrap"),
        ("linear", "bootstrap-t"),
        ("nonlinear", "t distribution"),
        ("nonlinear", "percentile bootstrap"),
        ("nonlinear", "bootstrap-t"),
    ]

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
        }
    )

    fig, ax = plt.subplots(figsize=(6.6, 3.25))
    fig.subplots_adjust(left=0.13, right=0.985, bottom=0.17, top=0.70)

    x = np.arange(len(REGION_ORDER))
    width = 0.115
    offsets = (np.arange(len(series)) - (len(series) - 1) / 2) * width

    for offset, (model, method) in zip(offsets, series):
        values = []
        for region in REGION_ORDER:
            row = summary[
                (summary["model"] == model)
                & (summary["method"] == method)
                & (summary["density_region"] == region)
            ]
            values.append(float(row["mean_rmse"].iloc[0]))
        ax.bar(
            x + offset,
            values,
            width=width,
            label=f"{model}, {METHOD_LABELS[method]}",
            color=COLORS[(model, method)],
            edgecolor="#333333",
            linewidth=0.45,
            hatch=HATCHES[model],
        )

    ax.set_xticks(x)
    ax.set_xticklabels([REGION_LABELS[r] for r in REGION_ORDER])
    ax.set_ylabel(r"mean RMSE over $R$")
    ax.set_ylim(0.0, 0.0148)
    ax.grid(False)
    handles, labels = ax.get_legend_handles_labels()
    legend_order = [0, 3, 1, 4, 2, 5]
    ax.legend(
        [handles[i] for i in legend_order],
        [labels[i] for i in legend_order],
        frameon=False,
        loc="lower center",
        bbox_to_anchor=(0.5, 1.04),
        ncol=3,
        fontsize=7.2,
        handlelength=1.9,
        columnspacing=1.05,
    )

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT, bbox_inches="tight", pad_inches=0.025)
    plt.close(fig)


if __name__ == "__main__":
    main()