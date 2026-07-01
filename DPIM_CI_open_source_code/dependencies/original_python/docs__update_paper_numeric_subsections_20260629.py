from __future__ import annotations

import json
import shutil
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import norm, t


PAPER_ROOT = Path(r"C:\Users\Wangg\Desktop\期刊论文")
FIG_DIR = PAPER_ROOT / "DPIM_CI_full_integrated_figures" / "formal_current"
BACKUP_DIR = PAPER_ROOT / "build" / "codex_backups" / "numeric_subsections_20260629"
OUT_DIR = PAPER_ROOT / "build" / "numeric_subsections_20260629"
DATA_ROOT = Path(
    r"C:\Users\Wangg\Desktop\biyeshejimatlab_portable_e1e8\results"
    r"\linear_nonlinear_multiy_centered_paper_20260626_202142"
)

METHODS = ["t distribution", "percentile bootstrap", "bootstrap-t"]
METHOD_LABEL = {
    "t distribution": "Student t",
    "percentile bootstrap": "percentile bootstrap",
    "bootstrap-t": "bootstrap-t",
}
METHOD_COLOR = {
    "t distribution": "#1f5a9d",
    "percentile bootstrap": "#b43c35",
    "bootstrap-t": "#2d7f46",
}


def ensure_dirs() -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    BACKUP_DIR.mkdir(parents=True, exist_ok=True)
    OUT_DIR.mkdir(parents=True, exist_ok=True)


def backup_existing(names: list[str]) -> None:
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    for name in names:
        src = FIG_DIR / name
        if src.exists():
            shutil.copy2(src, BACKUP_DIR / f"{src.stem}.before_{stamp}{src.suffix}")


def set_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "SimSun", "DejaVu Serif"],
            "mathtext.fontset": "stix",
            "axes.linewidth": 0.7,
            "axes.labelsize": 8,
            "axes.titlesize": 7.5,
            "xtick.labelsize": 7,
            "ytick.labelsize": 7,
            "legend.fontsize": 7,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def gaussian_kernel_ordinates(response: np.ndarray, weights: np.ndarray, y_grid: np.ndarray, h: float) -> np.ndarray:
    scaled = (y_grid[None, None, :] - response[:, :, None]) / h
    kernels = np.exp(-0.5 * scaled**2) / (np.sqrt(2.0 * np.pi) * h)
    return np.sum(weights[:, :, None] * kernels, axis=1)


def bootstrap_mean_ci(q: np.ndarray, method: str, rng: np.random.Generator, b_count: int = 399) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    r_count, g_count = q.shape
    mean = q.mean(axis=0)
    sd = q.std(axis=0, ddof=1)
    se = sd / np.sqrt(r_count)
    if method == "t distribution":
        z = t.ppf(0.975, r_count - 1)
        return mean, mean - z * se, mean + z * se

    idx = rng.integers(0, r_count, size=(b_count, r_count))
    boot = q[idx, :].mean(axis=1)
    if method == "percentile bootstrap":
        return mean, np.quantile(boot, 0.025, axis=0), np.quantile(boot, 0.975, axis=0)

    boot_sd = q[idx, :].std(axis=1, ddof=1)
    boot_se = boot_sd / np.sqrt(r_count)
    with np.errstate(divide="ignore", invalid="ignore"):
        t_star = (boot - mean[None, :]) / boot_se
    t_star = np.where(np.isfinite(t_star), t_star, np.nan)
    q_lo = np.nanquantile(t_star, 0.025, axis=0)
    q_hi = np.nanquantile(t_star, 0.975, axis=0)
    return mean, mean - q_hi * se, mean - q_lo * se


def plot_pointwise_interval_panels() -> Path:
    raw = np.load(DATA_ROOT / "curve_pool_raw.npz")
    weights = raw["weights"][:16, :]
    response = raw["linear_response"][:16, :]
    y_grid = np.r_[
        np.linspace(-3.6, -1.8, 5),
        np.linspace(-1.5, -0.5, 6),
        np.linspace(-0.4, 0.4, 9),
        np.linspace(0.5, 1.5, 6),
        np.linspace(1.8, 3.6, 5),
    ]
    h = 0.021
    rng = np.random.default_rng(20260629)

    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.35), sharex=True, sharey=True)
    for row, method in enumerate(METHODS):
        ax = axes[row]
        q = gaussian_kernel_ordinates(response, weights, y_grid, h)
        mean, low, high = bootstrap_mean_ci(q, method, rng)
        truth = norm.pdf(y_grid, loc=0.0, scale=np.sqrt(1.0 + h**2))
        low_plot = np.where((low >= -0.02) & (low <= 0.64), low, np.nan)
        high_plot = np.where((high >= -0.02) & (high <= 0.64), high, np.nan)
        ax.fill_between(y_grid, low_plot, high_plot, color=METHOD_COLOR[method], alpha=0.18, linewidth=0)
        ax.plot(y_grid, low_plot, color=METHOD_COLOR[method], linestyle="--", linewidth=0.95, marker="o", markersize=1.7)
        ax.plot(y_grid, high_plot, color=METHOD_COLOR[method], linestyle="--", linewidth=0.95, marker="^", markersize=1.7)
        ax.plot(y_grid, truth, color="#222222", linewidth=1.15, marker="s", markersize=1.7)
        ax.set_title(METHOD_LABEL[method], pad=2)
        ax.set_xlim(-3.8, 3.8)
        ax.set_ylim(0.0, 0.62)
        if row == 0:
            ax.set_ylabel("Density")
        ax.set_xlabel(r"$y$")
        ax.tick_params(direction="in", length=2.5, width=0.6)
        for spine in ax.spines.values():
            spine.set_color("#333333")
    handles = [
        plt.Line2D([0], [0], color="#222222", marker="s", markersize=2, linewidth=1.1, label="reference density"),
        plt.Line2D([0], [0], color="#555555", linestyle="--", marker="o", markersize=2, linewidth=0.9, label="pointwise bounds"),
        plt.Rectangle((0, 0), 1, 1, facecolor="#777777", alpha=0.18, edgecolor="none", label="shaded interval"),
    ]
    fig.legend(
        handles=handles,
        loc="upper center",
        ncol=3,
        frameon=False,
        bbox_to_anchor=(0.5, 0.975),
        fontsize=8,
        borderaxespad=0.05,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.90), w_pad=0.75)
    out = FIG_DIR / "linear_dpim_pointwise_ci_panels.pdf"
    fig.savefig(out, bbox_inches="tight")
    fig.savefig(out.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    return out


def representative_y_for_r(df_r: pd.DataFrame, max_points: int = 23) -> np.ndarray:
    score = (
        df_r.groupby("y", as_index=False)
        .agg(score=("residual", lambda s: np.mean(np.abs(s))), density=("relative_density", "mean"))
        .sort_values("y")
    )
    bins = pd.cut(score["y"], bins=np.linspace(score["y"].min(), score["y"].max(), 8), include_lowest=True)
    chosen: list[float] = []
    for _, part in score.groupby(bins, observed=True):
        part = part.sort_values(["score", "density"], ascending=[True, False])
        chosen.extend(part.head(3)["y"].tolist())
    for val in [0.0, -0.75, 0.75]:
        closest = score.iloc[(score["y"] - val).abs().argsort()[:1]]["y"].iloc[0]
        chosen.append(float(closest))
    selected = np.array(sorted(set(np.round(chosen, 12))))
    if selected.size > max_points:
        sub = score[score["y"].isin(selected)].sort_values("score").head(max_points)
        selected = np.array(sorted(sub["y"].to_numpy()))
    return selected


def plot_formula_by_r_and_write_table() -> tuple[Path, Path]:
    df = pd.read_csv(DATA_ROOT / "multiy_formula_predictions.csv")
    df = df[(df["model"] == "linear") & (df["B"] == 399)].copy()

    selected_r = [32, 128]
    fig, axes = plt.subplots(2, 3, figsize=(7.2, 4.25), sharey=True)
    for i, r_value in enumerate(selected_r):
        df_r = df[df["R"] == r_value].copy()
        selected_y = representative_y_for_r(df_r, max_points=18)
        for j, method in enumerate(METHODS):
            ax = axes[i, j]
            sub = df_r[(df_r["method"] == method) & (df_r["y"].isin(selected_y))].sort_values("y")
            color = METHOD_COLOR[method]
            ax.plot(sub["y"], sub["predicted_coverage"], color=color, linewidth=1.15, label="formula")
            ax.scatter(sub["y"], sub["coverage"], s=15, facecolor="white", edgecolor=color, linewidth=0.8, label="centered experiment", zorder=3)
            ax.axhline(0.95, color="#666666", linestyle=":", linewidth=0.75)
            if i == 0:
                ax.set_title(METHOD_LABEL[method])
            if j == 0:
                ax.set_ylabel(f"R={r_value}\nCentered coverage")
            ax.set_xlabel(r"$y$")
            ax.set_ylim(0.895, 0.980)
            ax.tick_params(direction="in", length=2.5, width=0.6)
            if i == 0 and j == 0:
                ax.legend(frameon=False, loc="lower right")
    fig.tight_layout(w_pad=0.9, h_pad=0.75)
    fig_out = FIG_DIR / "linear_direct_formula_by_R_selected.pdf"
    fig.savefig(fig_out, bbox_inches="tight")
    fig.savefig(fig_out.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)

    summary = (
        df.groupby(["method", "R"], as_index=False)
        .agg(
            point_count=("y", "count"),
            mean_coverage=("coverage", "mean"),
            mean_formula=("predicted_coverage", "mean"),
            mean_abs_error=("residual", lambda s: np.mean(np.abs(s))),
            rmse=("residual", lambda s: np.sqrt(np.mean(np.square(s)))),
        )
        .sort_values(["method", "R"])
    )
    order = {m: i for i, m in enumerate(METHODS)}
    summary["method_order"] = summary["method"].map(order)
    summary = summary.sort_values(["method_order", "R"]).drop(columns=["method_order"])
    csv_out = OUT_DIR / "linear_direct_formula_by_R_summary.csv"
    summary.to_csv(csv_out, index=False)

    lines = [
        r"\begin{table}[!htbp]",
        r"\centering",
        r"\caption{Euler梁算例中不同 \(R\) 下直接计算覆盖率与中心化实验覆盖率}",
        r"\label{tab:linear-direct-formula}",
        r"\begingroup",
        r"\small",
        r"\setlength{\tabcolsep}{4pt}",
        r"\begin{tabular}{cccccc}",
        r"\toprule",
        r"方法 & \(R\) & 点数 & 平均实验值 & 平均计算值 & RMSE \\",
        r"\midrule",
    ]
    last_method = None
    for _, row in summary.iterrows():
        method = METHOD_LABEL[row["method"]]
        method_cell = method if method != last_method else ""
        if last_method is not None and method != last_method:
            lines.append(r"\addlinespace")
        lines.append(
            f"{method_cell} & {int(row['R'])} & {int(row['point_count'])} & "
            f"{row['mean_coverage']:.5f} & {row['mean_formula']:.5f} & {row['rmse']:.5f} \\\\"
        )
        last_method = method
    lines.extend(
        [
            r"\bottomrule",
            r"\end{tabular}",
            r"\endgroup",
            r"\end{table}",
            "",
        ]
    )
    tex_out = OUT_DIR / "linear_direct_formula_by_R_table.tex"
    tex_out.write_text("\n".join(lines), encoding="utf-8")
    return fig_out, tex_out


def plot_n_stability_and_write_table() -> tuple[Path, Path]:
    raw = np.load(DATA_ROOT / "curve_pool_raw.npz")
    response_all = raw["linear_response"]
    weights_all = raw["weights"]
    y_grid = np.linspace(-3.0, 3.0, 81)
    h = 0.08
    truth = norm.pdf(y_grid, loc=0.0, scale=np.sqrt(1.0 + h**2))
    n_values = [192, 384, 768]
    rows = []
    means = {}
    sds = {}
    for n in n_values:
        weights_n = weights_all[:, :n].copy()
        weights_n = weights_n / weights_n.sum(axis=1, keepdims=True)
        q = gaussian_kernel_ordinates(response_all[:, :n], weights_n, y_grid, h)
        mean_q = q.mean(axis=0)
        sd_q = q.std(axis=0, ddof=1)
        means[n] = mean_q
        sds[n] = sd_q
        rows.append(
            {
                "n": n,
                "mean_abs_bias": float(np.mean(np.abs(mean_q - truth))),
                "rmse": float(np.sqrt(np.mean((mean_q - truth) ** 2))),
                "max_abs_bias": float(np.max(np.abs(mean_q - truth))),
                "mean_curve_sd": float(np.mean(sd_q)),
            }
        )
    summary = pd.DataFrame(rows)
    summary.to_csv(OUT_DIR / "linear_n_stability_summary.csv", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(7.1, 2.75))
    ax = axes[0]
    ax.plot(y_grid, truth, color="#222222", linewidth=1.2, label="reference density")
    palette = ["#6b8fb8", "#cf7f4f", "#2d7f46"]
    for color, n in zip(palette, n_values):
        ax.plot(y_grid, means[n], color=color, linewidth=1.0, label=f"n={n}")
    ax.set_xlabel(r"$y$")
    ax.set_ylabel("Density")
    ax.tick_params(direction="in", length=2.5, width=0.6)
    ax.legend(frameon=False, loc="upper right")

    ax = axes[1]
    x = np.arange(len(n_values))
    width = 0.26
    ax.bar(x - width, summary["mean_abs_bias"], width, color="#6b8fb8", label="mean abs. bias")
    ax.bar(x, summary["rmse"], width, color="#cf7f4f", label="RMSE")
    ax.bar(x + width, summary["mean_curve_sd"], width, color="#2d7f46", label="mean curve SD")
    ax.set_xticks(x, [str(n) for n in n_values])
    ax.set_xlabel("Internal sample size n")
    ax.set_ylabel("Density-scale metric")
    ax.tick_params(direction="in", length=2.5, width=0.6)
    ax.legend(frameon=False, loc="upper right")
    fig.tight_layout(w_pad=1.0)
    fig_out = FIG_DIR / "linear_internal_sample_size_stability.pdf"
    fig.savefig(fig_out, bbox_inches="tight")
    fig.savefig(fig_out.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)

    lines = [
        r"\begin{table}[!htbp]",
        r"\centering",
        r"\caption{Euler梁算例中不同内部样本数下曲线池均值的稳定性}",
        r"\label{tab:linear-n-stability}",
        r"\begingroup",
        r"\small",
        r"\setlength{\tabcolsep}{6pt}",
        r"\begin{tabular}{ccccc}",
        r"\toprule",
        r"\(n\) & 平均绝对偏差 & RMSE & 最大绝对偏差 & 平均曲线标准差 \\",
        r"\midrule",
    ]
    for _, row in summary.iterrows():
        lines.append(
            f"{int(row['n'])} & {row['mean_abs_bias']:.5f} & {row['rmse']:.5f} & "
            f"{row['max_abs_bias']:.5f} & {row['mean_curve_sd']:.5f} \\\\"
        )
    lines.extend([r"\bottomrule", r"\end{tabular}", r"\endgroup", r"\end{table}", ""])
    tex_out = OUT_DIR / "linear_n_stability_table.tex"
    tex_out.write_text("\n".join(lines), encoding="utf-8")
    return fig_out, tex_out


def main() -> None:
    ensure_dirs()
    backup_existing(
        [
            "linear_dpim_pointwise_ci_panels.pdf",
            "linear_dpim_pointwise_ci_panels.png",
            "linear_direct_formula_by_R_selected.pdf",
            "linear_direct_formula_by_R_selected.png",
        ]
    )
    set_style()
    outputs = [plot_pointwise_interval_panels()]
    fig_out, table_out = plot_formula_by_r_and_write_table()
    outputs.extend([fig_out, table_out])
    n_fig, n_table = plot_n_stability_and_write_table()
    outputs.extend([n_fig, n_table])
    manifest = {
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "source_data": [
            str(DATA_ROOT / "curve_pool_raw.npz"),
            str(DATA_ROOT / "multiy_formula_predictions.csv"),
        ],
        "outputs": [str(p) for p in outputs],
        "note": "All coverage summaries use centered multiy_formula_predictions.csv for the linear Euler beam direct formula comparison.",
    }
    (OUT_DIR / "manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    for p in outputs:
        print(p)


if __name__ == "__main__":
    main()