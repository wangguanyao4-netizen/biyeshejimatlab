from __future__ import annotations

import hashlib
import json
import math
import shutil
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import beta, norm, t


PAPER_ROOT = Path(r"C:\Users\Wangg\Desktop\期刊论文")
PROJECT_ROOT = Path(r"C:\Users\Wangg\Desktop\biyeshejimatlab")
DATA_ROOT = Path(
    r"C:\Users\Wangg\Desktop\biyeshejimatlab_portable_e1e8\results"
    r"\linear_nonlinear_multiy_centered_paper_20260626_202142"
)
FIG_DIR = PAPER_ROOT / "DPIM_CI_full_integrated_figures" / "formal_current"
TEX_PATH = PAPER_ROOT / "DPIM_CI_full_integrated_weighted_RBn_natural.tex"
PDF_PATH = PAPER_ROOT / "DPIM_CI_full_integrated_weighted_RBn_natural.pdf"

METHODS = ["t distribution", "percentile bootstrap", "bootstrap-t"]
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
MODEL_LABEL = {
    "linear": "Euler beam",
    "nonlinear": "nonlinear response",
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

N_VALUES = [192, 384, 768]
R_VALUE = 128
M_REPLICATES = 400
B_BOOTSTRAP = 399
ALPHA = 0.05
SEED = 2026063001


def timestamp() -> str:
    return datetime.now().strftime("%Y%m%d_%H%M%S")


STAMP = timestamp()
OUT_DIR = PAPER_ROOT / "build" / f"n_coverage_revision_{STAMP}"
BACKUP_DIR = PAPER_ROOT / "build" / "codex_backups" / f"n_coverage_revision_{STAMP}"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest().upper()


def ensure_dirs() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    BACKUP_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)


def backup_existing() -> None:
    names = [
        TEX_PATH,
        PDF_PATH,
        FIG_DIR / "linear_nonlinear_multiy_region_rmse.pdf",
        FIG_DIR / "linear_nonlinear_multiy_region_rmse.png",
        FIG_DIR / "linear_internal_sample_size_stability.pdf",
        FIG_DIR / "linear_internal_sample_size_stability.png",
        FIG_DIR / "linear_n_centered_coverage_by_n.pdf",
        FIG_DIR / "linear_n_centered_coverage_by_n.png",
        FIG_DIR / "nonlinear_n_centered_coverage_by_n.pdf",
        FIG_DIR / "nonlinear_n_centered_coverage_by_n.png",
    ]
    for src in names:
        if src.exists():
            shutil.copy2(src, BACKUP_DIR / src.name)


def configure_plots() -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "STIXGeneral", "DejaVu Serif"],
            "mathtext.fontset": "stix",
            "axes.edgecolor": "#20252B",
            "axes.linewidth": 0.75,
            "axes.grid": False,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": False,
            "ytick.right": False,
            "legend.frameon": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "axes.unicode_minus": False,
        }
    )


def gaussian_kernel_matrix(y: np.ndarray, response: np.ndarray, h: float) -> np.ndarray:
    scaled = (y[:, None] - response[None, :]) / h
    return np.exp(-0.5 * scaled**2) / (math.sqrt(2.0 * math.pi) * h)


def finite_b_ranks(alpha: float, b: int) -> tuple[int, int]:
    k_minus = math.floor((alpha / 2.0) * (b + 1))
    k_plus = math.ceil((1.0 - alpha / 2.0) * (b + 1))
    return max(k_minus, 1), min(k_plus, b)


def exact_interval(hits: np.ndarray, trials: int) -> tuple[np.ndarray, np.ndarray]:
    lower = np.where(hits == 0, 0.0, beta.ppf(0.025, hits, trials - hits + 1))
    upper = np.where(hits == trials, 1.0, beta.ppf(0.975, hits + 1, trials - hits))
    return lower, upper


def weighted_pool_values(raw: np.lib.npyio.NpzFile, model: str, n: int, y_grid: np.ndarray, h: float) -> np.ndarray:
    weights = raw["weights"][:, :n].copy()
    weights_sum = weights.sum(axis=1, keepdims=True)
    weights = weights / np.maximum(weights_sum, np.finfo(float).tiny)
    response = raw[f"{model}_response"][:, :n]
    out = np.empty((response.shape[0], y_grid.size), dtype=np.float64)
    for i in range(response.shape[0]):
        out[i, :] = weights[i] @ gaussian_kernel_matrix(y_grid, response[i], h).T
    return out


def coverage_for_pool(pool_values: np.ndarray, truth: np.ndarray, model: str, n: int, h: float) -> list[dict]:
    rng = np.random.default_rng(SEED + (0 if model == "linear" else 5_000_000) + n)
    y_count = pool_values.shape[1]
    pool_mean = pool_values.mean(axis=0)
    shift = truth - pool_mean
    centered = pool_values + shift[None, :]

    idx = rng.integers(0, centered.shape[0], size=(M_REPLICATES, R_VALUE))
    sample = centered[idx, :]
    mu = sample.mean(axis=1)
    sd = sample.std(axis=1, ddof=1)
    se = sd / math.sqrt(R_VALUE)

    lower = np.empty((M_REPLICATES, y_count, len(METHODS)), dtype=np.float64)
    upper = np.empty_like(lower)
    tcrit = t.ppf(1.0 - ALPHA / 2.0, R_VALUE - 1)
    lower[:, :, 0] = mu - tcrit * se
    upper[:, :, 0] = mu + tcrit * se

    k_minus, k_plus = finite_b_ranks(ALPHA, B_BOOTSTRAP)
    bootstrap_t_diagnostic = np.zeros((M_REPLICATES, y_count), dtype=bool)
    batch_size = 20
    for start in range(0, M_REPLICATES, batch_size):
        stop = min(M_REPLICATES, start + batch_size)
        xb = sample[start:stop, :, :]
        mb = stop - start
        boot_idx = rng.integers(0, R_VALUE, size=(mb, B_BOOTSTRAP, R_VALUE))
        boot = np.take_along_axis(xb[:, None, :, :], boot_idx[:, :, :, None], axis=2)
        boot_mean = boot.mean(axis=2)
        sorted_mean = np.sort(boot_mean, axis=1)
        lower[start:stop, :, 1] = sorted_mean[:, k_minus - 1, :]
        upper[start:stop, :, 1] = sorted_mean[:, k_plus - 1, :]

        boot_sd = boot.std(axis=2, ddof=1)
        with np.errstate(divide="ignore", invalid="ignore"):
            pivot = math.sqrt(R_VALUE) * (boot_mean - mu[start:stop, None, :]) / boot_sd
        zero_sd = boot_sd <= 0
        pivot[zero_sd & (boot_mean == mu[start:stop, None, :])] = 0.0
        pivot[zero_sd & (boot_mean > mu[start:stop, None, :])] = np.inf
        pivot[zero_sd & (boot_mean < mu[start:stop, None, :])] = -np.inf
        pivot = np.sort(pivot, axis=1)
        lower[start:stop, :, 2] = mu[start:stop, :] - pivot[:, k_plus - 1, :] * se[start:stop, :]
        upper[start:stop, :, 2] = mu[start:stop, :] - pivot[:, k_minus - 1, :] * se[start:stop, :]
        length_bt = upper[start:stop, :, 2] - lower[start:stop, :, 2]
        length_p = upper[start:stop, :, 1] - lower[start:stop, :, 1]
        bootstrap_t_diagnostic[start:stop, :] = (~np.isfinite(length_bt)) | (
            np.isfinite(length_bt) & (length_bt / np.maximum(length_p, np.finfo(float).tiny) > 5.0)
        )

    rows: list[dict] = []
    for iy in range(y_count):
        for im, method in enumerate(METHODS):
            finite = np.isfinite(lower[:, iy, im]) & np.isfinite(upper[:, iy, im])
            hit = finite & (lower[:, iy, im] <= truth[iy]) & (truth[iy] <= upper[:, iy, im])
            coverage = float(hit.mean())
            rows.append(
                {
                    "model": model,
                    "n": n,
                    "y": float(coverage_for_pool.y_grid[iy]),
                    "h": h,
                    "method": method,
                    "R": R_VALUE,
                    "M": M_REPLICATES,
                    "B": B_BOOTSTRAP,
                    "truth": float(truth[iy]),
                    "pool_mean": float(pool_mean[iy]),
                    "centering_shift": float(shift[iy]),
                    "coverage": coverage,
                    "coverage_mcse": math.sqrt(max(coverage * (1.0 - coverage), 0.0) / M_REPLICATES),
                    "hit_count": int(round(coverage * M_REPLICATES)),
                    "mean_interval_length": float(np.nanmean(upper[finite, iy, im] - lower[finite, iy, im]))
                    if finite.any()
                    else np.nan,
                    "interval_inf_rate": float((~finite).mean()),
                    "bootstrap_t_diagnostic_rate": float(bootstrap_t_diagnostic[:, iy].mean())
                    if method == "bootstrap-t"
                    else 0.0,
                }
            )
    return rows


coverage_for_pool.y_grid = np.array([], dtype=float)  # type: ignore[attr-defined]


def build_n_coverage() -> pd.DataFrame:
    raw = np.load(DATA_ROOT / "curve_pool_raw.npz")
    moments = pd.read_csv(DATA_ROOT / "multiy_kernel_moments.csv")
    rows: list[dict] = []
    for model in ["linear", "nonlinear"]:
        model_moments = moments[moments["model"] == model].sort_values("y").reset_index(drop=True)
        y_grid = model_moments["y"].to_numpy(dtype=float)
        truth = model_moments["truth"].to_numpy(dtype=float)
        h = float(model_moments["h"].iloc[0])
        coverage_for_pool.y_grid = y_grid  # type: ignore[attr-defined]
        for n in N_VALUES:
            print(f"n coverage: model={model}, n={n}", flush=True)
            pool_values = weighted_pool_values(raw, model, n, y_grid, h)
            rows.extend(coverage_for_pool(pool_values, truth, model, n, h))
    out = pd.DataFrame(rows)
    hits = out["hit_count"].to_numpy()
    lower, upper = exact_interval(hits, M_REPLICATES)
    out["coverage_exact95_lower"] = lower
    out["coverage_exact95_upper"] = upper
    out.to_csv(OUT_DIR / "n_centered_coverage_results.csv", index=False)
    summary = (
        out.groupby(["model", "n", "method"], as_index=False)
        .agg(
            point_count=("y", "size"),
            mean_coverage=("coverage", "mean"),
            min_coverage=("coverage", "min"),
            max_coverage=("coverage", "max"),
            mean_mcse=("coverage_mcse", "mean"),
            mean_interval_length=("mean_interval_length", "mean"),
            max_inf_rate=("interval_inf_rate", "max"),
        )
        .sort_values(["model", "method", "n"])
    )
    summary.to_csv(OUT_DIR / "n_centered_coverage_summary.csv", index=False)
    return out


def plot_n_coverage(df: pd.DataFrame, model: str, output_name: str) -> None:
    frame = df[df["model"] == model].copy()
    fig, axes = plt.subplots(len(N_VALUES), len(METHODS), figsize=(8.9, 5.65), sharex=True, sharey=True)
    for i, n in enumerate(N_VALUES):
        for j, method in enumerate(METHODS):
            ax = axes[i, j]
            sub = frame[(frame["n"] == n) & (frame["method"] == method)].sort_values("y")
            color = METHOD_COLOR[method]
            ax.plot(sub["y"], sub["coverage"], color=color, linewidth=1.05)
            selected = np.unique(np.round(np.linspace(0, len(sub) - 1, 17)).astype(int))
            ax.scatter(
                sub["y"].to_numpy()[selected],
                sub["coverage"].to_numpy()[selected],
                s=13,
                facecolor="white",
                edgecolor=color,
                linewidth=0.75,
                zorder=3,
            )
            ax.axhline(0.95, color="#5C626B", linestyle=":", linewidth=0.8)
            if i == 0:
                ax.set_title(METHOD_LABEL[method], fontsize=10.0, pad=3)
            if j == 0:
                ax.set_ylabel(f"$n={n}$\ncoverage", fontsize=9.1)
            if i == len(N_VALUES) - 1:
                ax.set_xlabel(r"response coordinate $y$", fontsize=9.1)
            ax.tick_params(labelsize=8.4)
            ax.grid(True, axis="y", color="#D7DCE2", linewidth=0.35, alpha=0.70)
    axes[0, 0].set_ylim(0.88, 0.985)
    fig.tight_layout(w_pad=0.85, h_pad=0.70)
    for ext in ("pdf", "png"):
        fig.savefig(FIG_DIR / f"{output_name}.{ext}", dpi=320, bbox_inches="tight", pad_inches=0.035)
        fig.savefig(OUT_DIR / f"{output_name}.{ext}", dpi=320, bbox_inches="tight", pad_inches=0.035)
    plt.close(fig)


def plot_region_rmse() -> None:
    region = pd.read_csv(DATA_ROOT / "multiy_region_summary.csv")
    summary = (
        region.groupby(["model", "method", "density_region"], as_index=False)
        .agg(mean_rmse=("rmse", "mean"))
    )
    fig, axes = plt.subplots(2, 3, figsize=(8.7, 4.75), sharey=True)
    x = np.arange(len(REGION_ORDER))
    for i, model in enumerate(["linear", "nonlinear"]):
        for j, method in enumerate(METHODS):
            ax = axes[i, j]
            sub = summary[(summary["model"] == model) & (summary["method"] == method)]
            vals = [
                float(sub[sub["density_region"] == region_name]["mean_rmse"].iloc[0])
                for region_name in REGION_ORDER
            ]
            color = METHOD_COLOR[method]
            ax.plot(x, vals, color=color, linewidth=1.05, marker="o", markersize=4.2, markerfacecolor="white")
            if i == 0:
                ax.set_title(METHOD_LABEL[method], fontsize=10.0, pad=3)
            if j == 0:
                ax.set_ylabel(f"{MODEL_LABEL[model]}\nmean RMSE", fontsize=9.2)
            ax.set_xticks(x)
            ax.set_xticklabels([REGION_LABEL[r] for r in REGION_ORDER], fontsize=8.8)
            ax.tick_params(labelsize=8.4)
            ax.grid(True, axis="y", color="#D7DCE2", linewidth=0.35, alpha=0.70)
    fig.tight_layout(w_pad=0.85, h_pad=0.80)
    for ext in ("pdf", "png"):
        fig.savefig(
            FIG_DIR / f"linear_nonlinear_multiy_region_rmse.{ext}",
            dpi=320,
            bbox_inches="tight",
            pad_inches=0.035,
        )
        fig.savefig(
            OUT_DIR / f"linear_nonlinear_multiy_region_rmse.{ext}",
            dpi=320,
            bbox_inches="tight",
            pad_inches=0.035,
        )
    plt.close(fig)


def plot_nonlinear_direct_formula_no_shadow() -> None:
    pred = pd.read_csv(DATA_ROOT / "multiy_formula_predictions.csv")
    frame = pred[
        (pred["model"] == "nonlinear")
        & (pred["R"].isin([16, 128]))
        & (pred["B"] == B_BOOTSTRAP)
    ].copy()
    fig, axes = plt.subplots(2, 3, figsize=(8.7, 4.75), sharex=True, sharey=True)
    for i, r_value in enumerate([16, 128]):
        for j, method in enumerate(METHODS):
            ax = axes[i, j]
            sub = frame[(frame["R"] == r_value) & (frame["method"] == method)].sort_values("y")
            color = METHOD_COLOR[method]
            x = sub["y"].to_numpy()
            selected = np.unique(np.round(np.linspace(0, len(x) - 1, 17)).astype(int))
            ax.plot(x, sub["predicted_coverage"].to_numpy(), color=color, linewidth=1.05)
            ax.scatter(
                x[selected],
                sub["coverage"].to_numpy()[selected],
                s=13,
                facecolor="white",
                edgecolor=color,
                linewidth=0.75,
                zorder=3,
            )
            ax.axhline(0.95, color="#5C626B", linestyle=":", linewidth=0.8)
            if i == 0:
                ax.set_title(METHOD_LABEL[method], fontsize=10.0, pad=3)
            if j == 0:
                ax.set_ylabel(f"$R={r_value}$\ncoverage", fontsize=9.2)
            if i == 1:
                ax.set_xlabel(r"response coordinate $y$", fontsize=9.1)
            ax.set_ylim(0.88, 0.985)
            ax.tick_params(labelsize=8.4)
            ax.grid(True, axis="y", color="#D7DCE2", linewidth=0.35, alpha=0.70)
    fig.tight_layout(w_pad=0.85, h_pad=0.80)
    for ext in ("pdf", "png"):
        fig.savefig(
            FIG_DIR / f"nonlinear_direct_formula_R16_R128.{ext}",
            dpi=320,
            bbox_inches="tight",
            pad_inches=0.035,
        )
        fig.savefig(
            OUT_DIR / f"nonlinear_direct_formula_R16_R128.{ext}",
            dpi=320,
            bbox_inches="tight",
            pad_inches=0.035,
        )
    plt.close(fig)


def write_manifest(outputs: list[Path]) -> None:
    files = [TEX_PATH, PDF_PATH] + outputs
    manifest = {
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "source_data": [
            str(DATA_ROOT / "curve_pool_raw.npz"),
            str(DATA_ROOT / "multiy_kernel_moments.csv"),
            str(DATA_ROOT / "multiy_region_summary.csv"),
        ],
        "settings": {
            "n_values": N_VALUES,
            "R": R_VALUE,
            "M": M_REPLICATES,
            "B": B_BOOTSTRAP,
            "centered": True,
            "probability_weighting": "Voronoi probability weights from existing RQMC curve pool",
        },
        "backup_dir": str(BACKUP_DIR),
        "outputs": [str(p) for p in outputs],
        "hashes": {str(p): sha256(p) for p in files if p.exists()},
    }
    (OUT_DIR / "manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")


def main() -> None:
    ensure_dirs()
    backup_existing()
    configure_plots()
    coverage = build_n_coverage()
    plot_n_coverage(coverage, "linear", "linear_n_centered_coverage_by_n")
    plot_n_coverage(coverage, "nonlinear", "nonlinear_n_centered_coverage_by_n")
    plot_region_rmse()
    plot_nonlinear_direct_formula_no_shadow()
    outputs = [
        OUT_DIR / "n_centered_coverage_results.csv",
        OUT_DIR / "n_centered_coverage_summary.csv",
        FIG_DIR / "linear_n_centered_coverage_by_n.pdf",
        FIG_DIR / "nonlinear_n_centered_coverage_by_n.pdf",
        FIG_DIR / "linear_nonlinear_multiy_region_rmse.pdf",
        FIG_DIR / "nonlinear_direct_formula_R16_R128.pdf",
    ]
    write_manifest(outputs)
    print(f"output_dir={OUT_DIR}")
    print(f"backup_dir={BACKUP_DIR}")


if __name__ == "__main__":
    main()