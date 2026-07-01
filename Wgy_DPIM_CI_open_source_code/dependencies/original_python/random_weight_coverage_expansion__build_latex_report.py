from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


METHOD_LABEL = {
    "t distribution": r"$t$ distribution",
    "percentile bootstrap": "percentile bootstrap",
    "bootstrap-t": r"bootstrap-\textnormal{t}",
}
SOURCE_LABEL = {
    "complete outer cumulants": "complete outer cumulants",
    "fixed-weight factorization": "fixed-weight factorization",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("plate_root", type=Path)
    parser.add_argument("linear_nonlinear_root", type=Path)
    return parser.parse_args()


def fmt(value: float, digits: int = 5) -> str:
    return f"{float(value):.{digits}f}"


def coefficient_table(plate_root: Path) -> str:
    data = pd.read_csv(plate_root / "direct_coefficients.csv")
    lines = []
    for B in sorted(data["B"].unique()):
        for method in METHOD_LABEL:
            block = data[(data["B"] == B) & (data["method"] == method)].set_index(
                "coefficient"
            )
            lines.append(
                f"{B} & {METHOD_LABEL[method]} & "
                f"{fmt(block.loc['A0', 'value'], 6)} & "
                f"{fmt(block.loc['A4', 'value'], 6)} & "
                f"{fmt(block.loc['A33', 'value'], 6)} \\\\"
            )
    return "\n".join(lines)


def plate_overall_table(plate_root: Path) -> str:
    data = pd.read_csv(plate_root / "formula_validation_summary.csv")
    overall = (
        data.groupby(["method", "moment_source"], as_index=False)
        .agg(
            MAE=("mean_abs_error", "mean"),
            RMSE=("root_mean_square_error", "mean"),
            inside=("share_prediction_inside_exact95", "mean"),
            max_diag=("max_bootstrap_t_diagnostic_rate", "max"),
        )
    )
    lines = []
    for method in METHOD_LABEL:
        for source in SOURCE_LABEL:
            row = overall[
                (overall["method"] == method)
                & (overall["moment_source"] == source)
            ].iloc[0]
            lines.append(
                f"{METHOD_LABEL[method]} & {SOURCE_LABEL[source]} & "
                f"{fmt(row.MAE)} & {fmt(row.RMSE)} & "
                f"{fmt(row.inside, 4)} & {fmt(row.max_diag, 4)} \\\\"
            )
    return "\n".join(lines)


def plate_core_table(plate_root: Path) -> str:
    data = pd.read_csv(plate_root / "formula_validation_by_density_region.csv")
    data = data[
        (data["moment_source"] == "complete outer cumulants")
        & (data["relative_density_threshold"] == 0.05)
    ]
    overall = (
        data.groupby("method", as_index=False)
        .agg(
            MAE=("mean_abs_error", "mean"),
            RMSE=("root_mean_square_error", "mean"),
            inside=("share_prediction_inside_exact95", "mean"),
            max_error=("max_abs_error", "max"),
            max_l3=("max_abs_lambda3_over_sqrt_R", "max"),
            max_l4=("max_abs_lambda4_over_R", "max"),
        )
    ).set_index("method")
    lines = []
    for method in METHOD_LABEL:
        row = overall.loc[method]
        lines.append(
            f"{METHOD_LABEL[method]} & {fmt(row.MAE)} & {fmt(row.RMSE)} & "
            f"{fmt(row.max_error)} & {fmt(row.inside, 4)} & "
            f"{fmt(row.max_l3, 3)} & {fmt(row.max_l4, 3)} \\\\"
        )
    return "\n".join(lines)


def linear_nonlinear_table(root: Path) -> str:
    data = pd.read_csv(root / "formula_validation_summary.csv")
    pivot = data.pivot_table(
        index=["model", "method"],
        columns="moment_source",
        values="root_mean_square_error",
    )
    lines = []
    for model in ("linear", "nonlinear"):
        for method in METHOD_LABEL:
            row = pivot.loc[(model, method)]
            lines.append(
                f"{model} & {METHOD_LABEL[method]} & "
                f"{fmt(row['complete outer cumulants'])} & "
                f"{fmt(row['fixed-weight factorization'])} \\\\"
            )
    return "\n".join(lines)


def write_report(plate_root: Path, linear_root: Path) -> Path:
    tex = r"""\documentclass[UTF8,11pt]{ctexart}
\usepackage[a4paper,margin=2.35cm]{geometry}
\usepackage{amsmath,amssymb,mathtools,bm}
\usepackage{booktabs,longtable,array}
\usepackage{graphicx}
\usepackage{xcolor}
\usepackage{hyperref}
\usepackage{microtype}
\usepackage{enumitem}
\hypersetup{colorlinks=true,linkcolor=black,citecolor=blue,urlcolor=blue}
\setlength{\parindent}{2em}
\setlength{\parskip}{0.25em}
\newcommand{\E}{\mathbb E}
\newcommand{\Pp}{\mathbb P}
\newcommand{\cum}{\operatorname{cum}}
\newcommand{\normalpdf}{\phi}
\newcommand{\normalcdf}{\Phi}
\newcommand{\Lop}{\mathcal L}
\title{随机 Voronoi 权重下 DPIM 置信区间覆盖展开：\\外层累积量公式、直接计算与薄板验证}
\author{独立理论与数值审计稿}
\date{\today}

\begin{document}
\maketitle

\begin{abstract}
本文建立一个不依赖曲线内部独立性的覆盖率计算框架。对固定的内层点数
$n$、带宽 $h$ 与响应位置 $y$，一次完整的 randomized QMC、Voronoi
概率质量计算和 DPIM 核估计被视为一个外层随机变量 $X_h(y)$。只要不同
随机化曲线相互独立，三种双侧区间的二阶覆盖修正即可直接用
$X_h(y)$ 的标准化三、四阶累积量计算。本文给出 $t$ distribution、
percentile bootstrap 与 bootstrap-\textnormal{t} 的有限 $B$ 系数、
一维积分公式和数值实现；固定权重因子化只作为特殊情形。线性、非线性和
Kirchhoff 薄板结果表明：固定权重公式在线性/非线性标量算例中仍是良好近似，
但在薄板响应尾部会严重高估外层偏度和峰度；改用独立 pilot 曲线估计的完整
外层累积量后，覆盖预测误差明显下降。本文的定理是固定
$(n,h,y,B)$、$R\to\infty$ 的点态条件定理，不声称已经证明
$n\to\infty$、$h\to0$、$B\to\infty$ 或连续响应区间上的联合一致余项。
\end{abstract}

\section{问题与结论边界}

令第 $r$ 次独立随机化产生代表点、Voronoi 概率权重和响应值
\[
\{(\Theta_{r,i},W_{r,i},g(\Theta_{r,i}))\}_{i=1}^{n},
\]
并在固定响应位置 $y$ 定义整条算法的一次输出
\begin{equation}\label{eq:outer}
X_r(y)=\sum_{i=1}^{n}W_{r,i}
K_h\!\left(y-g(\Theta_{r,i})\right).
\end{equation}
本文的统计样本是 $X_1(y),\ldots,X_R(y)$，而不是式
\eqref{eq:outer} 内部的 $n$ 个核项。因此，单条 scrambled Sobol 曲线内部
可以相关，$W_{r,i}$ 也可以与 $g(\Theta_{r,i})$ 相关；这些结构均由
$X_r(y)$ 的真实分布吸收。本文只要求不同 $r$ 使用独立随机化。

记
\[
\mu=\E X_r,\qquad \sigma^2=\cum_2(X_r)>0,
\qquad
\lambda_3=\frac{\cum_3(X_r)}{\sigma^3},\qquad
\lambda_4=\frac{\cum_4(X_r)}{\sigma^4}.
\]
这里 $\lambda_4$ 是超额峰度。以下覆盖展开首先针对 $\mu$；若需要覆盖
MC 参考密度，则必须额外处理有限 $n$、有限 Voronoi 质量估计和核平滑偏差。
本次数值验证使用点态平移使验证池均值与 MC 参考曲线一致，目的只是隔离区间
覆盖机制，不能被解释为未知真值下可实施的偏差修正算法。

\section{有限 \texorpdfstring{$B$}{B} 的精确基准}

设 bootstrap 端点由 $B$ 个条件独立重采样统计量的第 $k_-$ 与第 $k_+$
个次序统计量给出。对 percentile bootstrap，令
\[
U_{p,R}=G_R^*(\mu),
\]
其中 $G_R^*$ 是给定原始 $R$ 条曲线时 bootstrap 均值的条件分布函数。
bootstrap-\textnormal{t} 同理使用其 studentized root 的条件分布函数，
记为 $U_{bt,R}$。忽略概率为零的并列事件后，给定 $U_{m,R}=u$，
区间覆盖等价于
\[
k_-\le N_B(u)\le k_+-1,\qquad N_B(u)\sim\operatorname{Bin}(B,u).
\]
故定义
\begin{equation}\label{eq:PiB}
\Pi_B(u;k_-,k_+)
=\Pp\{k_-\le\operatorname{Bin}(B,u)\le k_+-1\},
\end{equation}
即可得到精确恒等式
\begin{equation}\label{eq:finiteBexact}
C_m(R,B;y)=\E\{\Pi_B(U_{m,R};k_-,k_+)\}.
\end{equation}
若校准变量严格服从 $\operatorname{Unif}(0,1)$，Beta--binomial 积分给出
\begin{equation}\label{eq:C0B}
C_{0,B}=\int_0^1\Pi_B(u;k_-,k_+)\,du
=\frac{k_+-k_-}{B+1}.
\end{equation}
因此 percentile bootstrap 与 bootstrap-\textnormal{t} 的基准不是先验
固定的 0.95，而是由实际端点秩确定的 $C_{0,B}$。

\section{外层 Edgeworth 输入与三种直接公式}

\subsection{工作假设}

以下结果采用 Hall 的标准化均值与 studentized 均值 Edgeworth 规范
\cite{Hall1988,Hall1992}。对每个固定的 $(n,h,y)$，假设外层变量具有高于
六阶的矩、满足适当的 Cram\'er 非格点条件，且方差远离零；同时假设原统计量
和非参数 bootstrap 条件分布在覆盖所需紧区间上具有一致到
$o(R^{-1})$ 的二阶展开。该组条件是充分而非最弱条件。RQMC 的中心极限定理
\cite{Loh2003,NakayamaTuffin2024} 本身只支持一阶正态性，不能自动推出
本文使用的二阶 Edgeworth 余项。

在上述规范下，标准化均值与 studentized 均值的分布函数分别写为
\begin{align}
\Pp(Z_R\le z)
&=\normalcdf(z)+R^{-1/2}p_1(z)\normalpdf(z)
+R^{-1}p_2(z)\normalpdf(z)+o(R^{-1}),\\
\Pp(T_R\le z)
&=\normalcdf(z)+R^{-1/2}q_1(z)\normalpdf(z)
+R^{-1}q_2(z)\normalpdf(z)+o(R^{-1}),
\end{align}
其中
\begin{align}
p_1(z)&=-\frac{\lambda_3}{6}(z^2-1),\\
p_2(z)&=-z\left\{\frac{\lambda_4}{24}(z^2-3)
+\frac{\lambda_3^2}{72}(z^4-10z^2+15)\right\},\\
q_1(z)&=\frac{\lambda_3}{6}(2z^2+1),\\
q_2(z)&=z\left\{\frac{\lambda_4}{12}(z^2-3)
-\frac{\lambda_3^2}{18}(z^4+2z^2-3)
-\frac{z^2+3}{4}\right\}.
\end{align}

\subsection{校准变量的二阶多项式}

令 $Z_{m,R}=\normalcdf^{-1}(U_{m,R})$。对原分布展开、bootstrap 条件展开、
反函数展开和随机样本累积量代入项逐项组合，可写成
\begin{equation}\label{eq:calibration-cdf}
\Pp(Z_{m,R}\le z)
=\normalcdf(z)+R^{-1/2}r_{m,1}(z)\normalpdf(z)
+R^{-1}r_{m,2}(z)\normalpdf(z)+o(R^{-1}).
\end{equation}
对 percentile bootstrap，
\begin{align}
r_{p,1}(z)&=-\frac{\lambda_3}{2}z^2,\\
r_{p,2}(z)
&=-\frac{z(z^2+3)}{4}
+\lambda_4\frac{z(7z^2-13)}{24}
-\lambda_3^2\frac{z(3z^4+6z^2-11)}{24}.
\label{eq:rp2}
\end{align}
对 bootstrap-\textnormal{t}，
\begin{equation}\label{eq:rbt2}
r_{bt,1}(z)=0,\qquad
r_{bt,2}(z)
=-\lambda_4\frac{z(2z^2+1)}{6}
+\lambda_3^2\frac{z(2z^2+1)}{4}.
\end{equation}
式 \eqref{eq:rp2} 中的分布无关项
$-z(z^2+3)/4$ 不能删去；它解释了正态外层变量下 percentile bootstrap
仍存在 $R^{-1}$ 覆盖误差。本文实现使用 SymPy 重新完成上述代数组合，并以
独立数值积分验证有限 $B$ 基准、接受函数对称性和双侧
$R^{-1/2}$ 项抵消。

\subsection{有限 \texorpdfstring{$B$}{B} 一维积分}

由 \eqref{eq:calibration-cdf} 求导，
\[
f_{Z_{m,R}}(z)
=\normalpdf(z)\left[
1+R^{-1/2}\{r_{m,1}'(z)-zr_{m,1}(z)\}
+R^{-1}\{r_{m,2}'(z)-zr_{m,2}(z)\}
\right]+o(R^{-1}).
\]
将其代入 \eqref{eq:finiteBexact}，并记
\[
\widetilde\Pi_B(z)
=\Pi_B(\normalcdf(z);k_-,k_+),
\]
得到
\begin{equation}\label{eq:main-bootstrap}
C_m(R,B;y)
=C_{0,B}+\frac{
A_m^{(0)}(B)+A_m^{(4)}(B)\lambda_4(y)
+A_m^{(33)}(B)\lambda_3^2(y)
}{R}+o(R^{-1}),
\end{equation}
其中每个系数均可直接计算：
\begin{equation}\label{eq:coefficient-integral}
A_m^{(j)}(B)
=\int_{-\infty}^{\infty}
\widetilde\Pi_B(z)
\left\{\frac{d}{dz}r_m^{(j)}(z)-zr_m^{(j)}(z)\right\}
\normalpdf(z)\,dz.
\end{equation}
因为 $\widetilde\Pi_B$ 为偶函数，而 percentile 的一阶密度修正为奇函数，
等尾双侧总覆盖中的 $R^{-1/2}$ 项严格抵消。该抵消不适用于单侧覆盖或左右尾
误差分别讨论。

对使用 Student 临界值 $t_{R-1,1-\alpha/2}$ 的 $t$ distribution 区间，
Cornish--Fisher 临界值修正抵消 $q_2$ 中的分布无关项，得到
\begin{equation}\label{eq:t-direct}
C_t(R;y)
=1-\alpha+\frac{2z\normalpdf(z)}{R}
\left\{
\frac{z^2-3}{12}\lambda_4(y)
-\frac{z^4+2z^2-3}{18}\lambda_3^2(y)
\right\}
+o(R^{-1}),
\end{equation}
其中 $z=z_{1-\alpha/2}$。当 $\lambda_3=\lambda_4=0$ 时，
式 \eqref{eq:t-direct} 退化为正态样本下 Student 区间的精确覆盖。

\begin{table}[htbp]
\centering
\caption{本次程序直接积分得到的覆盖系数。}
\begin{tabular}{clrrr}
\toprule
$B$ & method & $A^{(0)}$ & $A^{(4)}$ & $A^{(33)}$\\
\midrule
@@COEFFICIENT_TABLE@@
\bottomrule
\end{tabular}
\end{table}

\section{随机权重与交叉累积量}

令 $V_i=W_iK_h(y-g(\Theta_i))$。联合累积量的多线性给出
\begin{equation}\label{eq:joint-cumulant}
\cum_q(X)=\sum_{i_1,\ldots,i_q}
\cum(V_{i_1},\ldots,V_{i_q}),\qquad q=2,3,4.
\end{equation}
式 \eqref{eq:joint-cumulant} 是精确恒等式；条件累积量进一步分解可由
Brillinger 公式给出 \cite{Brillinger1969}。只有当 $W_i=w_i$ 为确定常数且
内部核变量独立同分布时，非对角联合累积量才消失，并有
\[
\lambda_3=\gamma_h\rho_3(w),\qquad
\lambda_4=\kappa_h\rho_4(w),
\]
\[
\rho_3(w)=\frac{\sum_iw_i^3}{(\sum_iw_i^2)^{3/2}},
\qquad
\rho_4(w)=\frac{\sum_iw_i^4}{(\sum_iw_i^2)^2}.
\]
对于每条曲线权重均随机的情形，平均覆盖的固定权重比较量应写成
\[
\lambda_{3,\mathrm{fac}}^2
=\gamma_h^2\E\{\rho_3(W)^2\},\qquad
\lambda_{4,\mathrm{fac}}
=\kappa_h\E\{\rho_4(W)\},
\]
而不是把随机权重简单替换为一组平均权重。定义
\[
\Delta_{33}=\lambda_3^2-\lambda_{3,\mathrm{fac}}^2,\qquad
\Delta_4=\lambda_4-\lambda_{4,\mathrm{fac}},
\]
则完整外层公式相对固定权重因子化公式的二阶修正恰为
\[
\frac{A_m^{(4)}\Delta_4+A_m^{(33)}\Delta_{33}}{R}.
\]
这一定义把所有权重--响应相关项和曲线内部 RQMC 相关项压缩为两个可估计的
外层量，但没有声称已经逐项解析计算式 \eqref{eq:joint-cumulant} 中的每个
非对角联合累积量。

\section{数值设计}

\subsection{线性与非线性标量算例}

线性和非线性正式结果沿用已锁定的概率加权 scrambled Sobol 曲线池：
$n=768$，每个模型有 1200 条曲线，$R\in
\{16,24,32,48,64,96,128,192\}$，$B=999$，每个配置重复 $M=1200$。
覆盖目标通过点态平移与解析或高精度 MC 核平滑真值对齐。本文重新读取原始
覆盖率和矩数据，所有系数均由式 \eqref{eq:coefficient-integral} 直接计算，
没有用覆盖率拟合系数。

\begin{table}[htbp]
\centering
\caption{线性与非线性算例的覆盖预测 RMSE。}
\begin{tabular}{llrr}
\toprule
model & method & complete outer & fixed-weight factorization\\
\midrule
@@LINEAR_TABLE@@
\bottomrule
\end{tabular}
\end{table}

两种累积量输入在标量算例中均给出约 $0.004$--$0.006$ 的 RMSE，48 个
模型--方法--$R$ 直接预测点全部落入实验覆盖率的逐点精确二项 95\% 区间。
完整外层输入并非在每一行都优于因子化输入，这与有限曲线池的三、四阶矩估计
噪声一致；因此标量结果支持“因子化在该组模型中近似可用”，不能支持
“随机权重交叉项普遍为零”。

\subsection{Kirchhoff 薄板}

薄板为 $1\,\mathrm m\times1\,\mathrm m$ 简支 Kirchhoff 方板，厚度
$5\,\mathrm{mm}$，均布荷载 $q=-100\,\mathrm{Pa}$，中心挠度以米为单位。
内层样本数 $n=200$，带宽
$h=1.10742439019\times10^{-6}\,\mathrm m$。200 条独立 RQMC 曲线只用于
估计 $\lambda_3,\lambda_4$ 及其 bootstrap 不确定性；另一组 600 条曲线只
用于覆盖试验。验证池按
\[
X_r^c(y)=X_r(y)+p_{\mathrm{MC}}(y)-\bar X_{600}(y)
\]
平移。该变换保持方差、偏度、峰度和区间长度不变。正式配置为
\[
R=\{16,32,64,128\},\quad B=\{399,999\},\quad M=1000,
\]
并在 MC 参考密度不低于峰值 1\% 的 123 个网格点上计算。GPU 仅用于
bootstrap 计数矩阵与曲线矩阵乘法；所有公式系数和最终统计量使用 double
precision。

\begin{figure}[htbp]
\centering
\includegraphics[width=\textwidth]{plate_coverage_profiles.pdf}
\caption{薄板逐点实验覆盖率与直接公式预测。灰带为 $M=1000$ 次覆盖计数的
Clopper--Pearson 95\% 区间；蓝线使用独立 pilot 的完整外层累积量，橙色虚线
使用固定权重因子化累积量。}
\end{figure}

\begin{figure}[htbp]
\centering
\includegraphics[width=0.92\textwidth]{outer_vs_factorized_cumulants.pdf}
\caption{完整外层累积量与固定权重因子化结果。差异在响应区两端显著增大，
直接量化了随机权重和曲线内部相关结构的综合影响。}
\end{figure}

\begin{table}[htbp]
\centering
\caption{薄板 1\% 有效区内全部 $R,B$ 配置的平均预测误差。}
\resizebox{\textwidth}{!}{%
\begin{tabular}{llrrrr}
\toprule
method & moment input & MAE & RMSE & inside exact 95\% & max diagnostic\\
\midrule
@@PLATE_OVERALL_TABLE@@
\bottomrule
\end{tabular}
}
\end{table}

完整外层累积量对 $t$ distribution 和 percentile bootstrap 的改善最明显：
平均 RMSE 分别由约 0.0288 降至 0.0078、由约 0.0470 降至 0.0108；
bootstrap-\textnormal{t} 由约 0.01025 降至 0.00783。所有正式配置的区间
端点均有限，预先记录的 bootstrap-\textnormal{t} 超长诊断率为零。

\begin{table}[htbp]
\centering
\caption{薄板核心区（参考密度不低于峰值 5\%）的完整外层公式表现。}
\resizebox{\textwidth}{!}{%
\begin{tabular}{lrrrrrr}
\toprule
method & MAE & RMSE & max error & inside 95\% &
$\max|\lambda_3|/\sqrt R$ & $\max|\lambda_4|/R$\\
\midrule
@@PLATE_CORE_TABLE@@
\bottomrule
\end{tabular}
}
\end{table}

在 5\% 核心区，三种方法的平均 RMSE 均约为 0.007--0.008，约
93\%--95\% 的公式预测落入实验覆盖率的逐点精确 95\% 区间。1\% 有效区
两端的最差情形出现在 $R=16$ 的 percentile bootstrap：独立 pilot 的
$|\lambda_3|/\sqrt R$ 可达约 0.52，$|\lambda_4|/R$ 可达约 0.49，最大
逐点残差约 0.076。此时二阶修正已不再相对基准足够小，不能期待点态
$o(R^{-1})$ 近似在尾部保持同样精度。

\section{证明状态、未解决问题与可投稿表述}

\begin{enumerate}[leftmargin=2.2em]
\item 式 \eqref{eq:finiteBexact}--\eqref{eq:C0B} 是有限样本精确恒等式。
\item 式 \eqref{eq:coefficient-integral} 的变量变换、奇偶抵消和一维积分
已由独立代码完成符号与数值自检；积分误差在本次 $B=399,999$ 下不超过
$1.5\times10^{-12}$。
\item 式 \eqref{eq:main-bootstrap} 和 \eqref{eq:t-direct} 是以 Hall 型
二阶 Edgeworth 与 bootstrap 条件展开为前提的条件定理。本文没有重新证明
一般随机 Voronoi/RQMC 算法对所有 $(n,h,y)$ 满足统一 Cram\'er 条件。
\item 完整外层累积量公式不要求曲线内部独立；但将 pilot 样本矩代入理论
累积量会引入额外误差，
\[
O_p(P^{-1/2}),
\]
尤其四阶累积量在尾部不稳定。
\item 本次数值证据支持如下表述：在固定 $n,h$ 的点态外层重复框架下，
三种区间的二阶覆盖预测可由真实外层三、四阶累积量直接计算；随机权重
因子化误差在薄板尾部不可忽略。它不支持如下更强表述：连续响应区域上的
一致余项已经证明，或所有随机权重交叉累积量都已解析闭式计算。
\end{enumerate}

\section{可复现文件}

正式结果目录包含：
\begin{itemize}[leftmargin=2.4em]
\item \path{coverage_results.csv}、\path{outer_moments_active_grid.csv}；
\item \path{direct_coefficients.csv}、\path{formula_predictions.csv}；
\item \path{formula_validation_summary.csv}、
      \path{formula_validation_by_density_region.csv}；
\item \path{coefficient_self_checks.json} 以及 MATLAB 原始结果。
\end{itemize}
主入口为 \path{run_plate_outer_cumulant_validation.m}；系数入口为
\path{coverage_coefficients.py}。两者均不调用旧的经验拟合或旧 Edgeworth
系数脚本。

\begin{thebibliography}{99}
\bibitem{Hall1988}
P. Hall, ``Theoretical comparison of bootstrap confidence intervals,''
\emph{The Annals of Statistics}, 16(3), 927--953, 1988.
\url{https://doi.org/10.1214/aos/1176350933}

\bibitem{Hall1992}
P. Hall, \emph{The Bootstrap and Edgeworth Expansion},
Springer, New York, 1992.
\url{https://doi.org/10.1007/978-1-4612-4384-7}

\bibitem{Brillinger1969}
D. R. Brillinger, ``The calculation of cumulants via conditioning,''
\emph{Annals of the Institute of Statistical Mathematics}, 21, 215--218, 1969.
\url{https://doi.org/10.1007/BF02532246}

\bibitem{Owen1997}
A. B. Owen, ``Monte Carlo variance of scrambled net quadrature,''
\emph{SIAM Journal on Numerical Analysis}, 34(5), 1884--1910, 1997.
\url{https://doi.org/10.1137/S0036142994277468}

\bibitem{Loh2003}
W.-L. Loh, ``On the asymptotic distribution of scrambled net quadrature,''
\emph{The Annals of Statistics}, 31(4), 1282--1324, 2003.
\url{https://doi.org/10.1214/aos/1059655914}

\bibitem{NakayamaTuffin2024}
M. K. Nakayama and B. Tuffin, ``Sufficient conditions for central limit
theorems and confidence intervals for randomized quasi-Monte Carlo methods,''
\emph{ACM Transactions on Modeling and Computer Simulation}, 2024.
\url{https://doi.org/10.1145/3643847}
\end{thebibliography}

\end{document}
"""
    tex = tex.replace("@@COEFFICIENT_TABLE@@", coefficient_table(plate_root))
    tex = tex.replace("@@PLATE_OVERALL_TABLE@@", plate_overall_table(plate_root))
    tex = tex.replace("@@PLATE_CORE_TABLE@@", plate_core_table(plate_root))
    tex = tex.replace("@@LINEAR_TABLE@@", linear_nonlinear_table(linear_root))
    output = plate_root / "outer_cumulant_coverage_report.tex"
    output.write_text(tex, encoding="utf-8")
    return output


def main() -> None:
    args = parse_args()
    output = write_report(args.plate_root.resolve(), args.linear_nonlinear_root.resolve())
    print(output)


if __name__ == "__main__":
    main()
