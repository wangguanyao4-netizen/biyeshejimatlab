function outputs = build_formal_e1e8_manuscript_assets(resultsRoot, paperDir)
%BUILD_FORMAL_E1E8_MANUSCRIPT_ASSETS Build paper-ready E1--E8 assets.
%
% This script is intentionally downstream-only: it reads the locked formal
% campaign outputs, writes reproducible chart/table sources, exports figures,
% and drafts a LaTeX results section. It does not rerun experiments.

if nargin < 1 || strlength(string(resultsRoot)) == 0
    resultsRoot = localLatestFormalRoot();
end
if nargin < 2 || strlength(string(paperDir)) == 0
    paperDir = "C:\Users\Wangg\Desktop\期刊论文";
end

resultsRoot = char(string(resultsRoot));
paperDir = char(string(paperDir));
paperTex = fullfile(paperDir, "DPIM_CI_full_integrated_weighted_RBn_natural.tex");
paperFigureDir = fullfile(paperDir, "DPIM_CI_full_integrated_figures", "formal_e1e8");
assetRoot = fullfile(resultsRoot, "_manuscript_assets");
chartDataDir = fullfile(assetRoot, "chart_data");
figureDir = fullfile(assetRoot, "figures");
tableDir = fullfile(assetRoot, "tables");

localEnsureDir(assetRoot);
localEnsureDir(chartDataDir);
localEnsureDir(figureDir);
localEnsureDir(tableDir);
localEnsureDir(paperFigureDir);

fprintf("Building formal E1--E8 manuscript assets from:\n  %s\n", resultsRoot);

S = localLoadSources(resultsRoot);
localWriteChartData(S, chartDataDir);

figureInventory = localBuildFigures(S, figureDir, paperFigureDir);
latexTables = localBuildLatexTables(S, tableDir);
sectionPath = fullfile(paperDir, "formal_e1e8_results_section.tex");
localWriteLatexSection(S, sectionPath);
blueprintPath = fullfile(assetRoot, "manuscript_integration_blueprint.md");
localWriteBlueprint(S, blueprintPath, paperTex, sectionPath, figureInventory, latexTables);
reportPath = fullfile(assetRoot, "formal_e1e8_manuscript_report.html");
localWriteHtmlReport(S, reportPath, figureInventory);

outputs = struct();
outputs.results_root = resultsRoot;
outputs.asset_root = assetRoot;
outputs.chart_data_dir = chartDataDir;
outputs.figure_dir = figureDir;
outputs.paper_figure_dir = paperFigureDir;
outputs.latex_section = sectionPath;
outputs.blueprint = blueprintPath;
outputs.report_html = reportPath;
outputs.figure_inventory = figureInventory;
outputs.latex_tables = latexTables;

dpimnumeric.writeJson(fullfile(assetRoot, "asset_manifest.json"), outputs);
fprintf("Manuscript assets completed:\n  %s\n", assetRoot);
end

function resultsRoot = localLatestFormalRoot()
runs = dir(fullfile(pwd, "results", "weighted_paper_formal_formal_*"));
assert(~isempty(runs), "No weighted_paper_formal_formal_* result directory found.");
[~, idx] = max([runs.datenum]);
resultsRoot = fullfile(runs(idx).folder, runs(idx).name);
end

function S = localLoadSources(root)
S = struct();
S.root = root;
S.claimAudit = localReadCsv(fullfile(root, "claim_audit.csv"));
S.methodAudit = localReadCsv(fullfile(root, "method_audit.csv"));
S.preflight = localReadCsv(fullfile(root, "preflight.csv"));
S.paperMap = localReadCsv(fullfile(root, "paper_experiment_map.csv"));
S.e1 = localReadCsv(fullfile(root, "evidence_suite", "E1_finite_B_rank_audit", "summary.csv"));
S.e2 = localReadCsv(fullfile(root, "evidence_suite", "E2_scalar_normal_R_order_weighted", "confirmation_summary.csv"));
S.e3 = localReadCsv(fullfile(root, "evidence_suite", "E4_standard_normal_nonlinear_tail_weighted", "confirmation_summary.csv"));
S.e4 = localReadCsv(fullfile(root, "evidence_suite", "E5_standard_normal_bootstrap_t_instability_weighted", "confirmation_summary.csv"));
S.e5 = localReadCsv(fullfile(root, "formal_experiments", "E5_probability_weighted_density", "summary.csv"));
S.e6 = localReadCsv(fullfile(root, "formal_experiments", "E6_rqmc_effective_order", "effective_order_summary.csv"));
S.e6Detail = localReadCsv(fullfile(root, "formal_experiments", "E6_rqmc_effective_order", "summary.csv"));
S.e7 = localReadCsv(fullfile(root, "formal_experiments", "E7_finite_grid_band", "summary.csv"));
S.e8 = localReadCsv(fullfile(root, "formal_experiments", "E8_weighted_formula_holdout", "holdout_summary.csv"));
S.e8Pred = localReadCsv(fullfile(root, "formal_experiments", "E8_weighted_formula_holdout", "holdout_predictions.csv"));
S.weightMoments = localReadCsv(fullfile(root, "coverage_formula_validation", "weight_moment_summary.csv"));
end

function T = localReadCsv(path)
if ~isfile(path)
    warning("Missing CSV: %s", path);
    T = table();
    return;
end
T = readtable(path, TextType="string", Delimiter=",");
end

function localWriteChartData(S, outDir)
localWriteTable(S.claimAudit, fullfile(outDir, "formal_claim_audit_clean.csv"));
localWriteTable(S.methodAudit, fullfile(outDir, "formal_method_audit_clean.csv"));
e1 = S.e1;
e1.grid_error_vs_nominal = double(e1.C0B) - 0.95;
e1.simulation_error_vs_C0B_signed = double(e1.simulated_coverage) - double(e1.C0B);
localWriteTable(e1, fullfile(outDir, "E1_finite_B_clean.csv"));
localWriteTable(localCoverageAggregate(S.e2), fullfile(outDir, "E2_coverage_by_R_method.csv"));
localWriteTable(localCoverageAggregate(S.e3), fullfile(outDir, "E3_coverage_by_R_method.csv"));
localWriteTable(localCoverageAggregate(S.e4), fullfile(outDir, "E4_coverage_by_R_method.csv"));
localWriteTable(localCoverageAggregate(S.e5), fullfile(outDir, "E5_coverage_by_R_method.csv"));
localWriteTable(S.e6, fullfile(outDir, "E6_effective_order_clean.csv"));
localWriteTable(localCoverageAggregate(S.e7), fullfile(outDir, "E7_band_coverage_by_R_method.csv"));
localWriteTable(S.e8, fullfile(outDir, "E8_holdout_summary_clean.csv"));
localWriteTable(localWeightMomentAggregate(S.weightMoments), fullfile(outDir, "E8_weight_moment_ranges_clean.csv"));
end

function localWriteTable(T, path)
if isempty(T)
    writetable(table(), path);
else
    writetable(T, path);
end
end

function G = localCoverageAggregate(T)
if isempty(T) || ~all(ismember(["method","R","coverage"], string(T.Properties.VariableNames)))
    G = table();
    return;
end
T.R = double(T.R);
T.coverage = double(T.coverage);
if ismember("coverage_mc_se", string(T.Properties.VariableNames))
    T.coverage_mc_se = double(T.coverage_mc_se);
else
    T.coverage_mc_se = nan(height(T), 1);
end
if ismember("formula_baseline", string(T.Properties.VariableNames))
    T.formula_baseline = double(T.formula_baseline);
else
    T.formula_baseline = nan(height(T), 1);
end
if ismember("abs_coverage_error", string(T.Properties.VariableNames))
    T.abs_coverage_error = double(T.abs_coverage_error);
else
    T.abs_coverage_error = abs(T.coverage - T.formula_baseline);
end
G = groupsummary(T, ["method","R"], "mean", ...
    ["coverage","coverage_mc_se","formula_baseline","abs_coverage_error"]);
G.Properties.VariableNames = strrep(G.Properties.VariableNames, "mean_", "mean_");
end

function G = localWeightMomentAggregate(T)
if isempty(T)
    G = table();
    return;
end
vars = string(T.Properties.VariableNames);
numericVars = intersect(["mean_rho3_w","mean_rho4_w","mean_n2_eff_w","mean_n3_eff_w","mean_n4_eff_w","mean_ess_ratio","mean_max_over_equal"], vars, "stable");
for v = numericVars
    T.(v) = double(T.(v));
end
G = groupsummary(T, ["scheme","scheme_kind"], ["mean","min","max"], numericVars);
end

function inventory = localBuildFigures(S, figDir, paperFigDir)
rows = cell(0, 1);
rows{end+1,1} = localPlotE1(S.e1, figDir, paperFigDir); %#ok<AGROW>
rows{end+1,1} = localPlotClaimAudit(S.claimAudit, figDir, paperFigDir); %#ok<AGROW>
rows{end+1,1} = localPlotMethodAudit(S.methodAudit, figDir, paperFigDir); %#ok<AGROW>
rows{end+1,1} = localPlotCoverageByR(S.e2, "E2 Gaussian/linear coverage by outer repeats", "E2_formal_coverage_by_R", figDir, paperFigDir); %#ok<AGROW>
rows{end+1,1} = localPlotCoverageByR(S.e3, "E3 nonlinear-tail coverage by outer repeats", "E3_formal_coverage_by_R", figDir, paperFigDir); %#ok<AGROW>
rows{end+1,1} = localPlotE4(S.methodAudit, figDir, paperFigDir); %#ok<AGROW>
rows{end+1,1} = localPlotCoverageByR(S.e5, "E5 probability-weighted DPIM coverage by outer repeats", "E5_formal_weighted_coverage_by_R", figDir, paperFigDir); %#ok<AGROW>
rows{end+1,1} = localPlotE6(S.e6, figDir, paperFigDir); %#ok<AGROW>
rows{end+1,1} = localPlotCoverageByR(S.e7, "E7 finite-grid simultaneous band coverage", "E7_formal_band_coverage_by_R", figDir, paperFigDir); %#ok<AGROW>
rows{end+1,1} = localPlotE8Pred(S.e8Pred, figDir, paperFigDir); %#ok<AGROW>
rows{end+1,1} = localPlotE8Weights(S.weightMoments, figDir, paperFigDir); %#ok<AGROW>
inventory = vertcat(rows{:});
writetable(inventory, fullfile(figDir, "formal_e1e8_figure_inventory.csv"));
end

function row = localFigureRow(stem, pngPath, pdfPath, paperPath, note)
row = table(string(stem), string(pngPath), string(pdfPath), string(paperPath), string(note), ...
    'VariableNames', {'figure_id','png_path','pdf_path','paper_path','note'});
end

function row = localPlotE1(T, figDir, paperFigDir)
stem = "E1_finite_B_grid_error";
fig = localFigure();
tl = tiledlayout(fig, 1, 2, TileSpacing="compact", Padding="compact");
nexttile(tl);
hold on;
semilogx(T.B, T.C0B, "-o", LineWidth=1.4, MarkerSize=5);
semilogx(T.B, T.simulated_coverage, "--s", LineWidth=1.2, MarkerSize=5);
yline(0.95, ":", "nominal 0.95", Color=[0.20 0.20 0.20], LineWidth=1.0);
grid on;
xlabel("B");
ylabel("Coverage");
title("Finite-B target is C_{0,B}");
legend(["C_{0,B}","simulation","nominal"], Location="best");
nexttile(tl);
bar(categorical(string(T.B)), double(T.C0B) - 0.95, FaceColor=[0.62 0.74 0.97], EdgeColor=[0.18 0.28 0.50]);
yline(0, ":", Color=[0.20 0.20 0.20]);
grid on;
xlabel("B");
ylabel("C_{0,B} - 0.95");
title("Grid error is finite-order-statistic bias");
sgtitle("E1 finite-B order-statistic baseline", FontWeight="bold");
row = localExportFigure(fig, figDir, paperFigDir, stem, "Finite-B C0B baseline and simulation check.");
end

function row = localPlotClaimAudit(T, figDir, paperFigDir)
stem = "formal_claim_status_overview";
fig = localFigure([100 100 980 660]);
statuses = ["main_text_ready","supplement_only","failed_boundary"];
counts = zeros(numel(statuses), 1);
for i = 1:numel(statuses)
    counts(i) = sum(T.claim_status == statuses(i));
end
bar(categorical(statuses), counts, FaceColor=[0.62 0.74 0.97], EdgeColor=[0.18 0.28 0.50], LineWidth=1.0);
grid on;
ylabel("Number of experiments");
title("Formal claim audit separates main evidence from boundary diagnostics");
ylim([0 max(counts)+1]);
row = localExportFigure(fig, figDir, paperFigDir, stem, "Claim status overview.");
end

function row = localPlotMethodAudit(T, figDir, paperFigDir)
stem = "formal_method_audit_error";
fig = localFigure([100 100 1180 780]);
T = T(~isnan(double(T.mean_abs_error)), :);
cats = categorical(T.paper_experiment + " / " + T.method);
bar(cats, double(T.mean_abs_error), FaceColor=[0.95 0.60 0.43], EdgeColor=[0.50 0.25 0.15], LineWidth=0.8);
grid on;
ylabel("Mean absolute error to baseline");
title("Method-level audit uses C_{0,B}/formula baselines, not a blanket 0.95 target");
xtickangle(45);
row = localExportFigure(fig, figDir, paperFigDir, stem, "Method-level mean absolute coverage error.");
end

function row = localPlotCoverageByR(T, plotTitle, stem, figDir, paperFigDir)
fig = localFigure();
G = localCoverageAggregate(T);
if isempty(G)
    text(0.1, 0.5, "Missing coverage data", FontSize=14);
else
    methods = unique(G.method, "stable");
    colors = localPalette(numel(methods));
    hold on;
    for i = 1:numel(methods)
        part = G(G.method == methods(i), :);
        [~, order] = sort(part.R);
        part = part(order, :);
        plot(part.R, part.mean_coverage, "-o", Color=colors(i,:), LineWidth=1.35, MarkerSize=5, DisplayName=methods(i));
    end
    baseline = localMeanBaseline(T);
    yline(baseline, ":", sprintf("baseline %.3f", baseline), Color=[0.20 0.20 0.20], LineWidth=1.0);
    grid on;
    xlabel("R");
    ylabel("Mean coverage");
    title(plotTitle);
    legend(Location="bestoutside");
    set(gca, XScale="log", XTick=unique(G.R));
end
row = localExportFigure(fig, figDir, paperFigDir, stem, "Mean coverage by R and method.");
end

function row = localPlotE4(T, figDir, paperFigDir)
stem = "E4_bootstrap_t_instability_audit";
T = T(T.paper_experiment == "E4", :);
fig = localFigure();
tl = tiledlayout(fig, 1, 2, TileSpacing="compact", Padding="compact");
nexttile(tl);
bar(categorical(T.method), double(T.min_coverage), FaceColor=[0.95 0.60 0.43], EdgeColor=[0.50 0.25 0.15]);
yline(0.95, ":", Color=[0.20 0.20 0.20]);
grid on;
ylabel("Minimum coverage");
title("Coverage collapses at unstable points");
xtickangle(35);
nexttile(tl);
x = categorical(T.method);
b = bar(x, [double(T.max_fallback_rate), double(T.max_inf_rate)], "grouped");
b(1).FaceColor = [0.62 0.74 0.97];
b(2).FaceColor = [0.95 0.60 0.43];
grid on;
ylabel("Event rate");
title("Failure mechanism is denominator instability");
legend(["fallback","Inf interval"], Location="best");
xtickangle(35);
sgtitle("E4 bootstrap-t small-denominator boundary", FontWeight="bold");
row = localExportFigure(fig, figDir, paperFigDir, stem, "E4 instability and failure-boundary audit.");
end

function row = localPlotE6(T, figDir, paperFigDir)
stem = "E6_rqmc_effective_order_diagnostic";
fig = localFigure();
tl = tiledlayout(fig, 1, 2, TileSpacing="compact", Padding="compact");
nexttile(tl);
bar(categorical(T.method), double(T.effective_order_variance), FaceColor=[0.62 0.74 0.97], EdgeColor=[0.18 0.28 0.50]);
grid on;
ylabel("Estimated variance order");
title("Scrambled RQMC order is weak in this formal run");
xtickangle(25);
nexttile(tl);
bar(categorical(T.method), double(T.variance_ratio_at_max_n_vs_mc), FaceColor=[0.64 0.84 0.46], EdgeColor=[0.20 0.36 0.08]);
yline(1, ":", Color=[0.20 0.20 0.20]);
grid on;
ylabel("Variance ratio at max n vs MC");
title("Use as construction/weight diagnostic only");
xtickangle(25);
sgtitle("E6 RQMC effective-order diagnostic", FontWeight="bold");
row = localExportFigure(fig, figDir, paperFigDir, stem, "RQMC construction valid but formal effective-order gain is weak.");
end

function row = localPlotE8Pred(T, figDir, paperFigDir)
stem = "E8_holdout_observed_vs_predicted";
fig = localFigure([100 100 900 760]);
if isempty(T)
    text(0.1, 0.5, "Missing holdout predictions", FontSize=14);
else
    x = double(T.predicted_error_to_C0B);
    y = double(T.coverage_error_to_C0B);
    isTrain = logical(double(T.is_train));
    hold on;
    scatter(x(isTrain), y(isTrain), 18, [0.70 0.72 0.78], "filled", MarkerFaceAlpha=0.25, DisplayName="train");
    scatter(x(~isTrain), y(~isTrain), 20, [0.33 0.47 0.77], "filled", MarkerFaceAlpha=0.45, DisplayName="holdout");
    lim = max(abs([x; y]), [], "omitnan");
    if isempty(lim) || isnan(lim) || lim == 0
        lim = 0.05;
    end
    plot([-lim lim], [-lim lim], ":", Color=[0.20 0.20 0.20], LineWidth=1.0, DisplayName="ideal");
    xlim([-lim lim]); ylim([-lim lim]);
    axis square;
    grid on;
    xlabel("Predicted error to C_{0,B}");
    ylabel("Observed error to C_{0,B}");
    legend(Location="best");
    title("E8 formula-structure diagnostic, not coefficient proof");
end
row = localExportFigure(fig, figDir, paperFigDir, stem, "Observed vs predicted weighted formula coverage error.");
end

function row = localPlotE8Weights(T, figDir, paperFigDir)
stem = "E8_weight_moment_effective_order";
fig = localFigure([100 100 1180 780]);
G = localWeightMomentAggregate(T);
if isempty(G)
    text(0.1, 0.5, "Missing weight moments", FontSize=14);
else
    labels = categorical(G.scheme);
    vals = [double(G.mean_mean_rho3_w), double(G.mean_mean_rho4_w), double(G.mean_mean_ess_ratio)];
    bar(labels, vals, "grouped");
    grid on;
    ylabel("Moment / effective sample diagnostic");
    title("Weight concentration controls the scale of coverage-error terms");
    legend(["mean rho_3(w)","mean rho_4(w)","mean ESS ratio"], Location="best");
    xtickangle(35);
end
row = localExportFigure(fig, figDir, paperFigDir, stem, "Weight moment and effective-order diagnostic.");
end

function baseline = localMeanBaseline(T)
if ismember("formula_baseline", string(T.Properties.VariableNames))
    baseline = mean(double(T.formula_baseline), "omitnan");
elseif ismember("C0B", string(T.Properties.VariableNames))
    baseline = mean(double(T.C0B), "omitnan");
else
    baseline = 0.95;
end
if isnan(baseline)
    baseline = 0.95;
end
end

function fig = localFigure(pos)
if nargin < 1
    pos = [100 100 1120 760];
end
fig = figure(Visible="off", Color="w", Position=pos);
set(fig, DefaultAxesFontName="Microsoft YaHei");
set(fig, DefaultTextFontName="Microsoft YaHei");
set(fig, DefaultAxesFontSize=10);
end

function colors = localPalette(n)
base = [
    0.32 0.47 0.77
    0.94 0.60 0.43
    0.64 0.84 0.46
    0.72 0.63 0.22
    0.74 0.34 0.61
    0.46 0.51 0.56];
idx = mod(0:n-1, size(base,1)) + 1;
colors = base(idx,:);
end

function row = localExportFigure(fig, figDir, paperFigDir, stem, note)
pngPath = fullfile(figDir, stem + ".png");
pdfPath = fullfile(figDir, stem + ".pdf");
paperPath = fullfile(paperFigDir, stem + ".png");
try
    exportgraphics(fig, pngPath, Resolution=450);
    exportgraphics(fig, pdfPath, ContentType="vector");
catch ME
    warning("exportgraphics failed for %s: %s. Falling back to saveas.", stem, ME.message);
    saveas(fig, pngPath);
    saveas(fig, pdfPath);
end
copyfile(pngPath, paperPath, "f");
close(fig);
row = localFigureRow(stem, pngPath, pdfPath, paperPath, note);
end

function latexTables = localBuildLatexTables(S, tableDir)
claimPath = fullfile(tableDir, "formal_claim_audit_compact.tex");
methodPath = fullfile(tableDir, "formal_method_audit_compact.tex");
e1Path = fullfile(tableDir, "E1_finite_B_table.tex");
localWriteText(claimPath, localClaimAuditTable(S.claimAudit));
localWriteText(methodPath, localMethodAuditTable(S.methodAudit));
localWriteText(e1Path, localE1Table(S.e1));
latexTables = [string(claimPath); string(methodPath); string(e1Path)];
end

function txt = localClaimAuditTable(T)
lines = strings(0,1);
lines(end+1) = "\begin{table}[H]";
lines(end+1) = "\centering";
lines(end+1) = "\small";
lines(end+1) = "\caption{正式 E1--E8 数值证据的主文使用边界。}";
lines(end+1) = "\label{tab:formal-e1e8-claim-audit}";
lines(end+1) = "\begin{tabularx}{\textwidth}{p{0.08\textwidth}p{0.23\textwidth}p{0.18\textwidth}X}";
lines(end+1) = "\toprule";
lines(end+1) = "实验 & 证据对象 & 审计状态 & 正文使用边界 \\";
lines(end+1) = "\midrule";
for i = 1:height(T)
    lines(end+1) = sprintf("%s & %s & \\texttt{%s} & %s \\\\", ...
        localTex(T.paper_experiment(i)), localChineseExperimentTitle(T.paper_experiment(i)), ...
        localTex(T.claim_status(i)), localChineseDecision(T.paper_experiment(i), T.claim_status(i), T.main_text_decision(i)));
end
lines(end+1) = "\bottomrule";
lines(end+1) = "\end{tabularx}";
lines(end+1) = "\end{table}";
txt = strjoin(lines, newline);
end

function txt = localMethodAuditTable(T)
coverageExperiments = ["E2","E3","E4","E5","E7"];
keep = T(ismember(T.paper_experiment, coverageExperiments) & ~isnan(double(T.mean_coverage)), :);
lines = strings(0,1);
lines(end+1) = "\begin{table}[H]";
lines(end+1) = "\centering";
lines(end+1) = "\footnotesize";
lines(end+1) = "\caption{方法层覆盖率误差审计摘要。误差基准为有限 $B$ 的 $C_{0,B}$ 或相应公式基准。}";
lines(end+1) = "\label{tab:formal-method-audit}";
lines(end+1) = "\begin{tabularx}{\textwidth}{p{0.07\textwidth}p{0.25\textwidth}rrrrX}";
lines(end+1) = "\toprule";
lines(end+1) = "实验 & 方法 & 均值覆盖 & 最小覆盖 & 平均绝对误差 & 最大 MCSE & 状态 \\";
lines(end+1) = "\midrule";
for i = 1:height(keep)
    lines(end+1) = sprintf("%s & %s & %s & %s & %s & %s & \\texttt{%s} \\\\", ...
        localTex(keep.paper_experiment(i)), localTex(keep.method(i)), ...
        localNum(double(keep.mean_coverage(i))), localNum(double(keep.min_coverage(i))), ...
        localNum(double(keep.mean_abs_error(i))), localNum(double(keep.max_mc_se(i))), ...
        localTex(keep.method_status(i)));
end
lines(end+1) = "\bottomrule";
lines(end+1) = "\end{tabularx}";
lines(end+1) = "\end{table}";
txt = strjoin(lines, newline);
end

function txt = localE1Table(T)
lines = strings(0,1);
lines(end+1) = "\begin{table}[H]";
lines(end+1) = "\centering";
lines(end+1) = "\caption{有限 $B$ 次序统计量基准覆盖率。}";
lines(end+1) = "\label{tab:formal-e1-finite-B}";
lines(end+1) = "\begin{tabular}{rrrrr}";
lines(end+1) = "\toprule";
lines(end+1) = "$B$ & $k_-$ & $k_+$ & $C_{0,B}$ & $C_{0,B}-0.95$ \\";
lines(end+1) = "\midrule";
for i = 1:height(T)
    lines(end+1) = sprintf("%d & %d & %d & %.5f & %.5f \\\\", ...
        round(double(T.B(i))), round(double(T.k_minus(i))), round(double(T.k_plus(i))), ...
        double(T.C0B(i)), double(T.C0B(i)) - 0.95);
end
lines(end+1) = "\bottomrule";
lines(end+1) = "\end{tabular}";
lines(end+1) = "\end{table}";
txt = strjoin(lines, newline);
end

function localWriteLatexSection(S, sectionPath)
audit = S.claimAudit;
method = S.methodAudit;
e1 = S.e1;
e6 = S.e6;
e8 = S.e8;
texRho = char(92) + "rho";

e2Student = localMethodRow(method, "E2", "Student-t");
e3Student = localMethodRow(method, "E3", "Student-t");
e4Bt = localMethodRow(method, "E4", "Bootstrap-t");
e5Student = localMethodRow(method, "E5", "Student-t");
e7Bt = localMethodRow(method, "E7", "Bootstrap-t max-stat band");

lines = strings(0,1);
lines(end+1) = "% Auto-generated by build_formal_e1e8_manuscript_assets.m";
lines(end+1) = "% Source root: " + string(S.root);
lines(end+1) = "\section{正式数值实验结果与覆盖率审计}";
lines(end+1) = "";
lines(end+1) = "本节把前述 E1--E8 实验协议落实为可复现的 formal 数值证据。所有主证据均采用随机化 RQMC 样本池与 Voronoi 概率权重；等权结果仅作为固定权重特例或对照，不进入主结论。覆盖率比较的基准不是统一的 $0.95$，而是有限 $B$ 次序统计量基准 $C_{0,B}$ 或对应的公式基准。";
lines(end+1) = "";
lines(end+1) = "本轮 formal campaign 的预检显示：$M=1000$、$B=999$ 的规模门槛通过，样本池为 scrambled Sobol RQMC，CPU/GPU Voronoi 权重 parity smoke test 通过。图表和表格均来自锁定的 formal 结果目录，正文只采用经审计允许进入主文的正面结果；失效边界和公式诊断按其真实状态写入。";
lines(end+1) = "";
lines(end+1) = localClaimAuditTable(audit);
lines(end+1) = "";
lines(end+1) = "\begin{figure}[H]";
lines(end+1) = "\centering";
lines(end+1) = "\SafeIncludeGraphics[width=0.82\textwidth]{formal_e1e8/formal_claim_status_overview.png}";
lines(end+1) = "\caption{E1--E8 formal 审计状态总览。该图说明本文只把 E1/E2/E3/E5/E7 写作主文正面证据，E4 写作失效边界，E6/E8 写作补充或诊断证据。}";
lines(end+1) = "\label{fig:formal-claim-status}";
lines(end+1) = "\end{figure}";
lines(end+1) = "";
lines(end+1) = "\subsection{E1：有限 $B$ 次序统计量基准}";
lines(end+1) = "表 \ref{tab:formal-e1-finite-B} 和图 \ref{fig:formal-e1-finite-B} 直接检验有限 $B$ 分位数端点的基准覆盖率。结果显示，模拟覆盖率贴近 $C_{0,B}$，而 $C_{0,B}$ 本身随 $B$ 和端点整数规则产生确定网格偏移。这证明后续 bootstrap 覆盖率不能只按名义 $0.95$ 判断。";
lines(end+1) = localE1Table(e1);
lines(end+1) = "\begin{figure}[H]";
lines(end+1) = "\centering";
lines(end+1) = "\SafeIncludeGraphics[width=0.90\textwidth]{formal_e1e8/E1_finite_B_grid_error.png}";
lines(end+1) = "\caption{有限 $B$ 基准 $C_{0,B}$ 与模拟覆盖率对比。左图检查程序端点约定，右图给出 $C_{0,B}-0.95$ 的确定网格误差。}";
lines(end+1) = "\label{fig:formal-e1-finite-B}";
lines(end+1) = "\end{figure}";
lines(end+1) = "";
lines(end+1) = "\subsection{E2：高斯线性闭式基准}";
lines(end+1) = sprintf("E2 作为闭式或高精度基准校验，Student-t 方法的平均覆盖率为 %s，最小覆盖率为 %s，平均绝对误差为 %s，未出现 fallback 或无穷端点。该结果支持程序端点、RQMC 样本池和覆盖率统计流程在对称线性基准下是可用的；但它只证明数值阶和端点实现合理，不等价于一般非线性响应的完整公式证明。", localNum(double(e2Student.mean_coverage)), localNum(double(e2Student.min_coverage)), localNum(double(e2Student.mean_abs_error)));
lines(end+1) = "\begin{figure}[H]";
lines(end+1) = "\centering";
lines(end+1) = "\SafeIncludeGraphics[width=0.90\textwidth]{formal_e1e8/E2_formal_coverage_by_R.png}";
lines(end+1) = "\caption{E2 高斯线性基准下各方法平均覆盖率随外层重复数 $R$ 的变化。基准线为相应 $C_{0,B}$ 或 nominal/formula baseline。}";
lines(end+1) = "\label{fig:formal-e2-coverage}";
lines(end+1) = "\end{figure}";
lines(end+1) = "";
lines(end+1) = "\subsection{E3：非线性尾部与偏度效应}";
lines(end+1) = sprintf("E3 的 Student-t 平均覆盖率为 %s，最小覆盖率为 %s，平均绝对误差为 %s；bootstrap-$t$ 平均覆盖率较高，但个别尾部点仍有欠覆盖。该结果说明偏态或非线性响应下覆盖误差与局部偏度、峰度和核平滑位置有关，不能简单归结为 $M$ 或 $B$ 不足。", localNum(double(e3Student.mean_coverage)), localNum(double(e3Student.min_coverage)), localNum(double(e3Student.mean_abs_error)));
lines(end+1) = "\begin{figure}[H]";
lines(end+1) = "\centering";
lines(end+1) = "\SafeIncludeGraphics[width=0.90\textwidth]{formal_e1e8/E3_formal_coverage_by_R.png}";
lines(end+1) = "\caption{E3 非线性尾部算例的覆盖率随 $R$ 变化。尾部位置的偏度和峰度会改变 percentile 与 bootstrap-$t$ 的相对误差。}";
lines(end+1) = "\label{fig:formal-e3-coverage}";
lines(end+1) = "\end{figure}";
lines(end+1) = "";
lines(end+1) = "\subsection{E4：Bootstrap-$t$ 小分母失效边界}";
lines(end+1) = sprintf("E4 不作为正面成功结果使用。Bootstrap-$t$ 行的平均覆盖率为 %s，但最小覆盖率降至 %s，最大 fallback 频率为 %s，最大无穷端点比例为 %s。该实验的作用是证明小分母事件会导致区间长度和覆盖行为失稳，因此混合替换规则或正则化规则必须作为边界机制讨论。", localNum(double(e4Bt.mean_coverage)), localNum(double(e4Bt.min_coverage)), localNum(double(e4Bt.max_fallback_rate)), localNum(double(e4Bt.max_inf_rate)));
lines(end+1) = "\begin{figure}[H]";
lines(end+1) = "\centering";
lines(end+1) = "\SafeIncludeGraphics[width=0.90\textwidth]{formal_e1e8/E4_bootstrap_t_instability_audit.png}";
lines(end+1) = "\caption{E4 小分母失效边界。左图显示覆盖率在不稳定位置退化，右图显示 fallback 与无穷端点事件是主要诊断信号。}";
lines(end+1) = "\label{fig:formal-e4-instability}";
lines(end+1) = "\end{figure}";
lines(end+1) = "";
lines(end+1) = "\subsection{E5：概率加权 DPIM 密度估计主实验}";
lines(end+1) = sprintf("E5 是本文正式概率加权 DPIM 主实验。Student-t 方法的平均覆盖率为 %s，最小覆盖率为 %s，平均绝对误差为 %s，最大 MCSE 为 %s；所有主结果均使用 scrambled Sobol 样本池和 Voronoi 概率权重。由于最小覆盖率仍低于平均水平，正文表述应强调位置相关的局部欠覆盖，而不是宣称所有位置都达到相同覆盖水平。", localNum(double(e5Student.mean_coverage)), localNum(double(e5Student.min_coverage)), localNum(double(e5Student.mean_abs_error)), localNum(double(e5Student.max_mc_se)));
lines(end+1) = "\begin{figure}[H]";
lines(end+1) = "\centering";
lines(end+1) = "\SafeIncludeGraphics[width=0.90\textwidth]{formal_e1e8/E5_formal_weighted_coverage_by_R.png}";
lines(end+1) = "\caption{E5 概率加权 DPIM 密度估计的平均覆盖率。该图只使用概率加权 RQMC 主证据，等权结果不作为主文结论。}";
lines(end+1) = "\label{fig:formal-e5-weighted}";
lines(end+1) = "\end{figure}";
lines(end+1) = "";
lines(end+1) = "\subsection{E6：RQMC 构造与有效阶诊断}";
lines(end+1) = sprintf("E6 证明 RQMC 样本池和 Voronoi 权重构造可运行且无 fallback，但不能支持强有效阶加速结论。scrambled Sobol 的估计方差阶为 %s，在最大 $n$ 下相对 MC 的方差比为 %s。因此本文只能把 E6 写作构造与权重矩诊断，不能把它写成 RQMC 覆盖率自动改善的证明。", localNum(double(e6.effective_order_variance(e6.method=="sobol_scrambled"))), localNum(double(e6.variance_ratio_at_max_n_vs_mc(e6.method=="sobol_scrambled"))));
lines(end+1) = "\begin{figure}[H]";
lines(end+1) = "\centering";
lines(end+1) = "\SafeIncludeGraphics[width=0.90\textwidth]{formal_e1e8/E6_rqmc_effective_order_diagnostic.png}";
lines(end+1) = "\caption{E6 RQMC 有效阶诊断。plain QMC 是确定性对照，scrambled RQMC 在本 formal 设置下方差阶增益较弱。}";
lines(end+1) = "\label{fig:formal-e6-rqmc}";
lines(end+1) = "\end{figure}";
lines(end+1) = "";
lines(end+1) = "\subsection{E7：有限网格同时置信带}";
lines(end+1) = sprintf("E7 的 bootstrap-$t$ 最大统计量置信带平均覆盖率为 %s，最小覆盖率为 %s，平均绝对误差为 %s，未出现无穷端点。该结论只针对有限网格 simultaneous band，不能外推为连续响应曲线覆盖。", localNum(double(e7Bt.mean_coverage)), localNum(double(e7Bt.min_coverage)), localNum(double(e7Bt.mean_abs_error)));
lines(end+1) = "\begin{figure}[H]";
lines(end+1) = "\centering";
lines(end+1) = "\SafeIncludeGraphics[width=0.90\textwidth]{formal_e1e8/E7_formal_band_coverage_by_R.png}";
lines(end+1) = "\caption{E7 有限网格同时置信带覆盖率。图中覆盖率只对应离散响应网格，不声明连续曲线 simultaneous coverage。}";
lines(end+1) = "\label{fig:formal-e7-band}";
lines(end+1) = "\end{figure}";
lines(end+1) = "";
lines(end+1) = "\subsection{E8：概率权重覆盖表达式的数值诊断}";
lines(end+1) = sprintf("E8 对概率权重覆盖表达式进行经验结构诊断。Percentile bootstrap 的 held-out MAE 为 %s，bootstrap-$t$ 的 held-out MAE 为 %s，最大 fallback 频率为 %s。该结果说明 $%s_3(w)$、$%s_4(w)$ 等权重矩项对覆盖误差有数值解释力，但不能声称未知系数 $A$ 已从理论上逐项证明。", localNum(double(e8.mae(e8.method=="Percentile bootstrap"))), localNum(double(e8.mae(e8.method=="Bootstrap-t"))), localNum(max(double(e8.max_fallback_rate), [], "omitnan")), texRho, texRho);
lines(end+1) = "\begin{figure}[H]";
lines(end+1) = "\centering";
lines(end+1) = "\SafeIncludeGraphics[width=0.82\textwidth]{formal_e1e8/E8_holdout_observed_vs_predicted.png}";
lines(end+1) = "\caption{E8 加权覆盖表达式的 observed--predicted 诊断。散点贴近对角线说明结构项有解释力；离散残差和 fallback 事件说明该图不是未知系数的理论证明。}";
lines(end+1) = "\label{fig:formal-e8-prediction}";
lines(end+1) = "\end{figure}";
lines(end+1) = "\begin{figure}[H]";
lines(end+1) = "\centering";
lines(end+1) = "\SafeIncludeGraphics[width=0.90\textwidth]{formal_e1e8/E8_weight_moment_effective_order.png}";
lines(end+1) = "\caption{E8 权重矩与有效样本量诊断。权重集中会增大 $" + texRho + "_3(w)$ 与 $" + texRho + "_4(w)$，从而放大覆盖误差展开中的高阶项。}";
lines(end+1) = "\label{fig:formal-e8-weight-moments}";
lines(end+1) = "\end{figure}";
lines(end+1) = "";
lines(end+1) = localMethodAuditTable(method);
lines(end+1) = "";
lines(end+1) = "\paragraph{综合判断.}";
lines(end+1) = "正式数值结果已经足以支撑本文的有限 $B$ 基准、线性闭式校验、非线性偏态机制、概率加权 DPIM 主实验和有限网格 simultaneous band 结论。仍需降级处理的部分有三点：E4 必须写成失效边界，E6 只能写成 RQMC 构造和权重矩诊断，E8 只能写成固定权重公式结构的数值校准与随机权重边界诊断。";

localWriteText(sectionPath, strjoin(lines, newline));
end

function row = localMethodRow(T, expId, methodName)
idx = find(T.paper_experiment == expId & T.method == methodName, 1);
assert(~isempty(idx), "Missing method audit row: %s / %s", expId, methodName);
row = T(idx, :);
end

function localWriteBlueprint(S, path, paperTex, sectionPath, figureInventory, latexTables)
lines = strings(0,1);
lines(end+1) = "# Formal E1--E8 Manuscript Integration Blueprint";
lines(end+1) = "";
lines(end+1) = "## Scope";
lines(end+1) = "";
lines(end+1) = "- Core paper TeX: `" + string(paperTex) + "`";
lines(end+1) = "- Generated LaTeX section: `" + string(sectionPath) + "`";
lines(end+1) = "- Formal results root: `" + string(S.root) + "`";
lines(end+1) = "- Evidence policy: main text uses probability-weighted RQMC/Voronoi outputs; equal-weight rows are controls only.";
lines(end+1) = "- Coverage policy: compare to finite-B `C0B` or formula baselines, not a blanket `0.95` target.";
lines(end+1) = "";
lines(end+1) = "## TeX Insertion";
lines(end+1) = "";
lines(end+1) = "- Insert `\input{formal_e1e8_results_section.tex}` immediately before `\section{解析性总表与论文使用边界}`.";
lines(end+1) = "- Do not rewrite the theory sections during this integration pass.";
lines(end+1) = "- After compilation, review whether the older fourth-chapter legacy examples should move to appendix or be explicitly marked as historical diagnostics.";
lines(end+1) = "";
lines(end+1) = "## Figure Inventory";
for i = 1:height(figureInventory)
    lines(end+1) = sprintf("- `%s`: `%s`", figureInventory.figure_id(i), figureInventory.paper_path(i));
end
lines(end+1) = "";
lines(end+1) = "## Table Inventory";
for i = 1:numel(latexTables)
    lines(end+1) = "- `" + latexTables(i) + "`";
end
lines(end+1) = "";
lines(end+1) = "## Claim Decisions";
for i = 1:height(S.claimAudit)
    lines(end+1) = sprintf("- `%s`: `%s`; decision `%s`; note: %s", ...
        S.claimAudit.paper_experiment(i), S.claimAudit.claim_status(i), ...
        S.claimAudit.main_text_decision(i), S.claimAudit.note(i));
end
lines(end+1) = "";
lines(end+1) = "## Required QA";
lines(end+1) = "";
lines(end+1) = "- Compile with XeLaTeX/MiKTeX and inspect the rendered PDF pages containing the inserted section.";
lines(end+1) = "- Confirm no missing figures, no severe overfull boxes, and no claim stronger than `claim_audit.csv`.";
lines(end+1) = "- Render final PDF pages to PNG and visually compare spacing, captions, and table density against the journal reference PDFs.";
localWriteText(path, strjoin(lines, newline));
end

function localWriteHtmlReport(S, path, figureInventory)
lines = strings(0,1);
lines(end+1) = "<!doctype html><html><head><meta charset=""utf-8"">";
lines(end+1) = "<title>Formal E1-E8 Manuscript Report</title>";
lines(end+1) = "<style>body{font-family:'Segoe UI','Microsoft YaHei',sans-serif;max-width:1120px;margin:32px auto;line-height:1.65;color:#1f2430}h1,h2{line-height:1.25}table{border-collapse:collapse;width:100%;font-size:13px}th,td{border-bottom:1px solid #e6e8f0;padding:7px;text-align:left}th{background:#f4f5f7}.status{font-family:Consolas,monospace}img{max-width:100%;border:1px solid #e6e8f0;margin:8px 0 24px}.note{color:#6f768a}</style></head><body>";
lines(end+1) = "<h1>Formal E1-E8 Manuscript Report</h1>";
lines(end+1) = "<p class=""note"">Source root: <code>" + string(S.root) + "</code></p>";
lines(end+1) = "<h2>Technical Summary</h2>";
lines(end+1) = "<p>The formal campaign supports E1/E2/E3/E5/E7 as main-text candidates, E4 as a failure-boundary result, and E6/E8 as supplement or diagnostic evidence. Coverage is interpreted against finite-B C0B or formula baselines.</p>";
lines(end+1) = "<h2>Claim Audit</h2><table><tr><th>Experiment</th><th>Claim</th><th>Status</th><th>Decision</th></tr>";
for i = 1:height(S.claimAudit)
    lines(end+1) = "<tr><td>" + S.claimAudit.paper_experiment(i) + "</td><td>" + S.claimAudit.claim_type(i) + "</td><td class=""status"">" + S.claimAudit.claim_status(i) + "</td><td>" + S.claimAudit.main_text_decision(i) + "</td></tr>";
end
lines(end+1) = "</table><h2>Figures</h2>";
for i = 1:height(figureInventory)
    rel = erase(figureInventory.png_path(i), fileparts(path) + filesep);
    lines(end+1) = "<h3>" + figureInventory.figure_id(i) + "</h3>";
    lines(end+1) = "<p>" + figureInventory.note(i) + "</p>";
    lines(end+1) = "<img src=""figures/" + figureInventory.figure_id(i) + ".png"" alt=""" + figureInventory.figure_id(i) + """>";
end
lines(end+1) = "<h2>Limitations</h2>";
lines(end+1) = "<p>E6 does not establish strong RQMC effective-order acceleration. E8 is a numerical structure diagnostic, not a theorem proof for unknown expansion coefficients. E4 must be written as denominator-instability boundary evidence.</p>";
lines(end+1) = "</body></html>";
localWriteText(path, strjoin(lines, newline));
end

function s = localChineseDecision(expId, status, decision)
switch string(expId)
    case "E1"
        s = "可作为有限 $B$ 基准进入主文。";
    case "E2"
        s = "可作为线性闭式或高精度基准进入主文。";
    case "E3"
        s = "可作为偏态/非线性机制证据进入主文，但不得写成完整系数证明。";
    case "E4"
        s = "只作为 bootstrap-$t$ 小分母失效边界。";
    case "E5"
        s = "可作为概率加权 DPIM 主实验，但需说明局部欠覆盖。";
    case "E6"
        s = "仅作 RQMC 构造与权重矩诊断。";
    case "E7"
        s = "可作为有限网格 simultaneous band 证据。";
    case "E8"
        s = "仅作加权覆盖表达式结构诊断。";
    otherwise
        s = sprintf("\\texttt{%s}/\\texttt{%s}", localTex(status), localTex(decision));
end
end

function s = localChineseExperimentTitle(expId)
switch string(expId)
    case "E1"
        s = "有限 $B$ 次序统计量网格误差";
    case "E2"
        s = "高斯/线性 DPIM 闭式基准";
    case "E3"
        s = "偏态/非线性响应覆盖修正";
    case "E4"
        s = "Bootstrap-$t$ 小分母失稳";
    case "E5"
        s = "概率加权 DPIM 密度估计";
    case "E6"
        s = "RQMC 样本池有效阶诊断";
    case "E7"
        s = "有限网格同时置信带";
    case "E8"
        s = "加权 $R$--$B$--$w$ 覆盖表达式诊断";
    otherwise
        s = localTex(expId);
end
end

function s = localTex(x)
s = char(string(x));
s = strrep(s, "\", "\textbackslash{}");
s = strrep(s, "_", "\_");
s = strrep(s, "%", "\%");
s = strrep(s, "&", "\&");
s = strrep(s, "#", "\#");
end

function s = localNum(x)
if isempty(x) || isnan(x)
    s = "--";
elseif abs(x) >= 1000 || (abs(x) > 0 && abs(x) < 1e-3)
    s = sprintf("%.3g", x);
else
    s = sprintf("%.4f", x);
end
end

function localEnsureDir(path)
if ~isfolder(path)
    mkdir(path);
end
end

function localWriteText(path, txt)
fid = fopen(path, "w", "n", "UTF-8");
assert(fid > 0, "Cannot write %s", path);
c = onCleanup(@() fclose(fid));
fprintf(fid, "%s", txt);
end
