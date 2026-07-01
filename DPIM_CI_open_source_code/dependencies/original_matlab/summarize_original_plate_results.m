function summaryRoot = summarize_original_plate_results(resultRoot)
%summarize_original_plate_results Summarize original thin-plate E7/E8 outputs.
%
% This post-processor does not rerun GF/KL/FEM. It recomputes curve pools
% from the saved wc responses and Voronoi weights, then reports all-grid and
% active-grid coverage diagnostics for the regenerated original plate run.

if nargin < 1 || strlength(string(resultRoot)) == 0
    resultRoot = localLatestOriginalPlateRoot();
end
resultRoot = char(string(resultRoot));
summaryRoot = fullfile(resultRoot, "_summary");
dpimnumeric.ensureDir(summaryRoot);

config = jsondecode(fileread(fullfile(resultRoot, "config.json")));
inputs = localLoadOriginalPlateInputs(resultRoot);

E7 = readtable(fullfile(resultRoot, "E7_plate_SFEM", "summary.csv"));
E8 = readtable(fullfile(resultRoot, "E8_simultaneous_band", "summary.csv"));
P7 = readtable(fullfile(resultRoot, "E7_plate_SFEM", "pointwise_by_h.csv"));
P8 = readtable(fullfile(resultRoot, "E8_simultaneous_band", "pointwise_band_by_h.csv"));
D8 = load(fullfile(resultRoot, "E8_simultaneous_band", ...
    "curve_coverage_diagnostics.mat"), "diagnosticBlocks");

activeSummary = [
    localActiveSummaryForExperiment("E7_plate_SFEM", P7, inputs, config, false, {})
    localActiveSummaryForExperiment("E8_simultaneous_band", P8, inputs, config, ...
    true, D8.diagnosticBlocks)
    ];
writetable(activeSummary, fullfile(summaryRoot, "active_grid_summary.csv"));

overall = [
    localOverallRows("E7_plate_SFEM", E7)
    localOverallRows("E8_simultaneous_band", E8)
    ];
writetable(overall, fullfile(summaryRoot, "overall_method_summary.csv"));

bestH = localBestH(activeSummary);
writetable(bestH, fullfile(summaryRoot, "best_h_active_grid.csv"));

localWriteReport(summaryRoot, resultRoot, config, overall, activeSummary, bestH);
fprintf("Original plate summary written: %s\n", summaryRoot);
end

function resultRoot = localLatestOriginalPlateRoot()
files = dir(fullfile(pwd, "results", "original_plate_full_*"));
files = files([files.isdir]);
if isempty(files)
    error("No results/original_plate_full_* directory found.");
end
[~, idx] = max([files.datenum]);
resultRoot = fullfile(files(idx).folder, files(idx).name);
end

function inputs = localLoadOriginalPlateInputs(resultRoot)
refData = load(fullfile(resultRoot, "reference", "thin_reference.mat"), "reference");
poolData = load(fullfile(resultRoot, "rqmc_pool", "thin_rqmc_pool.mat"), "pool");
inputs = struct();
inputs.reference = refData.reference;
inputs.pool = poolData.pool;
inputs.w_list = refData.reference.w_list;
inputs.h_list = refData.reference.h_list;
end

function tbl = localOverallRows(experimentName, summaryTable)
methods = unique(string(summaryTable.method), "stable");
rows = cell(numel(methods), 1);
for iMethod = 1:numel(methods)
    mask = string(summaryTable.method) == methods(iMethod);
    sub = summaryTable(mask, :);
    rows{iMethod} = table(string(experimentName), methods(iMethod), height(sub), ...
        mean(sub.coverage, "omitnan"), mean(sub.abs_coverage_error, "omitnan"), ...
        mean(sub.fallback_rate, "omitnan"), mean(sub.bootstrap_t_inf_rate, "omitnan"), ...
        mean(sub.mean_interval_length, "omitnan"), ...
        'VariableNames', {'experiment', 'method', 'h_count', 'mean_coverage', ...
        'mean_abs_coverage_error', 'mean_fallback_rate', 'mean_bootstrap_t_inf_rate', ...
        'mean_interval_length'});
end
tbl = vertcat(rows{:});
end

function tbl = localActiveSummaryForExperiment(experimentName, pointwise, inputs, ...
        config, hasBand, diagnosticBlocks)
floors = [0, 1e-12, 1e-10, 1e-8];
hList = inputs.h_list(:).';
rows = cell(numel(hList) * numel(floors), 1);
rowId = 0;
for ih = 1:numel(hList)
    h = hList(ih);
    [referenceCurve, curvePool] = localCurvePoolAtH(inputs, h);
    [degenerateGrid, activeGrid] = localActiveGrid(curvePool, referenceCurve);
    pointRows = pointwise(pointwise.h_index == ih, :);
    for iFloor = 1:numel(floors)
        floorValue = floors(iFloor);
        if floorValue == 0
            mask = activeGrid;
            activeName = "active_legacy";
        else
            mask = activeGrid & (abs(referenceCurve) >= floorValue);
            activeName = "active_ref_abs_ge_" + string(floorValue);
        end
        rowId = rowId + 1;
        bandDiagnostic = [];
        if hasBand
            bandDiagnostic = diagnosticBlocks{ih}.band;
        end
        rows{rowId} = localOneActiveRow(experimentName, ih, h, activeName, floorValue, ...
            mask, degenerateGrid, pointRows, config, hasBand, bandDiagnostic);
    end
end
tbl = vertcat(rows{:});
end

function row = localOneActiveRow(experimentName, hIndex, h, activeName, floorValue, ...
        mask, degenerateGrid, pointRows, config, hasBand, bandDiagnostic)
activeCount = sum(mask);
if activeCount > 0
    btCoverage = mean(pointRows.bootstrap_t_coverage(mask), "omitnan");
    fbCoverage = mean(pointRows.bootstrap_t_fallback_rule_coverage(mask), "omitnan");
    fbRate = mean(pointRows.fallback_count(mask) / config.M, "omitnan");
    infRate = mean(pointRows.bootstrap_t_inf_count(mask) / config.M, "omitnan");
else
    btCoverage = NaN;
    fbCoverage = NaN;
    fbRate = NaN;
    infRate = NaN;
end
if hasBand && activeCount > 0
    bandPointwiseMean = mean(pointRows.simultaneous_band_coverage(mask), "omitnan");
    bandCoverage = mean(all(bandDiagnostic.band_point_hit(:, mask), 2));
else
    bandCoverage = NaN;
    bandPointwiseMean = NaN;
end
row = table(string(experimentName), hIndex, h, string(activeName), floorValue, ...
    activeCount, sum(degenerateGrid), numel(mask), btCoverage, fbCoverage, ...
    fbRate, infRate, bandCoverage, bandPointwiseMean, ...
    'VariableNames', {'experiment', 'h_index', 'h', 'active_definition', ...
    'reference_abs_floor', 'active_grid_count', 'degenerate_grid_count', ...
    'grid_count', 'bootstrap_t_coverage', 'bootstrap_t_fallback_rule_coverage', ...
    'fallback_rate', 'bootstrap_t_inf_rate', 'simultaneous_band_coverage', ...
    'simultaneous_band_mean_pointwise_coverage'});
end

function [referenceCurve, curvePool] = localCurvePoolAtH(inputs, h)
hIdx = find(abs(inputs.h_list - h) <= max(1e-12, 1e-12 * abs(h)), 1);
if isempty(hIdx)
    error("h %.16g is not present in the saved h list.", h);
end
referenceCurve = inputs.reference.curves_by_h(:, hIdx).';
nPool = inputs.pool.sample_pool_count;
nGrid = numel(inputs.w_list);
curvePool = zeros(nPool, nGrid);
for iRun = 1:nPool
    curvePool(iRun, :) = compute_probability_estimate_curve( ...
        inputs.pool.wc_runs_thin{iRun}, inputs.pool.weight_runs{iRun}, h, inputs.w_list);
end
end

function [degenerateGrid, activeGrid] = localActiveGrid(curvePool, referenceCurve)
poolMean = mean(curvePool, 1);
poolSd = std(curvePool, 0, 1);
toleranceValue = localElementTolerance(poolMean, referenceCurve);
degenerateGrid = (poolSd <= toleranceValue) & (abs(poolMean - referenceCurve) <= toleranceValue);
activeGrid = ~degenerateGrid;
end

function toleranceValue = localElementTolerance(varargin)
scale = zeros(size(varargin{1}));
for iInput = 1:nargin
    scale = max(scale, abs(varargin{iInput}));
end
toleranceValue = max(realmin, 100 * eps(scale));
end

function bestH = localBestH(activeSummary)
methods = ["bootstrap_t_coverage"; "bootstrap_t_fallback_rule_coverage"; ...
    "simultaneous_band_coverage"];
defs = unique(string(activeSummary.active_definition), "stable");
experiments = unique(string(activeSummary.experiment), "stable");
rows = {};
for iExp = 1:numel(experiments)
    for iDef = 1:numel(defs)
        mask = string(activeSummary.experiment) == experiments(iExp) ...
            & string(activeSummary.active_definition) == defs(iDef);
        sub = activeSummary(mask, :);
        for iMethod = 1:numel(methods)
            values = sub.(methods(iMethod));
            valid = isfinite(values);
            if any(valid)
                [~, localIdx] = min(abs(values(valid) - 0.95));
                validRows = find(valid);
                idx = validRows(localIdx);
                rows{end + 1, 1} = table(experiments(iExp), defs(iDef), methods(iMethod), ...
                    sub.h_index(idx), sub.h(idx), values(idx), abs(values(idx) - 0.95), ...
                    sub.active_grid_count(idx), sub.fallback_rate(idx), sub.bootstrap_t_inf_rate(idx), ...
                    'VariableNames', {'experiment', 'active_definition', 'metric', ...
                    'h_index', 'h', 'coverage', 'abs_error_to_095', 'active_grid_count', ...
                    'fallback_rate', 'bootstrap_t_inf_rate'}); %#ok<AGROW>
            end
        end
    end
end
bestH = vertcat(rows{:});
end

function localWriteReport(summaryRoot, resultRoot, config, overall, activeSummary, bestH)
reportPath = fullfile(summaryRoot, "original_plate_summary_report.md");
overallBt = overall(contains(string(overall.method), "Bootstrap-t"), :);
activeLegacy = activeSummary(string(activeSummary.active_definition) == "active_legacy", :);
lines = [
    "# Original Thin-Plate FEM E7/E8 Summary"
    ""
    "Result root: `" + string(resultRoot) + "`"
    ""
    "This is the regenerated original GF/Voronoi/KL2D/20x20 thin-plate FEM workflow, not the standalone diagnostic RFEM approximation used in the h/fallback sweep."
    ""
    "Key parameters:"
    "- run_mode = `" + string(config.run_mode) + "`"
    "- R = " + string(config.R) + ", M = " + string(config.M) + ", B = " + string(config.B)
    "- reference_sample_count = " + string(config.reference_sample_count)
    "- sample_pool_count = " + string(config.sample_pool_count) + ", inner_sample_count = " + string(config.inner_sample_count)
    "- w_grid_count = " + string(config.w_grid_count) + ", h_count = " + string(numel(config.h_list))
    ""
    "Overall Bootstrap-t rows:"
    localFormatOverall(overallBt)
    ""
    "Active-grid legacy means:"
    localFormatActive(activeLegacy)
    ""
    "Best h by active-grid coverage closeness to 0.95:"
    localFormatBest(bestH)
    ""
    "Caution: active-grid definitions are diagnostics. The legacy definition only removes numerically degenerate points; the reference-floor definitions additionally remove near-zero reference-curve locations."
    ];
dpimnumeric.writeText(reportPath, strjoin(lines, newline));
end

function lines = localFormatOverall(tbl)
lines = strings(height(tbl), 1);
for i = 1:height(tbl)
    lines(i) = sprintf("- %s / %s: coverage=%.6f, fallback_rate=%.6f, inf_rate=%.6f", ...
        tbl.experiment(i), tbl.method(i), tbl.mean_coverage(i), ...
        tbl.mean_fallback_rate(i), tbl.mean_bootstrap_t_inf_rate(i));
end
end

function lines = localFormatActive(tbl)
g = groupsummary(tbl, "experiment", "mean", ...
    ["bootstrap_t_coverage", "bootstrap_t_fallback_rule_coverage", ...
    "fallback_rate", "bootstrap_t_inf_rate", "active_grid_count"]);
lines = strings(height(g), 1);
for i = 1:height(g)
    lines(i) = sprintf("- %s: raw BT=%.6f, fallback BT=%.6f, fallback_rate=%.6f, inf_rate=%.6f, active_count_mean=%.1f", ...
        g.experiment(i), g.mean_bootstrap_t_coverage(i), ...
        g.mean_bootstrap_t_fallback_rule_coverage(i), g.mean_fallback_rate(i), ...
        g.mean_bootstrap_t_inf_rate(i), g.mean_active_grid_count(i));
end
end

function lines = localFormatBest(tbl)
if isempty(tbl)
    lines = "- no valid best-h rows";
    return;
end
lines = strings(height(tbl), 1);
for i = 1:height(tbl)
    lines(i) = sprintf("- %s / %s / %s: h=%.16g, coverage=%.6f, abs_error=%.6f", ...
        tbl.experiment(i), tbl.active_definition(i), tbl.metric(i), ...
        tbl.h(i), tbl.coverage(i), tbl.abs_error_to_095(i));
end
end
