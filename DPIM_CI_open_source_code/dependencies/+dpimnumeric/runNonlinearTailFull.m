function artifact = runNonlinearTailFull(projectRoot, outDir, runMode)
%runNonlinearTailFull Run locked-mean nonlinear custom-grid experiments.

runMode = lower(string(runMode));
nonlinearDir = localFindNonlinearDir(projectRoot);
cases = localCaseConfigs(runMode);

startDir = pwd;
addpath(nonlinearDir);
cleanupObj = onCleanup(@() localRestorePathAndCwd(nonlinearDir, startDir)); %#ok<NASGU>

summaryRows = table();
coverageRows = table();
artifactRows = table();

for iCase = 1:numel(cases)
    caseCfg = cases(iCase);
    gridCfg = localLockedMeanGrid(caseCfg);
    expCfg = localBuildExpCfg(caseCfg, gridCfg);
    runOutputs = run_nonlinear_lognormal_custom_range_spacingh_experiment(expCfg);

    summaryTable = readtable(runOutputs.summary_csv_path);
    summaryTable.source_case = repmat(string(caseCfg.case_id), height(summaryTable), 1);
    summaryTable.model_z = repmat("mean(theta)", height(summaryTable), 1);
    summaryTable.h_rule = repmat("locked_model_quantile_grid_spacing_dy", height(summaryTable), 1);
    summaryTable.grid_q_low = repmat(gridCfg.q_low, height(summaryTable), 1);
    summaryTable.grid_q_high = repmat(gridCfg.q_high, height(summaryTable), 1);
    summaryRows = [summaryRows; summaryTable]; %#ok<AGROW>

    coverageTable = localReadCoverage(runOutputs.coverage_paths.Bootstrap_t, caseCfg.case_id, caseCfg.d);
    coverageRows = [coverageRows; coverageTable]; %#ok<AGROW>

    artifactRows = [artifactRows; table( ...
        string(caseCfg.case_id), string(runOutputs.summary_csv_path), string(runOutputs.y_grid_csv_path), ...
        string(runOutputs.coverage_paths.Bootstrap_t), caseCfg.d, caseCfg.R, caseCfg.M, caseCfg.B, ...
        gridCfg.y_min, gridCfg.y_max, runOutputs.h_value, gridCfg.sample_count, ...
        'VariableNames', {'case_id', 'summary_csv', 'y_grid_csv', 'bootstrap_t_coverage_csv', ...
        'd', 'R', 'M', 'B', 'y_min', 'y_max', 'h', 'grid_pilot_sample_count'})]; %#ok<AGROW>
end

writetable(summaryRows, fullfile(outDir, "summary.csv"));
writetable(coverageRows, fullfile(outDir, "bootstrap_t_coverage_by_y_combined.csv"));
writetable(artifactRows, fullfile(outDir, "artifact_map.csv"));

raw = struct();
raw.summary = summaryRows;
raw.coverage_by_y = coverageRows;
raw.artifact_map = artifactRows;
raw.config = struct( ...
    "experiment", "E4_nonlinear_tail", ...
    "run_mode", char(runMode), ...
    "seed", 20260430, ...
    "studentize_std_floor_rel", 0.02, ...
    "model_z", "mean(theta)", ...
    "h_rule", "locked_model_quantile_grid_spacing_dy", ...
    "grid_quantiles", [0.001, 0.999], ...
    "cases", cases);
save(fullfile(outDir, "raw_results.mat"), "-struct", "raw", "-v7.3");

fig = figure("Visible", "off", "Color", "w");
hold on;
for iCase = 1:numel(cases)
    mask = coverageRows.case_id == string(cases(iCase).case_id);
    plot(coverageRows.y_value(mask), coverageRows.mean_bootstrap_t_coverage(mask), ...
        "LineWidth", 1.4, "DisplayName", string(cases(iCase).case_id));
end
hold off;
xlabel("y");
ylabel("mean Bootstrap-t coverage");
title("Locked-mean nonlinear coverage by y", "Interpreter", "none");
grid on;
legend("Location", "best");
try
    exportgraphics(fig, fullfile(outDir, "figures", "nonlinear_tail_density.png"), "Resolution", 180);
catch
    saveas(fig, fullfile(outDir, "figures", "nonlinear_tail_density.png"));
end
close(fig);

artifact = struct();
artifact.config = raw.config;
artifact.dependencies = [ ...
    "run_nonlinear_lognormal_custom_range_spacingh_experiment.m"; ...
    "nonlinear_lognormal_pilotfixed_build_opts.m"; ...
    "main_rqmc_ci_experiment_lognormal_nonlinear.m"; ...
    "ci_methods_lognormal_nonlinear.m"];
artifact.note = "Full nonlinear tail rerun with h=dy on a locked-model quantile grid and explicit mean(theta) coefficients.";
artifact.status = "completed_" + runMode + "_lockedmean_quantile_grid";
end

function cases = localCaseConfigs(runMode)
switch runMode
    case "medium"
        cases = struct( ...
            'case_id', {"d10_lockedgrid_medium", "d5_lockedgrid_medium"}, ...
            'd', {10, 5}, ...
            'k_list', {6:11, 5:10}, ...
            'R', {10, 10}, ...
            'M', {200, 200}, ...
            'B', {1800, 1800}, ...
            'N_pool', {1200, 1200}, ...
            'point_count', {800, 800}, ...
            'grid_sample_count', {120000, 120000}, ...
            'tag', {"codex_lockedmean_d10_lockedgrid800_spacingh_B1800_M200_R10_medium", ...
                    "codex_lockedmean_d5_lockedgrid800_spacingh_B1800_M200_R10_medium"});
    otherwise
        cases = struct( ...
            'case_id', {"d10_lockedgrid_full", "d5_lockedgrid_full"}, ...
            'd', {10, 5}, ...
            'k_list', {6:10, 6:10}, ...
            'R', {20, 20}, ...
            'M', {400, 300}, ...
            'B', {999, 399}, ...
            'N_pool', {1000, 800}, ...
            'point_count', {800, 800}, ...
            'grid_sample_count', {200000, 200000}, ...
            'tag', {"codex_lockedmean_d10_lockedgrid800_spacingh_B999_M400_R20_N1000_k6_10_full", ...
                    "codex_lockedmean_d5_lockedgrid800_spacingh_B399_M300_R20_N800_k6_10_full"});
end
end

function expCfg = localBuildExpCfg(caseCfg, gridCfg)
expCfg = struct();
expCfg.d = caseCfg.d;
expCfg.k_list = caseCfg.k_list;
expCfg.R = caseCfg.R;
expCfg.M = caseCfg.M;
expCfg.B = caseCfg.B;
expCfg.N_pool = caseCfg.N_pool;
expCfg.point_count = caseCfg.point_count;
expCfg.y_min = gridCfg.y_min;
expCfg.y_max = gridCfg.y_max;
expCfg.base_seed = 20260430;
expCfg.run_mode = "full";
expCfg.tag = string(caseCfg.tag);
expCfg.range_source_label = sprintf("locked mean(theta) pilot quantiles %.4g--%.4g with %d samples", ...
    gridCfg.q_low, gridCfg.q_high, gridCfg.sample_count);
expCfg.range_note = "Custom y-range built from the locked nonlinear model; kernel sigma equals the custom y-grid spacing.";
expCfg.save_endpoint_details = false;
expCfg.save_endpoint_csv = false;
expCfg.execution_stage = "scan";
expCfg.endpoint_export_mode = "aggregate";
expCfg.bootstrap = struct("studentize_std_floor_rel", 0.02);
expCfg.integrand_params = struct( ...
    "g_vector_mode", "explicit_vector", ...
    "g_linear_coeffs", ones(caseCfg.d, 1) / caseCfg.d);
end

function gridCfg = localLockedMeanGrid(caseCfg)
qLow = 0.001;
qHigh = 0.999;
sampleCount = caseCfg.grid_sample_count;
seed = 880000 + 101 * caseCfg.d + sampleCount;

previousState = rng;
rng(seed, "twister");
cleanupObj = onCleanup(@() rng(previousState)); %#ok<NASGU>

theta = exp(randn(sampleCount, caseCfg.d));
z = mean(theta, 2);
e = mean(theta.^2, 2) - 1;
response = z + 0.4 * z.^3 + 0.8 * z .* e;
q = localQuantiles(response, [qLow, qHigh]);
width = max(q(2) - q(1), eps);

gridCfg = struct();
gridCfg.q_low = qLow;
gridCfg.q_high = qHigh;
gridCfg.sample_count = sampleCount;
gridCfg.y_min = max(0, q(1) - 0.02 * width);
gridCfg.y_max = q(2) + 0.02 * width;
gridCfg.seed = seed;
end

function q = localQuantiles(x, probs)
x = sort(x(:));
n = numel(x);
q = zeros(size(probs));
for i = 1:numel(probs)
    pos = 1 + (n - 1) * probs(i);
    lo = floor(pos);
    hi = ceil(pos);
    if lo == hi
        q(i) = x(lo);
    else
        q(i) = x(lo) + (pos - lo) * (x(hi) - x(lo));
    end
end
end

function coverageTable = localReadCoverage(csvPath, caseId, dValue)
tbl = readtable(csvPath);
coverageMatrix = tbl{:, 3:width(tbl)};
meanCoverage = mean(coverageMatrix, 2);
minCoverage = min(coverageMatrix, [], 2);
maxCoverage = max(coverageMatrix, [], 2);
coverageTable = table( ...
    repmat(string(caseId), height(tbl), 1), ...
    repmat(dValue, height(tbl), 1), ...
    tbl.y_index, tbl.y_value, meanCoverage, minCoverage, maxCoverage, ...
    'VariableNames', {'case_id', 'd', 'y_index', 'y_value', ...
    'mean_bootstrap_t_coverage', 'min_bootstrap_t_coverage', 'max_bootstrap_t_coverage'});
end

function nonlinearDir = localFindNonlinearDir(projectRoot)
matches = dir(fullfile(projectRoot, "**", "run_nonlinear_lognormal_custom_range_spacingh_experiment.m"));
if isempty(matches)
    error("Cannot locate run_nonlinear_lognormal_custom_range_spacingh_experiment.m under %s.", projectRoot);
end

requiredHelpers = [
    "nonlinear_lognormal_pilotfixed_build_opts.m"
    "main_rqmc_ci_experiment_lognormal_nonlinear.m"
    "ci_methods_lognormal_nonlinear.m"
    "export_coverage_by_y_tables.m"];

for iMatch = 1:numel(matches)
    candidateDir = string(matches(iMatch).folder);
    hasAllHelpers = true;
    for iHelper = 1:numel(requiredHelpers)
        if ~isfile(fullfile(candidateDir, requiredHelpers(iHelper)))
            hasAllHelpers = false;
            break;
        end
    end
    if hasAllHelpers
        nonlinearDir = char(candidateDir);
        return;
    end
end

error("Located nonlinear wrappers under %s, but none contains all required pilotfixed helper files.", projectRoot);
end

function localRestorePathAndCwd(addedPath, startDir)
if strlength(string(addedPath)) > 0 && contains(path, char(addedPath))
    rmpath(addedPath);
end
if isfolder(startDir)
    cd(startDir);
end
end
