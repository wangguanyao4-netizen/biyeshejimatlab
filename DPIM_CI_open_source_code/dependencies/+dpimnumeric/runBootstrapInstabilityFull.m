function artifact = runBootstrapInstabilityFull(projectRoot, outDir, runMode)
%runBootstrapInstabilityFull Run locked-mean targeted bootstrap-t diagnostics.

runMode = lower(string(runMode));
nonlinearDir = localFindNonlinearDir(projectRoot);
baseCfg = localBaseConfig(runMode);

startDir = pwd;
addpath(nonlinearDir);
cleanupObj = onCleanup(@() localRestorePathAndCwd(nonlinearDir, startDir)); %#ok<NASGU>

cases = localCaseConfigs(baseCfg);

aggregateSummary = table();
replacementOverview = table();
coverageRows = table();
artifactRows = table();

for iCase = 1:numel(cases)
    caseCfg = cases(iCase);
    expCfg = localBuildExpCfg(baseCfg, caseCfg);
    runOutputs = run_nonlinear_lognormal_custom_range_spacingh_experiment(expCfg);
    endpointFallback = localApplyEndpointFallbackIfAvailable(runOutputs, expCfg.tag);
    summaryFallback = apply_bootstrapt_summary_yfallback_5x_to_percentile(runOutputs.summary_csv_path, expCfg.tag + "_summaryyfb5xpb");

    rawSummary = readtable(runOutputs.summary_csv_path);
    summaryYFallback = readtable(summaryFallback.adjusted_summary_csv_path);

    aggregateSummary = [aggregateSummary; ...
        localAggregateSummary(rawSummary, caseCfg.case_id, "raw", runOutputs.h_value)]; %#ok<AGROW>
    if endpointFallback.available
        endpointSummary = readtable(endpointFallback.adjusted_summary_csv_path);
        aggregateSummary = [aggregateSummary; ...
            localAggregateSummary(endpointSummary, caseCfg.case_id, "endpoint_fallback", runOutputs.h_value)]; %#ok<AGROW>
    end
    aggregateSummary = [aggregateSummary; ...
        localAggregateSummary(summaryYFallback, caseCfg.case_id, "summary_y_fallback", runOutputs.h_value)]; %#ok<AGROW>

    replacementOverview = [replacementOverview; table( ...
        string(caseCfg.case_id), runOutputs.h_value, ...
        endpointFallback.replacement_count, endpointFallback.replacement_y_point_count, ...
        summaryFallback.replacement_count, summaryFallback.replacement_y_point_count, ...
        caseCfg.y_min, caseCfg.y_max, caseCfg.q_low, caseCfg.q_high, ...
        'VariableNames', {'case_id', 'h', 'endpoint_replacement_rows', 'endpoint_replacement_y_points', ...
        'summary_replacement_rows', 'summary_replacement_y_points', 'y_min', 'y_max', 'q_low', 'q_high'})]; %#ok<AGROW>

    coverageRows = [coverageRows; ...
        localReadCoverage(runOutputs.coverage_paths.Bootstrap_t, caseCfg.case_id, "raw")]; %#ok<AGROW>
    if endpointFallback.available
        coverageRows = [coverageRows; ...
            localReadCoverage(endpointFallback.coverage_paths.Bootstrap_t, caseCfg.case_id, "endpoint_fallback")]; %#ok<AGROW>
    end
    coverageRows = [coverageRows; ...
        localReadCoverage(summaryFallback.coverage_paths.Bootstrap_t, caseCfg.case_id, "summary_y_fallback")]; %#ok<AGROW>

    artifactRows = [artifactRows; table( ...
        string(caseCfg.case_id), string(runOutputs.summary_csv_path), ...
        string(endpointFallback.adjusted_summary_csv_path), string(summaryFallback.adjusted_summary_csv_path), ...
        string(runOutputs.endpoint_csv_path), string(endpointFallback.replacement_log_csv_path), ...
        string(summaryFallback.replacement_log_csv_path), ...
        'VariableNames', {'case_id', 'raw_summary_csv', 'endpoint_fallback_summary_csv', ...
        'summary_y_fallback_csv', 'endpoint_csv', 'endpoint_replacement_csv', ...
        'summary_replacement_csv'})]; %#ok<AGROW>
end

writetable(aggregateSummary, fullfile(outDir, "summary.csv"));
writetable(replacementOverview, fullfile(outDir, "replacement_overview.csv"));
writetable(coverageRows, fullfile(outDir, "bootstrap_t_coverage_compare.csv"));
writetable(artifactRows, fullfile(outDir, "artifact_map.csv"));

raw = struct();
raw.summary = aggregateSummary;
raw.replacement_overview = replacementOverview;
raw.coverage_compare = coverageRows;
raw.artifact_map = artifactRows;
raw.config = struct( ...
    "experiment", "E5_bootstrap_t_instability", ...
    "run_mode", char(runMode), ...
    "seed", 20260430, ...
    "studentize_std_floor_rel", 0.02, ...
    "model_z", "mean(theta)", ...
    "h_rule", "targeted_locked_model_quantile_grid_spacing_dy", ...
    "cases", cases);
save(fullfile(outDir, "raw_results.mat"), "-struct", "raw", "-v7.3");

fig = figure("Visible", "off", "Color", "w");
tiledlayout(numel(cases), 1, "TileSpacing", "compact");
for iCase = 1:numel(cases)
    nexttile;
    mask = coverageRows.case_id == string(cases(iCase).case_id);
    stageOrder = ["raw", "endpoint_fallback", "summary_y_fallback"];
    hold on;
    for iStage = 1:numel(stageOrder)
        stageMask = mask & coverageRows.stage == stageOrder(iStage);
        plot(coverageRows.y_value(stageMask), coverageRows.mean_bootstrap_t_coverage(stageMask), ...
            "LineWidth", 1.2, "DisplayName", stageOrder(iStage));
    end
    hold off;
    title(string(cases(iCase).case_id), "Interpreter", "none");
    xlabel("y");
    ylabel("mean Bootstrap-t coverage");
    grid on;
    legend("Location", "best");
end
try
    exportgraphics(fig, fullfile(outDir, "figures", "bootstrap_t_diagnostics.png"), "Resolution", 180);
catch
    saveas(fig, fullfile(outDir, "figures", "bootstrap_t_diagnostics.png"));
end
close(fig);

artifact = struct();
artifact.config = raw.config;
artifact.dependencies = [ ...
    "run_nonlinear_lognormal_custom_range_spacingh_experiment.m"; ...
    "apply_bootstrapt_fallback_5x_to_percentile.m"; ...
    "apply_bootstrapt_summary_yfallback_5x_to_percentile.m"; ...
    "ci_methods_lognormal_nonlinear.m"];
artifact.note = "Full bootstrap-t instability rerun on targeted d5 locked-model quantile intervals with h=dy and fallback summaries.";
artifact.status = "completed_" + runMode + "_lockedmean_quantile_grid";
end

function baseCfg = localBaseConfig(runMode)
switch runMode
    case "medium"
        baseCfg = struct( ...
            "R", 10, "M", 200, "B", 1800, "N_pool", 1200, ...
            "grid_sample_count", 120000, ...
            "tag_prefix", "codex_lockedmean_d5_targeted_spacingh_B1800_M200_R10_medium");
    otherwise
        baseCfg = struct( ...
            "R", 20, "M", 300, "B", 399, "N_pool", 800, ...
            "grid_sample_count", 200000, ...
            "tag_prefix", "codex_lockedmean_d5_targeted_spacingh_B399_M300_R20_N800_k6_10_full");
end
end

function cases = localCaseConfigs(baseCfg)
targetDimension = 5;
quantiles = localLockedMeanQuantiles(targetDimension, baseCfg.grid_sample_count, [0.50, 0.95, 0.999]);

cases = struct( ...
    'case_id', {"core_q50_q95_400", "righttail_q95_q999_400"}, ...
    'tag', {baseCfg.tag_prefix + "_core_q50_q95_400", baseCfg.tag_prefix + "_righttail_q95_q999_400"}, ...
    'q_low', {0.50, 0.95}, ...
    'q_high', {0.95, 0.999}, ...
    'y_min', {quantiles(1), quantiles(2)}, ...
    'y_max', {quantiles(2), quantiles(3)}, ...
    'point_count', {400, 400}, ...
    'd', {targetDimension, targetDimension}, ...
    'k_list', {6:10, 6:10});
end

function expCfg = localBuildExpCfg(baseCfg, caseCfg)
expCfg = struct();
expCfg.d = caseCfg.d;
expCfg.k_list = caseCfg.k_list;
expCfg.R = baseCfg.R;
expCfg.M = baseCfg.M;
expCfg.B = baseCfg.B;
expCfg.N_pool = baseCfg.N_pool;
expCfg.point_count = caseCfg.point_count;
expCfg.y_min = caseCfg.y_min;
expCfg.y_max = caseCfg.y_max;
expCfg.base_seed = 20260430;
expCfg.run_mode = "full";
expCfg.tag = string(caseCfg.tag);
expCfg.range_source_label = "locked mean(theta) quantiles q=" + string(caseCfg.q_low) + "--" + string(caseCfg.q_high);
expCfg.range_note = "Locked-mean targeted d5 interval; source response cache is aligned with E4 d5 N_pool=800, k=6:10.";
expCfg.save_endpoint_details = true;
expCfg.save_endpoint_csv = false;
expCfg.execution_stage = "endpoint";
expCfg.endpoint_export_mode = "stream_per_r";
expCfg.bootstrap = struct("studentize_std_floor_rel", 0.02);
expCfg.integrand_params = struct("g_vector_mode", "explicit_vector", "g_linear_coeffs", ones(expCfg.d, 1) / expCfg.d);
end

function endpointFallback = localApplyEndpointFallbackIfAvailable(runOutputs, tag)
endpointCsvPath = string(runOutputs.endpoint_csv_path);

endpointFallback = struct();
endpointFallback.available = false;
endpointFallback.adjusted_summary_csv_path = "";
endpointFallback.adjusted_endpoint_csv_path = "";
endpointFallback.replacement_log_csv_path = "";
endpointFallback.coverage_paths = struct();
endpointFallback.replacement_count = NaN;
endpointFallback.replacement_y_point_count = NaN;

if strlength(endpointCsvPath) == 0 || ~isfile(endpointCsvPath)
    fprintf("Endpoint fallback skipped for %s: aggregate endpoint CSV is unavailable; using summary-level fallback only.\n", string(tag));
    return;
end

endpointFallback = apply_bootstrapt_fallback_5x_to_percentile(tag);
endpointFallback.available = true;
end

function q = localLockedMeanQuantiles(d, sampleCount, probs)
seed = 990000 + 101 * d + sampleCount;
previousState = rng;
rng(seed, "twister");
cleanupObj = onCleanup(@() rng(previousState)); %#ok<NASGU>

theta = exp(randn(sampleCount, d));
z = mean(theta, 2);
e = mean(theta.^2, 2) - 1;
response = z + 0.4 * z.^3 + 0.8 * z .* e;
q = localQuantiles(response, probs);
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

function summaryOut = localAggregateSummary(tbl, caseId, stage, hValue)
methods = unique(string(tbl.CI_method), "stable");
summaryOut = table();
for iMethod = 1:numel(methods)
    mask = string(tbl.CI_method) == methods(iMethod);
    block = tbl(mask, :);
    summaryOut = [summaryOut; table( ...
        string(caseId), string(stage), methods(iMethod), hValue, height(block), ...
        mean(block.coverage), min(block.coverage), max(block.coverage), ...
        mean(block.mean_interval_length), median(block.mean_interval_length), max(block.mean_interval_length), ...
        'VariableNames', {'case_id', 'stage', 'CI_method', 'h', 'row_count', ...
        'mean_coverage', 'min_coverage', 'max_coverage', ...
        'mean_interval_length', 'median_interval_length', 'max_interval_length'})]; %#ok<AGROW>
end
end

function coverageTable = localReadCoverage(csvPath, caseId, stage)
tbl = readtable(csvPath);
coverageMatrix = tbl{:, 3:width(tbl)};
coverageTable = table( ...
    repmat(string(caseId), height(tbl), 1), ...
    repmat(string(stage), height(tbl), 1), ...
    tbl.y_index, tbl.y_value, mean(coverageMatrix, 2), ...
    'VariableNames', {'case_id', 'stage', 'y_index', 'y_value', 'mean_bootstrap_t_coverage'});
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
    "export_coverage_by_y_tables.m"
    "apply_bootstrapt_summary_yfallback_5x_to_percentile.m"];

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
