function outputs = runOriginalPlateFullSuite(projectRoot, resultsRoot, runMode, startAt)
%runOriginalPlateFullSuite Source-driven original thin-plate FEM E7/E8 suite.
%
% The older E7/E8 wrappers only read an existing curve pool. This suite
% rebuilds the curve pool from the original project scripts in an isolated
% output directory, preserving the raw FE responses and probability weights.

if nargin < 3 || strlength(string(runMode)) == 0
    runMode = "smoke";
end
if nargin < 4 || strlength(string(startAt)) == 0
    startAt = "";
end

runMode = lower(string(runMode));
startAt = string(startAt);
if ~any(runMode == ["smoke", "medium", "full", "journal"])
    error("runMode must be 'smoke', 'medium', 'full', or 'journal'.");
end

config = localConfig(runMode);
dpimnumeric.ensureDir(resultsRoot);
dpimnumeric.ensureDir(fullfile(resultsRoot, "reference"));
dpimnumeric.ensureDir(fullfile(resultsRoot, "rqmc_pool"));
dpimnumeric.ensureDir(fullfile(resultsRoot, "E7_plate_SFEM"));
dpimnumeric.ensureDir(fullfile(resultsRoot, "E8_simultaneous_band"));

manifest = struct();
manifest.project_root = projectRoot;
manifest.results_root = resultsRoot;
manifest.run_mode = char(runMode);
manifest.start_at = char(startAt);
manifest.created_at = char(datetime("now"));
manifest.note = "Original thin-plate GF/KL/FEM rebuild; E6 cancelled; no manuscript edits.";
dpimnumeric.writeJson(fullfile(resultsRoot, "original_plate_manifest.json"), manifest);
dpimnumeric.writeJson(fullfile(resultsRoot, "config.json"), config);

oldPath = path;
oldDir = pwd;
cleanupObj = onCleanup(@() localRestore(oldPath, oldDir)); %#ok<NASGU>
addpath(projectRoot);

steps = ["reference"; "rqmc_pool"; "E7_plate_SFEM"; "E8_simultaneous_band"];
if strlength(startAt) > 0
    idx = find(steps == startAt, 1);
    if isempty(idx)
        error("Unknown startAt value %s.", startAt);
    end
    steps = steps(idx:end);
end

outputs = struct();
if any(steps == "reference")
    outputs.reference = localBuildReference(projectRoot, resultsRoot, config);
else
    outputs.reference = fullfile(resultsRoot, "reference", "thin_reference.mat");
end

if any(steps == "rqmc_pool")
    outputs.rqmc_pool = localBuildRqmcPool(projectRoot, resultsRoot, config);
else
    outputs.rqmc_pool = fullfile(resultsRoot, "rqmc_pool", "thin_rqmc_pool.mat");
end

inputs = localLoadInputs(resultsRoot);
if any(steps == "E7_plate_SFEM")
    outputs.E7_plate_SFEM = localRunE7(resultsRoot, config, inputs);
end
if any(steps == "E8_simultaneous_band")
    outputs.E8_simultaneous_band = localRunE8(resultsRoot, config, inputs);
end

localWriteExperimentStatus(resultsRoot, config);
end

function config = localConfig(runMode)
dy = 1010 / 799;
hBase = logspace(log10(0.005), log10(5.0), 18);
hList = sort([hBase(:); 1; dy]);

switch runMode
    case "journal"
        config = struct("run_mode", "journal", "alpha", 0.05, "lambda", 5, ...
            "R", 20, "M", 1000, "B", 399, "seed", 20260605 + 810, ...
            "dimension", 10, "mesh_nx", 20, "mesh_ny", 20, ...
            "sample_pool_count", 600, "inner_sample_count", 200, ...
            "reference_sample_count", 10000, "w_grid_count", 800, ...
            "gf_aux_sample_count", 50000, "gf_block_size", 2000, ...
            "voronoi_aux_sample_count", 50000, "voronoi_train_aux_count", 5000, ...
            "voronoi_exact_block_size", 1000, "h_list", hList(:).');
    case "full"
        config = struct("run_mode", "full", "alpha", 0.05, "lambda", 5, ...
            "R", 20, "M", 300, "B", 399, "seed", 20260605 + 810, ...
            "dimension", 10, "mesh_nx", 20, "mesh_ny", 20, ...
            "sample_pool_count", 200, "inner_sample_count", 200, ...
            "reference_sample_count", 10000, "w_grid_count", 800, ...
            "gf_aux_sample_count", 50000, "gf_block_size", 2000, ...
            "voronoi_aux_sample_count", 50000, "voronoi_train_aux_count", 5000, ...
            "voronoi_exact_block_size", 1000, "h_list", hList(:).');
    case "medium"
        config = struct("run_mode", "medium", "alpha", 0.05, "lambda", 5, ...
            "R", 10, "M", 120, "B", 199, "seed", 20260605 + 810, ...
            "dimension", 10, "mesh_nx", 20, "mesh_ny", 20, ...
            "sample_pool_count", 40, "inner_sample_count", 80, ...
            "reference_sample_count", 1000, "w_grid_count", 160, ...
            "gf_aux_sample_count", 8000, "gf_block_size", 1000, ...
            "voronoi_aux_sample_count", 8000, "voronoi_train_aux_count", 1500, ...
            "voronoi_exact_block_size", 1000, "h_list", hList(:).');
    otherwise
        config = struct("run_mode", "smoke", "alpha", 0.05, "lambda", 5, ...
            "R", 4, "M", 12, "B", 49, "seed", 20260605 + 810, ...
            "dimension", 10, "mesh_nx", 20, "mesh_ny", 20, ...
            "sample_pool_count", 8, "inner_sample_count", 16, ...
            "reference_sample_count", 80, "w_grid_count", 40, ...
            "gf_aux_sample_count", 1000, "gf_block_size", 500, ...
            "voronoi_aux_sample_count", 1000, "voronoi_train_aux_count", 200, ...
            "voronoi_exact_block_size", 500, "h_list", hList(:).');
end

config.confidence_level = 1 - config.alpha;
config.sigma_reference = 1;
config.paper_text_sigma_dy = dy;
config.w_list = linspace(-10, 1000, config.w_grid_count);
config.source_workflow = "GF -> Voronoi -> KL2D -> 20x20 Morley/Kirchhoff thin plate FEM";
config.scope_note = "E6 cancelled; E7/E8 only; this does not modify manuscript files.";
end

function referencePath = localBuildReference(projectRoot, resultsRoot, config)
outDir = fullfile(resultsRoot, "reference");
dpimnumeric.ensureDir(outDir);
cd(outDir);

localEnsureEvalCenters("eval_centers_mc.txt", config.mesh_nx, config.mesh_ny);
mc_pipeline_config = localPipelineConfig(config, config.reference_sample_count, config.seed + 11);
save("mc_pipeline_config.mat", "mc_pipeline_config");

logPath = fullfile(outDir, "reference_build.log");
diary(logPath);
fprintf("Original thin-plate reference build started: %s\n", char(datetime("now")));
fprintf("projectRoot = %s\n", projectRoot);
fprintf("sample_count = %d\n", config.reference_sample_count);

evalin("base", "gf_normal_final_mc;");
evalin("base", "voronoi_pca_toolbox_final_blockwise_mc;");
evalin("base", "KL2D_final_complete_3D_mc;");
evalin("base", "jisuanpro_mc_quad_paperstyle_cloud;");

wcData = load("wc_all_mc.mat", "wc_all_mc");
PV_mc = load("pexact_mc.txt");
reference = struct();
reference.sample_count = config.reference_sample_count;
reference.w_list = config.w_list;
reference.wc_all_mc_thin = wcData.wc_all_mc;
reference.PV_mc = PV_mc;
reference.h_list = config.h_list;
reference.curves_by_h = localCurvesForH(reference.wc_all_mc_thin, reference.PV_mc, config.h_list, config.w_list);
reference.source = config.source_workflow;

referencePath = fullfile(outDir, "thin_reference.mat");
save(referencePath, "reference", "config", "-v7.3");
writematrix([config.w_list(:), reference.curves_by_h(:, 1)], ...
    fullfile(outDir, "thin_reference_first_h_curve.txt"), "Delimiter", "tab");
fprintf("Reference saved: %s\n", referencePath);
diary off;

if exist("mc_pipeline_config.mat", "file") == 2
    delete("mc_pipeline_config.mat");
end
end

function poolPath = localBuildRqmcPool(projectRoot, resultsRoot, config)
outDir = fullfile(resultsRoot, "rqmc_pool");
dpimnumeric.ensureDir(outDir);
cd(outDir);

localEnsureEvalCenters("eval_centers_rqmc.txt", config.mesh_nx, config.mesh_ny);
curveCount = config.sample_pool_count;
wcRuns = cell(curveCount, 1);
weightRuns = cell(curveCount, 1);
seedList = config.seed + 1000 + (0:(curveCount - 1));

logPath = fullfile(outDir, "rqmc_pool_build.log");
diary(logPath);
fprintf("Original thin-plate RQMC pool build started: %s\n", char(datetime("now")));
fprintf("projectRoot = %s\n", projectRoot);
fprintf("outer samples = %d, inner sample count = %d\n", curveCount, config.inner_sample_count);

for iRun = 1:curveCount
    fprintf("RQMC outer sample %d / %d\n", iRun, curveCount);
    runPath = fullfile(outDir, sprintf("rqmc_outer_sample_thin_original_%03d.mat", iRun));
    if isfile(runPath)
        existing = load(runPath, "wc_all_rqmc_thin", "PV_rqmc", "seed");
        if isfield(existing, "seed") && existing.seed == seedList(iRun)
            wcRuns{iRun} = existing.wc_all_rqmc_thin;
            weightRuns{iRun} = existing.PV_rqmc;
            fprintf("Reused completed outer sample %d.\n", iRun);
            continue;
        end
    end
    rqmc_pipeline_config = localPipelineConfig(config, config.inner_sample_count, seedList(iRun));
    save("rqmc_pipeline_config.mat", "rqmc_pipeline_config");

    evalin("base", "gf_normal_final_rqmc;");
    evalin("base", "voronoi_pca_toolbox_final_blockwise_rqmc;");
    evalin("base", "KL2D_final_complete_3D_rqmc;");
    evalin("base", "jisuanpro_rqmc_quad_paperstyle_cloud;");

    wcData = load("wc_all_rqmc.mat", "wc_all_rqmc");
    PV_rqmc = load("pexact_rqmc.txt");
    wcRuns{iRun} = wcData.wc_all_rqmc;
    weightRuns{iRun} = PV_rqmc;

    single = struct();
    single.run_id = iRun;
    single.seed = seedList(iRun);
    single.w_list = config.w_list;
    single.wc_all_rqmc_thin = wcRuns{iRun};
    single.PV_rqmc = weightRuns{iRun};
    save(runPath, "-struct", "single");
end

pool = struct();
pool.sample_pool_count = curveCount;
pool.inner_sample_count = config.inner_sample_count;
pool.seed_list = seedList;
pool.w_list = config.w_list;
pool.h_list = config.h_list;
pool.wc_runs_thin = wcRuns;
pool.weight_runs = weightRuns;
pool.source = config.source_workflow;

poolPath = fullfile(outDir, "thin_rqmc_pool.mat");
save(poolPath, "pool", "config", "-v7.3");
fprintf("RQMCPool saved: %s\n", poolPath);
diary off;

if exist("rqmc_pipeline_config.mat", "file") == 2
    delete("rqmc_pipeline_config.mat");
end
end

function result = localRunE7(resultsRoot, config, inputs)
outDir = fullfile(resultsRoot, "E7_plate_SFEM");
dpimnumeric.ensureDir(outDir);
[summary, pointwise, diagnosticBlocks] = localRunPointwiseSweep(config, inputs);
writetable(summary, fullfile(outDir, "summary.csv"));
writetable(pointwise, fullfile(outDir, "pointwise_by_h.csv"));
save(fullfile(outDir, "curve_coverage_diagnostics.mat"), "diagnosticBlocks", "config", "-v7.3");
localWriteReadme(outDir, "Original thin-plate FEM pointwise CI coverage over regenerated RQMC curve pool.", config);
result = localResultPaths(outDir);
end

function result = localRunE8(resultsRoot, config, inputs)
outDir = fullfile(resultsRoot, "E8_simultaneous_band");
dpimnumeric.ensureDir(outDir);
[summary, pointwise, diagnosticBlocks] = localRunBandSweep(config, inputs);
writetable(summary, fullfile(outDir, "summary.csv"));
writetable(pointwise, fullfile(outDir, "pointwise_band_by_h.csv"));
save(fullfile(outDir, "curve_coverage_diagnostics.mat"), "diagnosticBlocks", "config", "-v7.3");
localWriteReadme(outDir, "Original thin-plate FEM finite-grid simultaneous-band coverage over regenerated RQMC curve pool.", config);
result = localResultPaths(outDir);
end

function inputs = localLoadInputs(resultsRoot)
refData = load(fullfile(resultsRoot, "reference", "thin_reference.mat"), "reference");
poolData = load(fullfile(resultsRoot, "rqmc_pool", "thin_rqmc_pool.mat"), "pool");
inputs = struct();
inputs.reference = refData.reference;
inputs.pool = poolData.pool;
inputs.h_list = refData.reference.h_list;
inputs.w_list = refData.reference.w_list;
end

function [summaryAll, pointwiseAll, diagnosticBlocks] = localRunPointwiseSweep(config, inputs)
hList = inputs.h_list(:).';
summaryAll = table();
pointwiseAll = table();
diagnosticBlocks = cell(numel(hList), 1);
for ih = 1:numel(hList)
    [referenceCurve, curvePool] = localCurvePoolAtH(inputs, hList(ih));
    [summary, pointwise, diagnostics] = localCurvePoolCoverageWithFallback( ...
        curvePool, referenceCurve, inputs.w_list, config.alpha, config.B, config.M, config.R, ...
        config.seed + 30000 + ih * 100, config.lambda);
    summary.h_index = repmat(ih, height(summary), 1);
    summary.h = repmat(hList(ih), height(summary), 1);
    pointwise.h_index = repmat(ih, height(pointwise), 1);
    pointwise.h = repmat(hList(ih), height(pointwise), 1);
    summaryAll = [summaryAll; summary]; %#ok<AGROW>
    pointwiseAll = [pointwiseAll; pointwise]; %#ok<AGROW>
    diagnosticBlocks{ih} = diagnostics;
end
summaryAll = movevars(summaryAll, ["h_index", "h"], "Before", 1);
pointwiseAll = movevars(pointwiseAll, ["h_index", "h"], "Before", 1);
end

function [summaryAll, pointwiseAll, diagnosticBlocks] = localRunBandSweep(config, inputs)
hList = inputs.h_list(:).';
summaryAll = table();
pointwiseAll = table();
diagnosticBlocks = cell(numel(hList), 1);
for ih = 1:numel(hList)
    [referenceCurve, curvePool] = localCurvePoolAtH(inputs, hList(ih));
    [summary, pointwise, diagnostics] = localSimultaneousCoverageWithFallback( ...
        curvePool, referenceCurve, inputs.w_list, config.alpha, config.B, config.M, config.R, ...
        config.seed + 50000 + ih * 100, config.lambda);
    summary.h_index = repmat(ih, height(summary), 1);
    summary.h = repmat(hList(ih), height(summary), 1);
    pointwise.h_index = repmat(ih, height(pointwise), 1);
    pointwise.h = repmat(hList(ih), height(pointwise), 1);
    summaryAll = [summaryAll; summary]; %#ok<AGROW>
    pointwiseAll = [pointwiseAll; pointwise]; %#ok<AGROW>
    diagnosticBlocks{ih} = diagnostics;
end
summaryAll = movevars(summaryAll, ["h_index", "h"], "Before", 1);
pointwiseAll = movevars(pointwiseAll, ["h_index", "h"], "Before", 1);
end

function [referenceCurve, curvePool] = localCurvePoolAtH(inputs, h)
hIdx = find(abs(inputs.h_list - h) <= max(1e-12, 1e-12 * abs(h)), 1);
if isempty(hIdx)
    error("h %.16g is not present in inputs.h_list.", h);
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

function curves = localCurvesForH(wcAll, weights, hList, wList)
curves = zeros(numel(wList), numel(hList));
for ih = 1:numel(hList)
    curves(:, ih) = compute_probability_estimate_curve(wcAll, weights, hList(ih), wList).';
end
end

function cfg = localPipelineConfig(config, sampleCount, seed)
cfg = struct();
cfg.sample_count = sampleCount;
cfg.confidence_level = config.confidence_level;
cfg.sigma = config.sigma_reference;
cfg.w_list = config.w_list;
cfg.rng_seed = seed;
cfg.gf_aux_sample_count = config.gf_aux_sample_count;
cfg.gf_block_size = config.gf_block_size;
cfg.voronoi_aux_sample_count = config.voronoi_aux_sample_count;
cfg.voronoi_train_aux_count = config.voronoi_train_aux_count;
cfg.voronoi_exact_block_size = config.voronoi_exact_block_size;
cfg.save_figures = false;
cfg.sample_fig_on = false;
cfg.draw_figures = false;
cfg.save_png_figures = false;
cfg.plot_gf_figures = false;
end

function localEnsureEvalCenters(fileName, nx, ny)
a = 1.0;
b = 1.0;
alpha = pi / 2;
nnode = (nx + 1) * (ny + 1);
nodes = zeros(nnode, 2);
id = 0;
for j = 0:ny
    eta = j / ny;
    for i = 0:nx
        xi = i / nx;
        x = a * xi + b * eta * cos(alpha);
        y = b * eta * sin(alpha);
        id = id + 1;
        nodes(id, :) = [x, y];
    end
end
writematrix(nodes, fileName, "Delimiter", "tab");
end

function [summary, pointwise, diagnostics] = localCurvePoolCoverageWithFallback(curvePool, referenceCurve, wGrid, alpha, B, M, R, seed, lambda)
nGrid = numel(wGrid);
coverageHits = zeros(4, nGrid);
lengthSums = zeros(4, nGrid);
fallbackCounts = zeros(4, nGrid);
btInfCounts = zeros(4, nGrid);
diagnostics = localEmptyCurveDiagnostics(M, R, nGrid, localScalarMethods());
previousState = rng;
rng(double(seed), "twister");
cleanupObj = onCleanup(@() rng(previousState)); %#ok<NASGU>
for trial = 1:M
    idx = randperm(size(curvePool, 1), R);
    diagnostics.sample_indices(trial, :) = idx;
    ci = localCiMethodsMatrix(curvePool(idx, :), alpha, B, seed + trial, referenceCurve, lambda);
    for iMethod = 1:4
        coverageHits(iMethod, :) = coverageHits(iMethod, :) + ci(iMethod).contains;
        finiteLength = isfinite(ci(iMethod).length);
        lengthSums(iMethod, finiteLength) = lengthSums(iMethod, finiteLength) + ci(iMethod).length(finiteLength);
        fallbackCounts(iMethod, :) = fallbackCounts(iMethod, :) + ci(iMethod).fallback_trigger;
        btInfCounts(iMethod, :) = btInfCounts(iMethod, :) + ci(iMethod).source_bootstrap_t_infinite;
        diagnostics.lower(trial, iMethod, :) = reshape(ci(iMethod).lower, 1, 1, nGrid);
        diagnostics.upper(trial, iMethod, :) = reshape(ci(iMethod).upper, 1, 1, nGrid);
        diagnostics.length(trial, iMethod, :) = reshape(ci(iMethod).length, 1, 1, nGrid);
        diagnostics.contains(trial, iMethod, :) = reshape(ci(iMethod).contains, 1, 1, nGrid);
        diagnostics.fallback_trigger(trial, iMethod, :) = reshape(ci(iMethod).fallback_trigger, 1, 1, nGrid);
        diagnostics.source_bootstrap_t_infinite(trial, iMethod, :) = reshape(ci(iMethod).source_bootstrap_t_infinite, 1, 1, nGrid);
    end
end
method = localScalarMethods();
coverage = mean(coverageHits / M, 2);
meanLength = mean(lengthSums / M, 2, "omitnan");
fallbackCount = sum(fallbackCounts, 2);
btInfCount = sum(btInfCounts, 2);
nominal = 1 - alpha;
coverageObservationCount = M * nGrid;
coverageMcSe = sqrt(max(coverage .* (1 - coverage), 0) / coverageObservationCount);
summary = table(method, coverage, repmat(nominal, 4, 1), coverage - nominal, ...
    abs(coverage - nominal), coverageMcSe, 1.96 * coverageMcSe, meanLength, ...
    fallbackCount, fallbackCount / coverageObservationCount, ...
    btInfCount, btInfCount / coverageObservationCount, repmat(nGrid, 4, 1), ...
    repmat(coverageObservationCount, 4, 1), ...
    'VariableNames', {'method', 'coverage', 'nominal_coverage', ...
    'coverage_error', 'abs_coverage_error', 'coverage_mc_se', 'coverage_ci_half_width_95', ...
    'mean_interval_length', 'fallback_count', 'fallback_rate', ...
    'bootstrap_t_inf_count', 'bootstrap_t_inf_rate', 'grid_count', 'coverage_observation_count'});
pointwise = table(wGrid(:), referenceCurve(:), ...
    (coverageHits(3, :)' / M), (coverageHits(4, :)' / M), ...
    fallbackCounts(4, :)', btInfCounts(4, :)', ...
    'VariableNames', {'w', 'reference_curve', 'bootstrap_t_coverage', ...
    'bootstrap_t_fallback_rule_coverage', 'fallback_count', 'bootstrap_t_inf_count'});
diagnostics.reference_curve = referenceCurve;
diagnostics.w_grid = wGrid;
diagnostics.coverage_hits = coverageHits;
diagnostics.length_sums = lengthSums;
diagnostics.fallback_counts = fallbackCounts;
diagnostics.bootstrap_t_inf_counts = btInfCounts;
end

function [summary, pointwise, diagnostics] = localSimultaneousCoverageWithFallback(curvePool, referenceCurve, wGrid, alpha, B, M, R, seed, lambda)
[pointSummary, pointwise, pointDiagnostics] = localCurvePoolCoverageWithFallback(curvePool, referenceCurve, wGrid, alpha, B, M, R, seed, lambda);
nGrid = numel(wGrid);
fullHits = 0;
bandPointHits = zeros(1, nGrid);
bandLengthSum = zeros(1, nGrid);
bandInfCount = 0;
bandDiagnostics = struct();
bandDiagnostics.sample_indices = zeros(M, R);
bandDiagnostics.critical_value = nan(M, 1);
bandDiagnostics.full_band_hit = false(M, 1);
bandDiagnostics.band_point_hit = false(M, nGrid);
bandDiagnostics.lower = nan(M, nGrid);
bandDiagnostics.upper = nan(M, nGrid);
bandDiagnostics.length = nan(M, nGrid);
bandDiagnostics.bootstrap_max_stats = nan(M, B);
previousState = rng;
rng(double(seed + 900000), "twister");
cleanupObj = onCleanup(@() rng(previousState)); %#ok<NASGU>
for trial = 1:M
    idx = randperm(size(curvePool, 1), R);
    bandDiagnostics.sample_indices(trial, :) = idx;
    Y = curvePool(idx, :);
    mu = mean(Y, 1);
    sd = std(Y, 0, 1);
    se = sd / sqrt(R);
    bootIdx = randi(R, R, B);
    maxStats = zeros(B, 1);
    statGrid = se > 100 * eps(max(abs(mu), abs(referenceCurve)));
    for b = 1:B
        bootY = Y(bootIdx(:, b), :);
        bootMean = mean(bootY, 1);
        bootSe = std(bootY, 0, 1) / sqrt(R);
        if any(statGrid)
            maxStats(b) = max(abs((bootMean(statGrid) - mu(statGrid)) ./ bootSe(statGrid)));
        else
            maxStats(b) = 0;
        end
    end
    maxStats(~isfinite(maxStats)) = Inf;
    bandDiagnostics.bootstrap_max_stats(trial, :) = maxStats(:).';
    crit = sort(maxStats);
    crit = crit(max(1, ceil((1 - alpha) * B)));
    lower = mu - crit * se;
    upper = mu + crit * se;
    bandInfCount = bandInfCount + double(~isfinite(crit) || any(~isfinite(lower) | ~isfinite(upper)));
    hit = (lower <= referenceCurve) & (referenceCurve <= upper);
    fullHits = fullHits + all(hit);
    bandPointHits = bandPointHits + hit;
    bandLengthSum = bandLengthSum + (upper - lower);
    bandDiagnostics.critical_value(trial) = crit;
    bandDiagnostics.full_band_hit(trial) = all(hit);
    bandDiagnostics.band_point_hit(trial, :) = hit;
    bandDiagnostics.lower(trial, :) = lower;
    bandDiagnostics.upper(trial, :) = upper;
    bandDiagnostics.length(trial, :) = upper - lower;
end
bandCoverage = fullHits / M;
bandMcSe = sqrt(max(bandCoverage * (1 - bandCoverage), 0) / M);
bandRow = table("Bootstrap-t simultaneous-band", bandCoverage, 1 - alpha, ...
    bandCoverage - (1 - alpha), abs(bandCoverage - (1 - alpha)), bandMcSe, ...
    1.96 * bandMcSe, mean(bandLengthSum / M, "omitnan"), 0, 0, ...
    bandInfCount, bandInfCount / M, nGrid, M, ...
    'VariableNames', {'method', 'coverage', 'nominal_coverage', ...
    'coverage_error', 'abs_coverage_error', 'coverage_mc_se', 'coverage_ci_half_width_95', ...
    'mean_interval_length', 'fallback_count', 'fallback_rate', ...
    'bootstrap_t_inf_count', 'bootstrap_t_inf_rate', 'grid_count', 'coverage_observation_count'});
summary = [pointSummary; bandRow];
pointwise.simultaneous_band_coverage = bandPointHits(:) / M;
pointwise.simultaneous_band_mean_length = bandLengthSum(:) / M;
diagnostics = struct();
diagnostics.pointwise = pointDiagnostics;
diagnostics.band = bandDiagnostics;
end

function ci = localCiMethodsMatrix(Y, alpha, B, seed, truth, lambda)
previousState = rng;
rng(double(seed), "twister");
cleanupObj = onCleanup(@() rng(previousState)); %#ok<NASGU>
[R, nGrid] = size(Y);
mu = mean(Y, 1);
sd = std(Y, 0, 1);
se = sd / sqrt(R);
tCrit = localStudentTCrit(1 - alpha / 2, R - 1);
ci = repmat(localEmptyMatrixCi(nGrid), 4, 1);
studentLower = mu - tCrit * se;
studentUpper = mu + tCrit * se;
ci(1) = localPackMatrixCi("Student-t", studentLower, studentUpper, truth, false(1, nGrid), false(1, nGrid));
bootIdx = randi(R, R, B);
bootMeans = zeros(B, nGrid);
bootStd = zeros(B, nGrid);
for b = 1:B
    bootY = Y(bootIdx(:, b), :);
    bootMeans(b, :) = mean(bootY, 1);
    bootStd(b, :) = std(bootY, 0, 1);
end
[pctLower, pctUpper] = localQuantileRows(bootMeans, alpha);
ci(2) = localPackMatrixCi("Percentile bootstrap", pctLower, pctUpper, truth, false(1, nGrid), false(1, nGrid));
Tstar = sqrt(R) * (bootMeans - mu) ./ bootStd;
zeroStd = bootStd == 0;
if any(zeroStd, "all")
    delta = bootMeans - mu;
    signDelta = sign(delta);
    signDelta(signDelta == 0) = 1;
    Tstar(zeroStd) = signDelta(zeroStd) .* Inf;
end
[tLower, tUpper] = localQuantileRows(Tstar, alpha);
btLower = mu - sd .* tUpper / sqrt(R);
btUpper = mu - sd .* tLower / sqrt(R);
ci(3) = localPackMatrixCi("Bootstrap-t", btLower, btUpper, truth, false(1, nGrid), false(1, nGrid));
ci(3).source_bootstrap_t_infinite = ci(3).infinite;
trigger = ci(3).infinite | (ci(3).length > lambda * max(ci(2).length, realmin));
fbLower = btLower;
fbUpper = btUpper;
fbLower(trigger) = pctLower(trigger);
fbUpper(trigger) = pctUpper(trigger);
ci(4) = localPackMatrixCi("Bootstrap-t fallback-rule", fbLower, fbUpper, truth, trigger, ci(3).infinite);
end

function ci = localEmptyMatrixCi(nGrid)
ci = struct("name", "", "lower", nan(1, nGrid), "upper", nan(1, nGrid), ...
    "contains", false(1, nGrid), "length", nan(1, nGrid), ...
    "infinite", false(1, nGrid), "fallback_trigger", false(1, nGrid), ...
    "source_bootstrap_t_infinite", false(1, nGrid));
end

function ci = localPackMatrixCi(name, lower, upper, truth, fallbackTrigger, sourceBtInf)
swap = lower > upper;
tmp = lower(swap);
lower(swap) = upper(swap);
upper(swap) = tmp;
ci = localEmptyMatrixCi(numel(lower));
ci.name = string(name);
ci.lower = lower;
ci.upper = upper;
ci.infinite = ~(isfinite(lower) & isfinite(upper));
ci.length = upper - lower;
ci.length(ci.infinite) = Inf;
ci.contains = (lower <= truth) & (truth <= upper);
ci.fallback_trigger = fallbackTrigger;
ci.source_bootstrap_t_infinite = sourceBtInf;
end

function [qLower, qUpper] = localQuantileRows(x, alpha)
x = sort(x, 1);
n = size(x, 1);
qLower = x(max(1, floor(n * alpha / 2)), :);
qUpper = x(min(n, ceil(n * (1 - alpha / 2))), :);
end

function q = localStudentTCrit(p, nu)
tailProb = 2 * min(p, 1 - p);
x = betaincinv(tailProb, nu / 2, 0.5);
q = abs(sqrt(nu * (1 / x - 1)));
end

function methods = localScalarMethods()
methods = ["Student-t"; "Percentile bootstrap"; "Bootstrap-t"; "Bootstrap-t fallback-rule"];
end

function diagnostics = localEmptyCurveDiagnostics(M, R, nGrid, methods)
nMethods = numel(methods);
diagnostics = struct();
diagnostics.methods = methods;
diagnostics.sample_indices = zeros(M, R);
diagnostics.lower = nan(M, nMethods, nGrid);
diagnostics.upper = nan(M, nMethods, nGrid);
diagnostics.length = nan(M, nMethods, nGrid);
diagnostics.contains = false(M, nMethods, nGrid);
diagnostics.fallback_trigger = false(M, nMethods, nGrid);
diagnostics.source_bootstrap_t_infinite = false(M, nMethods, nGrid);
end

function localWriteReadme(outDir, description, config)
lines = [
    "# Original Plate FEM Experiment"
    ""
    string(description)
    ""
    "This run rebuilds the thin-plate curve pool from original GF/KL/FEM scripts."
    "E6 is cancelled and is not part of this suite."
    "Raw Bootstrap-t and fallback-rule Bootstrap-t are reported separately."
    ""
    "Key parameters:"
    "- R = " + config.R
    "- M = " + config.M
    "- B = " + config.B
    "- sample_pool_count = " + config.sample_pool_count
    "- inner_sample_count = " + config.inner_sample_count
    "- reference_sample_count = " + config.reference_sample_count
    "- w_grid_count = " + config.w_grid_count
    "- h_count = " + numel(config.h_list)
    ];
dpimnumeric.writeText(fullfile(outDir, "README.md"), strjoin(lines, newline));
dpimnumeric.writeJson(fullfile(outDir, "config.json"), config);
end

function result = localResultPaths(outDir)
result = struct();
result.summary = fullfile(outDir, "summary.csv");
result.config = fullfile(outDir, "config.json");
result.readme = fullfile(outDir, "README.md");
end

function localWriteExperimentStatus(resultsRoot, config)
status = [
    "# Original Plate FEM Experimental Status"
    ""
    "Completed scope: E7/E8 original thin-plate FEM rebuild runner."
    ""
    "Skipped by user instruction:"
    "- GitHub push/authentication."
    "- Manuscript TeX/PDF replacement or compilation."
    "- E6 GF discrepancy experiment; E6 is cancelled."
    ""
    "Important distinction:"
    "- `results/hfallback_sweep_20260605_140246/E7_plate_SFEM` and E8 are standalone diagnostic RFEM approximations."
    "- This `original_plate_full_*` runner uses the original GF/KL/20x20 thin-plate FEM scripts and stores raw `wc` and probability weights."
    ""
    "Run mode: " + string(config.run_mode)
    "R = " + string(config.R) + ", M = " + string(config.M) + ", B = " + string(config.B)
    "sample_pool_count = " + string(config.sample_pool_count) + ", inner_sample_count = " + string(config.inner_sample_count)
    "reference_sample_count = " + string(config.reference_sample_count)
    ];
dpimnumeric.writeText(fullfile(resultsRoot, "EXPERIMENT_STATUS.md"), strjoin(status, newline));
end

function localRestore(oldPath, oldDir)
path(oldPath);
if isfolder(oldDir)
    cd(oldDir);
end
end
