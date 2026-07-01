function resultsRoot = run_kirchhoff_plate_small_deflection_reanalysis(runMode)
%RUN_KIRCHHOFF_PLATE_SMALL_DEFLECTION_REANALYSIS
% Rebuild the formal Kirchhoff-plate density and coverage experiment in SI units.
%
% The stored finite-element responses were generated with q=-1e5 Pa. Because
% the plate solve is linear in q, multiplying those raw responses by 1e-3 is
% exactly equivalent to rerunning the same systems with q=-100 Pa.

if nargin < 1 || strlength(string(runMode)) == 0
    runMode = "journal";
end
runMode = lower(string(runMode));

projectRoot = fileparts(mfilename("fullpath"));
sourceRoot = localFindSourceRoot(projectRoot);
config = localConfig(runMode);
resultsRoot = fullfile(projectRoot, "results", ...
    "kirchhoff_plate_small_deflection_" + runMode + "_" + ...
    string(datetime("now", "Format", "yyyyMMdd_HHmmss")));
dpimnumeric.ensureDir(resultsRoot);

refData = load(fullfile(sourceRoot, "reference", "thin_reference.mat"), "reference");
poolData = load(fullfile(sourceRoot, "rqmc_pool", "thin_rqmc_pool.mat"), "pool");
reference = refData.reference;
pool = poolData.pool;

if pool.sample_pool_count < config.P
    error("Source pool only has %d curves; %d requested.", pool.sample_pool_count, config.P);
end
if numel(reference.wc_all_mc_thin) ~= config.N_ref
    error("Reference sample count is %d; expected %d.", ...
        numel(reference.wc_all_mc_thin), config.N_ref);
end

scale = config.q_new / config.q_source;
referenceResponse = reference.wc_all_mc_thin(:) * scale;
referenceWeights = reference.PV_mc(:);
referenceWeights = referenceWeights / sum(referenceWeights);

poolResponse = cell(config.P, 1);
poolWeights = cell(config.P, 1);
for i = 1:config.P
    poolResponse{i} = pool.wc_runs_thin{i}(:) * scale;
    poolWeights{i} = pool.weight_runs{i}(:);
    poolWeights{i} = poolWeights{i} / sum(poolWeights{i});
end

[wGrid, responseAudit] = localResponseGrid(referenceResponse, config.grid_count, ...
    config.quantile_limits, config.grid_margin_fraction);
[hList, bandwidthAudit] = localBandwidthGrid(referenceResponse, ...
    config.n, config.h_multiplier_range, config.h_count);

config.source_results_root = sourceRoot;
config.response_scale_from_source = scale;
config.response_unit = "m";
config.w_grid = wGrid;
config.h_list = hList;
config.active_relative_density_threshold = config.active_relative_threshold;
config.created_at = char(datetime("now", "Format", "yyyy-MM-dd HH:mm:ss"));
dpimnumeric.writeJson(fullfile(resultsRoot, "config.json"), config);

analytic = localAnalyticalAudit(config, referenceResponse);
writetable(struct2table(analytic), fullfile(resultsRoot, "analytical_small_deflection_audit.csv"));
writetable(struct2table(responseAudit), fullfile(resultsRoot, "response_grid_audit.csv"));
writetable(struct2table(bandwidthAudit), fullfile(resultsRoot, "bandwidth_grid_audit.csv"));

fprintf("Building probability-weighted KDE curve pools: P=%d, G=%d, H=%d\n", ...
    config.P, numel(wGrid), numel(hList));
[referenceCurves, curvePools] = localBuildCurves( ...
    referenceResponse, referenceWeights, poolResponse, poolWeights, wGrid, hList);

tuningRows = cell(numel(hList), 1);
for ih = 1:numel(hList)
    fprintf("Tuning h %d/%d: %.12g m\n", ih, numel(hList), hList(ih));
    refCurve = referenceCurves(:, ih).';
    active = refCurve >= config.active_relative_threshold * max(refCurve);
    tuningRows{ih} = localCoverageSummary(curvePools(:, :, ih), refCurve, active, ...
        config.alpha, config.R, config.B_tuning, config.M_tuning, ...
        config.seed_tuning + ih * 10000, config.fallback_lambda);
    tuningRows{ih}.h_index(:) = ih;
    tuningRows{ih}.h(:) = hList(ih);
    tuningRows{ih}.active_grid_count(:) = nnz(active);
end
tuning = vertcat(tuningRows{:});
writetable(tuning, fullfile(resultsRoot, "tuning_summary.csv"));

selected = localSelectBandwidth(tuning);
btTuning = tuning(tuning.method == "bootstrap-t" & ...
    tuning.bootstrap_t_fallback_rate == 0 & tuning.bootstrap_t_inf_rate == 0, :);
btTuning = sortrows(btTuning, "h", "ascend");
candidatePosition = find(btTuning.h >= selected.h(1), 1, "first");
if isempty(candidatePosition)
    candidatePosition = height(btTuning);
end

while true
    selected = btTuning(candidatePosition, :);
    selectedIndex = selected.h_index(1);
    selectedH = selected.h(1);
    referenceCurve = referenceCurves(:, selectedIndex).';
    curvePool = curvePools(:, :, selectedIndex);
    active = referenceCurve >= config.active_relative_threshold * max(referenceCurve);

    fprintf("Formal validation at h=%.12g m with %d active grid points.\n", ...
        selectedH, nnz(active));
    [formalSummary, pointwise, diagnostics] = localFormalCoverage( ...
        curvePool, referenceCurve, active, wGrid, config);
    btFormal = formalSummary(formalSummary.method == "bootstrap-t", :);
    if btFormal.bootstrap_t_fallback_rate == 0 && btFormal.bootstrap_t_inf_rate == 0
        break;
    end
    if candidatePosition >= height(btTuning)
        warning("No formally validated bandwidth achieved zero fallback.");
        break;
    end
    candidatePosition = candidatePosition + 1;
end
selected.selection_rule = ...
    "smallest tuning-zero-fallback h at or above the best coverage candidate that also has zero fallback in formal validation";
selected.formal_bootstrap_t_coverage = btFormal.mean_pointwise_coverage;
selected.formal_bootstrap_t_fallback_rate = btFormal.bootstrap_t_fallback_rate;
selected.formal_bootstrap_t_inf_rate = btFormal.bootstrap_t_inf_rate;
writetable(selected, fullfile(resultsRoot, "selected_h.csv"));
writetable(formalSummary, fullfile(resultsRoot, "formal_summary.csv"));
writetable(pointwise, fullfile(resultsRoot, "formal_pointwise.csv"));

curvePoolMean = mean(curvePool, 1);
curvePoolSd = std(curvePool, 0, 1);
curveTable = table(wGrid(:), referenceCurve(:), curvePoolMean(:), curvePoolSd(:), active(:), ...
    'VariableNames', {'response_m', 'reference_density_per_m', ...
    'pool_mean_density_per_m', 'pool_sd_density_per_m', 'active_grid'});
writetable(curveTable, fullfile(resultsRoot, "selected_density_curve.csv"));

[centerTarget, centerPool, centerMeta] = localIndependentCenterTarget( ...
    projectRoot, wGrid, selectedH, scale, config.P);
validationMean = mean(curvePool, 1);
centeringShift = centerTarget - validationMean;
centeredCurvePool = curvePool + centeringShift;

centerConfig = config;
centerConfig.seed_formal = config.seed_formal + 700000;
fprintf("Independent-target coverage before centering.\n");
[independentTargetSummary, independentTargetPointwise, independentTargetDiagnostics] = ...
    localFormalCoverage(curvePool, centerTarget, active, wGrid, centerConfig);
fprintf("Centered coverage with independent target and translated validation pool.\n");
[centeredSummary, centeredPointwise, centeredDiagnostics] = localFormalCoverage( ...
    centeredCurvePool, centerTarget, active, wGrid, centerConfig);

independentTargetSummary.coverage_mode = repmat("independent target, unshifted validation pool", ...
    height(independentTargetSummary), 1);
centeredSummary.coverage_mode = repmat("independent target, centered validation pool", ...
    height(centeredSummary), 1);
formalSummary.coverage_mode = repmat("MC reference, end-to-end", height(formalSummary), 1);
coverageComparison = [formalSummary; independentTargetSummary; centeredSummary];
writetable(independentTargetSummary, ...
    fullfile(resultsRoot, "independent_target_unshifted_summary.csv"));
writetable(independentTargetPointwise, ...
    fullfile(resultsRoot, "independent_target_unshifted_pointwise.csv"));
writetable(centeredSummary, fullfile(resultsRoot, "centered_summary.csv"));
writetable(centeredPointwise, fullfile(resultsRoot, "centered_pointwise.csv"));
writetable(coverageComparison, fullfile(resultsRoot, "coverage_mode_comparison.csv"));

centeringAudit = table(wGrid(:), referenceCurve(:), centerTarget(:), validationMean(:), ...
    centeringShift(:), std(curvePool, 0, 1).', active(:), ...
    'VariableNames', {'response_m', 'mc_reference_density_per_m', ...
    'independent_center_target_per_m', 'validation_pool_mean_per_m', ...
    'centering_shift_per_m', 'validation_pool_sd_per_m', 'active_grid'});
centeringAudit.standardized_shift_for_R = centeringAudit.centering_shift_per_m ./ ...
    (centeringAudit.validation_pool_sd_per_m / sqrt(config.R));
writetable(centeringAudit, fullfile(resultsRoot, "centering_audit.csv"));
writetable(struct2table(centerMeta), fullfile(resultsRoot, "centering_pool_manifest.csv"));

originalVariance = var(curvePool, 0, 1);
centeredVariance = var(centeredCurvePool, 0, 1);
originalSkewness = skewness(curvePool, 0, 1);
centeredSkewness = skewness(centeredCurvePool, 0, 1);
originalExcessKurtosis = kurtosis(curvePool, 0, 1) - 3;
centeredExcessKurtosis = kurtosis(centeredCurvePool, 0, 1) - 3;
invarianceAudit = struct( ...
    "max_abs_mean_alignment_error", max(abs(mean(centeredCurvePool, 1) - centerTarget)), ...
    "max_abs_variance_difference", max(abs(centeredVariance - originalVariance)), ...
    "max_abs_skewness_difference", max(abs(centeredSkewness - originalSkewness)), ...
    "max_abs_excess_kurtosis_difference", ...
        max(abs(centeredExcessKurtosis - originalExcessKurtosis)), ...
    "max_abs_interval_length_difference", ...
        max(abs(centeredSummary.mean_interval_length - ...
        independentTargetSummary.mean_interval_length)));
writetable(struct2table(invarianceAudit), ...
    fullfile(resultsRoot, "centering_invariance_audit.csv"));

save(fullfile(resultsRoot, "formal_results.mat"), "config", "analytic", ...
    "responseAudit", "bandwidthAudit", "selected", "formalSummary", "pointwise", ...
    "diagnostics", "referenceCurve", "curvePoolMean", "curvePoolSd", ...
    "active", "wGrid", "hList", "centerTarget", "centerPool", ...
    "centeringShift", "independentTargetSummary", "independentTargetPointwise", ...
    "independentTargetDiagnostics", "centeredSummary", "centeredPointwise", ...
    "centeredDiagnostics", "centeringAudit", "centerMeta", ...
    "invarianceAudit", "-v7.3");

localMakeFigures(resultsRoot, tuning, selectedH, curveTable, pointwise, ...
    formalSummary, centeredPointwise, centeredSummary);
localWriteReadme(resultsRoot, config, analytic, selected, formalSummary, ...
    independentTargetSummary, centeredSummary, centerMeta);

fprintf("Kirchhoff plate small-deflection reanalysis completed:\n%s\n", resultsRoot);
end

function config = localConfig(runMode)
config = struct();
config.run_mode = char(runMode);
config.alpha = 0.05;
config.q_source = 1.0e5;
config.q_new = 100;
config.a = 1;
config.b = 1;
config.thickness = 0.005;
config.poisson = 0.30;
config.mean_E = 2.184e11;
config.n = 200;
config.P = 600;
config.R = 20;
config.N_ref = 10000;
config.grid_count = 201;
config.quantile_limits = [0.001, 0.999];
config.grid_margin_fraction = 0.15;
config.active_relative_threshold = 1e-2;
config.h_multiplier_range = [0.01, 25.0];
config.h_count = 21;
config.fallback_lambda = 5;
config.seed_tuning = 2026062201;
config.seed_formal = 2026062202;
switch runMode
    case "journal"
        config.M_tuning = 200;
        config.B_tuning = 199;
        config.M_formal = 1000;
        config.B_formal = 399;
    case "medium"
        config.M_tuning = 80;
        config.B_tuning = 99;
        config.M_formal = 300;
        config.B_formal = 199;
    otherwise
        config.P = 80;
        config.M_tuning = 20;
        config.B_tuning = 49;
        config.M_formal = 60;
        config.B_formal = 99;
        config.grid_count = 101;
        config.h_count = 9;
end
end

function sourceRoot = localFindSourceRoot(projectRoot)
candidates = dir(fullfile(projectRoot, "results", "original_plate_journal_*"));
candidates = candidates([candidates.isdir]);
for i = numel(candidates):-1:1
    candidate = fullfile(candidates(i).folder, candidates(i).name);
    if isfile(fullfile(candidate, "reference", "thin_reference.mat")) && ...
            isfile(fullfile(candidate, "rqmc_pool", "thin_rqmc_pool.mat"))
        sourceRoot = candidate;
        return;
    end
end
error("No completed original_plate_journal source result was found.");
end

function [target, curvePool, meta] = localIndependentCenterTarget( ...
    projectRoot, wGrid, h, scale, validationPoolCount)
files = dir(fullfile(projectRoot, "rqmc_outer_sample_dual_plate_*.mat"));
if isempty(files)
    error("No independent dual-plate RQMC samples were found for centering.");
end
[~, order] = sort({files.name});
files = files(order);
curvePool = zeros(numel(files), numel(wGrid));
seeds = nan(numel(files), 1);
for i = 1:numel(files)
    data = load(fullfile(files(i).folder, files(i).name), "thin", "PV_rqmc", "seed");
    if ~isfield(data, "thin") || ~isfield(data.thin, "wc_all_rqmc")
        error("%s does not contain thin.wc_all_rqmc.", files(i).name);
    end
    response = data.thin.wc_all_rqmc(:) * scale;
    weights = data.PV_rqmc(:);
    weights = weights / sum(weights);
    curvePool(i, :) = localWeightedKde(response, weights, wGrid, h);
    if isfield(data, "seed")
        seeds(i) = data.seed;
    end
end
target = mean(curvePool, 1);
meta = struct( ...
    "center_pool_source", "rqmc_outer_sample_dual_plate_*.mat", ...
    "center_pool_count", numel(files), ...
    "center_pool_seed_min", min(seeds, [], "omitnan"), ...
    "center_pool_seed_max", max(seeds, [], "omitnan"), ...
    "validation_pool_count", validationPoolCount, ...
    "target_definition", "mean probability-weighted KDE ordinate over independent center pool", ...
    "centering_definition", "add target minus validation-pool mean at each response grid point");
end

function [wGrid, audit] = localResponseGrid(response, gridCount, probs, marginFraction)
q = quantile(response, probs);
span = q(2) - q(1);
if ~(isfinite(span) && span > 0)
    error("Reference response quantile span is not positive.");
end
wMin = q(1) - marginFraction * span;
wMax = q(2) + marginFraction * span;
wGrid = linspace(wMin, wMax, gridCount);
audit = struct( ...
    "quantile_low_probability", probs(1), ...
    "quantile_high_probability", probs(2), ...
    "quantile_low_m", q(1), ...
    "quantile_high_m", q(2), ...
    "margin_fraction", marginFraction, ...
    "grid_min_m", wMin, ...
    "grid_max_m", wMax, ...
    "grid_count", gridCount, ...
    "grid_spacing_m", wGrid(2) - wGrid(1));
end

function [hList, audit] = localBandwidthGrid(response, n, multiplierRange, hCount)
sigmaRobust = min(std(response), iqr(response) / 1.34);
hSilverman = 1.06 * sigmaRobust * n^(-1 / 5);
hList = hSilverman * logspace(log10(multiplierRange(1)), ...
    log10(multiplierRange(2)), hCount);
audit = struct( ...
    "response_sd_m", std(response), ...
    "response_iqr_m", iqr(response), ...
    "robust_sigma_m", sigmaRobust, ...
    "inner_sample_count", n, ...
    "silverman_h_m", hSilverman, ...
    "h_min_m", min(hList), ...
    "h_max_m", max(hList), ...
    "h_count", hCount);
end

function analytic = localAnalyticalAudit(config, referenceResponse)
q = config.q_new;
a = config.a;
b = config.b;
E = config.mean_E;
nu = config.poisson;
t = config.thickness;
D = E * t^3 / (12 * (1 - nu^2));
w = 0;
for m = 1:2:501
    for n = 1:2:501
        signCenter = (-1)^((m - 1) / 2 + (n - 1) / 2);
        qmn = 16 * q / (pi^2 * m * n);
        wave2 = (m * pi / a)^2 + (n * pi / b)^2;
        w = w + signCenter * qmn / (D * wave2^2);
    end
end
analytic = struct( ...
    "model", "Kirchhoff simply supported square plate", ...
    "load_pa", -q, ...
    "flexural_rigidity_Nm", D, ...
    "navier_center_deflection_m", -w, ...
    "navier_abs_w_over_t", abs(w) / t, ...
    "small_deflection_limit_w_over_t", 0.1, ...
    "reference_mean_deflection_m", mean(referenceResponse), ...
    "reference_sd_deflection_m", std(referenceResponse), ...
    "reference_max_abs_deflection_m", max(abs(referenceResponse)), ...
    "reference_max_abs_w_over_t", max(abs(referenceResponse)) / t, ...
    "mean_relative_difference_from_navier", ...
    abs(abs(mean(referenceResponse)) - abs(w)) / abs(w), ...
    "navier_series_max_odd_index", 501);
end

function [referenceCurves, curvePools] = localBuildCurves( ...
    referenceResponse, referenceWeights, poolResponse, poolWeights, wGrid, hList)
H = numel(hList);
P = numel(poolResponse);
G = numel(wGrid);
referenceCurves = zeros(G, H);
curvePools = zeros(P, G, H);
for ih = 1:H
    h = hList(ih);
    referenceCurves(:, ih) = localWeightedKde(referenceResponse, referenceWeights, wGrid, h);
    for ip = 1:P
        curvePools(ip, :, ih) = localWeightedKde( ...
            poolResponse{ip}, poolWeights{ip}, wGrid, h);
    end
end
end

function curve = localWeightedKde(response, weights, wGrid, h)
z = (wGrid(:) - response(:).') / h;
curve = (exp(-0.5 * z.^2) * weights(:)) / (sqrt(2 * pi) * h);
end

function summary = localCoverageSummary(curvePool, truth, active, alpha, R, B, M, seed, lambda)
[methodCoverage, meanLength, fallbackRate, infRate] = ...
    localPointwiseCoverage(curvePool, truth, active, alpha, R, B, M, seed, lambda, false);
methods = ["t distribution"; "percentile bootstrap"; "bootstrap-t"];
summary = table(methods(:), methodCoverage(:), meanLength(:), ...
    repmat(fallbackRate, 3, 1), repmat(infRate, 3, 1), ...
    abs(methodCoverage(:) - (1 - alpha)), ...
    zeros(3, 1), zeros(3, 1), zeros(3, 1), ...
    'VariableNames', {'method', 'mean_pointwise_coverage', 'mean_interval_length', ...
    'bootstrap_t_fallback_rate', 'bootstrap_t_inf_rate', ...
    'abs_error_to_nominal', 'h_index', 'h', 'active_grid_count'});
end

function selected = localSelectBandwidth(tuning)
bt = tuning(tuning.method == "bootstrap-t", :);
eligible = bt.bootstrap_t_fallback_rate == 0 & bt.bootstrap_t_inf_rate == 0;
if any(eligible)
    candidates = bt(eligible, :);
else
    minFallback = min(bt.bootstrap_t_fallback_rate);
    candidates = bt(bt.bootstrap_t_fallback_rate == minFallback, :);
end
withinTolerance = candidates.abs_error_to_nominal <= 0.015;
if any(withinTolerance)
    candidates = candidates(withinTolerance, :);
    [~, order] = sort(candidates.h, "ascend");
    selected = candidates(order(1), :);
    selected.selection_rule = "smallest h with |coverage-0.95| <= 0.015, zero fallback, and zero infinite bootstrap-t intervals";
else
    [~, order] = sortrows([candidates.abs_error_to_nominal, candidates.h], [1, 2]);
    selected = candidates(order(1), :);
    selected.selection_rule = "minimum |coverage-0.95| because no zero-fallback candidate met the 0.015 tolerance";
end
end

function [summary, pointwise, diagnostics] = localFormalCoverage(curvePool, truth, active, wGrid, config)
[coverage, meanLength, fallbackRate, infRate, details] = localPointwiseCoverage( ...
    curvePool, truth, active, config.alpha, config.R, config.B_formal, ...
    config.M_formal, config.seed_formal, config.fallback_lambda, true);

methods = ["t distribution"; "percentile bootstrap"; "bootstrap-t"];
summary = table(methods(:), coverage(:), meanLength(:), ...
    repmat(fallbackRate, 3, 1), repmat(infRate, 3, 1), ...
    abs(coverage(:) - (1 - config.alpha)), ...
    repmat(nnz(active), 3, 1), repmat(config.M_formal, 3, 1), ...
    repmat(config.B_formal, 3, 1), ...
    'VariableNames', {'method', 'mean_pointwise_coverage', 'mean_interval_length', ...
    'bootstrap_t_fallback_rate', 'bootstrap_t_inf_rate', ...
    'abs_error_to_nominal', 'active_grid_count', 'M', 'B'});

bandMcse = sqrt(details.simultaneous_coverage * ...
    (1 - details.simultaneous_coverage) / config.M_formal);
bandRow = table("bootstrap-t simultaneous band", details.simultaneous_coverage, ...
    details.mean_band_length, fallbackRate, details.band_inf_rate, ...
    abs(details.simultaneous_coverage - (1 - config.alpha)), nnz(active), ...
    config.M_formal, config.B_formal, ...
    'VariableNames', summary.Properties.VariableNames);
summary = [summary; bandRow];
summary.mcse = [nan(3, 1); bandMcse];
summary.mcse_95_half_width = [nan(3, 1); 1.96 * bandMcse];

pointwise = table(wGrid(:), truth(:), active(:), ...
    details.coverage_by_method(1, :).', details.coverage_by_method(2, :).', ...
    details.coverage_by_method(3, :).', details.band_pointwise_coverage(:), ...
    'VariableNames', {'response_m', 'reference_density_per_m', 'active_grid', ...
    't_distribution_coverage', 'percentile_bootstrap_coverage', ...
    'bootstrap_t_coverage', 'simultaneous_band_pointwise_coverage'});
diagnostics = details;
end

function [coverage, meanLength, fallbackRate, infRate, details] = localPointwiseCoverage( ...
    curvePool, truth, active, alpha, R, B, M, seed, lambda, retainDetails)
G = size(curvePool, 2);
coverageHits = zeros(3, G);
lengthSums = zeros(3, G);
fallbackCount = 0;
infCount = 0;
activeCount = nnz(active);
simHits = 0;
bandPointHits = zeros(1, G);
bandLengthSum = zeros(1, G);
bandInfCount = 0;

if retainDetails
    sampleIndices = zeros(M, R);
    bandCritical = nan(M, 1);
else
    sampleIndices = [];
    bandCritical = [];
end

previousState = rng;
rng(double(seed), "twister");
cleanupObj = onCleanup(@() rng(previousState)); %#ok<NASGU>
tCrit = localStudentTCrit(1 - alpha / 2, R - 1);
kLower = max(1, floor((B + 1) * alpha / 2));
kUpper = min(B, ceil((B + 1) * (1 - alpha / 2)));
kBand = min(B, ceil((B + 1) * (1 - alpha)));

for trial = 1:M
    idx = randperm(size(curvePool, 1), R);
    if retainDetails
        sampleIndices(trial, :) = idx;
    end
    Y = curvePool(idx, :);
    mu = mean(Y, 1);
    sd = std(Y, 0, 1);
    se = sd / sqrt(R);

    lowerT = mu - tCrit * se;
    upperT = mu + tCrit * se;

    bootIdx = randi(R, R, B);
    bootCube = reshape(Y(bootIdx(:), :), R, B, G);
    bootMeans = reshape(mean(bootCube, 1), B, G);
    bootStd = reshape(std(bootCube, 0, 1), B, G);
    sortedMeans = sort(bootMeans, 1);
    lowerPct = sortedMeans(kLower, :);
    upperPct = sortedMeans(kUpper, :);

    Tstar = sqrt(R) * (bootMeans - mu) ./ bootStd;
    bad = bootStd == 0 | ~isfinite(Tstar);
    if any(bad, "all")
        delta = bootMeans - mu;
        s = sign(delta);
        s(s == 0) = 1;
        Tstar(bad) = s(bad) .* Inf;
    end
    sortedT = sort(Tstar, 1);
    qLower = sortedT(kLower, :);
    qUpper = sortedT(kUpper, :);
    lowerBt = mu - sd .* qUpper / sqrt(R);
    upperBt = mu - sd .* qLower / sqrt(R);

    lower = [lowerT; lowerPct; lowerBt];
    upper = [upperT; upperPct; upperBt];
    for im = 1:3
        coverageHits(im, :) = coverageHits(im, :) + ...
            (lower(im, :) <= truth & truth <= upper(im, :));
        len = upper(im, :) - lower(im, :);
        len(~isfinite(len)) = 0;
        lengthSums(im, :) = lengthSums(im, :) + len;
    end

    btLength = upperBt - lowerBt;
    pctLength = upperPct - lowerPct;
    btInf = ~(isfinite(lowerBt) & isfinite(upperBt));
    trigger = btInf | btLength > lambda * max(pctLength, realmin);
    fallbackCount = fallbackCount + nnz(trigger(active));
    infCount = infCount + nnz(btInf(active));

    if retainDetails
        bootSe = bootStd / sqrt(R);
        z = abs((bootMeans(:, active) - mu(active)) ./ bootSe(:, active));
        z(~isfinite(z)) = Inf;
        maxStats = max(z, [], 2);
        sortedMax = sort(maxStats);
        critical = sortedMax(kBand);
        bandCritical(trial) = critical;
        lowerBand = mu - critical * se;
        upperBand = mu + critical * se;
        bandHit = lowerBand <= truth & truth <= upperBand;
        simHits = simHits + all(bandHit(active));
        bandPointHits = bandPointHits + bandHit;
        finiteLength = upperBand - lowerBand;
        finiteLength(~isfinite(finiteLength)) = 0;
        bandLengthSum = bandLengthSum + finiteLength;
        bandInfCount = bandInfCount + double(~isfinite(critical));
    end
end

coverageByMethod = coverageHits / M;
coverage = mean(coverageByMethod(:, active), 2);
meanLengthByMethod = lengthSums / M;
meanLength = mean(meanLengthByMethod(:, active), 2);
fallbackRate = fallbackCount / (M * activeCount);
infRate = infCount / (M * activeCount);

details = struct();
details.coverage_by_method = coverageByMethod;
details.mean_length_by_method = meanLengthByMethod;
details.sample_indices = sampleIndices;
details.band_critical_value = bandCritical;
details.simultaneous_coverage = simHits / max(M, 1);
details.band_pointwise_coverage = bandPointHits / max(M, 1);
meanBandLengthByGrid = bandLengthSum / max(M, 1);
details.mean_band_length = mean(meanBandLengthByGrid(active), "omitnan");
details.band_inf_rate = bandInfCount / max(M, 1);
end

function q = localStudentTCrit(p, nu)
tailProb = 2 * min(p, 1 - p);
x = betaincinv(tailProb, nu / 2, 0.5);
q = abs(sqrt(nu * (1 / x - 1)));
end

function localMakeFigures(resultsRoot, tuning, selectedH, curveTable, pointwise, ...
    summary, centeredPointwise, centeredSummary)
figDir = fullfile(resultsRoot, "figures");
dpimnumeric.ensureDir(figDir);

f1 = figure("Visible", "off", "Color", "w", "Position", [100, 100, 760, 500]);
hold on;
methods = ["t distribution", "percentile bootstrap", "bootstrap-t"];
colors = [0.122, 0.306, 0.475; 0.180, 0.490, 0.196; 0.698, 0.133, 0.133];
for i = 1:3
    rows = tuning.method == methods(i);
    semilogx(tuning.h(rows), tuning.mean_pointwise_coverage(rows), ...
        "o-", "LineWidth", 1.4, "Color", colors(i, :), "DisplayName", methods(i));
end
yline(0.95, "k--", "LineWidth", 1.1, "DisplayName", "nominal");
xline(selectedH, "Color", [0.2, 0.2, 0.2], "LineStyle", ":", ...
    "LineWidth", 1.2, "DisplayName", "selected h");
xlabel("Bandwidth h (m)");
ylabel("Mean pointwise coverage");
legend("Location", "best");
grid on;
exportgraphics(f1, fullfile(figDir, "plate_bandwidth_tuning.pdf"), "ContentType", "vector");
close(f1);

f2 = figure("Visible", "off", "Color", "w", "Position", [100, 100, 760, 500]);
plot(curveTable.response_m, curveTable.reference_density_per_m, ...
    "k-", "LineWidth", 1.7, "DisplayName", "MC reference");
hold on;
plot(curveTable.response_m, curveTable.pool_mean_density_per_m, ...
    "-", "LineWidth", 1.4, "Color", [0.10, 0.45, 0.72], "DisplayName", "RQMC pool mean");
xlabel("Center deflection w_c (m)");
ylabel("Probability density (m^{-1})");
legend("Location", "best");
grid on;
exportgraphics(f2, fullfile(figDir, "plate_density_reference_vs_pool.pdf"), "ContentType", "vector");
close(f2);

active = pointwise.active_grid;
f3 = figure("Visible", "off", "Color", "w", "Position", [100, 100, 760, 500]);
plot(pointwise.response_m(active), pointwise.t_distribution_coverage(active), ...
    "-", "LineWidth", 1.3, "Color", colors(1, :), "DisplayName", "t distribution");
hold on;
plot(pointwise.response_m(active), pointwise.percentile_bootstrap_coverage(active), ...
    "--", "LineWidth", 1.3, "Color", colors(2, :), "DisplayName", "percentile bootstrap");
plot(pointwise.response_m(active), pointwise.bootstrap_t_coverage(active), ...
    "-.", "LineWidth", 1.3, "Color", colors(3, :), "DisplayName", "bootstrap-t");
yline(0.95, "k--", "LineWidth", 1.1, "DisplayName", "nominal");
xlabel("Center deflection w_c (m)");
ylabel("Pointwise coverage");
legend("Location", "best");
grid on;
exportgraphics(f3, fullfile(figDir, "plate_formal_pointwise_coverage.pdf"), "ContentType", "vector");
close(f3);

f4 = figure("Visible", "off", "Color", "w", "Position", [100, 100, 660, 440]);
bar(1:height(summary), ...
    [summary.mean_pointwise_coverage, centeredSummary.mean_pointwise_coverage], 0.75);
yline(0.95, "k--", "LineWidth", 1.1);
xticks(1:height(summary));
xticklabels({"t distribution", "percentile bootstrap", "bootstrap-t", "simultaneous band"});
xtickangle(12);
ylabel("Coverage");
ylim([max(0, min(summary.mean_pointwise_coverage) - 0.08), 1]);
legend({"MC reference", "centered target", "nominal"}, "Location", "southoutside", ...
    "Orientation", "horizontal");
grid on;
exportgraphics(f4, fullfile(figDir, "plate_formal_coverage_summary.pdf"), "ContentType", "vector");
close(f4);

centerActive = centeredPointwise.active_grid;
f5 = figure("Visible", "off", "Color", "w", "Position", [100, 100, 760, 500]);
plot(centeredPointwise.response_m(centerActive), ...
    centeredPointwise.t_distribution_coverage(centerActive), ...
    "-", "LineWidth", 1.3, "Color", colors(1, :), "DisplayName", "t distribution");
hold on;
plot(centeredPointwise.response_m(centerActive), ...
    centeredPointwise.percentile_bootstrap_coverage(centerActive), ...
    "--", "LineWidth", 1.3, "Color", colors(2, :), "DisplayName", "percentile bootstrap");
plot(centeredPointwise.response_m(centerActive), ...
    centeredPointwise.bootstrap_t_coverage(centerActive), ...
    "-.", "LineWidth", 1.3, "Color", colors(3, :), "DisplayName", "bootstrap-t");
yline(0.95, "k--", "LineWidth", 1.1, "DisplayName", "nominal");
xlabel("Center deflection w_c (m)");
ylabel("Centered pointwise coverage");
legend("Location", "best");
grid on;
exportgraphics(f5, fullfile(figDir, "plate_centered_pointwise_coverage.pdf"), ...
    "ContentType", "vector");
close(f5);
end

function localWriteReadme(resultsRoot, config, analytic, selected, summary, ...
    independentTargetSummary, centeredSummary, centerMeta)
lines = [
    "# Kirchhoff plate small-deflection reanalysis"
    ""
    "- Mechanical model: simply supported Kirchhoff square plate."
    "- Load: q = -" + config.q_new + " Pa."
    "- Response unit: metre; no hidden /1000 conversion is used in the KDE."
    "- Stored q=-1e5 Pa responses were multiplied by 1e-3 using exact load linearity."
    "- Navier mean-field center deflection: " + sprintf("%.8e m", analytic.navier_center_deflection_m)
    "- Navier |w|/t: " + sprintf("%.6f", analytic.navier_abs_w_over_t)
    "- Maximum sampled |w|/t: " + sprintf("%.6f", analytic.reference_max_abs_w_over_t)
    "- Selected h: " + sprintf("%.8e m", selected.h)
    "- Active grid count: " + selected.active_grid_count
    "- Formal M: " + config.M_formal + ", B: " + config.B_formal + ", R: " + config.R
    "- Independent center-pool curves: " + centerMeta.center_pool_count
    ""
    "The simultaneous result covers only the finite active response grid."
    ""
    "## End-to-end summary"
];
for i = 1:height(summary)
    lines(end + 1) = "- " + summary.method(i) + ": " + ...
        sprintf("%.6f", summary.mean_pointwise_coverage(i)); %#ok<AGROW>
end
lines(end + 1) = "";
lines(end + 1) = "## Independent target without centering";
for i = 1:height(independentTargetSummary)
    lines(end + 1) = "- " + independentTargetSummary.method(i) + ": " + ...
        sprintf("%.6f", independentTargetSummary.mean_pointwise_coverage(i)); %#ok<AGROW>
end
lines(end + 1) = "";
lines(end + 1) = "## Centered summary";
for i = 1:height(centeredSummary)
    lines(end + 1) = "- " + centeredSummary.method(i) + ": " + ...
        sprintf("%.6f", centeredSummary.mean_pointwise_coverage(i)); %#ok<AGROW>
end
dpimnumeric.writeText(fullfile(resultsRoot, "README.md"), strjoin(lines, newline));
end
