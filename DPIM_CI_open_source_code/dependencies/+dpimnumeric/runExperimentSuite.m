function outputs = runExperimentSuite(projectRoot, resultsRoot, runMode, startAt)
%runExperimentSuite Generate standardized experiment outputs for E1--E8.

if nargin < 3 || strlength(string(runMode)) == 0
    runMode = "small";
end
if nargin < 4 || strlength(string(startAt)) == 0
    startAt = "";
end

runMode = lower(string(runMode));
if ~any(runMode == ["small", "medium", "full"])
    error("runMode must be 'small', 'medium', or 'full'.");
end

dpimnumeric.ensureDir(resultsRoot);

experimentNames = [ ...
    "E1_finite_B"; ...
    "E2_R_order"; ...
    "E3_linear_beam_n_h"; ...
    "E4_nonlinear_tail"; ...
    "E5_bootstrap_t_instability"; ...
    "E7_plate_SFEM"; ...
    "E8_simultaneous_band"];

outputs = struct();
if strlength(string(startAt)) == 0
    startIndex = 1;
else
    startIndex = find(experimentNames == string(startAt), 1, "first");
    if isempty(startIndex)
        error("Unknown startAt experiment %s.", string(startAt));
    end
end

for iExp = startIndex:numel(experimentNames)
    expName = experimentNames(iExp);
    outDir = fullfile(resultsRoot, char(expName));
    dpimnumeric.ensureDir(outDir);
    dpimnumeric.ensureDir(fullfile(outDir, "figures"));

    switch expName
        case "E1_finite_B"
            result = localRunE1(projectRoot, outDir, runMode);
        case "E2_R_order"
            result = localRunE2(projectRoot, outDir, runMode);
        case "E3_linear_beam_n_h"
            result = localRunE3(projectRoot, outDir, runMode);
        case "E4_nonlinear_tail"
            result = localRunE4(projectRoot, outDir, runMode);
        case "E5_bootstrap_t_instability"
            result = localRunE5(projectRoot, outDir, runMode);
        case "E7_plate_SFEM"
            result = localRunE7(projectRoot, outDir, runMode);
        case "E8_simultaneous_band"
            result = localRunE8(projectRoot, outDir, runMode);
        otherwise
            error("Unknown experiment %s.", expName);
    end

    outputs.(matlab.lang.makeValidName(char(expName))) = result;
end
end

function result = localRunE1(projectRoot, outDir, runMode)
switch runMode
    case "full"
        config = struct("experiment", "E1_finite_B", "alpha", 0.05, ...
            "B_list", [399, 999, 4999, 9999], "simulation_count", 200000, ...
            "seed", 20260604, "source_reuse", "run_exp1_finite_B_baseline.m and smoke_finite_B_exact_baseline.m");
    case "medium"
        config = struct("experiment", "E1_finite_B", "alpha", 0.05, ...
            "B_list", [99, 199, 399, 999], "simulation_count", 80000, ...
            "seed", 20260604, "source_reuse", "run_exp1_finite_B_baseline.m and smoke_finite_B_exact_baseline.m");
    otherwise
        config = struct("experiment", "E1_finite_B", "alpha", 0.05, ...
            "B_list", [19, 39, 99, 199], "simulation_count", 20000, ...
            "seed", 20260604, "source_reuse", "run_exp1_finite_B_baseline.m and smoke_finite_B_exact_baseline.m");
end
localWriteConfig(outDir, config);
rng(config.seed, "twister");

nominal = 1 - config.alpha;
BList = config.B_list(:);
nB = numel(BList);
rankCount = BList + 1;
kMinus = zeros(nB, 1);
kPlus = zeros(nB, 1);
coverageRankCount = zeros(nB, 1);
c0b = zeros(nB, 1);
gridError = zeros(nB, 1);
absGridError = zeros(nB, 1);
lowerTailProbability = zeros(nB, 1);
upperTailProbability = zeros(nB, 1);
tailImbalance = zeros(nB, 1);
simCoverage = zeros(nB, 1);
simAbsError = zeros(nB, 1);
simStdError = zeros(nB, 1);
simCiHalfWidth95 = zeros(nB, 1);
simZScoreVsExact = zeros(nB, 1);

for iB = 1:nB
    B = BList(iB);
    kMinus(iB) = floor((config.alpha / 2) * (B + 1));
    kPlus(iB) = ceil((1 - config.alpha / 2) * (B + 1));
    coverageRankCount(iB) = kPlus(iB) - kMinus(iB);
    c0b(iB) = coverageRankCount(iB) / rankCount(iB);
    gridError(iB) = c0b(iB) - nominal;
    absGridError(iB) = abs(gridError(iB));
    lowerTailProbability(iB) = kMinus(iB) / rankCount(iB);
    upperTailProbability(iB) = (rankCount(iB) - kPlus(iB)) / rankCount(iB);
    tailImbalance(iB) = lowerTailProbability(iB) - upperTailProbability(iB);

    ranks = randi(B + 1, config.simulation_count, 1);
    simCoverage(iB) = mean((ranks > kMinus(iB)) & (ranks <= kPlus(iB)));
    simAbsError(iB) = abs(simCoverage(iB) - c0b(iB));
    simStdError(iB) = sqrt(c0b(iB) * (1 - c0b(iB)) / config.simulation_count);
    simCiHalfWidth95(iB) = 1.96 * sqrt(simCoverage(iB) * (1 - simCoverage(iB)) / config.simulation_count);
    simZScoreVsExact(iB) = (simCoverage(iB) - c0b(iB)) / max(simStdError(iB), realmin);
end

summary = table(BList, rankCount, kMinus, kPlus, coverageRankCount, c0b, ...
    repmat(nominal, nB, 1), gridError, absGridError, ...
    repmat(config.alpha / 2, nB, 1), lowerTailProbability, upperTailProbability, tailImbalance, ...
    simCoverage, simAbsError, simStdError, simCiHalfWidth95, simZScoreVsExact, ...
    repmat(config.simulation_count, nB, 1), ...
    'VariableNames', {'B', 'rank_count', 'k_minus', 'k_plus', 'coverage_rank_count', ...
    'C0B', 'nominal', 'grid_error', 'abs_grid_error', ...
    'nominal_one_sided_tail', 'lower_tail_probability', 'upper_tail_probability', 'tail_imbalance', ...
    'simulated_coverage', 'simulation_abs_error', 'simulation_std_error_vs_exact', ...
    'simulation_ci_half_width_95', 'simulation_z_score_vs_exact', 'simulation_count'});
writetable(summary, fullfile(outDir, "summary.csv"));

finiteBMethodSummary = localFiniteBMethodSummary(summary);
writetable(finiteBMethodSummary, fullfile(outDir, "finite_B_ci_method_summary.csv"));

raw = struct("summary", summary, "finite_B_method_summary", finiteBMethodSummary, "config", config);
save(fullfile(outDir, "raw_results.mat"), "-struct", "raw");

fig = figure("Visible", "off", "Color", "w");
tiledlayout(1, 2, "TileSpacing", "compact");
nexttile;
plot(BList, gridError, "-o", "LineWidth", 1.5);
xlabel("B");
ylabel("C0B - nominal");
title("Finite-B grid coverage error", "Interpreter", "none");
grid on;
nexttile;
errorbar(BList, simCoverage, simCiHalfWidth95, "-o", "LineWidth", 1.2);
hold on;
plot(BList, c0b, "--", "LineWidth", 1.2);
hold off;
xlabel("B");
ylabel("coverage");
title("Rank simulation vs exact C0B", "Interpreter", "none");
grid on;
localSaveFigure(fig, fullfile(outDir, "figures", "finite_B_grid_error.png"));

localWriteStandardFiles(projectRoot, outDir, config, ...
    ["run_exp1_finite_B_baseline.m", "smoke_finite_B_exact_baseline.m"], ...
    "Exact finite-B integerization output.", "completed_" + runMode, runMode);
result = localResultPaths(outDir);
end

function methodSummary = localFiniteBMethodSummary(summary)
methods = ["Percentile bootstrap finite-B rank"; "Bootstrap-t finite-B rank"];
endpointMapping = ["percentile endpoints use bootstrap order statistics directly"; ...
    "bootstrap-t endpoints use studentized pivot quantiles with lower/upper signs reversed"];
methodSummary = table();
for iMethod = 1:numel(methods)
    block = table( ...
        repmat(methods(iMethod), height(summary), 1), ...
        summary.B, summary.rank_count, summary.k_minus, summary.k_plus, ...
        summary.coverage_rank_count, summary.C0B, summary.nominal, ...
        summary.grid_error, summary.abs_grid_error, ...
        summary.lower_tail_probability, summary.upper_tail_probability, summary.tail_imbalance, ...
        repmat("k_minus < rank <= k_plus", height(summary), 1), ...
        repmat(endpointMapping(iMethod), height(summary), 1), ...
        repmat("C0B=(k_plus-k_minus)/(B+1)", height(summary), 1), ...
        'VariableNames', {'CI_method', 'B', 'rank_count', 'k_minus', 'k_plus', ...
        'coverage_rank_count', 'C0B', 'nominal_coverage', 'coverage_error', ...
        'abs_coverage_error', 'lower_tail_probability', 'upper_tail_probability', ...
        'tail_imbalance', 'rank_interval_expression', 'endpoint_mapping', ...
        'finite_B_expression'});
    methodSummary = [methodSummary; block]; %#ok<AGROW>
end
end

function result = localRunE2(projectRoot, outDir, runMode)
switch runMode
    case "full"
        config = struct("experiment", "E2_R_order", "alpha", 0.05, "R_list", [10, 20, 40, 80, 160, 320], ...
            "B", 4999, "M", 2000, "seed", 20260604 + 2, "lambda", 5, ...
            "data_model", "centered exponential X=-log(U)-1; true mean 0", ...
            "source_reuse", "ci_methods.m and compute_hybrid_coverage.m");
    case "medium"
        config = struct("experiment", "E2_R_order", "alpha", 0.05, "R_list", [10, 20, 40, 80, 160], ...
            "B", 999, "M", 800, "seed", 20260604 + 2, "lambda", 5, ...
            "data_model", "centered exponential X=-log(U)-1; true mean 0", ...
            "source_reuse", "ci_methods.m and compute_hybrid_coverage.m");
    otherwise
        config = struct("experiment", "E2_R_order", "alpha", 0.05, "R_list", [8, 16, 32, 64], ...
            "B", 199, "M", 160, "seed", 20260604 + 2, "lambda", 5, ...
            "data_model", "centered exponential X=-log(U)-1; true mean 0", ...
            "source_reuse", "ci_methods.m and compute_hybrid_coverage.m");
end
localWriteConfig(outDir, config);
rng(config.seed, "twister");

summary = table();
for iR = 1:numel(config.R_list)
    R = config.R_list(iR);
    accumulator = localEmptyAccumulator(["Student-t", "Percentile bootstrap", "Bootstrap-t", "Hybrid"]);

    for trial = 1:config.M
        sample = -log(rand(R, 1)) - 1;
        ci = localCiMethods(sample, config.alpha, config.B, config.seed + 1000 * trial + R, 0);
        ci = localAppendHybrid(ci, config.lambda);
        accumulator = localAccumulate(accumulator, ci, 0);
    end

    block = localAccumulatorTable(accumulator, config.M);
    block.R = repmat(R, height(block), 1);
    summary = [summary; block]; %#ok<AGROW>
end

summary = movevars(summary, "R", "Before", 1);
summary.loglog_slope = localAttachSlopes(summary, config.R_list, 1 - config.alpha);
writetable(summary, fullfile(outDir, "summary.csv"));

raw = struct("summary", summary, "config", config);
save(fullfile(outDir, "raw_results.mat"), "-struct", "raw");

fig = figure("Visible", "off", "Color", "w");
methodList = unique(summary.method, "stable");
hold on;
for iMethod = 1:numel(methodList)
    mask = summary.method == methodList(iMethod);
    err = abs(summary.coverage(mask) - (1 - config.alpha));
    loglog(summary.R(mask), max(err, eps), "-o", "LineWidth", 1.5, "DisplayName", methodList(iMethod));
end
hold off;
xlabel("R");
ylabel("|coverage - nominal|");
title("Small R-order coverage diagnostic", "Interpreter", "none");
grid on;
legend("Location", "best");
localSaveFigure(fig, fullfile(outDir, "figures", "R_order_coverage_error.png"));

localWriteStandardFiles(projectRoot, outDir, config, ["ci_methods.m", "compute_hybrid_coverage.m"], ...
    "Centered-exponential coverage diagnostic for R-order calibration.", "completed_" + runMode, runMode);
result = localResultPaths(outDir);
end

function result = localRunE3(projectRoot, outDir, runMode)
switch runMode
    case "full"
        config = struct("experiment", "E3_linear_beam_n_h", "alpha", 0.05, "d", 10, ...
            "x0", 0.5, "n_list", [64, 128, 256, 512, 1024], "h_list", [0.05, 0.10, 0.15, 0.30, 0.50], ...
            "R_list", [10, 20, 40], "B", 399, "M", 300, "seed", 20260604 + 3, ...
            "dy_dimension", 1, "source_reuse", "smoke_linear_gaussian_closed_form.m and run_exp4_beam_benchmark.m");
    case "medium"
        config = struct("experiment", "E3_linear_beam_n_h", "alpha", 0.05, "d", 10, ...
            "x0", 0.5, "n_list", [32, 64, 128, 256], "h_list", [0.15, 0.30, 0.50], ...
            "R_list", [10, 20], "B", 199, "M", 160, "seed", 20260604 + 3, ...
            "dy_dimension", 1, "source_reuse", "smoke_linear_gaussian_closed_form.m and run_exp4_beam_benchmark.m");
    otherwise
        config = struct("experiment", "E3_linear_beam_n_h", "alpha", 0.05, "d", 10, ...
            "x0", 0.5, "n_list", [16, 32], "h_list", [0.35, 0.70], ...
            "R_list", [8, 16], "B", 99, "M", 80, "seed", 20260604 + 3, ...
            "dy_dimension", 1, "source_reuse", "smoke_linear_gaussian_closed_form.m and run_exp4_beam_benchmark.m");
end
localWriteConfig(outDir, config);
rng(config.seed, "twister");

coeffs = localGreenCoefficients(config.x0, config.d);
responseStd = norm(coeffs, 2);
yTarget = 0;
summary = table();

for iN = 1:numel(config.n_list)
    nInner = config.n_list(iN);
    for iH = 1:numel(config.h_list)
        h = config.h_list(iH);
        truth = localNormalPdf(yTarget, 0, sqrt(responseStd^2 + h^2));
        for iR = 1:numel(config.R_list)
            R = config.R_list(iR);
            accumulator = localEmptyAccumulator(["Student-t", "Percentile bootstrap", "Bootstrap-t", "Hybrid"]);
            for trial = 1:config.M
                estimates = localLinearKernelEstimates(R, nInner, coeffs, yTarget, h);
                ci = localCiMethods(estimates, config.alpha, config.B, config.seed + 1000 * trial + R, truth);
                ci = localAppendHybrid(ci, 5);
                accumulator = localAccumulate(accumulator, ci, truth);
            end
            block = localAccumulatorTable(accumulator, config.M);
            block.n = repmat(nInner, height(block), 1);
            block.h = repmat(h, height(block), 1);
            block.R = repmat(R, height(block), 1);
            block.Rnh = repmat(R * nInner * h^config.dy_dimension, height(block), 1);
            block.truth_smoothed_density = repmat(truth, height(block), 1);
            summary = [summary; block]; %#ok<AGROW>
        end
    end
end

summary = movevars(summary, ["n", "h", "R", "Rnh", "truth_smoothed_density"], "Before", 1);
writetable(summary, fullfile(outDir, "summary.csv"));

yGrid = linspace(-3 * responseStd, 3 * responseStd, 121);
pdfTable = table(yGrid(:), 'VariableNames', {'y'});
for iH = 1:numel(config.h_list)
    h = config.h_list(iH);
    pdfTable.("pdf_h_" + string(strrep(num2str(h), ".", "p"))) = localNormalPdf(yGrid(:), 0, sqrt(responseStd^2 + h^2));
end
writetable(pdfTable, fullfile(outDir, "linear_beam_pdf_reference.csv"));

raw = struct("summary", summary, "pdf_reference", pdfTable, "coefficients", coeffs, "config", config);
save(fullfile(outDir, "raw_results.mat"), "-struct", "raw");

fig = figure("Visible", "off", "Color", "w");
hold on;
for iH = 1:numel(config.h_list)
    h = config.h_list(iH);
    plot(yGrid, localNormalPdf(yGrid, 0, sqrt(responseStd^2 + h^2)), "LineWidth", 1.5, ...
        "DisplayName", sprintf("h=%.2f", h));
end
hold off;
xlabel("response y");
ylabel("smoothed density");
title("Linear beam closed-form smoke reference", "Interpreter", "none");
grid on;
legend("Location", "best");
localSaveFigure(fig, fullfile(outDir, "figures", "linear_beam_pdf_reference.png"));

localWriteStandardFiles(projectRoot, outDir, config, ["smoke_linear_gaussian_closed_form.m", "run_exp4_beam_benchmark.m"], ...
    "Closed-form linear beam output using Green coefficients.", "completed_" + runMode, runMode);
result = localResultPaths(outDir);
end

function result = localRunE4(projectRoot, outDir, runMode)
if runMode ~= "small"
    artifact = dpimnumeric.runNonlinearTailFull(projectRoot, outDir, runMode);
    config = artifact.config;
    localWriteConfig(outDir, config);
    localWriteStandardFiles(projectRoot, outDir, config, artifact.dependencies, artifact.note, artifact.status, runMode);
    result = localResultPaths(outDir);
    return;
end

config = struct("experiment", "E4_nonlinear_tail", "alpha", 0.05, "d", 10, ...
    "b", 0, "c_t", 0, "alpha_model", 0.40, "beta_model", 0.80, ...
    "theta_lognormal_mu", 0, "theta_lognormal_sigma", 1, ...
    "distributions", ["normal", "lognormal"], "quantile_grid", [0.05, 0.50, 0.95], ...
    "truth_sample_count", 25000, "n_inner", 32, "R", 10, "B", 99, "M", 80, ...
    "seed", 20260604 + 4, "source_reuse", "nonlinear pilotfixed family; new wrapper uses z=mean(theta)");
localWriteConfig(outDir, config);
rng(config.seed, "twister");

summary = table();
densityRows = table();

for iDist = 1:numel(config.distributions)
    distName = config.distributions(iDist);
    thetaTruth = localDrawTheta(distName, config.truth_sample_count, config.d, config);
    responseTruth = localLockedNonlinearResponse(thetaTruth, config);
    yValues = localQuantiles(responseTruth, config.quantile_grid);
    h = 0.20 * max(localIqr(responseTruth), std(responseTruth));
    h = max(h, eps);
    skewValue = localSkewness(responseTruth);
    kurtValue = localExcessKurtosis(responseTruth);

    yPlot = linspace(min(yValues), max(yValues), 121);
    density = zeros(size(yPlot));
    for iy = 1:numel(yPlot)
        density(iy) = mean(localNormalPdf(yPlot(iy), responseTruth, h));
    end
    densityRows = [densityRows; table(repmat(distName, numel(yPlot), 1), yPlot(:), density(:), ...
        'VariableNames', {'distribution', 'y', 'smoothed_density'})]; %#ok<AGROW>

    for iy = 1:numel(yValues)
        yTarget = yValues(iy);
        truth = mean(localNormalPdf(yTarget, responseTruth, h));
        tailEffN = localEffectiveCount(localNormalPdf(yTarget, responseTruth, h));
        accumulator = localEmptyAccumulator(["Student-t", "Percentile bootstrap", "Bootstrap-t", "Hybrid"]);
        for trial = 1:config.M
            estimates = localNonlinearKernelEstimates(config.R, config.n_inner, distName, yTarget, h, config);
            ci = localCiMethods(estimates, config.alpha, config.B, config.seed + 1000 * trial + iy, truth);
            ci = localAppendHybrid(ci, 5);
            accumulator = localAccumulate(accumulator, ci, truth);
        end
        block = localAccumulatorTable(accumulator, config.M);
        block.distribution = repmat(distName, height(block), 1);
        block.y_quantile = repmat(config.quantile_grid(iy), height(block), 1);
        block.y_value = repmat(yTarget, height(block), 1);
        block.h = repmat(h, height(block), 1);
        block.truth_smoothed_density = repmat(truth, height(block), 1);
        block.response_skewness = repmat(skewValue, height(block), 1);
        block.response_excess_kurtosis = repmat(kurtValue, height(block), 1);
        block.tail_effective_sample_size = repmat(tailEffN, height(block), 1);
        summary = [summary; block]; %#ok<AGROW>
    end
end

summary = movevars(summary, ["distribution", "y_quantile", "y_value", "h", ...
    "truth_smoothed_density", "response_skewness", "response_excess_kurtosis", ...
    "tail_effective_sample_size"], "Before", 1);
writetable(summary, fullfile(outDir, "summary.csv"));
writetable(densityRows, fullfile(outDir, "nonlinear_tail_density_curves.csv"));

raw = struct("summary", summary, "density_curves", densityRows, "config", config);
save(fullfile(outDir, "raw_results.mat"), "-struct", "raw");

fig = figure("Visible", "off", "Color", "w");
hold on;
for iDist = 1:numel(config.distributions)
    mask = densityRows.distribution == config.distributions(iDist);
    plot(densityRows.y(mask), densityRows.smoothed_density(mask), "LineWidth", 1.5, ...
        "DisplayName", config.distributions(iDist));
end
hold off;
xlabel("response y");
ylabel("smoothed density");
title("Locked nonlinear tail smoke densities", "Interpreter", "none");
grid on;
legend("Location", "best");
localSaveFigure(fig, fullfile(outDir, "figures", "nonlinear_tail_density.png"));

localWriteStandardFiles(projectRoot, outDir, config, ["数值实验_对数正态非线性/integrands_rqmc_lognormal_nonlinear.m"], ...
    "Small nonlinear tail output using locked z=mean(theta) model.", "completed_small", runMode);
result = localResultPaths(outDir);
end

function result = localRunE5(projectRoot, outDir, runMode)
if runMode ~= "small"
    artifact = dpimnumeric.runBootstrapInstabilityFull(projectRoot, outDir, runMode);
    config = artifact.config;
    localWriteConfig(outDir, config);
    localWriteStandardFiles(projectRoot, outDir, config, artifact.dependencies, artifact.note, artifact.status, runMode);
    result = localResultPaths(outDir);
    return;
end

config = struct("experiment", "E5_bootstrap_t_instability", "alpha", 0.05, ...
    "d", 10, "b", 0, "c_t", 0, "alpha_model", 0.40, "beta_model", 0.80, ...
    "theta_lognormal_mu", 0, "theta_lognormal_sigma", 1, "distribution", "lognormal", ...
    "truth_sample_count", 25000, "n_inner", 24, "R", 10, "B", 199, "M", 120, ...
    "lambda", 5, "seed", 20260604 + 5, ...
    "source_reuse", "ci_methods_lognormal_nonlinear.m and compute_hybrid_coverage.m concepts");
localWriteConfig(outDir, config);
rng(config.seed, "twister");

thetaTruth = localDrawTheta(config.distribution, config.truth_sample_count, config.d, config);
responseTruth = localLockedNonlinearResponse(thetaTruth, config);
yTarget = localQuantiles(responseTruth, 0.95);
h = 0.20 * max(localIqr(responseTruth), std(responseTruth));
truth = mean(localNormalPdf(yTarget, responseTruth, h));

btContains = false(config.M, 1);
pContains = false(config.M, 1);
hybridContains = false(config.M, 1);
btLength = inf(config.M, 1);
pLength = inf(config.M, 1);
lengthRatio = inf(config.M, 1);
minStdRatio = nan(config.M, 1);
eventA = false(config.M, 1);

for trial = 1:config.M
    estimates = localNonlinearKernelEstimates(config.R, config.n_inner, config.distribution, yTarget, h, config);
    ci = localCiMethods(estimates, config.alpha, config.B, config.seed + 1000 * trial, truth);
    p = ci(2);
    bt = ci(3);
    pContains(trial) = p.contains;
    btContains(trial) = bt.contains;
    pLength(trial) = p.length;
    btLength(trial) = bt.length;
    lengthRatio(trial) = bt.length / max(p.length, realmin);
    minStdRatio(trial) = bt.min_std_ratio;
    eventA(trial) = bt.infinite || (bt.length > config.lambda * p.length);
    hybridContains(trial) = (~eventA(trial) && bt.contains) || (eventA(trial) && p.contains);
end

summary = table( ...
    string("Bootstrap-t"), mean(btContains), mean(btLength, "omitnan"), NaN, NaN, ...
    mean(lengthRatio(isfinite(lengthRatio)), "omitnan"), mean(minStdRatio, "omitnan"), mean(eventA), ...
    yTarget, h, truth, ...
    'VariableNames', {'method', 'coverage', 'mean_interval_length', 'left_miss', 'right_miss', ...
    'mean_d_bt_over_d_p', 'mean_min_boot_std_ratio', 'P_A_lambda', ...
    'y_value', 'h', 'truth_smoothed_density'});
summary = [summary; table(string("Percentile bootstrap"), mean(pContains), mean(pLength, "omitnan"), NaN, NaN, ...
    NaN, NaN, NaN, yTarget, h, truth, 'VariableNames', summary.Properties.VariableNames)]; %#ok<AGROW>
summary = [summary; table(string("Hybrid"), mean(hybridContains), ...
    mean((~eventA).*btLength + eventA.*pLength, "omitnan"), NaN, NaN, ...
    NaN, NaN, mean(eventA), yTarget, h, truth, 'VariableNames', summary.Properties.VariableNames)]; %#ok<AGROW>
writetable(summary, fullfile(outDir, "summary.csv"));

diagnostics = table((1:config.M)', pContains, btContains, hybridContains, pLength, btLength, ...
    lengthRatio, minStdRatio, eventA, 'VariableNames', {'trial', 'percentile_contains', ...
    'bootstrap_t_contains', 'hybrid_contains', 'percentile_length', 'bootstrap_t_length', ...
    'd_bt_over_d_p', 'min_boot_std_ratio', 'event_A_lambda'});
writetable(diagnostics, fullfile(outDir, "bootstrap_t_diagnostics.csv"));

raw = struct("summary", summary, "diagnostics", diagnostics, "config", config);
save(fullfile(outDir, "raw_results.mat"), "-struct", "raw");

fig = figure("Visible", "off", "Color", "w");
tiledlayout(1, 2);
nexttile;
histogram(min(lengthRatio, 20), 20);
xlabel("min(d_bt/d_p, 20)");
ylabel("count");
title("Length ratio", "Interpreter", "none");
nexttile;
histogram(minStdRatio, 20);
xlabel("min bootstrap std / sample std");
ylabel("count");
title("Small denominator", "Interpreter", "none");
localSaveFigure(fig, fullfile(outDir, "figures", "bootstrap_t_diagnostics.png"));

localWriteStandardFiles(projectRoot, outDir, config, ["ci_methods_lognormal.m", "compute_hybrid_coverage.m"], ...
    "Small bootstrap-t denominator diagnostic with hybrid fallback.", "completed_small", runMode);
result = localResultPaths(outDir);
end

function result = localRunE7(projectRoot, outDir, runMode)
switch runMode
    case "full"
        config = struct("experiment", "E7_plate_SFEM", "alpha", 0.05, "R", 20, ...
            "B", 399, "M", 300, "seed", 20260604 + 7, "model", "dual_plate_thin", ...
            "studentize_std_floor_rel", 0.02, "max_curve_files", 800, ...
            "reference_pool_sigma", 1.0, "paper_text_sigma_dy", 1010 / 799, ...
            "matched_e4_parameter_block", "d5_M300_B399_N800_R20_k6_10", ...
            "source_reuse", "rqmc_outer_sample_dual_plate_*.mat and mc_reference_truth_dual_plate_thin_n10000_curve.txt");
    case "medium"
        config = struct("experiment", "E7_plate_SFEM", "alpha", 0.05, "R", 10, ...
            "B", 400, "M", 200, "seed", 20260604 + 7, "model", "dual_plate_thin", ...
            "studentize_std_floor_rel", 0.02, "max_curve_files", 120, ...
            "reference_pool_sigma", 1.0, "paper_text_sigma_dy", 1010 / 799, ...
            "source_reuse", "rqmc_outer_sample_dual_plate_*.mat and mc_reference_truth_dual_plate_thin_n10000_curve.txt");
    otherwise
        config = struct("experiment", "E7_plate_SFEM", "alpha", 0.05, "R", 10, ...
            "B", 99, "M", 120, "seed", 20260604 + 7, "model", "dual_plate_thin", ...
            "studentize_std_floor_rel", 0.02, "max_curve_files", 80, ...
            "reference_pool_sigma", 1.0, "paper_text_sigma_dy", 1010 / 799, ...
            "source_reuse", "rqmc_outer_sample_dual_plate_*.mat and mc_reference_truth_dual_plate_thin_n10000_curve.txt");
end
localWriteConfig(outDir, config);
rng(config.seed, "twister");

[wGrid, referenceCurve, curvePool, statusNote] = localLoadDualPlatePool(projectRoot, "thin", config.max_curve_files);

if isempty(curvePool)
    summary = table(string("TODO"), string(statusNote), 'VariableNames', {'status', 'reason'});
    raw = struct("summary", summary, "config", config);
else
    [summary, pointwise] = localCurvePoolCoverage(curvePool, referenceCurve, wGrid, config);
    raw = struct("summary", summary, "pointwise", pointwise, "config", config);
    writetable(pointwise, fullfile(outDir, "pointwise_summary.csv"));
end
writetable(summary, fullfile(outDir, "summary.csv"));
save(fullfile(outDir, "raw_results.mat"), "-struct", "raw");

fig = figure("Visible", "off", "Color", "w");
if exist("pointwise", "var") && height(pointwise) > 0
    plot(pointwise.w, pointwise.student_t_coverage, "-o", pointwise.w, pointwise.bootstrap_t_coverage, "-s");
    legend("Student-t", "Bootstrap-t", "Location", "best");
    xlabel("w");
    ylabel("pointwise coverage");
else
    text(0.1, 0.5, "Plate files missing; interface only.");
end
title("Dual-plate SFEM coverage over existing pool", "Interpreter", "none");
grid on;
localSaveFigure(fig, fullfile(outDir, "figures", "plate_sfem_coverage.png"));

localWriteStandardFiles(projectRoot, outDir, config, ["run_exp5_plate_benchmark.m", "run_dual_plate_rqmc_sample_pool_ci_driver.m"], ...
    "Dual-plate SFEM output over the existing sigma=1 pool; the thesis-text sigma=dy pool is not rebuilt here.", ...
    "completed_" + runMode + "_existing_pool", runMode);
result = localResultPaths(outDir);
end

function result = localRunE8(projectRoot, outDir, runMode)
switch runMode
    case "full"
        config = struct("experiment", "E8_simultaneous_band", "alpha", 0.05, ...
            "R", 20, "B", 399, "M", 300, "seed", 20260604 + 8, ...
            "studentize_se_floor_rel", 0.02, "max_curve_files", 800, ...
            "reference_pool_sigma", 1.0, "paper_text_sigma_dy", 1010 / 799, ...
            "matched_e4_parameter_block", "d5_M300_B399_N800_R20_k6_10", ...
            "source_reuse", "finite-grid simultaneous-band wrapper over the existing dual-plate pool");
    case "medium"
        config = struct("experiment", "E8_simultaneous_band", "alpha", 0.05, ...
            "R", 10, "B", 400, "M", 200, "seed", 20260604 + 8, ...
            "studentize_se_floor_rel", 0.02, "max_curve_files", 120, ...
            "reference_pool_sigma", 1.0, "paper_text_sigma_dy", 1010 / 799, ...
            "source_reuse", "finite-grid simultaneous-band wrapper over the existing dual-plate pool");
    otherwise
        config = struct("experiment", "E8_simultaneous_band", "alpha", 0.05, ...
            "R", 10, "B", 99, "M", 120, "seed", 20260604 + 8, ...
            "studentize_se_floor_rel", 0.02, "max_curve_files", 80, ...
            "reference_pool_sigma", 1.0, "paper_text_sigma_dy", 1010 / 799, ...
            "source_reuse", "finite-grid simultaneous-band wrapper over the existing dual-plate pool");
end
localWriteConfig(outDir, config);
rng(config.seed, "twister");

[wGrid, referenceCurve, curvePool, statusNote] = localLoadDualPlatePool(projectRoot, "thin", config.max_curve_files);
if isempty(curvePool)
    summary = table(string("TODO"), string(statusNote), 'VariableNames', {'status', 'reason'});
    raw = struct("summary", summary, "config", config);
else
    [summary, pointwise] = localSimultaneousBandCoverage(curvePool, referenceCurve, wGrid, config);
    raw = struct("summary", summary, "pointwise", pointwise, "config", config);
    writetable(pointwise, fullfile(outDir, "pointwise_band_summary.csv"));
end
writetable(summary, fullfile(outDir, "summary.csv"));
save(fullfile(outDir, "raw_results.mat"), "-struct", "raw");

fig = figure("Visible", "off", "Color", "w");
if exist("pointwise", "var") && height(pointwise) > 0
    plot(pointwise.w, pointwise.pointwise_coverage, "-o", pointwise.w, pointwise.simultaneous_band_pointwise_coverage, "-s");
    xlabel("w");
    ylabel("coverage");
    legend("Pointwise interval", "Simultaneous band", "Location", "best");
else
    text(0.1, 0.5, "Plate files missing; interface only.");
end
title("Finite-grid simultaneous band over existing pool", "Interpreter", "none");
grid on;
localSaveFigure(fig, fullfile(outDir, "figures", "simultaneous_band_coverage.png"));

localWriteStandardFiles(projectRoot, outDir, config, ["run_dual_plate_rqmc_sample_pool_ci_driver.m"], ...
    "Finite-grid simultaneous-band output over the existing sigma=1 pool; the thesis-text sigma=dy pool is not rebuilt here.", ...
    "completed_" + runMode + "_existing_pool", runMode);
result = localResultPaths(outDir);
end

function localWriteConfig(outDir, config)
dpimnumeric.writeJson(fullfile(outDir, "config.json"), config);
end

function localWriteStandardFiles(projectRoot, outDir, config, dependencies, note, status, runMode)
if nargin < 7 || strlength(string(runMode)) == 0
    runMode = "small";
end

runMode = lower(string(runMode));
metadata = struct();
metadata.run_time = char(datetime("now"));
metadata.matlab_version = version;
if isfield(config, "seed")
    metadata.random_seed = config.seed;
else
    metadata.random_seed = NaN;
end
metadata.git_commit = "";
metadata.git_status = "not_a_git_repository_or_not_checked";
metadata.file_hashes = localWrapperHashes(projectRoot);
metadata.model_formula = "g(theta,t)=b+c_t*t+z+alpha*z^3+beta*z*e";
metadata.model_z = "mean(theta)";
metadata.model_e = "mean(theta.^2)-1";
metadata.parameters = config;
metadata.script_entry = which("run_dpim_ci_numeric_rebuild");
metadata.dependencies = dependencies;
metadata.status = status;
metadata.note = note;
metadata.run_mode = char(runMode);
dpimnumeric.writeJson(fullfile(outDir, "metadata.json"), metadata);

logLines = [
    "DPIM CI rebuild log"
    "status: " + string(status)
    "run_mode: " + string(runMode)
    "time: " + string(datetime("now"))
    "note: " + string(note)
    "entry: run_dpim_ci_numeric_rebuild('" + string(runMode) + "')"
    ""
];
dpimnumeric.writeText(fullfile(outDir, "run_log.txt"), strjoin(logLines, newline));

modeLabel = upper(char(runMode));
smallNote = "Small outputs are not final full-scale paper evidence.";
if runMode == "full"
    smallNote = "Full outputs are intended as paper-scale numerical artifacts, subject to the remaining risk log.";
elseif runMode == "medium"
    smallNote = "Medium outputs are trend-check artifacts, not final paper-scale evidence.";
end

if isfield(config, "experiment")
    experimentName = string(config.experiment);
else
    experimentName = string(getfield(config, "run_mode")); %#ok<GFLD>
end

readmeLines = [
    "# " + experimentName
    ""
    "This directory was generated by the DPIM CI rebuild wrapper."
    ""
    "- `config.json`: experiment configuration."
    "- `metadata.json`: MATLAB version, seed, wrapper hashes, model formula, and dependencies."
    "- `raw_results.mat`: raw MATLAB data used to make the summary."
    "- `summary.csv`: main tabular output."
    "- `figures/*.png`: figures regenerated from data."
    ""
    "Run mode: " + string(modeLabel)
    ""
    "Status: " + string(status)
    ""
    "Note: " + string(note)
    ""
    smallNote
    ""
];
dpimnumeric.writeText(fullfile(outDir, "README.md"), strjoin(readmeLines, newline));
end

function result = localResultPaths(outDir)
result = struct();
result.output_dir = outDir;
result.config = fullfile(outDir, "config.json");
result.run_log = fullfile(outDir, "run_log.txt");
result.raw_results = fullfile(outDir, "raw_results.mat");
result.summary = fullfile(outDir, "summary.csv");
result.metadata = fullfile(outDir, "metadata.json");
result.readme = fullfile(outDir, "README.md");
end

function hashes = localWrapperHashes(projectRoot)
files = ["run_dpim_ci_numeric_rebuild.m", "run_dpim_ci_numeric_rebuild_small.m", ...
    "+dpimnumeric/runExperimentSuite.m", "+dpimnumeric/runNonlinearTailFull.m", ...
    "+dpimnumeric/runBootstrapInstabilityFull.m", "+dpimnumeric/buildInventory.m", ...
    "+dpimnumeric/writeModelLock.m", "+dpimnumeric/writePaperReports.m"];
hashes = struct();
for iFile = 1:numel(files)
    absPath = fullfile(projectRoot, char(files(iFile)));
    fieldToken = char(files(iFile));
    fieldToken = strrep(fieldToken, "/", "_");
    fieldToken = strrep(fieldToken, "+", "_");
    fieldToken = strrep(fieldToken, ".", "_");
    fieldName = matlab.lang.makeValidName(fieldToken);
    if isfile(absPath)
        hashes.(fieldName) = localMd5(absPath);
    else
        hashes.(fieldName) = "";
    end
end
end

function hashText = localMd5(filePath)
try
    bytes = uint8(fileread(filePath));
    md = java.security.MessageDigest.getInstance("MD5");
    md.update(bytes);
    digest = typecast(md.digest(), "uint8");
    hashText = lower(reshape(dec2hex(digest)', 1, []));
catch
    hashText = "";
end
end

function accumulator = localEmptyAccumulator(methods)
n = numel(methods);
accumulator.methods = methods(:);
accumulator.coverageHits = zeros(n, 1);
accumulator.leftMiss = zeros(n, 1);
accumulator.rightMiss = zeros(n, 1);
accumulator.lengthSum = zeros(n, 1);
accumulator.finiteLengthCount = zeros(n, 1);
end

function accumulator = localAccumulate(accumulator, ci, truth)
for i = 1:numel(accumulator.methods)
    name = accumulator.methods(i);
    idx = find([ci.name] == name, 1);
    if isempty(idx)
        continue;
    end
    c = ci(idx);
    accumulator.coverageHits(i) = accumulator.coverageHits(i) + double(c.contains);
    accumulator.leftMiss(i) = accumulator.leftMiss(i) + double(c.upper < truth);
    accumulator.rightMiss(i) = accumulator.rightMiss(i) + double(c.lower > truth);
    if isfinite(c.length)
        accumulator.lengthSum(i) = accumulator.lengthSum(i) + c.length;
        accumulator.finiteLengthCount(i) = accumulator.finiteLengthCount(i) + 1;
    end
end
end

function tbl = localAccumulatorTable(accumulator, M)
n = numel(accumulator.methods);
meanLength = nan(n, 1);
mask = accumulator.finiteLengthCount > 0;
meanLength(mask) = accumulator.lengthSum(mask) ./ accumulator.finiteLengthCount(mask);
tbl = table(accumulator.methods, accumulator.coverageHits / M, meanLength, ...
    accumulator.leftMiss / M, accumulator.rightMiss / M, ...
    'VariableNames', {'method', 'coverage', 'mean_interval_length', 'left_miss', 'right_miss'});
end

function slopes = localAttachSlopes(summary, RList, nominal)
slopes = nan(height(summary), 1);
methods = unique(summary.method, "stable");
for iMethod = 1:numel(methods)
    mask = summary.method == methods(iMethod);
    rows = summary(mask, :);
    [~, idx] = ismember(RList(:), rows.R);
    idx = idx(idx > 0);
    if numel(idx) < 2
        continue;
    end
    err = abs(rows.coverage(idx) - nominal);
    valid = err > 0 & isfinite(err);
    if sum(valid) >= 2
        p = polyfit(log(rows.R(idx(valid))), log(err(valid)), 1);
        slopes(mask) = p(1);
    end
end
end

function ci = localCiMethods(Y, alpha, B, seed, truth)
previousState = rng;
rng(double(seed), "twister");
cleanupObj = onCleanup(@() rng(previousState));

Y = double(Y);
if isvector(Y)
    Y = Y(:);
end
R = size(Y, 1);
Ybar = mean(Y, 1);
SR = std(Y, 0, 1);
sqrtR = sqrt(R);

tCrit = localStudentTCrit(1 - alpha / 2, R - 1);
studentLower = Ybar - tCrit .* SR ./ sqrtR;
studentUpper = Ybar + tCrit .* SR ./ sqrtR;

bootIndices = randi(R, R, B);
Ystar = reshape(Y(bootIndices(:), :), R, B, []);
bootMeans = reshape(mean(Ystar, 1), B, []);
[pctLower, pctUpper] = localQuantileRows(bootMeans, alpha);

bootStd = reshape(std(Ystar, 0, 1), B, []);
sampleStd = max(SR, realmin);
Tstar = sqrtR .* (bootMeans - Ybar) ./ bootStd;
zeroStd = bootStd == 0;
if any(zeroStd, "all")
    delta = bootMeans - Ybar;
    signDelta = sign(delta);
    signDelta(signDelta == 0) = 1;
    Tstar(zeroStd) = signDelta(zeroStd) .* Inf;
end
[tLower, tUpper] = localQuantileRows(Tstar, alpha);
btLower = Ybar - SR .* tUpper ./ sqrtR;
btUpper = Ybar - SR .* tLower ./ sqrtR;
minStdRatio = min(bootStd ./ sampleStd, [], 1);

ci = repmat(struct("name", "", "lower", NaN, "upper", NaN, "length", NaN, ...
    "contains", false, "infinite", false, "min_std_ratio", NaN), 3, 1);
ci(1) = localPackCi("Student-t", studentLower, studentUpper, truth, NaN);
ci(2) = localPackCi("Percentile bootstrap", pctLower, pctUpper, truth, NaN);
ci(3) = localPackCi("Bootstrap-t", btLower, btUpper, truth, minStdRatio);
end

function ci = localAppendHybrid(ci, lambda)
p = ci(2);
bt = ci(3);
trigger = bt.infinite || (bt.length > lambda * max(p.length, realmin));
if trigger
    hybrid = p;
else
    hybrid = bt;
end
hybrid.name = "Hybrid";
ci(end + 1) = hybrid;
end

function ci = localPackCi(name, lower, upper, truth, minStdRatio)
lower = lower(1);
upper = upper(1);
if lower > upper
    tmp = lower;
    lower = upper;
    upper = tmp;
end
ci = struct();
ci.name = string(name);
ci.lower = lower;
ci.upper = upper;
ci.length = upper - lower;
ci.infinite = ~(isfinite(lower) && isfinite(upper));
if ci.infinite
    ci.length = Inf;
end
ci.contains = (lower <= truth) && (truth <= upper);
ci.min_std_ratio = minStdRatio(1);
end

function [qLower, qUpper] = localQuantileRows(x, alpha)
xSorted = sort(x, 1);
n = size(xSorted, 1);
lowerIndex = max(1, floor((alpha / 2) * n));
upperIndex = min(n, ceil((1 - alpha / 2) * n));
qLower = xSorted(lowerIndex, :);
qUpper = xSorted(upperIndex, :);
end

function q = localStudentTCrit(p, nu)
if p == 0.5
    q = 0;
    return;
end
tailProb = 2 * min(p, 1 - p);
x = betaincinv(tailProb, nu / 2, 0.5);
q = sign(p - 0.5) * sqrt(nu * (1 / x - 1));
end

function coeffs = localGreenCoefficients(x0, d)
xi = (1:d)' / (d + 1);
coeffs = zeros(d, 1);
for i = 1:d
    if x0 <= xi(i)
        coeffs(i) = x0 * (1 - xi(i));
    else
        coeffs(i) = xi(i) * (1 - x0);
    end
end
end

function estimates = localLinearKernelEstimates(R, nInner, coeffs, yTarget, h)
d = numel(coeffs);
estimates = zeros(R, 1);
for r = 1:R
    theta = randn(nInner, d);
    response = theta * coeffs;
    estimates(r) = mean(localNormalPdf(yTarget, response, h));
end
end

function theta = localDrawTheta(distribution, n, d, config)
switch lower(string(distribution))
    case "normal"
        theta = randn(n, d);
    case "lognormal"
        theta = exp(config.theta_lognormal_mu + config.theta_lognormal_sigma * randn(n, d));
    otherwise
        error("Unknown theta distribution %s.", distribution);
end
end

function response = localLockedNonlinearResponse(theta, config)
z = mean(theta, 2);
e = mean(theta.^2, 2) - 1;
response = config.b + config.c_t + z + config.alpha_model * z.^3 + config.beta_model * z .* e;
end

function estimates = localNonlinearKernelEstimates(R, nInner, distribution, yTarget, h, config)
estimates = zeros(R, 1);
for r = 1:R
    theta = localDrawTheta(distribution, nInner, config.d, config);
    response = localLockedNonlinearResponse(theta, config);
    estimates(r) = mean(localNormalPdf(yTarget, response, h));
end
end

function y = localNormalPdf(x, mu, sigma)
y = exp(-0.5 * ((x - mu) ./ sigma).^2) ./ (sqrt(2 * pi) .* sigma);
end

function q = localQuantiles(x, p)
x = sort(x(:));
n = numel(x);
q = zeros(size(p));
for i = 1:numel(p)
    pos = 1 + (n - 1) * p(i);
    lo = floor(pos);
    hi = ceil(pos);
    if lo == hi
        q(i) = x(lo);
    else
        q(i) = x(lo) + (pos - lo) * (x(hi) - x(lo));
    end
end
end

function value = localIqr(x)
q = localQuantiles(x, [0.25, 0.75]);
value = q(2) - q(1);
end

function s = localSkewness(x)
x = x(:);
xc = x - mean(x);
sd = std(x);
s = mean(xc.^3) / max(sd^3, realmin);
end

function k = localExcessKurtosis(x)
x = x(:);
xc = x - mean(x);
sd = std(x);
k = mean(xc.^4) / max(sd^4, realmin) - 3;
end

function neff = localEffectiveCount(weights)
w = max(weights(:), 0);
if sum(w) <= 0
    neff = 0;
else
    w = w / sum(w);
    neff = 1 / sum(w.^2);
end
end

function [wGrid, referenceCurve, curvePool, statusNote] = localLoadDualPlatePool(projectRoot, modelName, maxFiles)
if nargin < 3 || isempty(maxFiles)
    maxFiles = inf;
end

referencePath = fullfile(projectRoot, "mc_reference_truth_dual_plate_" + modelName + "_n10000_curve.txt");
sampleFiles = dir(fullfile(projectRoot, "rqmc_outer_sample_dual_plate_*.mat"));
if ~isempty(sampleFiles)
    [~, sortIdx] = sort({sampleFiles.name});
    sampleFiles = sampleFiles(sortIdx);
end
wGrid = [];
referenceCurve = [];
curvePool = [];
statusNote = "";

if ~isfile(referencePath)
    statusNote = "Missing reference curve text file: " + string(referencePath);
    return;
end
if isempty(sampleFiles)
    statusNote = "Missing rqmc_outer_sample_dual_plate_*.mat files.";
    return;
end

ref = readmatrix(referencePath);
wGrid = ref(:, 1).';
referenceCurve = ref(:, 2).';
maxFiles = min(double(maxFiles), numel(sampleFiles));
curvePool = nan(maxFiles, numel(wGrid));
loaded = 0;
for iFile = 1:maxFiles
    s = load(fullfile(sampleFiles(iFile).folder, sampleFiles(iFile).name));
    if isfield(s, modelName) && isfield(s.(modelName), "p_rqmc")
        curve = s.(modelName).p_rqmc(:).';
        if numel(curve) == numel(wGrid)
            loaded = loaded + 1;
            curvePool(loaded, :) = curve;
        end
    end
end
curvePool = curvePool(1:loaded, :);
if loaded == 0
    curvePool = [];
    statusNote = "No matching p_rqmc curves in sample MAT files.";
else
    statusNote = sprintf("Loaded %d dual-plate %s curves.", loaded, modelName);
end
end

function [summary, pointwise] = localCurvePoolCoverage(curvePool, referenceCurve, wGrid, config)
nGrid = numel(wGrid);
coverageHits = zeros(3, nGrid);
lengthSums = zeros(3, nGrid);
lengthCounts = zeros(3, nGrid);
degenerateGrid = localDegenerateGridMask(curvePool, referenceCurve);
activeGrid = ~degenerateGrid;

for trial = 1:config.M
    idx = randperm(size(curvePool, 1), config.R);
    Y = curvePool(idx, :);
    ci = localCiMethodsMatrix(Y, config.alpha, config.B, config.seed + trial, referenceCurve, ...
        config.studentize_std_floor_rel);
    for iMethod = 1:3
        coverageHits(iMethod, :) = coverageHits(iMethod, :) + ci(iMethod).contains;
        finiteLength = isfinite(ci(iMethod).length);
        lengthSums(iMethod, finiteLength) = lengthSums(iMethod, finiteLength) + ci(iMethod).length(finiteLength);
        lengthCounts(iMethod, finiteLength) = lengthCounts(iMethod, finiteLength) + 1;
    end
end

meanLengths = nan(3, nGrid);
finiteMask = lengthCounts > 0;
meanLengths(finiteMask) = lengthSums(finiteMask) ./ lengthCounts(finiteMask);

pointwise = table(wGrid(:), referenceCurve(:), ...
    degenerateGrid(:), activeGrid(:), ...
    coverageHits(1,:)' / config.M, meanLengths(1,:)', ...
    coverageHits(2,:)' / config.M, meanLengths(2,:)', ...
    coverageHits(3,:)' / config.M, meanLengths(3,:)', ...
    'VariableNames', {'w', 'reference_curve', 'is_degenerate_grid', 'is_active_grid', ...
    'student_t_coverage', 'student_t_mean_length', ...
    'percentile_bootstrap_coverage', 'percentile_bootstrap_mean_length', ...
    'bootstrap_t_coverage', 'bootstrap_t_mean_length'});

method = ["Student-t"; "Percentile bootstrap"; "Bootstrap-t"];
coverageAll = [mean(pointwise.student_t_coverage); ...
    mean(pointwise.percentile_bootstrap_coverage); ...
    mean(pointwise.bootstrap_t_coverage)];
lengthAll = [mean(pointwise.student_t_mean_length, "omitnan"); ...
    mean(pointwise.percentile_bootstrap_mean_length, "omitnan"); ...
    mean(pointwise.bootstrap_t_mean_length, "omitnan")];
coverageActive = localMeanByMask(pointwise, activeGrid, ...
    ["student_t_coverage", "percentile_bootstrap_coverage", "bootstrap_t_coverage"]);
lengthActive = localMeanByMask(pointwise, activeGrid, ...
    ["student_t_mean_length", "percentile_bootstrap_mean_length", "bootstrap_t_mean_length"]);

summary = table(method, coverageAll, lengthAll, coverageActive, lengthActive, ...
    repmat(sum(activeGrid), 3, 1), repmat(sum(degenerateGrid), 3, 1), ...
    repmat(config.studentize_std_floor_rel, 3, 1), ...
    'VariableNames', {'method', 'mean_pointwise_coverage', 'mean_interval_length', ...
    'active_grid_mean_coverage', 'active_grid_mean_interval_length', ...
    'active_grid_count', 'degenerate_grid_count', 'studentize_std_floor_rel'});
end

function [summary, pointwise] = localSimultaneousBandCoverage(curvePool, referenceCurve, wGrid, config)
nGrid = numel(wGrid);
pointwiseHits = zeros(1, nGrid);
bandHitsByPoint = zeros(1, nGrid);
fullBandHits = 0;
pointwiseLength = zeros(1, nGrid);
bandLength = zeros(1, nGrid);
degenerateGrid = localDegenerateGridMask(curvePool, referenceCurve);
activeGrid = ~degenerateGrid;

for trial = 1:config.M
    idx = randperm(size(curvePool, 1), config.R);
    Y = curvePool(idx, :);
    mu = mean(Y, 1);
    sd = std(Y, 0, 1);
    se = sd / sqrt(config.R);
    seTolerance = localElementTolerance(mu, referenceCurve);
    trialStatGrid = activeGrid & (se > seTolerance);
    tCrit = localStudentTCrit(1 - config.alpha / 2, config.R - 1);
    pwLower = mu - tCrit * se;
    pwUpper = mu + tCrit * se;

    bootIdx = randi(config.R, config.R, config.B);
    maxStats = zeros(config.B, 1);
    for b = 1:config.B
        bootY = Y(bootIdx(:, b), :);
        bootMean = mean(bootY, 1);
        bootSd = std(bootY, 0, 1);
        if any(trialStatGrid)
            seFloor = max(config.studentize_se_floor_rel * se(trialStatGrid), realmin);
            bootSe = max(bootSd(trialStatGrid) / sqrt(config.R), seFloor);
            maxStats(b) = max(abs((bootMean(trialStatGrid) - mu(trialStatGrid)) ./ bootSe));
        else
            maxStats(b) = 0;
        end
    end
    crit = sort(maxStats);
    crit = crit(max(1, ceil((1 - config.alpha) * config.B)));
    bandLower = mu - crit * se;
    bandUpper = mu + crit * se;

    pointwiseHits = pointwiseHits + ((pwLower <= referenceCurve) & (referenceCurve <= pwUpper));
    bandPointHit = (bandLower <= referenceCurve) & (referenceCurve <= bandUpper);
    bandHitsByPoint = bandHitsByPoint + bandPointHit;
    fullBandHits = fullBandHits + all(bandPointHit);
    pointwiseLength = pointwiseLength + (pwUpper - pwLower);
    bandLength = bandLength + (bandUpper - bandLower);
end

pointwise = table(wGrid(:), referenceCurve(:), degenerateGrid(:), activeGrid(:), ...
    pointwiseHits(:) / config.M, bandHitsByPoint(:) / config.M, ...
    pointwiseLength(:) / config.M, bandLength(:) / config.M, ...
    'VariableNames', {'w', 'reference_curve', 'is_degenerate_grid', 'is_active_grid', 'pointwise_coverage', ...
    'simultaneous_band_pointwise_coverage', 'pointwise_mean_length', 'simultaneous_band_mean_length'});

summary = table(string("finite_grid_max_bootstrap_t"), fullBandHits / config.M, ...
    mean(pointwise.pointwise_coverage), mean(pointwise.simultaneous_band_pointwise_coverage), ...
    mean(pointwise.pointwise_coverage(activeGrid), "omitnan"), ...
    mean(pointwise.simultaneous_band_pointwise_coverage(activeGrid), "omitnan"), ...
    mean(pointwise.pointwise_mean_length, "omitnan"), mean(pointwise.simultaneous_band_mean_length, "omitnan"), ...
    mean(pointwise.pointwise_mean_length(activeGrid), "omitnan"), ...
    mean(pointwise.simultaneous_band_mean_length(activeGrid), "omitnan"), ...
    sum(activeGrid), sum(degenerateGrid), config.studentize_se_floor_rel, ...
    'VariableNames', {'method', 'full_grid_simultaneous_coverage', ...
    'mean_pointwise_interval_coverage', 'mean_band_pointwise_coverage', ...
    'active_grid_pointwise_interval_coverage', 'active_grid_band_pointwise_coverage', ...
    'mean_pointwise_interval_length', 'mean_simultaneous_band_length', ...
    'active_grid_pointwise_interval_length', 'active_grid_simultaneous_band_length', ...
    'active_grid_count', 'degenerate_grid_count', 'studentize_se_floor_rel'});
end

function ci = localCiMethodsMatrix(Y, alpha, B, seed, truth, studentizeStdFloorRel)
if nargin < 6 || isempty(studentizeStdFloorRel)
    studentizeStdFloorRel = 0.0;
end

previousState = rng;
rng(double(seed), "twister");
cleanupObj = onCleanup(@() rng(previousState));

[R, nGrid] = size(Y);
mu = mean(Y, 1);
sd = std(Y, 0, 1);
se = sd / sqrt(R);
tCrit = localStudentTCrit(1 - alpha / 2, R - 1);

ci = repmat(struct("contains", false(1, nGrid), "length", nan(1, nGrid)), 3, 1);
lower = mu - tCrit * se;
upper = mu + tCrit * se;
ci(1).contains = (lower <= truth) & (truth <= upper);
ci(1).length = upper - lower;

bootIdx = randi(R, R, B);
bootMeans = zeros(B, nGrid);
bootStd = zeros(B, nGrid);
for b = 1:B
    bootY = Y(bootIdx(:, b), :);
    bootMeans(b, :) = mean(bootY, 1);
    bootStd(b, :) = std(bootY, 0, 1);
end
[pctLower, pctUpper] = localQuantileRows(bootMeans, alpha);
ci(2).contains = (pctLower <= truth) & (truth <= pctUpper);
ci(2).length = pctUpper - pctLower;

stdFloor = max(studentizeStdFloorRel * sd, realmin);
bootStdSafe = max(bootStd, stdFloor);
Tstar = sqrt(R) * (bootMeans - mu) ./ bootStdSafe;
degenerateTolerance = localElementTolerance(mu, truth);
degenerateMask = (sd <= degenerateTolerance) & (abs(mu - truth) <= degenerateTolerance);
Tstar(:, degenerateMask) = 0;
[tLower, tUpper] = localQuantileRows(Tstar, alpha);
btLower = mu - sd .* tUpper / sqrt(R);
btUpper = mu - sd .* tLower / sqrt(R);
btLower(degenerateMask) = mu(degenerateMask);
btUpper(degenerateMask) = mu(degenerateMask);
ci(3).contains = (btLower <= truth) & (truth <= btUpper);
ci(3).length = btUpper - btLower;
end

function mask = localDegenerateGridMask(curvePool, referenceCurve)
poolMean = mean(curvePool, 1);
poolSd = std(curvePool, 0, 1);
toleranceValue = localElementTolerance(poolMean, referenceCurve);
mask = (poolSd <= toleranceValue) & (abs(poolMean - referenceCurve) <= toleranceValue);
end

function values = localMeanByMask(tbl, mask, variableNames)
values = nan(numel(variableNames), 1);
for iVar = 1:numel(variableNames)
    column = tbl.(variableNames(iVar));
    values(iVar) = mean(column(mask), "omitnan");
end
end

function toleranceValue = localElementTolerance(varargin)
scale = zeros(size(varargin{1}));
for iInput = 1:nargin
    scale = max(scale, abs(varargin{iInput}));
end
toleranceValue = max(realmin, 100 * eps(scale));
end

function localSaveFigure(fig, filePath)
try
    exportgraphics(fig, filePath, "Resolution", 180);
catch
    saveas(fig, filePath);
end
close(fig);
end
