function outputs = runSmallSuite(projectRoot, resultsRoot)
%runSmallSuite Generate standardized small outputs for E1--E8.

dpimnumeric.ensureDir(resultsRoot);

experimentNames = [ ...
    "E1_finite_B"; ...
    "E2_R_order"; ...
    "E3_linear_beam_n_h"; ...
    "E4_nonlinear_tail"; ...
    "E5_bootstrap_t_instability"; ...
    "E6_weights_GF_RQMC"; ...
    "E7_plate_SFEM"; ...
    "E8_simultaneous_band"];

outputs = struct();
for iExp = 1:numel(experimentNames)
    expName = experimentNames(iExp);
    outDir = fullfile(resultsRoot, char(expName));
    dpimnumeric.ensureDir(outDir);
    dpimnumeric.ensureDir(fullfile(outDir, "figures"));

    switch expName
        case "E1_finite_B"
            result = localRunE1(projectRoot, outDir);
        case "E2_R_order"
            result = localRunE2(projectRoot, outDir);
        case "E3_linear_beam_n_h"
            result = localRunE3(projectRoot, outDir);
        case "E4_nonlinear_tail"
            result = localRunE4(projectRoot, outDir);
        case "E5_bootstrap_t_instability"
            result = localRunE5(projectRoot, outDir);
        case "E6_weights_GF_RQMC"
            result = localRunE6(projectRoot, outDir);
        case "E7_plate_SFEM"
            result = localRunE7(projectRoot, outDir);
        case "E8_simultaneous_band"
            result = localRunE8(projectRoot, outDir);
        otherwise
            error("Unknown experiment %s.", expName);
    end

    outputs.(matlab.lang.makeValidName(char(expName))) = result;
end
end

function result = localRunE1(projectRoot, outDir)
config = struct("experiment", "E1_finite_B", "alpha", 0.05, ...
    "B_list", [19, 39, 99, 199], "simulation_count", 20000, ...
    "seed", 20260604, "source_reuse", "run_exp1_finite_B_baseline.m and smoke_finite_B_exact_baseline.m");
localWriteConfig(outDir, config);
rng(config.seed, "twister");

nominal = 1 - config.alpha;
BList = config.B_list(:);
nB = numel(BList);
kMinus = zeros(nB, 1);
kPlus = zeros(nB, 1);
c0b = zeros(nB, 1);
gridError = zeros(nB, 1);
simCoverage = zeros(nB, 1);
simAbsError = zeros(nB, 1);

for iB = 1:nB
    B = BList(iB);
    kMinus(iB) = floor((config.alpha / 2) * (B + 1));
    kPlus(iB) = ceil((1 - config.alpha / 2) * (B + 1));
    c0b(iB) = (kPlus(iB) - kMinus(iB)) / (B + 1);
    gridError(iB) = c0b(iB) - nominal;

    ranks = randi(B + 1, config.simulation_count, 1);
    simCoverage(iB) = mean((ranks > kMinus(iB)) & (ranks <= kPlus(iB)));
    simAbsError(iB) = abs(simCoverage(iB) - c0b(iB));
end

summary = table(BList, kMinus, kPlus, c0b, gridError, simCoverage, simAbsError, ...
    repmat(nominal, nB, 1), 'VariableNames', ...
    {'B', 'k_minus', 'k_plus', 'C0B', 'grid_error', ...
    'simulated_coverage', 'simulation_abs_error', 'nominal'});
writetable(summary, fullfile(outDir, "summary.csv"));

raw = struct("summary", summary, "config", config);
save(fullfile(outDir, "raw_results.mat"), "-struct", "raw");

fig = figure("Visible", "off", "Color", "w");
plot(BList, gridError, "-o", "LineWidth", 1.5);
xlabel("B");
ylabel("C0B - nominal");
title("Finite-B grid coverage error", "Interpreter", "none");
grid on;
localSaveFigure(fig, fullfile(outDir, "figures", "finite_B_grid_error.png"));

localWriteStandardFiles(projectRoot, outDir, config, ...
    ["run_exp1_finite_B_baseline.m", "smoke_finite_B_exact_baseline.m"], ...
    "Exact finite-B integerization smoke output.", "completed_small");
result = localResultPaths(outDir);
end

function result = localRunE2(projectRoot, outDir)
config = struct("experiment", "E2_R_order", "alpha", 0.05, "R_list", [8, 16, 32, 64], ...
    "B", 199, "M", 160, "seed", 20260604 + 2, "lambda", 5, ...
    "data_model", "centered exponential X=-log(U)-1; true mean 0", ...
    "source_reuse", "ci_methods.m and compute_hybrid_coverage.m");
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
    "Small centered-exponential coverage diagnostic for R-order interface.", "completed_small");
result = localResultPaths(outDir);
end

function result = localRunE3(projectRoot, outDir)
config = struct("experiment", "E3_linear_beam_n_h", "alpha", 0.05, "d", 10, ...
    "x0", 0.5, "n_list", [16, 32], "h_list", [0.35, 0.70], ...
    "R_list", [8, 16], "B", 99, "M", 80, "seed", 20260604 + 3, ...
    "dy_dimension", 1, "source_reuse", "smoke_linear_gaussian_closed_form.m and run_exp4_beam_benchmark.m");
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
    "Closed-form linear beam smoke output using Green coefficients.", "completed_small");
result = localResultPaths(outDir);
end

function result = localRunE4(projectRoot, outDir)
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
    "Small nonlinear tail output using locked z=mean(theta) model.", "completed_small");
result = localResultPaths(outDir);
end

function result = localRunE5(projectRoot, outDir)
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
    "Small bootstrap-t denominator diagnostic with hybrid fallback.", "completed_small");
result = localResultPaths(outDir);
end

function result = localRunE6(projectRoot, outDir)
config = struct("experiment", "E6_weights_GF_RQMC", "alpha", 0.05, "N", 64, ...
    "d", 3, "R", 10, "B", 99, "M", 80, "seed", 20260604 + 6, ...
    "methods", ["mc_equal", "rqmc_equal", "gf_weight_interface"], ...
    "integrand", "exp(-sum(u.^2,2))", ...
    "source_reuse", "Voronoi_gf.m, generate_rqmc_points.m, voronoi_ci_probability_weights_provider.m");
localWriteConfig(outDir, config);
rng(config.seed, "twister");

truthPoints = rand(50000, config.d);
truth = mean(exp(-sum(truthPoints.^2, 2)));
summary = table();

for iMethod = 1:numel(config.methods)
    methodName = config.methods(iMethod);
    accumulator = localEmptyAccumulator(["Student-t", "Percentile bootstrap", "Bootstrap-t", "Hybrid"]);
    momentRows = table();
    tic;
    for trial = 1:config.M
        estimates = zeros(config.R, 1);
        wLast = [];
        for r = 1:config.R
            [points, weights] = localPointSetWithWeights(methodName, config.N, config.d, config.seed + trial * 100 + r);
            values = exp(-sum(points.^2, 2));
            estimates(r) = sum(weights .* values);
            wLast = weights;
        end
        ci = localCiMethods(estimates, config.alpha, config.B, config.seed + 1000 * trial + iMethod, truth);
        ci = localAppendHybrid(ci, 5);
        accumulator = localAccumulate(accumulator, ci, truth);
        if trial == 1
            momentRows = localWeightMomentTable(methodName, wLast);
        end
    end
    elapsedSeconds = toc;
    block = localAccumulatorTable(accumulator, config.M);
    block.pointset_method = repmat(methodName, height(block), 1);
    block.truth_integral = repmat(truth, height(block), 1);
    block.elapsed_seconds = repmat(elapsedSeconds, height(block), 1);
    block.weight_sum_w2 = repmat(momentRows.sum_w2, height(block), 1);
    block.weight_sum_w3 = repmat(momentRows.sum_w3, height(block), 1);
    block.weight_sum_w4 = repmat(momentRows.sum_w4, height(block), 1);
    block.effective_order2 = repmat(momentRows.effective_order2, height(block), 1);
    block.effective_order3 = repmat(momentRows.effective_order3, height(block), 1);
    block.effective_order4 = repmat(momentRows.effective_order4, height(block), 1);
    summary = [summary; block]; %#ok<AGROW>
end

summary = movevars(summary, ["pointset_method", "truth_integral", "elapsed_seconds", ...
    "weight_sum_w2", "weight_sum_w3", "weight_sum_w4", ...
    "effective_order2", "effective_order3", "effective_order4"], "Before", 1);
writetable(summary, fullfile(outDir, "summary.csv"));

raw = struct("summary", summary, "config", config);
save(fullfile(outDir, "raw_results.mat"), "-struct", "raw");

fig = figure("Visible", "off", "Color", "w");
bar(categorical(summary.pointset_method(summary.method == "Student-t")), ...
    summary.effective_order2(summary.method == "Student-t"));
ylabel("effective order 2");
title("Small probability-weight diagnostics", "Interpreter", "none");
grid on;
localSaveFigure(fig, fullfile(outDir, "figures", "weight_effective_orders.png"));

localWriteStandardFiles(projectRoot, outDir, config, ["Voronoi_gf.m", "generate_rqmc_points.m", "voronoi_ci_probability_weights_provider.m"], ...
    "Small pointset/weight interface. GF row is a stress-test surrogate unless replaced by project GF optimizer.", "interface_small");
result = localResultPaths(outDir);
end

function result = localRunE7(projectRoot, outDir)
config = struct("experiment", "E7_plate_SFEM", "alpha", 0.05, "R", 10, ...
    "B", 99, "M", 120, "seed", 20260604 + 7, "model", "dual_plate_thin", ...
    "studentize_std_floor_rel", 0.02, ...
    "source_reuse", "rqmc_outer_sample_dual_plate_*.mat and mc_reference_truth_dual_plate_thin_n10000_curve.txt");
localWriteConfig(outDir, config);
rng(config.seed, "twister");

[wGrid, referenceCurve, curvePool, statusNote] = localLoadDualPlatePool(projectRoot, "thin");

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
title("Dual-plate small SFEM coverage", "Interpreter", "none");
grid on;
localSaveFigure(fig, fullfile(outDir, "figures", "plate_sfem_coverage.png"));

localWriteStandardFiles(projectRoot, outDir, config, ["run_exp5_plate_benchmark.m", "run_dual_plate_rqmc_sample_pool_ci_driver.m"], ...
    "Small dual-plate SFEM output when precomputed pool exists.", "completed_small_or_interface");
result = localResultPaths(outDir);
end

function result = localRunE8(projectRoot, outDir)
config = struct("experiment", "E8_simultaneous_band", "alpha", 0.05, ...
    "R", 10, "B", 99, "M", 120, "seed", 20260604 + 8, ...
    "studentize_se_floor_rel", 0.02, ...
    "source_reuse", "new finite-grid simultaneous-band wrapper over existing dual-plate pool");
localWriteConfig(outDir, config);
rng(config.seed, "twister");

[wGrid, referenceCurve, curvePool, statusNote] = localLoadDualPlatePool(projectRoot, "thin");
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
title("Finite-grid simultaneous band smoke output", "Interpreter", "none");
grid on;
localSaveFigure(fig, fullfile(outDir, "figures", "simultaneous_band_coverage.png"));

localWriteStandardFiles(projectRoot, outDir, config, ["run_rqmc_sample_pool_ci_driver.m"], ...
    "Small finite-grid simultaneous-band interface over curve pool.", "completed_small_or_interface");
result = localResultPaths(outDir);
end

function localWriteConfig(outDir, config)
dpimnumeric.writeJson(fullfile(outDir, "config.json"), config);
end

function localWriteStandardFiles(projectRoot, outDir, config, dependencies, note, status)
metadata = struct();
metadata.run_time = char(datetime("now"));
metadata.matlab_version = version;
metadata.random_seed = config.seed;
metadata.git_commit = "";
metadata.git_status = "not_a_git_repository_or_not_checked";
metadata.file_hashes = localWrapperHashes(projectRoot);
metadata.model_formula = "g(theta,t)=b+c_t*t+z+alpha*z^3+beta*z*e";
metadata.model_z = "mean(theta)";
metadata.model_e = "mean(theta.^2)-1";
metadata.parameters = config;
metadata.script_entry = which("run_dpim_ci_numeric_rebuild_small");
metadata.dependencies = dependencies;
metadata.status = status;
metadata.note = note;
dpimnumeric.writeJson(fullfile(outDir, "metadata.json"), metadata);

logLines = [
    "DPIM CI small rebuild log"
    "status: " + string(status)
    "time: " + string(datetime("now"))
    "note: " + string(note)
    "entry: run_dpim_ci_numeric_rebuild_small('small')"
    ""
];
dpimnumeric.writeText(fullfile(outDir, "run_log.txt"), strjoin(logLines, newline));

readmeLines = [
    "# " + string(config.experiment)
    ""
    "This directory was generated by the small DPIM CI rebuild wrapper."
    ""
    "- `config.json`: small-run configuration."
    "- `metadata.json`: MATLAB version, seed, wrapper hashes, model formula, and dependencies."
    "- `raw_results.mat`: raw MATLAB data used to make the summary."
    "- `summary.csv`: main tabular output."
    "- `figures/*.png`: figures regenerated from data."
    ""
    "Status: " + string(status)
    ""
    "Note: " + string(note)
    ""
    "Small outputs are not final full-scale paper evidence."
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
files = ["run_dpim_ci_numeric_rebuild_small.m", "+dpimnumeric/runSmallSuite.m", ...
    "+dpimnumeric/buildInventory.m", "+dpimnumeric/writeModelLock.m"];
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

function [points, weights] = localPointSetWithWeights(methodName, N, d, seed)
rng(seed, "twister");
switch lower(string(methodName))
    case "mc_equal"
        points = rand(N, d);
        weights = ones(N, 1) / N;
    case "rqmc_equal"
        points = localSobolPoints(N, d, seed);
        weights = ones(N, 1) / N;
    case "gf_weight_interface"
        points = localSobolPoints(N, d, seed);
        rawWeights = exp(-5 * points(:, 1));
        weights = rawWeights / sum(rawWeights);
    otherwise
        error("Unknown pointset method %s.", methodName);
end
end

function points = localSobolPoints(N, d, seed)
if exist("sobolset", "file") == 2
    p = sobolset(d);
    p = scramble(p, "MatousekAffineOwen");
    points = net(p, N);
    shift = rand(1, d);
    points = mod(points + shift, 1);
else
    points = localLatinPoints(N, d, seed);
end
end

function points = localLatinPoints(N, d, seed)
rng(seed, "twister");
points = zeros(N, d);
for j = 1:d
    points(:, j) = ((randperm(N)' - rand(N, 1)) / N);
end
end

function tbl = localWeightMomentTable(methodName, weights)
w = weights(:);
sumW2 = sum(w.^2);
sumW3 = sum(w.^3);
sumW4 = sum(w.^4);
tbl = table(string(methodName), sumW2, sumW3, sumW4, ...
    sumW2^(-1), sumW3^(-1/2), sumW4^(-1/3), ...
    'VariableNames', {'pointset_method', 'sum_w2', 'sum_w3', 'sum_w4', ...
    'effective_order2', 'effective_order3', 'effective_order4'});
end

function [wGrid, referenceCurve, curvePool, statusNote] = localLoadDualPlatePool(projectRoot, modelName)
referencePath = fullfile(projectRoot, "mc_reference_truth_dual_plate_" + modelName + "_n10000_curve.txt");
sampleFiles = dir(fullfile(projectRoot, "rqmc_outer_sample_dual_plate_*.mat"));
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
maxFiles = min(80, numel(sampleFiles));
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
