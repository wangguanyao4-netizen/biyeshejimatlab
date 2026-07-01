function outputs = runHFallbackSweepSuite(projectRoot, resultsRoot, runMode, startAt)
%runHFallbackSweepSuite Run h-sweep diagnostics for Bootstrap-t fallback rules.

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

experimentNames = [ ...
    "E1_finite_B"; ...
    "E2_R_order"; ...
    "E3_linear_beam_n_h"; ...
    "E4_nonlinear_tail"; ...
    "E5_bootstrap_t_instability"; ...
    "E7_plate_SFEM"; ...
    "E8_simultaneous_band"];

if strlength(string(startAt)) == 0
    startIndex = 1;
else
    startIndex = find(experimentNames == string(startAt), 1, "first");
    if isempty(startIndex)
        error("Unknown startAt experiment %s.", string(startAt));
    end
end

outputs = struct();
for iExp = startIndex:numel(experimentNames)
    expName = experimentNames(iExp);
    outDir = fullfile(resultsRoot, char(expName));
    dpimnumeric.ensureDir(outDir);
    dpimnumeric.ensureDir(fullfile(outDir, "figures"));

    switch expName
        case "E1_finite_B"
            result = localRunE1(outDir, runMode);
        case "E2_R_order"
            result = localRunE2(outDir, runMode);
        case "E3_linear_beam_n_h"
            result = localRunE3(outDir, runMode);
        case "E4_nonlinear_tail"
            result = localRunE4(outDir, runMode);
        case "E5_bootstrap_t_instability"
            result = localRunE5(outDir, runMode);
        case "E7_plate_SFEM"
            result = localRunE7(outDir, runMode);
        case "E8_simultaneous_band"
            result = localRunE8(outDir, runMode);
    end

    outputs.(matlab.lang.makeValidName(char(expName))) = result;
end

manifest = struct();
manifest.project_root = projectRoot;
manifest.results_root = resultsRoot;
manifest.run_mode = char(runMode);
manifest.experiments = cellstr(experimentNames);
manifest.note = "E1 finite-B baseline plus 20-point h sweep for E2/E3/E4/E5/E7/E8 with raw Bootstrap-t and Bootstrap-t fallback-rule reported separately.";
dpimnumeric.writeJson(fullfile(resultsRoot, "hfallback_sweep_manifest.json"), manifest);
end

function result = localRunE1(outDir, runMode)
switch runMode
    case "full"
        config = struct("experiment", "E1_finite_B", "alpha", 0.05, ...
            "B_list", [399, 999, 4999, 9999], "simulation_count", 200000, ...
            "seed", 20260604, "note", "Finite-B integerization baseline; no h or fallback applies.");
    case "medium"
        config = struct("experiment", "E1_finite_B", "alpha", 0.05, ...
            "B_list", [99, 199, 399, 999], "simulation_count", 80000, ...
            "seed", 20260604, "note", "Finite-B integerization baseline; no h or fallback applies.");
    otherwise
        config = struct("experiment", "E1_finite_B", "alpha", 0.05, ...
            "B_list", [19, 39, 99, 199], "simulation_count", 20000, ...
            "seed", 20260604, "note", "Finite-B integerization baseline; no h or fallback applies.");
end
localWriteConfig(outDir, config);
rng(config.seed, "twister");

[summary, finiteBSimulationBlocks] = localFiniteBSummary(config);
localWriteSummary(outDir, summary, config, "E1 finite-B baseline; h-sweep and Bootstrap-t fallback are not defined here.");
save(fullfile(outDir, "finite_B_rank_simulations.mat"), "finiteBSimulationBlocks", "config", "-v7.3");
plotSummary = summary(summary.CI_method == "Percentile bootstrap finite-B rank", :);

fig = figure("Visible", "off", "Color", "w");
tiledlayout(1, 2, "TileSpacing", "compact");
nexttile;
plot(plotSummary.B, plotSummary.coverage_error, "-o", "LineWidth", 1.5);
xlabel("B");
ylabel("C0B - nominal");
title("Finite-B grid error", "Interpreter", "none");
grid on;
nexttile;
errorbar(plotSummary.B, plotSummary.simulated_coverage, plotSummary.simulation_ci_half_width_95, "-o", "LineWidth", 1.2);
hold on;
plot(plotSummary.B, plotSummary.C0B, "--", "LineWidth", 1.2);
hold off;
xlabel("B");
ylabel("coverage");
title("Simulation vs exact C0B", "Interpreter", "none");
grid on;
try
    exportgraphics(fig, fullfile(outDir, "figures", "E1_finite_B_baseline.png"), "Resolution", 180);
catch
    saveas(fig, fullfile(outDir, "figures", "E1_finite_B_baseline.png"));
end
close(fig);

result = localResultPaths(outDir);
end

function result = localRunE2(outDir, runMode)
base = localBaseConfig(runMode);
config = struct( ...
    "experiment", "E2_R_order", ...
    "alpha", 0.05, ...
    "lambda", 5, ...
    "R_list", base.e2_R_list, ...
    "B", base.B, ...
    "M", base.M, ...
    "seed", 20260604 + 20, ...
    "h_list", localLogHList(0.005, 1.0, 20), ...
    "y_target", 0.0, ...
    "data_model", "centered exponential X=-log(U)-1; smoothed density at y=0");
localWriteConfig(outDir, config);
rng(config.seed, "twister");

summary = table();
detailBlocks = {};
detailBlockIndex = 0;
for iH = 1:numel(config.h_list)
    h = config.h_list(iH);
    truth = localCenteredExponentialSmoothedPdf(config.y_target, h);
    for iR = 1:numel(config.R_list)
        R = config.R_list(iR);
        acc = localEmptyAccumulator(localScalarMethods());
        trialDetails = localEmptyScalarTrialDetails(config.M, R);
        for trial = 1:config.M
            x = -log(rand(R, 1)) - 1;
            estimates = localNormalPdf(config.y_target, x, h);
            ciSeed = config.seed + 1000 * trial + iH + R;
            ci = localCiMethods(estimates, config.alpha, config.B, ciSeed, truth);
            ci = localAppendFallbackRule(ci, config.lambda);
            acc = localAccumulate(acc, ci, truth);
            trialDetails = localFillScalarTrialDetails(trialDetails, trial, ciSeed, truth, estimates, ci);
        end
        block = localAccumulatorTable(acc, config.M, config.alpha);
        block.h_index = repmat(iH, height(block), 1);
        block.h = repmat(h, height(block), 1);
        block.R = repmat(R, height(block), 1);
        block.truth_smoothed_density = repmat(truth, height(block), 1);
        summary = [summary; block]; %#ok<AGROW>
        detailBlockIndex = detailBlockIndex + 1;
        detailBlocks{detailBlockIndex, 1} = localPackScalarDetailBlock( ...
            struct("experiment", "E2_R_order", "h_index", iH, "h", h, "R", R, ...
            "truth_smoothed_density", truth), trialDetails); %#ok<AGROW>
    end
end
summary = movevars(summary, ["h_index", "h", "R", "truth_smoothed_density"], "Before", 1);
localWriteSummary(outDir, summary, config, "E2 centered-exponential h/fallback sweep.");
localSaveDetailBlocks(outDir, detailBlocks, config);
localPlotCoverageByH(outDir, summary, "E2 coverage by h", "E2_hfallback_coverage.png");
result = localResultPaths(outDir);
end

function result = localRunE3(outDir, runMode)
base = localBaseConfig(runMode);
config = struct( ...
    "experiment", "E3_linear_beam_n_h", ...
    "alpha", 0.05, ...
    "lambda", 5, ...
    "d", 10, ...
    "x0", 0.5, ...
    "n_list", base.e3_n_list, ...
    "R_list", base.e3_R_list, ...
    "B", base.B, ...
    "M", base.M, ...
    "seed", 20260604 + 30, ...
    "h_list", localLogHList(0.005, 0.75, 20), ...
    "y_target", 0.0);
localWriteConfig(outDir, config);
rng(config.seed, "twister");

coeffs = localGreenCoefficients(config.x0, config.d);
responseStd = norm(coeffs, 2);
summary = table();
detailBlocks = {};
detailBlockIndex = 0;
for iH = 1:numel(config.h_list)
    h = config.h_list(iH);
    truth = localNormalPdf(config.y_target, 0, sqrt(responseStd ^ 2 + h ^ 2));
    for iN = 1:numel(config.n_list)
        nInner = config.n_list(iN);
        for iR = 1:numel(config.R_list)
            R = config.R_list(iR);
            acc = localEmptyAccumulator(localScalarMethods());
            trialDetails = localEmptyScalarTrialDetails(config.M, R);
            for trial = 1:config.M
                estimates = localLinearKernelEstimates(R, nInner, coeffs, config.y_target, h);
                ciSeed = config.seed + 1000 * trial + iH + R;
                ci = localCiMethods(estimates, config.alpha, config.B, ciSeed, truth);
                ci = localAppendFallbackRule(ci, config.lambda);
                acc = localAccumulate(acc, ci, truth);
                trialDetails = localFillScalarTrialDetails(trialDetails, trial, ciSeed, truth, estimates, ci);
            end
            block = localAccumulatorTable(acc, config.M, config.alpha);
            block.h_index = repmat(iH, height(block), 1);
            block.h = repmat(h, height(block), 1);
            block.n = repmat(nInner, height(block), 1);
            block.R = repmat(R, height(block), 1);
            block.Rnh = repmat(R * nInner * h, height(block), 1);
            block.truth_smoothed_density = repmat(truth, height(block), 1);
            summary = [summary; block]; %#ok<AGROW>
            detailBlockIndex = detailBlockIndex + 1;
            detailBlocks{detailBlockIndex, 1} = localPackScalarDetailBlock( ...
                struct("experiment", "E3_linear_beam_n_h", "h_index", iH, "h", h, ...
                "n", nInner, "R", R, "Rnh", R * nInner * h, ...
                "truth_smoothed_density", truth), trialDetails); %#ok<AGROW>
        end
    end
end
summary = movevars(summary, ["h_index", "h", "n", "R", "Rnh", "truth_smoothed_density"], "Before", 1);
localWriteSummary(outDir, summary, config, "E3 linear beam h/fallback sweep.");
localSaveDetailBlocks(outDir, detailBlocks, config);
localPlotCoverageByH(outDir, summary, "E3 coverage by h", "E3_hfallback_coverage.png");
result = localResultPaths(outDir);
end

function result = localRunE4(outDir, runMode)
base = localBaseConfig(runMode);
cases = localE4Cases(runMode, base);
summary = table();
detailBlocks = {};
detailBlockIndex = 0;
for iCase = 1:numel(cases)
    cfg = cases(iCase);
    rng(cfg.seed, "twister");
    responseTruth = localLockedLognormalTruth(cfg.d, cfg.truth_sample_count, cfg.seed + 17);
    yValues = localQuantiles(responseTruth, linspace(cfg.q_low, cfg.q_high, cfg.y_point_count));
    hList = localLogHList(cfg.h_min, cfg.h_max, 20);
    for iH = 1:numel(hList)
        h = hList(iH);
        for iN = 1:numel(cfg.n_list)
            nInner = cfg.n_list(iN);
            for iy = 1:numel(yValues)
                yTarget = yValues(iy);
                truth = mean(localNormalPdf(yTarget, responseTruth, h));
                acc = localEmptyAccumulator(localScalarMethods());
                trialDetails = localEmptyScalarTrialDetails(cfg.M, cfg.R);
                for trial = 1:cfg.M
                    estimates = localNonlinearKernelEstimates(cfg.R, nInner, cfg.d, yTarget, h);
                    ciSeed = cfg.seed + 1000 * trial + 10 * iH + iy;
                    ci = localCiMethods(estimates, cfg.alpha, cfg.B, ciSeed, truth);
                    ci = localAppendFallbackRule(ci, cfg.lambda);
                    acc = localAccumulate(acc, ci, truth);
                    trialDetails = localFillScalarTrialDetails(trialDetails, trial, ciSeed, truth, estimates, ci);
                end
                block = localAccumulatorTable(acc, cfg.M, cfg.alpha);
                block.case_id = repmat(string(cfg.case_id), height(block), 1);
                block.d = repmat(cfg.d, height(block), 1);
                block.h_index = repmat(iH, height(block), 1);
                block.h = repmat(h, height(block), 1);
                block.n = repmat(nInner, height(block), 1);
                block.y_index = repmat(iy, height(block), 1);
                block.y_value = repmat(yTarget, height(block), 1);
                block.truth_smoothed_density = repmat(truth, height(block), 1);
                summary = [summary; block]; %#ok<AGROW>
                detailBlockIndex = detailBlockIndex + 1;
                detailBlocks{detailBlockIndex, 1} = localPackScalarDetailBlock( ...
                    struct("experiment", "E4_nonlinear_tail", "case_id", string(cfg.case_id), ...
                    "d", cfg.d, "h_index", iH, "h", h, "n", nInner, ...
                    "y_index", iy, "y_value", yTarget, "truth_smoothed_density", truth), ...
                    trialDetails); %#ok<AGROW>
            end
        end
    end
end
summary = movevars(summary, ["case_id", "d", "h_index", "h", "n", "y_index", "y_value", ...
    "truth_smoothed_density"], "Before", 1);
config = struct("experiment", "E4_nonlinear_tail", "run_mode", char(runMode), ...
    "cases", cases, "note", "Direct h/fallback sweep; h is varied explicitly instead of using h=dy.");
localWriteSummary(outDir, summary, config, "E4 nonlinear h/fallback sweep.");
localSaveDetailBlocks(outDir, detailBlocks, config);
localPlotCoverageByH(outDir, summary, "E4 coverage by h", "E4_hfallback_coverage.png");
result = localResultPaths(outDir);
end

function result = localRunE5(outDir, runMode)
base = localBaseConfig(runMode);
targetDimension = 5;
responseTruth = localLockedLognormalTruth(targetDimension, base.truth_sample_count, 20260604 + 51);
q = localQuantiles(responseTruth, [0.50, 0.95, 0.999]);
cases = struct( ...
    "case_id", {"core_q50_q95", "righttail_q95_q999"}, ...
    "q_low", {0.50, 0.95}, ...
    "q_high", {0.95, 0.999}, ...
    "y_min", {q(1), q(2)}, ...
    "y_max", {q(2), q(3)});

summary = table();
detailBlocks = {};
detailBlockIndex = 0;
for iCase = 1:numel(cases)
    cfg = cases(iCase);
    yValues = linspace(cfg.y_min, cfg.y_max, base.e5_y_point_count);
    hList = localLogHList(base.e5_h_min, base.e5_h_max, 20);
    for iH = 1:numel(hList)
        h = hList(iH);
        for iN = 1:numel(base.e5_n_list)
            nInner = base.e5_n_list(iN);
            for iy = 1:numel(yValues)
                yTarget = yValues(iy);
                truth = mean(localNormalPdf(yTarget, responseTruth, h));
                acc = localEmptyAccumulator(localScalarMethods());
                trialDetails = localEmptyScalarTrialDetails(base.M, base.e5_R);
                for trial = 1:base.M
                    estimates = localNonlinearKernelEstimates(base.e5_R, nInner, targetDimension, yTarget, h);
                    ciSeed = 20260604 + 5000 + 1000 * trial + iH + iy;
                    ci = localCiMethods(estimates, base.alpha, base.B, ciSeed, truth);
                    ci = localAppendFallbackRule(ci, base.lambda);
                    acc = localAccumulate(acc, ci, truth);
                    trialDetails = localFillScalarTrialDetails(trialDetails, trial, ciSeed, truth, estimates, ci);
                end
                block = localAccumulatorTable(acc, base.M, base.alpha);
                block.case_id = repmat(string(cfg.case_id), height(block), 1);
                block.d = repmat(targetDimension, height(block), 1);
                block.h_index = repmat(iH, height(block), 1);
                block.h = repmat(h, height(block), 1);
                block.n = repmat(nInner, height(block), 1);
                block.y_index = repmat(iy, height(block), 1);
                block.y_value = repmat(yTarget, height(block), 1);
                block.truth_smoothed_density = repmat(truth, height(block), 1);
                summary = [summary; block]; %#ok<AGROW>
                detailBlockIndex = detailBlockIndex + 1;
                detailBlocks{detailBlockIndex, 1} = localPackScalarDetailBlock( ...
                    struct("experiment", "E5_bootstrap_t_instability", "case_id", string(cfg.case_id), ...
                    "d", targetDimension, "h_index", iH, "h", h, "n", nInner, ...
                    "y_index", iy, "y_value", yTarget, "truth_smoothed_density", truth), ...
                    trialDetails); %#ok<AGROW>
            end
        end
    end
end
summary = movevars(summary, ["case_id", "d", "h_index", "h", "n", "y_index", "y_value", ...
    "truth_smoothed_density"], "Before", 1);
config = struct("experiment", "E5_bootstrap_t_instability", "run_mode", char(runMode), ...
    "d", targetDimension, "lambda", base.lambda, "cases", cases, ...
    "note", "E5 is retained because it remains the targeted core/right-tail instability diagnostic.");
localWriteSummary(outDir, summary, config, "E5 targeted nonlinear h/fallback sweep.");
localSaveDetailBlocks(outDir, detailBlocks, config);
localPlotCoverageByH(outDir, summary, "E5 coverage by h", "E5_hfallback_coverage.png");
result = localResultPaths(outDir);
end

function result = localRunE7(outDir, runMode)
base = localBaseConfig(runMode);
config = localPlateConfig(runMode, "E7_plate_SFEM");
localWriteConfig(outDir, config);
rng(config.seed, "twister");

[wGrid, referenceResponses, curveResponses] = localBuildIndependentPlateResponses(config);
hList = localLogHList(config.h_min, config.h_max, 20);
summary = table();
pointwise = table();
diagnosticBlocks = {};
for iH = 1:numel(hList)
    h = hList(iH);
    referenceCurve = localKernelCurve(referenceResponses, wGrid, h);
    curvePool = localKernelCurvePool(curveResponses, wGrid, h);
    [hSummary, hPointwise, hDiagnostics] = localCurvePoolCoverageWithFallback(curvePool, referenceCurve, wGrid, ...
        config.alpha, config.B, config.M, config.R, config.seed + iH, base.lambda);
    hSummary.h_index = repmat(iH, height(hSummary), 1);
    hSummary.h = repmat(h, height(hSummary), 1);
    hPointwise.h_index = repmat(iH, height(hPointwise), 1);
    hPointwise.h = repmat(h, height(hPointwise), 1);
    summary = [summary; hSummary]; %#ok<AGROW>
    pointwise = [pointwise; hPointwise]; %#ok<AGROW>
    hDiagnostics.h_index = iH;
    hDiagnostics.h = h;
    diagnosticBlocks{iH, 1} = hDiagnostics; %#ok<AGROW>
end
summary = movevars(summary, ["h_index", "h"], "Before", 1);
pointwise = movevars(pointwise, ["h_index", "h"], "Before", 1);
writetable(pointwise, fullfile(outDir, "pointwise_by_h.csv"));
save(fullfile(outDir, "plate_generated_inputs.mat"), "wGrid", "referenceResponses", ...
    "curveResponses", "hList", "config", "-v7.3");
save(fullfile(outDir, "curve_coverage_diagnostics.mat"), "diagnosticBlocks", "config", "-v7.3");
localWriteSummary(outDir, summary, config, "E7 independent RFEM pointwise h/fallback sweep.");
localPlotCoverageByH(outDir, summary, "E7 coverage by h", "E7_hfallback_coverage.png");
result = localResultPaths(outDir);
end

function result = localRunE8(outDir, runMode)
base = localBaseConfig(runMode);
config = localPlateConfig(runMode, "E8_simultaneous_band");
localWriteConfig(outDir, config);
rng(config.seed + 100, "twister");

[wGrid, referenceResponses, curveResponses] = localBuildIndependentPlateResponses(config);
hList = localLogHList(config.h_min, config.h_max, 20);
summary = table();
pointwise = table();
diagnosticBlocks = {};
for iH = 1:numel(hList)
    h = hList(iH);
    referenceCurve = localKernelCurve(referenceResponses, wGrid, h);
    curvePool = localKernelCurvePool(curveResponses, wGrid, h);
    [hSummary, hPointwise, hDiagnostics] = localSimultaneousCoverageWithFallback(curvePool, referenceCurve, wGrid, ...
        config.alpha, config.B, config.M, config.R, config.seed + iH, base.lambda);
    hSummary.h_index = repmat(iH, height(hSummary), 1);
    hSummary.h = repmat(h, height(hSummary), 1);
    hPointwise.h_index = repmat(iH, height(hPointwise), 1);
    hPointwise.h = repmat(h, height(hPointwise), 1);
    summary = [summary; hSummary]; %#ok<AGROW>
    pointwise = [pointwise; hPointwise]; %#ok<AGROW>
    hDiagnostics.h_index = iH;
    hDiagnostics.h = h;
    diagnosticBlocks{iH, 1} = hDiagnostics; %#ok<AGROW>
end
summary = movevars(summary, ["h_index", "h"], "Before", 1);
pointwise = movevars(pointwise, ["h_index", "h"], "Before", 1);
writetable(pointwise, fullfile(outDir, "pointwise_band_by_h.csv"));
save(fullfile(outDir, "plate_generated_inputs.mat"), "wGrid", "referenceResponses", ...
    "curveResponses", "hList", "config", "-v7.3");
save(fullfile(outDir, "curve_coverage_diagnostics.mat"), "diagnosticBlocks", "config", "-v7.3");
localWriteSummary(outDir, summary, config, "E8 independent RFEM simultaneous-band h/fallback sweep.");
localPlotCoverageByH(outDir, summary, "E8 coverage by h", "E8_hfallback_coverage.png");
result = localResultPaths(outDir);
end

function base = localBaseConfig(runMode)
base = struct();
base.alpha = 0.05;
base.lambda = 5;
base.truth_sample_count = 50000;
switch runMode
    case "full"
        base.B = 399;
        base.M = 300;
        base.e2_R_list = [10, 20, 40, 80, 160, 320];
        base.e3_n_list = [64, 128, 256, 512, 1024];
        base.e3_R_list = [10, 20, 40];
        base.e4_y_point_count = 20;
        base.e5_y_point_count = 20;
        base.e5_n_list = 2 .^ (6:10);
        base.e5_R = 20;
        base.e5_h_min = 0.02;
        base.e5_h_max = 50;
    case "medium"
        base.B = 199;
        base.M = 120;
        base.e2_R_list = [10, 20, 40, 80, 160];
        base.e3_n_list = [32, 64, 128, 256];
        base.e3_R_list = [10, 20];
        base.e4_y_point_count = 12;
        base.e5_y_point_count = 12;
        base.e5_n_list = 2 .^ (6:9);
        base.e5_R = 12;
        base.e5_h_min = 0.02;
        base.e5_h_max = 50;
    otherwise
        base.B = 79;
        base.M = 40;
        base.e2_R_list = [8, 16, 32, 64];
        base.e3_n_list = [16, 32];
        base.e3_R_list = [8, 16];
        base.e4_y_point_count = 6;
        base.e5_y_point_count = 6;
        base.e5_n_list = 2 .^ (5:7);
        base.e5_R = 8;
        base.e5_h_min = 0.02;
        base.e5_h_max = 50;
end
end

function cases = localE4Cases(runMode, base)
switch runMode
    case "full"
        cases = struct( ...
            "case_id", {"d5_q001_q999", "d10_q001_q999"}, ...
            "d", {5, 10}, ...
            "n_list", {2 .^ (6:10), 2 .^ (6:10)}, ...
            "R", {20, 20}, ...
            "M", {300, 400}, ...
            "B", {399, 999}, ...
            "alpha", {base.alpha, base.alpha}, ...
            "lambda", {base.lambda, base.lambda}, ...
            "q_low", {0.001, 0.001}, ...
            "q_high", {0.999, 0.999}, ...
            "truth_sample_count", {base.truth_sample_count, base.truth_sample_count}, ...
            "y_point_count", {base.e4_y_point_count, base.e4_y_point_count}, ...
            "h_min", {0.02, 0.02}, ...
            "h_max", {50, 50}, ...
            "seed", {20260604 + 41, 20260604 + 42});
    otherwise
        cases = struct( ...
            "case_id", {"d5_q001_q999"}, ...
            "d", {5}, ...
            "n_list", {base.e5_n_list}, ...
            "R", {base.e5_R}, ...
            "M", {base.M}, ...
            "B", {base.B}, ...
            "alpha", {base.alpha}, ...
            "lambda", {base.lambda}, ...
            "q_low", {0.001}, ...
            "q_high", {0.999}, ...
            "truth_sample_count", {base.truth_sample_count}, ...
            "y_point_count", {base.e4_y_point_count}, ...
            "h_min", {0.02}, ...
            "h_max", {50}, ...
            "seed", {20260604 + 41});
end
end

function config = localPlateConfig(runMode, experimentName)
switch runMode
    case "full"
        config = struct("experiment", experimentName, "alpha", 0.05, "lambda", 5, ...
            "R", 20, "M", 300, "B", 399, "seed", 20260604 + 70, ...
            "mesh_nx", 4, "mesh_ny", 4, "d", 10, "curve_pool_count", 320, ...
            "inner_sample_count", 80, "reference_sample_count", 1200, "w_grid_count", 80, ...
            "h_min", 0.005, "h_max", 5.0, "E0", 2.184e11, "sigma_E", 0.10, ...
            "nu", 0.30, "t", 0.005, "q0", -1.0e5, "source_reuse", ...
            "standalone randomE -> Morley/Kirchhoff FE wc -> kernel probability curves");
    case "medium"
        config = struct("experiment", experimentName, "alpha", 0.05, "lambda", 5, ...
            "R", 12, "M", 120, "B", 199, "seed", 20260604 + 70, ...
            "mesh_nx", 4, "mesh_ny", 4, "d", 10, "curve_pool_count", 180, ...
            "inner_sample_count", 50, "reference_sample_count", 700, "w_grid_count", 60, ...
            "h_min", 0.005, "h_max", 5.0, "E0", 2.184e11, "sigma_E", 0.10, ...
            "nu", 0.30, "t", 0.005, "q0", -1.0e5, "source_reuse", ...
            "standalone randomE -> Morley/Kirchhoff FE wc -> kernel probability curves");
    otherwise
        config = struct("experiment", experimentName, "alpha", 0.05, "lambda", 5, ...
            "R", 8, "M", 40, "B", 79, "seed", 20260604 + 70, ...
            "mesh_nx", 3, "mesh_ny", 3, "d", 10, "curve_pool_count", 80, ...
            "inner_sample_count", 24, "reference_sample_count", 250, "w_grid_count", 40, ...
            "h_min", 0.005, "h_max", 5.0, "E0", 2.184e11, "sigma_E", 0.10, ...
            "nu", 0.30, "t", 0.005, "q0", -1.0e5, "source_reuse", ...
            "standalone randomE -> Morley/Kirchhoff FE wc -> kernel probability curves");
end
end

function methods = localScalarMethods()
methods = ["Student-t"; "Percentile bootstrap"; "Bootstrap-t"; "Bootstrap-t fallback-rule"];
end

function hList = localLogHList(hMin, hMax, count)
hList = logspace(log10(hMin), log10(hMax), count);
hList = sort(hList(:).');
end

function localWriteConfig(outDir, config)
dpimnumeric.writeJson(fullfile(outDir, "config.json"), config);
end

function localWriteSummary(outDir, summary, config, note)
writetable(summary, fullfile(outDir, "summary.csv"));
raw = struct("summary", summary, "config", config);
save(fullfile(outDir, "raw_results.mat"), "-struct", "raw", "-v7.3");
readmeLines = [
    "# " + string(config.experiment)
    ""
    "H/fallback sweep diagnostic or baseline output."
    ""
    "- `summary.csv`: main tabular output."
    "- `raw_results.mat`: summary and configuration."
    ""
    "Note: " + string(note)
    ];
if ismember("fallback_count", string(summary.Properties.VariableNames))
    readmeLines = [
        readmeLines(1:4)
        "- raw Bootstrap-t and Bootstrap-t fallback-rule are separate rows."
        "- `trial_details.mat` or `curve_coverage_diagnostics.mat`: trial-level data for post-processing when available."
        "- `fallback_count`: number of times the replacement rule selected Percentile bootstrap."
        "- `bootstrap_t_inf_count`: number of raw Bootstrap-t infinite intervals."
        readmeLines(5:end)
        ];
end
dpimnumeric.writeText(fullfile(outDir, "README.md"), strjoin(readmeLines, newline));
end

function result = localResultPaths(outDir)
result = struct();
result.output_dir = outDir;
result.config = fullfile(outDir, "config.json");
result.raw_results = fullfile(outDir, "raw_results.mat");
result.summary = fullfile(outDir, "summary.csv");
result.readme = fullfile(outDir, "README.md");
end

function localSaveDetailBlocks(outDir, detailBlocks, config)
save(fullfile(outDir, "trial_details.mat"), "detailBlocks", "config", "-v7.3");
end

function details = localEmptyScalarTrialDetails(M, R)
methodNames = localScalarMethods();
details = struct();
details.method_names = methodNames;
details.trial = (1:M).';
details.ci_seed = nan(M, 1);
details.truth = nan(M, 1);
details.estimates = nan(M, R);
for iMethod = 1:numel(methodNames)
    field = localMethodField(methodNames(iMethod));
    details.(field + "_lower") = nan(M, 1);
    details.(field + "_upper") = nan(M, 1);
    details.(field + "_length") = nan(M, 1);
    details.(field + "_contains") = false(M, 1);
    details.(field + "_infinite") = false(M, 1);
    details.(field + "_fallback_trigger") = false(M, 1);
    details.(field + "_source_bootstrap_t_infinite") = false(M, 1);
end
details.bootstrap_t_min_std_ratio = nan(M, 1);
end

function details = localFillScalarTrialDetails(details, trial, ciSeed, truth, estimates, ci)
details.ci_seed(trial) = ciSeed;
details.truth(trial) = truth;
details.estimates(trial, 1:numel(estimates)) = estimates(:).';
for iMethod = 1:numel(details.method_names)
    methodName = details.method_names(iMethod);
    idx = find([ci.name] == methodName, 1);
    if isempty(idx)
        continue;
    end
    field = localMethodField(methodName);
    details.(field + "_lower")(trial) = ci(idx).lower;
    details.(field + "_upper")(trial) = ci(idx).upper;
    details.(field + "_length")(trial) = ci(idx).length;
    details.(field + "_contains")(trial) = ci(idx).contains;
    details.(field + "_infinite")(trial) = ci(idx).infinite;
    details.(field + "_fallback_trigger")(trial) = ci(idx).fallback_trigger;
    details.(field + "_source_bootstrap_t_infinite")(trial) = ci(idx).source_bootstrap_t_infinite;
    if methodName == "Bootstrap-t"
        details.bootstrap_t_min_std_ratio(trial) = ci(idx).min_std_ratio;
    end
end
end

function block = localPackScalarDetailBlock(context, details)
block = struct();
block.context = context;
block.details = details;
end

function field = localMethodField(methodName)
field = lower(regexprep(char(methodName), '[^A-Za-z0-9]+', '_'));
field = regexprep(field, '_$', '');
field = string(matlab.lang.makeValidName(field));
end

function [summary, simulationBlocks] = localFiniteBSummary(config)
nominal = 1 - config.alpha;
BList = config.B_list(:);
nB = numel(BList);
simulationBlocks = cell(nB, 1);
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
    rankHits = (ranks > kMinus(iB)) & (ranks <= kPlus(iB));
    simCoverage(iB) = mean(rankHits);
    simAbsError(iB) = abs(simCoverage(iB) - c0b(iB));
    simStdError(iB) = sqrt(c0b(iB) * (1 - c0b(iB)) / config.simulation_count);
    simCiHalfWidth95(iB) = 1.96 * sqrt(simCoverage(iB) * (1 - simCoverage(iB)) / config.simulation_count);
    simZScoreVsExact(iB) = (simCoverage(iB) - c0b(iB)) / max(simStdError(iB), realmin);
    simulationBlocks{iB} = struct("B", B, "rank_count", rankCount(iB), ...
        "k_minus", kMinus(iB), "k_plus", kPlus(iB), "ranks", ranks, ...
        "rank_hits", rankHits, "C0B", c0b(iB), "simulated_coverage", simCoverage(iB));
end

baseTable = table(BList, rankCount, kMinus, kPlus, coverageRankCount, c0b, ...
    repmat(nominal, nB, 1), gridError, absGridError, ...
    repmat(config.alpha / 2, nB, 1), lowerTailProbability, upperTailProbability, tailImbalance, ...
    simCoverage, simAbsError, simStdError, simCiHalfWidth95, simZScoreVsExact, ...
    repmat(config.simulation_count, nB, 1), ...
    'VariableNames', {'B', 'rank_count', 'k_minus', 'k_plus', 'coverage_rank_count', ...
    'C0B', 'nominal_coverage', 'coverage_error', 'abs_coverage_error', ...
    'nominal_one_sided_tail', 'lower_tail_probability', 'upper_tail_probability', 'tail_imbalance', ...
    'simulated_coverage', 'simulation_abs_error', 'simulation_std_error_vs_exact', ...
    'simulation_ci_half_width_95', 'simulation_z_score_vs_exact', 'simulation_count'});

methods = ["Percentile bootstrap finite-B rank"; "Bootstrap-t finite-B rank"];
endpointMapping = ["percentile endpoints use bootstrap order statistics directly"; ...
    "bootstrap-t endpoints use studentized pivot quantiles with lower/upper signs reversed"];
summary = table();
for iMethod = 1:numel(methods)
    block = baseTable;
    block.CI_method = repmat(methods(iMethod), height(block), 1);
    block.rank_interval_expression = repmat("k_minus < rank <= k_plus", height(block), 1);
    block.endpoint_mapping = repmat(endpointMapping(iMethod), height(block), 1);
    block.finite_B_expression = repmat("C0B=(k_plus-k_minus)/(B+1)", height(block), 1);
    summary = [summary; block]; %#ok<AGROW>
end
summary = movevars(summary, ["CI_method", "rank_interval_expression", ...
    "endpoint_mapping", "finite_B_expression"], "Before", 1);
end

function localPlotCoverageByH(outDir, summary, titleText, fileName)
fig = figure("Visible", "off", "Color", "w");
if any(summary.method == "Bootstrap-t") || any(summary.method == "Bootstrap-t fallback-rule")
    hold on;
    methodList = ["Bootstrap-t", "Bootstrap-t fallback-rule"];
    for iMethod = 1:numel(methodList)
        mask = summary.method == methodList(iMethod);
        if any(mask)
            [hVals, covVals] = localMeanCoverageByH(summary(mask, :));
            semilogx(hVals, covVals, "-o", "LineWidth", 1.3, "DisplayName", methodList(iMethod));
        end
    end
    hold off;
    xlabel("h");
    ylabel("mean coverage");
    title(titleText, "Interpreter", "none");
    grid on;
    legend("Location", "best");
else
    text(0.1, 0.5, "No Bootstrap-t rows.");
end
try
    exportgraphics(fig, fullfile(outDir, "figures", fileName), "Resolution", 180);
catch
    saveas(fig, fullfile(outDir, "figures", fileName));
end
close(fig);
end

function [hVals, covVals] = localMeanCoverageByH(tbl)
hVals = unique(tbl.h, "stable");
covVals = nan(size(hVals));
for iH = 1:numel(hVals)
    covVals(iH) = mean(tbl.coverage(tbl.h == hVals(iH)), "omitnan");
end
end

function accumulator = localEmptyAccumulator(methods)
n = numel(methods);
accumulator.methods = string(methods(:));
accumulator.coverageHits = zeros(n, 1);
accumulator.leftMiss = zeros(n, 1);
accumulator.rightMiss = zeros(n, 1);
accumulator.lengthSum = zeros(n, 1);
accumulator.finiteLengthCount = zeros(n, 1);
accumulator.fallbackCount = zeros(n, 1);
accumulator.bootstrapTInfCount = zeros(n, 1);
end

function accumulator = localAccumulate(accumulator, ci, truth)
for i = 1:numel(accumulator.methods)
    idx = find([ci.name] == accumulator.methods(i), 1);
    if isempty(idx)
        continue;
    end
    c = ci(idx);
    accumulator.coverageHits(i) = accumulator.coverageHits(i) + double(c.contains);
    accumulator.leftMiss(i) = accumulator.leftMiss(i) + double(c.upper < truth);
    accumulator.rightMiss(i) = accumulator.rightMiss(i) + double(c.lower > truth);
    accumulator.fallbackCount(i) = accumulator.fallbackCount(i) + double(c.fallback_trigger);
    accumulator.bootstrapTInfCount(i) = accumulator.bootstrapTInfCount(i) + double(c.source_bootstrap_t_infinite);
    if isfinite(c.length)
        accumulator.lengthSum(i) = accumulator.lengthSum(i) + c.length;
        accumulator.finiteLengthCount(i) = accumulator.finiteLengthCount(i) + 1;
    end
end
end

function tbl = localAccumulatorTable(accumulator, M, alpha)
if nargin < 3 || isempty(alpha)
    alpha = NaN;
end
n = numel(accumulator.methods);
meanLength = nan(n, 1);
mask = accumulator.finiteLengthCount > 0;
meanLength(mask) = accumulator.lengthSum(mask) ./ accumulator.finiteLengthCount(mask);
coverage = accumulator.coverageHits / M;
nominal = 1 - alpha;
coverageMcSe = sqrt(max(coverage .* (1 - coverage), 0) / M);
tbl = table(accumulator.methods, coverage, repmat(nominal, n, 1), ...
    coverage - nominal, abs(coverage - nominal), coverageMcSe, 1.96 * coverageMcSe, meanLength, ...
    accumulator.leftMiss / M, accumulator.rightMiss / M, ...
    accumulator.fallbackCount, accumulator.fallbackCount / M, ...
    accumulator.bootstrapTInfCount, accumulator.bootstrapTInfCount / M, ...
    accumulator.finiteLengthCount, accumulator.finiteLengthCount / M, ...
    'VariableNames', {'method', 'coverage', 'nominal_coverage', ...
    'coverage_error', 'abs_coverage_error', 'coverage_mc_se', 'coverage_ci_half_width_95', ...
    'mean_interval_length', 'left_miss', 'right_miss', 'fallback_count', 'fallback_rate', ...
    'bootstrap_t_inf_count', 'bootstrap_t_inf_rate', 'finite_length_count', 'finite_length_rate'});
end

function ci = localCiMethods(Y, alpha, B, seed, truth)
previousState = rng;
rng(double(seed), "twister");
cleanupObj = onCleanup(@() rng(previousState)); %#ok<NASGU>

Y = double(Y(:));
R = numel(Y);
Ybar = mean(Y);
SR = std(Y);
sqrtR = sqrt(R);

tCrit = localStudentTCrit(1 - alpha / 2, R - 1);
studentLower = Ybar - tCrit * SR / sqrtR;
studentUpper = Ybar + tCrit * SR / sqrtR;

bootIndices = randi(R, R, B);
Ystar = reshape(Y(bootIndices(:)), R, B);
bootMeans = mean(Ystar, 1).';
bootStd = std(Ystar, 0, 1).';
[pctLower, pctUpper] = localQuantileRows(bootMeans, alpha);

Tstar = sqrtR * (bootMeans - Ybar) ./ bootStd;
zeroStd = bootStd == 0;
if any(zeroStd)
    delta = bootMeans - Ybar;
    signDelta = sign(delta);
    signDelta(signDelta == 0) = 1;
    Tstar(zeroStd) = signDelta(zeroStd) .* Inf;
end
[tLower, tUpper] = localQuantileRows(Tstar, alpha);
btLower = Ybar - SR * tUpper / sqrtR;
btUpper = Ybar - SR * tLower / sqrtR;
minStdRatio = min(bootStd ./ max(SR, realmin));

ci = repmat(localEmptyCi(), 3, 1);
ci(1) = localPackCi("Student-t", studentLower, studentUpper, truth, NaN, false, false);
ci(2) = localPackCi("Percentile bootstrap", pctLower, pctUpper, truth, NaN, false, false);
ci(3) = localPackCi("Bootstrap-t", btLower, btUpper, truth, minStdRatio, false, false);
ci(3).source_bootstrap_t_infinite = ci(3).infinite;
end

function ci = localAppendFallbackRule(ci, lambda)
p = ci(2);
bt = ci(3);
trigger = bt.infinite || (bt.length > lambda * max(p.length, realmin));
if trigger
    hybrid = p;
else
    hybrid = bt;
end
hybrid.name = "Bootstrap-t fallback-rule";
hybrid.fallback_trigger = trigger;
hybrid.source_bootstrap_t_infinite = bt.infinite;
ci(end + 1) = hybrid;
end

function ci = localEmptyCi()
ci = struct("name", "", "lower", NaN, "upper", NaN, "length", NaN, ...
    "contains", false, "infinite", false, "min_std_ratio", NaN, ...
    "fallback_trigger", false, "source_bootstrap_t_infinite", false);
end

function ci = localPackCi(name, lower, upper, truth, minStdRatio, fallbackTrigger, sourceBtInf)
lower = lower(1);
upper = upper(1);
if lower > upper
    tmp = lower;
    lower = upper;
    upper = tmp;
end
ci = localEmptyCi();
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
ci.fallback_trigger = fallbackTrigger;
ci.source_bootstrap_t_infinite = sourceBtInf;
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

function y = localCenteredExponentialSmoothedPdf(x, h)
y = exp(-1 - x + 0.5 * h ^ 2) .* localNormalCdf((x + 1 - h ^ 2) ./ h);
end

function F = localNormalCdf(x)
F = 0.5 * (1 + erf(x ./ sqrt(2)));
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

function response = localLockedLognormalTruth(d, sampleCount, seed)
previousState = rng;
rng(seed, "twister");
cleanupObj = onCleanup(@() rng(previousState)); %#ok<NASGU>
theta = exp(randn(sampleCount, d));
response = localLockedResponse(theta);
end

function estimates = localNonlinearKernelEstimates(R, nInner, d, yTarget, h)
estimates = zeros(R, 1);
for r = 1:R
    theta = exp(randn(nInner, d));
    response = localLockedResponse(theta);
    estimates(r) = mean(localNormalPdf(yTarget, response, h));
end
end

function response = localLockedResponse(theta)
z = mean(theta, 2);
e = mean(theta .^ 2, 2) - 1;
response = z + 0.4 * z .^ 3 + 0.8 * z .* e;
end

function y = localNormalPdf(x, mu, sigma)
y = exp(-0.5 * ((x - mu) ./ sigma) .^ 2) ./ (sqrt(2 * pi) .* sigma);
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

function [wGrid, referenceResponses, curveResponses] = localBuildIndependentPlateResponses(config)
shared = localPlateShared(config);
referenceResponses = abs(localPlateResponseSamples(shared, config.reference_sample_count, ...
    config.d, config.sigma_E, config.seed + 11));
curveResponses = zeros(config.curve_pool_count, config.inner_sample_count);
for iPool = 1:config.curve_pool_count
    curveResponses(iPool, :) = abs(localPlateResponseSamples(shared, config.inner_sample_count, ...
        config.d, config.sigma_E, config.seed + 1000 + iPool));
end
q = localQuantiles(referenceResponses, [0.001, 0.999]);
width = max(q(2) - q(1), eps);
wGrid = linspace(max(0, q(1) - 0.05 * width), q(2) + 0.05 * width, config.w_grid_count);
end

function shared = localPlateShared(config)
[nodes, tris, edges, tri2edge, edgeNormals, boundary] = localGenerateSkewTriMesh(config.mesh_nx, config.mesh_ny);
shared = struct();
shared.nodes = nodes;
shared.tris = tris;
shared.edges = edges;
shared.tri2edge = tri2edge;
shared.edgeNormals = edgeNormals;
shared.boundary = boundary;
shared.nnode = size(nodes, 1);
shared.nelem = size(tris, 1);
shared.ndof = shared.nnode + size(edges, 1);
shared.E0 = config.E0;
shared.nu = config.nu;
shared.t = config.t;
shared.q0 = config.q0;
shared.xc = 0.5;
shared.yc = 0.5;
end

function responses = localPlateResponseSamples(shared, sampleCount, d, sigmaE, seed)
previousState = rng;
rng(seed, "twister");
cleanupObj = onCleanup(@() rng(previousState)); %#ok<NASGU>
Z = randn(sampleCount, d);
randomE = localPlateRandomField(shared.nodes, Z, shared.E0, sigmaE);
responses = zeros(sampleCount, 1);
for iSample = 1:sampleCount
    responses(iSample) = localSolvePlateSample(shared, randomE, iSample);
end
end

function randomE = localPlateRandomField(nodes, Z, E0, sigmaE)
[~, d] = size(Z);
nnode = size(nodes, 1);
modePairs = localModePairs(d);
Phi = zeros(nnode, d);
for j = 1:d
    px = modePairs(j, 1);
    py = modePairs(j, 2);
    Phi(:, j) = sin(pi * px * nodes(:, 1)) .* sin(pi * py * nodes(:, 2)) / sqrt(px ^ 2 + py ^ 2);
end
field = Phi * Z.';
field = field ./ max(std(field, 0, 1), realmin);
randomE = E0 * exp(sigmaE * field);
end

function pairs = localModePairs(d)
pairs = zeros(d, 2);
idx = 0;
level = 1;
while idx < d
    for px = 1:level
        py = level + 1 - px;
        idx = idx + 1;
        pairs(idx, :) = [px, py];
        if idx >= d
            return;
        end
    end
    level = level + 1;
end
end

function response = localSolvePlateSample(shared, randomE, sampleIndex)
K = sparse(shared.ndof, shared.ndof);
F = zeros(shared.ndof, 1);
for elem = 1:shared.nelem
    vids = shared.tris(elem, :);
    eids = shared.tri2edge(elem, :);
    xe = shared.nodes(vids, 1);
    ye = shared.nodes(vids, 2);
    [Ke, Fe] = localMorleyPlateElement(xe, ye, eids, shared.edgeNormals, shared.nu, shared.t, ...
        shared.q0, vids, randomE, sampleIndex, shared.E0);
    edofs = [vids(:); shared.nnode + eids(:)];
    K(edofs, edofs) = K(edofs, edofs) + Ke;
    F(edofs) = F(edofs) + Fe;
end
fixedDofs = unique(shared.boundary.nodes(:));
freeDofs = setdiff((1:shared.ndof).', fixedDofs);
U = zeros(shared.ndof, 1);
U(freeDofs) = K(freeDofs, freeDofs) \ F(freeDofs);
response = localEvalMorleyPoint(shared.xc, shared.yc, shared.nodes, shared.tris, shared.tri2edge, ...
    shared.edgeNormals, U, shared.nnode);
end

function [nodes, tris, edges, tri2edge, edgeNormals, boundary] = localGenerateSkewTriMesh(nx, ny)
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
tris = zeros(2 * nx * ny, 3);
elem = 0;
for j = 1:ny
    for i = 1:nx
        n1 = (j - 1) * (nx + 1) + i;
        n2 = n1 + 1;
        n4 = j * (nx + 1) + i;
        n3 = n4 + 1;
        elem = elem + 1;
        tris(elem, :) = [n1, n2, n3];
        elem = elem + 1;
        tris(elem, :) = [n1, n3, n4];
    end
end
allLocalEdges = zeros(3 * size(tris, 1), 2);
for elem = 1:size(tris, 1)
    n1 = tris(elem, 1);
    n2 = tris(elem, 2);
    n3 = tris(elem, 3);
    row0 = 3 * (elem - 1);
    allLocalEdges(row0 + 1, :) = sort([n2, n3]);
    allLocalEdges(row0 + 2, :) = sort([n3, n1]);
    allLocalEdges(row0 + 3, :) = sort([n1, n2]);
end
[edges, ~, ic] = unique(allLocalEdges, "rows");
tri2edge = reshape(ic, 3, size(tris, 1)).';
edgeNormals = zeros(size(edges, 1), 2);
for idx = 1:size(edges, 1)
    i = edges(idx, 1);
    j = edges(idx, 2);
    dx = nodes(j, 1) - nodes(i, 1);
    dy = nodes(j, 2) - nodes(i, 2);
    L = sqrt(dx ^ 2 + dy ^ 2);
    edgeNormals(idx, :) = [dy, -dx] / L;
end
boundary.bottom = 1:(nx + 1);
boundary.top = ny * (nx + 1) + (1:(nx + 1));
boundary.left = 1:(nx + 1):(ny * (nx + 1) + 1);
boundary.right = (nx + 1):(nx + 1):((ny + 1) * (nx + 1));
boundary.nodes = unique([boundary.bottom(:); boundary.top(:); boundary.left(:); boundary.right(:)]);
end

function [Ke, Fe] = localMorleyPlateElement(xe, ye, eids, edgeNormals, nu, t, q0, vids, randomE, sampleIndex, E0)
xe = xe(:);
ye = ye(:);
x1 = xe(1); y1 = ye(1);
x2 = xe(2); y2 = ye(2);
x3 = xe(3); y3 = ye(3);
A2 = det([1, x1, y1; 1, x2, y2; 1, x3, y3]);
A = abs(A2) / 2;
E = E0;
if size(randomE, 1) >= max(vids) && size(randomE, 2) >= sampleIndex
    E = mean(randomE(vids, sampleIndex));
end
D0 = E * t ^ 3 / (12 * (1 - nu ^ 2));
Db = D0 * [1, nu, 0; nu, 1, 0; 0, 0, (1 - nu) / 2];
M = zeros(6, 6);
M(1, :) = localPolyRow(x1, y1);
M(2, :) = localPolyRow(x2, y2);
M(3, :) = localPolyRow(x3, y3);
m23 = 0.5 * [x2 + x3, y2 + y3];
m31 = 0.5 * [x3 + x1, y3 + y1];
m12 = 0.5 * [x1 + x2, y1 + y2];
M(4, :) = localDNormalRow(m23(1), m23(2), edgeNormals(eids(1), 1), edgeNormals(eids(1), 2));
M(5, :) = localDNormalRow(m31(1), m31(2), edgeNormals(eids(2), 1), edgeNormals(eids(2), 2));
M(6, :) = localDNormalRow(m12(1), m12(2), edgeNormals(eids(3), 1), edgeNormals(eids(3), 2));
C = M \ eye(6);
B = zeros(3, 6);
for idx = 1:6
    ci = C(:, idx);
    B(:, idx) = [2 * ci(4); 2 * ci(6); 2 * ci(5)];
end
Ke = A * (B.' * Db * B);
Fe = zeros(6, 1);
Lq = [1 / 6, 1 / 6, 2 / 3; 1 / 6, 2 / 3, 1 / 6; 2 / 3, 1 / 6, 1 / 6];
for q = 1:3
    xq = Lq(q, 1) * x1 + Lq(q, 2) * x2 + Lq(q, 3) * x3;
    yq = Lq(q, 1) * y1 + Lq(q, 2) * y2 + Lq(q, 3) * y3;
    Nq = localPolyRow(xq, yq) * C;
    Fe = Fe + q0 * Nq.' * (A / 3);
end
end

function wVal = localEvalMorleyPoint(xp, yp, nodes, tris, tri2edge, edgeNormals, U, nnode)
for elem = 1:size(tris, 1)
    vids = tris(elem, :);
    pts = nodes(vids, :);
    x1 = pts(1, 1); y1 = pts(1, 2);
    x2 = pts(2, 1); y2 = pts(2, 2);
    x3 = pts(3, 1); y3 = pts(3, 2);
    A2 = det([x2 - x1, x3 - x1; y2 - y1, y3 - y1]);
    l1 = det([x2 - xp, x3 - xp; y2 - yp, y3 - yp]) / A2;
    l2 = det([xp - x1, x3 - x1; yp - y1, y3 - y1]) / A2;
    l3 = 1 - l1 - l2;
    if (l1 >= -1e-10) && (l2 >= -1e-10) && (l3 >= -1e-10)
        eids = tri2edge(elem, :);
        M = zeros(6, 6);
        M(1, :) = localPolyRow(x1, y1);
        M(2, :) = localPolyRow(x2, y2);
        M(3, :) = localPolyRow(x3, y3);
        m23 = 0.5 * [x2 + x3, y2 + y3];
        m31 = 0.5 * [x3 + x1, y3 + y1];
        m12 = 0.5 * [x1 + x2, y1 + y2];
        M(4, :) = localDNormalRow(m23(1), m23(2), edgeNormals(eids(1), 1), edgeNormals(eids(1), 2));
        M(5, :) = localDNormalRow(m31(1), m31(2), edgeNormals(eids(2), 1), edgeNormals(eids(2), 2));
        M(6, :) = localDNormalRow(m12(1), m12(2), edgeNormals(eids(3), 1), edgeNormals(eids(3), 2));
        C = M \ eye(6);
        edofs = [vids(:); nnode + eids(:)];
        Np = localPolyRow(xp, yp) * C;
        wVal = Np * U(edofs);
        return;
    end
end
dist2 = (nodes(:, 1) - xp) .^ 2 + (nodes(:, 2) - yp) .^ 2;
[~, id] = min(dist2);
wVal = U(id);
end

function row = localPolyRow(x, y)
row = [1, x, y, x ^ 2, x * y, y ^ 2];
end

function row = localDNormalRow(x, y, nx, ny)
row = [0, nx, ny, 2 * x * nx, y * nx + x * ny, 2 * y * ny];
end

function curve = localKernelCurve(responses, wGrid, h)
curve = zeros(1, numel(wGrid));
for iw = 1:numel(wGrid)
    curve(iw) = mean(localNormalPdf(wGrid(iw), responses(:), h));
end
end

function curvePool = localKernelCurvePool(responseBatches, wGrid, h)
curvePool = zeros(size(responseBatches, 1), numel(wGrid));
for i = 1:size(responseBatches, 1)
    curvePool(i, :) = localKernelCurve(responseBatches(i, :), wGrid, h);
end
end

function [summary, pointwise, diagnostics] = localCurvePoolCoverageWithFallback(curvePool, referenceCurve, wGrid, alpha, B, M, R, seed, lambda)
nGrid = numel(wGrid);
coverageHits = zeros(4, nGrid);
lengthSums = zeros(4, nGrid);
fallbackCounts = zeros(4, nGrid);
btInfCounts = zeros(4, nGrid);
diagnostics = localEmptyCurveDiagnostics(M, R, nGrid, localScalarMethods());
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
bandDiagnostics.reference_curve = referenceCurve;
bandDiagnostics.w_grid = wGrid;
bandDiagnostics.band_point_hits = bandPointHits;
bandDiagnostics.band_length_sum = bandLengthSum;
bandDiagnostics.band_inf_count = bandInfCount;
diagnostics = struct();
diagnostics.pointwise = pointDiagnostics;
diagnostics.simultaneous_band = bandDiagnostics;
end

function diagnostics = localEmptyCurveDiagnostics(M, R, nGrid, methodNames)
diagnostics = struct();
diagnostics.method_names = methodNames;
diagnostics.sample_indices = zeros(M, R);
diagnostics.lower = nan(M, numel(methodNames), nGrid);
diagnostics.upper = nan(M, numel(methodNames), nGrid);
diagnostics.length = nan(M, numel(methodNames), nGrid);
diagnostics.contains = false(M, numel(methodNames), nGrid);
diagnostics.fallback_trigger = false(M, numel(methodNames), nGrid);
diagnostics.source_bootstrap_t_infinite = false(M, numel(methodNames), nGrid);
end
