function outputs = run_weighted_paper_evidence_suite(runMode, resultsRoot)
%RUN_WEIGHTED_PAPER_EVIDENCE_SUITE Probability-weighted paper evidence loop.
%
% Usage:
%   run_weighted_paper_evidence_suite("diagnostic")
%   run_weighted_paper_evidence_suite("small")
%   run_weighted_paper_evidence_suite("full")

if nargin < 1 || strlength(string(runMode)) == 0
    runMode = "diagnostic";
end

projectRoot = fileparts(mfilename("fullpath"));
if nargin < 2 || strlength(string(resultsRoot)) == 0
    cfg = weighted_paper_config(runMode, projectRoot);
else
    cfg = weighted_paper_config(runMode, projectRoot, resultsRoot);
end

dpimnumeric.ensureDir(cfg.results_root);
localAddPaths(cfg);
preflight = localPreflight(cfg);
dpimnumeric.writeJson(fullfile(cfg.results_root, "weighted_paper_config.json"), cfg);
writetable(preflight, fullfile(cfg.results_root, "preflight.csv"));

if any(preflight.status == "FAIL")
    error("Preflight failed. Inspect %s.", fullfile(cfg.results_root, "preflight.csv"));
end

manifest = struct();
manifest.schema_version = cfg.schema_version;
manifest.run_mode = cfg.run_mode;
manifest.project_root = cfg.project_root;
manifest.results_root = cfg.results_root;
manifest.paper_tex = cfg.paper_tex;
manifest.reference_dir = cfg.reference_dir;
manifest.thesis_pdf = cfg.thesis_pdf;
manifest.note = "Probability-weighted paper evidence suite. Equal-weight controls are excluded from main evidence.";
dpimnumeric.writeJson(fullfile(cfg.results_root, "manifest.json"), manifest);

fprintf("Weighted paper evidence suite: %s\n", cfg.run_mode);
fprintf("Output root: %s\n", cfg.results_root);

outputs = struct();
outputs.results_root = cfg.results_root;
outputs.E1 = localRunFiniteB(cfg);

pointExperiments = localPointExperiments();
for i = 1:numel(pointExperiments)
    expSpec = pointExperiments(i);
    fprintf("\n=== %s ===\n", expSpec.name);
    outputs.(matlab.lang.makeValidName(expSpec.name)) = localRunPointwiseExperiment(cfg, expSpec);
end

fprintf("\n=== E8_plate_RFEM_simultaneous_band_weighted_diagnostic ===\n");
outputs.E8_plate_RFEM_simultaneous_band_weighted_diagnostic = localRunE8BandDiagnostic(cfg);

outputs.original_plate_reference = localIndexOriginalPlate(cfg);

if exist("analyze_weighted_paper_evidence_suite", "file") == 2
    outputs.analysis = analyze_weighted_paper_evidence_suite(cfg.results_root);
end
if exist("build_weighted_paper_report", "file") == 2
    outputs.report = build_weighted_paper_report(cfg.results_root);
end

save(fullfile(cfg.results_root, "suite_outputs.mat"), "outputs", "cfg", "preflight", "-v7.3");
fprintf("\nWeighted paper evidence suite completed: %s\n", cfg.results_root);
end

function localAddPaths(cfg)
addpath(cfg.project_root, "-begin");
addpath(fullfile(cfg.weighted_root, "external_weight_providers"), "-begin");
addpath(fullfile(cfg.weighted_root, "common"), "-begin");
end

function T = localPreflight(cfg)
rows = cell(0, 1);
rows{end+1,1} = localCheckRow("paper_tex", isfile(cfg.paper_tex), cfg.paper_tex);
rows{end+1,1} = localCheckRow("reference_dir", isfolder(cfg.reference_dir), cfg.reference_dir);
rows{end+1,1} = localCheckRow("thesis_pdf", isfile(cfg.thesis_pdf), cfg.thesis_pdf);
rows{end+1,1} = localCheckRow("weighted_root", isfolder(cfg.weighted_root), cfg.weighted_root);
for i = 1:numel(cfg.required_functions)
    name = cfg.required_functions(i);
    rows{end+1,1} = localCheckRow("function:" + name, exist(name, "file") == 2, which(name)); %#ok<AGROW>
end

[gpuOk, gpuMsg] = localGpuStatus();
rows{end+1,1} = localCheckRow("gpu:" + cfg.assignment_backend, gpuOk || cfg.assignment_backend ~= "gpu", gpuMsg);
T = vertcat(rows{:});
end

function row = localCheckRow(name, ok, detail)
if ok
    status = "PASS";
else
    status = "FAIL";
end
row = table(string(name), string(status), string(detail), ...
    'VariableNames', {'check','status','detail'});
end

function [ok, msg] = localGpuStatus()
try
    dev = gpuDevice;
    ok = true;
    msg = sprintf("%s, available %.2f GiB", dev.Name, dev.AvailableMemory / 2^30);
catch ME
    ok = false;
    msg = ME.message;
end
end

function output = localRunFiniteB(cfg)
outDir = fullfile(cfg.results_root, "E1_finite_B_rank_audit");
dpimnumeric.ensureDir(outDir);
BList = unique([cfg.B, 99, 199, 399, 999, 1999]);
Msim = max(1000, cfg.M * 20);
rng(cfg.seed + 101, "twister");
rows = table();
for i = 1:numel(BList)
    B = BList(i);
    kMinus = floor((cfg.alpha / 2) * (B + 1));
    kPlus = ceil((1 - cfg.alpha / 2) * (B + 1));
    C0B = (kPlus - kMinus) / (B + 1);
    sim = localSimulateFiniteBOrderStatisticCoverage(B, kMinus, kPlus, Msim);
    coverage = mean(sim);
    se = sqrt(max(coverage * (1 - coverage), realmin) / Msim);
    rows = [rows; table("E1_finite_B_rank_audit", B, kMinus, kPlus, C0B, ...
        coverage, coverage - cfg.nominal, abs(coverage - C0B), 1.96 * se, Msim, ...
        'VariableNames', {'experiment','B','k_minus','k_plus','C0B', ...
        'simulated_coverage','coverage_error_vs_nominal','simulation_error_vs_C0B', ...
        'simulation_ci_half_width_95','simulation_count'})]; %#ok<AGROW>
end
writetable(rows, fullfile(outDir, "summary.csv"));
save(fullfile(outDir, "raw_results.mat"), "rows", "cfg", "-v7.3");
output = struct("output_dir", outDir, "summary_csv", fullfile(outDir, "summary.csv"));
end

function hit = localSimulateFiniteBOrderStatisticCoverage(B, kMinus, kPlus, Msim)
% Independently simulate a true Uniform pivot and B bootstrap pivots.
% The event U_(kMinus) <= U0 <= U_(kPlus) is equivalent to the rank of U0
% among the B+1 IID uniforms lying in {kMinus+1, ..., kPlus}.
hit = false(Msim, 1);
chunkSize = min(Msim, 1000);
for first = 1:chunkSize:Msim
    last = min(Msim, first + chunkSize - 1);
    count = last - first + 1;
    u0 = rand(count, 1);
    bootstrapPivots = rand(count, B);
    trueRank = 1 + sum(bootstrapPivots < u0, 2);
    hit(first:last) = trueRank >= (kMinus + 1) & trueRank <= kPlus;
end
end

function specs = localPointExperiments()
specs = [ ...
    struct("name", "E2_scalar_normal_R_order_weighted", ...
    "d", 1, "y_points", [-1, 0, 1], "anchor_h", 0.233572146909012, ...
    "paper_role", "main_candidate"), ...
    struct("name", "E3_linear_beam_weighted", ...
    "d", 5, "y_points", [-0.5, 0, 0.5], "anchor_h", 0.00183298071083244, ...
    "paper_role", "main_candidate"), ...
    struct("name", "E4_standard_normal_nonlinear_tail_weighted", ...
    "d", 5, "y_points", [1.6, 2.0, 2.4], "anchor_h", 4e-4, ...
    "paper_role", "main_candidate"), ...
    struct("name", "E5_standard_normal_bootstrap_t_instability_weighted", ...
    "d", 5, "y_points", [2.0, 2.5, 3.0], "anchor_h", 4e-4, ...
    "paper_role", "failed_boundary_candidate") ...
    ];
end

function output = localRunPointwiseExperiment(cfg, expSpec)
outDir = fullfile(cfg.results_root, expSpec.name);
dpimnumeric.ensureDir(outDir);
dpimnumeric.ensureDir(fullfile(outDir, "figures"));
problem = localProblem(expSpec);
expCfg = localExperimentCfg(cfg, outDir);
tuningCfg = expCfg;
tuningCfg.M = cfg.tuning_M;
tuningCfg.B = cfg.tuning_B;
dpimnumeric.writeJson(fullfile(outDir, "config.json"), localJsonSafe(expCfg, expSpec));

tuning = table();
confirmation = table();
mechanism = table();
weightRows = table();
candidateRows = table();
protocol = table();
momentValidation = table();
rawBlocks = cell(0, 1);
confirmBlocks = cell(0, 1);
protocolBlocks = cell(0, 1);

for im = 1:numel(cfg.methods)
    method = cfg.methods(im);
    fprintf("Building weighted curve pool: %s, method=%s, curves=%d, n=%d\n", ...
        expSpec.name, method, cfg.curve_pool_size, cfg.n);
    curves = dpim_build_weighted_curve_pool(problem, method, cfg.n, expCfg, cfg.curve_pool_size);
    methodWeights = localWeightTable(curves, expSpec.name, method);
    weightRows = [weightRows; methodWeights]; %#ok<AGROW>
    weightStats = localAggregateWeightStats(methodWeights);
    yPoints = double(expSpec.y_points(:)).';
    poolEst = localEvaluatePool(curves, yPoints, cfg.h_list);
    truth = localTruthGrid(problem, yPoints, cfg.h_list, expCfg);
    tuningIdx = 1:cfg.tuning_pool_size;
    confirmationIdx = (cfg.tuning_pool_size + 1):cfg.curve_pool_size;

    for iy = 1:numel(yPoints)
        for ih = 1:numel(cfg.h_list)
            for ir = 1:numel(cfg.R_list)
                R = cfg.R_list(ir);
                seed = localBlockSeed(cfg.seed, 1, im, iy, ih, ir);
                [rows, mech, raw] = localCoverageBlock( ...
                    squeeze(poolEst(tuningIdx, iy, ih)), truth(iy, ih), tuningCfg, ...
                    expSpec.name, method, yPoints(iy), cfg.h_list(ih), ih, R, ...
                    "tuning", seed, weightStats);
                tuning = [tuning; rows]; %#ok<AGROW>
                mechanism = [mechanism; mech]; %#ok<AGROW>
                rawBlocks{end+1,1} = raw; %#ok<AGROW>
            end
        end
    end

    candidates = localSelectCandidates(tuning(tuning.point_method == method, :), cfg, expSpec, method);
    candidateRows = [candidateRows; candidates]; %#ok<AGROW>
    for ic = 1:height(candidates)
        ih = candidates.h_index(ic);
        R = candidates.R(ic);
        ir = find(cfg.R_list == R, 1, "first");
        iy = find(abs(yPoints - candidates.y0(ic)) <= 10 * eps(max(1, abs(candidates.y0(ic)))), 1, "first");
        if isempty(iy)
            error("Candidate y0=%g is not in experiment y_points.", candidates.y0(ic));
        end
            seed = localBlockSeed(cfg.seed, 2 + ic, im, iy, ih, ir);
            [rows, mech, raw] = localCoverageBlock( ...
                squeeze(poolEst(confirmationIdx, iy, ih)), truth(iy, ih), expCfg, ...
                expSpec.name, method, yPoints(iy), cfg.h_list(ih), ih, R, ...
                "confirmation", seed, weightStats);
            rows.selection_reason = repmat(candidates.selection_reason(ic), height(rows), 1);
            confirmation = [confirmation; rows]; %#ok<AGROW>
            mechanism = [mechanism; mech]; %#ok<AGROW>
            confirmBlocks{end+1,1} = raw; %#ok<AGROW>
    end

    anchorIndex = localNearestIndex(cfg.h_list, expSpec.anchor_h);
    protocolBList = localProtocolBList(cfg, expSpec.name);
    for iy = 1:numel(yPoints)
        for ir = 1:numel(cfg.R_list)
            R = cfg.R_list(ir);
            for ib = 1:numel(protocolBList)
                protocolCfg = expCfg;
                protocolCfg.B = protocolBList(ib);
                seed = localBlockSeed(cfg.seed, 40 + ib, im, iy, anchorIndex, ir);
                [rows, mech, raw] = localCoverageBlock( ...
                    squeeze(poolEst(confirmationIdx, iy, anchorIndex)), ...
                    truth(iy, anchorIndex), protocolCfg, expSpec.name, method, ...
                    yPoints(iy), cfg.h_list(anchorIndex), anchorIndex, R, ...
                    "protocol_fixed_h", seed, weightStats);
                rows.protocol_role = repmat("fixed_h_factorial_validation", height(rows), 1);
                protocol = [protocol; rows]; %#ok<AGROW>
                mechanism = [mechanism; mech]; %#ok<AGROW>
                protocolBlocks{end+1,1} = raw; %#ok<AGROW>
            end
        end
    end
end

if string(expSpec.name) == "E2_scalar_normal_R_order_weighted"
    anchorIndex = localNearestIndex(cfg.h_list, expSpec.anchor_h);
    momentValidation = localGaussianMomentValidation( ...
        yPoints, cfg.h_list(anchorIndex), cfg, expSpec.name);
end

writetable(tuning, fullfile(outDir, "tuning_summary.csv"));
writetable(candidateRows, fullfile(outDir, "candidate_h.csv"));
writetable(confirmation, fullfile(outDir, "confirmation_summary.csv"));
writetable(protocol, fullfile(outDir, "protocol_summary.csv"));
writetable(mechanism, fullfile(outDir, "mechanism_summary.csv"));
writetable(weightRows, fullfile(outDir, "weight_diagnostics.csv"));
if ~isempty(momentValidation)
    writetable(momentValidation, fullfile(outDir, "kernel_moment_validation.csv"));
end
save(fullfile(outDir, "raw_results.mat"), "cfg", "expSpec", "tuning", ...
    "confirmation", "protocol", "momentValidation", "candidateRows", ...
    "mechanism", "weightRows", "rawBlocks", "confirmBlocks", ...
    "protocolBlocks", "-v7.3");
localPlotPointwise(outDir, tuning, confirmation, expSpec.name);
output = struct("output_dir", outDir, ...
    "tuning_csv", fullfile(outDir, "tuning_summary.csv"), ...
    "confirmation_csv", fullfile(outDir, "confirmation_summary.csv"), ...
    "protocol_csv", fullfile(outDir, "protocol_summary.csv"));
end

function expCfg = localExperimentCfg(cfg, outDir)
expCfg = cfg;
expCfg.results_root = outDir;
expCfg.project_root = cfg.project_root;
expCfg.weighting_cfg = cfg.weighting_cfg;
expCfg.weighting_cfg.voronoi_output_dir = fullfile(outDir, "ci_probability_weights");
end

function S = localJsonSafe(expCfg, expSpec)
S = expCfg;
S.experiment = expSpec.name;
S.y_points = expSpec.y_points;
S.anchor_h = expSpec.anchor_h;
S.paper_role = expSpec.paper_role;
end

function problem = localProblem(expSpec)
problem = struct();
problem.name = char(expSpec.name);
problem.short_name = char(extractBefore(string(expSpec.name) + "_", "_"));
problem.d = expSpec.d;
problem.target_distribution = "standard_normal";
problem.center_transform = "normal_icdf";
problem.provider = "voronoi_ci_probability_weights_provider";
switch string(expSpec.name)
    case "E2_scalar_normal_R_order_weighted"
        problem.response_fun = @(theta) theta(:,1);
    case "E3_linear_beam_weighted"
        coeff = localBeamCoeffs(problem.d);
        problem.response_fun = @(theta) theta * coeff(:);
    otherwise
        problem.response_fun = @(theta) localNormalNonlinearResponse(theta);
end
end

function coeff = localBeamCoeffs(d)
x0 = 0.5;
j = (1:d)';
coeff = sin(j * pi * x0) ./ ((j * pi).^2);
coeff = coeff / max(norm(coeff), realmin);
end

function y = localNormalNonlinearResponse(theta)
y = 2.0 + 0.12 * sum(theta, 2) ...
    + 0.25 * tanh(theta(:,1) .* theta(:,2)) ...
    + 0.08 * (theta(:,3).^2 - 1) ...
    + 0.06 * sin(theta(:,4) .* theta(:,5));
end

function poolEst = localEvaluatePool(curves, yPoints, hList)
poolEst = zeros(numel(curves), numel(yPoints), numel(hList));
for iy = 1:numel(yPoints)
    for ih = 1:numel(hList)
        poolEst(:, iy, ih) = dpim_curve_point_estimates(curves, yPoints(iy), hList(ih));
    end
end
end

function truth = localTruthGrid(problem, yPoints, hList, cfg)
truth = zeros(numel(yPoints), numel(hList));
for iy = 1:numel(yPoints)
    for ih = 1:numel(hList)
        truth(iy, ih) = localTruth(problem, yPoints(iy), hList(ih), cfg, 1000 * iy + ih);
    end
end
end

function truth = localTruth(problem, y0, h, cfg, seedOffset)
name = string(problem.name);
if name == "E2_scalar_normal_R_order_weighted"
    truth = dpim_gaussian_kernel(y0, sqrt(1 + h^2));
elseif name == "E3_linear_beam_weighted"
    coeff = localBeamCoeffs(problem.d);
    truth = dpim_gaussian_kernel(y0, sqrt(sum(coeff.^2) + h^2));
else
    truth = dpim_truth_smoothed_density(problem, y0, h, cfg, seedOffset);
end
end

function [rows, mechanismRow, raw] = localCoverageBlock(poolVals, truth, cfg, ...
    experiment, pointMethod, y0, h, hIndex, R, phase, blockSeed, weightStats)
poolVals = double(poolVals(:));
rng(blockSeed, "twister");
methodNames = ["Student-t"; "Percentile bootstrap"; "Bootstrap-t"; "Bootstrap-t fallback-rule"];
nMethods = numel(methodNames);
cover = false(cfg.M, nMethods);
leftMiss = false(cfg.M, nMethods);
rightMiss = false(cfg.M, nMethods);
isInf = false(cfg.M, nMethods);
lengths = nan(cfg.M, nMethods);
fallback = false(cfg.M, 1);
btInf = false(cfg.M, 1);
overlong = false(cfg.M, 1);
zeroSdCount = zeros(cfg.M, 1);
minPositiveSd = nan(cfg.M, 1);
lengthRatio = nan(cfg.M, 1);
sampleSd = nan(cfg.M, 1);
poolSize = numel(poolVals);
poolStats = localVectorStats(poolVals);
poolSe = poolStats.std / sqrt(R);
signedBiasSeRatio = (poolStats.mean - truth) / max(poolSe, realmin);

for m = 1:cfg.M
    idx = randi(poolSize, R, 1);
    x = poolVals(idx);
    [ci, diag] = localCiMethods(x, truth, cfg.alpha, cfg.B, cfg.lambda);
    for j = 1:nMethods
        cover(m,j) = ci(j).contains;
        leftMiss(m,j) = ci(j).left_miss;
        rightMiss(m,j) = ci(j).right_miss;
        isInf(m,j) = ci(j).infinite;
        lengths(m,j) = ci(j).length;
    end
    fallback(m) = diag.fallback_trigger;
    btInf(m) = diag.bootstrap_t_infinite;
    overlong(m) = diag.overlong_trigger;
    zeroSdCount(m) = diag.zero_sd_count;
    minPositiveSd(m) = diag.min_positive_bootstrap_sd;
    lengthRatio(m) = diag.bt_to_percentile_length_ratio;
    sampleSd(m) = diag.sample_sd;
end

rowCells = cell(nMethods, 1);
for j = 1:nMethods
    coverage = mean(cover(:,j));
    finiteLength = isfinite(lengths(:,j));
    [formulaBaseline, formulaBaselineType] = localFormulaBaseline(methodNames(j), cfg.alpha, cfg.B);
    firstOrderScaledError = sqrt(R) * (coverage - formulaBaseline);
    A_p1_hat = NaN;
    A_bt1_hat = NaN;
    if methodNames(j) == "Percentile bootstrap"
        A_p1_hat = firstOrderScaledError;
    elseif methodNames(j) == "Bootstrap-t"
        A_bt1_hat = firstOrderScaledError;
    end
    rowCells{j} = table(string(phase), string(experiment), string(pointMethod), ...
        methodNames(j), y0, hIndex, h, R, poolSize, cfg.M, cfg.B, truth, ...
        poolStats.mean, poolStats.mean - truth, poolStats.std, ...
        poolStats.skewness, poolStats.excess_kurtosis, signedBiasSeRatio, ...
        weightStats.mean_rho3_w, weightStats.mean_rho4_w, ...
        weightStats.mean_n2_eff_w, weightStats.mean_n3_eff_w, ...
        weightStats.mean_n4_eff_w, coverage, cfg.nominal, ...
        formulaBaseline, string(formulaBaselineType), ...
        coverage - formulaBaseline, abs(coverage - formulaBaseline), ...
        firstOrderScaledError, A_p1_hat, A_bt1_hat, ...
        sign(poolStats.skewness) * sign(coverage - formulaBaseline), ...
        sqrt(max(coverage * (1 - coverage), realmin) / cfg.M), ...
        mean(leftMiss(:,j)), mean(rightMiss(:,j)), sum(fallback), ...
        mean(fallback), sum(btInf), mean(btInf), sum(isInf(:,j)), ...
        mean(isInf(:,j)), mean(lengths(finiteLength,j), "omitnan"), ...
        median(lengths(finiteLength,j), "omitnan"), sum(finiteLength), ...
        sum(finiteLength) / cfg.M, blockSeed, ...
        'VariableNames', {'phase','experiment','point_method','method','y0', ...
        'h_index','h','R','pool_size','M','B','truth','pool_mean', ...
        'estimator_bias','pool_sd','pool_skewness','pool_excess_kurtosis', ...
        'bias_se_ratio','mean_rho3_w','mean_rho4_w','mean_n2_eff_w', ...
        'mean_n3_eff_w','mean_n4_eff_w','coverage','nominal_coverage','formula_baseline', ...
        'formula_baseline_type','coverage_error', ...
        'abs_coverage_error','first_order_scaled_error','A_p1_hat','A_bt1_hat', ...
        'skewness_coverage_sign_product','coverage_mc_se','left_miss','right_miss', ...
        'fallback_count','fallback_rate','bootstrap_t_inf_count', ...
        'bootstrap_t_inf_rate','interval_inf_count','interval_inf_rate', ...
        'mean_interval_length','median_interval_length','finite_length_count', ...
        'finite_length_rate','block_seed'});
end
rows = vertcat(rowCells{:});

mechanismRow = table(string(phase), string(experiment), string(pointMethod), ...
    y0, hIndex, h, R, poolSize, cfg.M, cfg.B, blockSeed, sum(fallback), ...
    mean(fallback), sum(btInf), mean(btInf), sum(overlong), mean(overlong), ...
    sum(zeroSdCount), sum(zeroSdCount) / max(cfg.M * cfg.B, 1), ...
    mean(minPositiveSd, "omitnan"), min(minPositiveSd, [], "omitnan"), ...
    median(lengthRatio, "omitnan"), localPercentile(lengthRatio, 0.95), ...
    mean(sampleSd), min(sampleSd), ...
    'VariableNames', {'phase','experiment','point_method','y0','h_index','h', ...
    'R','pool_size','M','B','block_seed','fallback_count','fallback_rate', ...
    'bootstrap_t_inf_count','bootstrap_t_inf_rate','overlong_count', ...
    'overlong_rate','zero_bootstrap_sd_count','zero_bootstrap_sd_rate', ...
    'mean_min_positive_bootstrap_sd','minimum_positive_bootstrap_sd', ...
    'median_bt_to_percentile_length_ratio','p95_bt_to_percentile_length_ratio', ...
    'mean_sample_sd','minimum_sample_sd'});

raw = struct("phase", string(phase), "experiment", string(experiment), ...
    "point_method", string(pointMethod), "y0", y0, "h", h, "h_index", hIndex, ...
    "R", R, "block_seed", blockSeed, "cover", cover, ...
    "interval_inf", isInf, "interval_lengths", single(lengths), ...
    "fallback", fallback, "bootstrap_t_inf", btInf, "overlong", overlong, ...
    "zero_bootstrap_sd_count", uint16(zeroSdCount), ...
    "min_positive_bootstrap_sd", single(minPositiveSd), ...
    "bt_to_percentile_length_ratio", single(lengthRatio), ...
    "sample_sd", single(sampleSd), "pool_stats", poolStats, ...
    "weight_stats", weightStats);
end

function [ci, diag] = localCiMethods(x, truth, alpha, B, lambda)
x = x(:);
R = numel(x);
mu = mean(x);
s = std(x, 0);
se = s / sqrt(R);
tCrit = tinv(1 - alpha / 2, max(R - 1, 1));
ci(1) = localPackCi("Student-t", mu - tCrit * se, mu + tCrit * se, truth, false);

bootIdx = randi(R, B, R);
boot = x(bootIdx);
bm = mean(boot, 2);
[kMinus, kPlus] = localFiniteBRanks(alpha, B);
bmSorted = sort(bm);
ci(2) = localPackCi("Percentile bootstrap", bmSorted(kMinus), bmSorted(kPlus), truth, false);

bs = std(boot, 0, 2);
T = sqrt(R) * (bm - mu) ./ bs;
T(bs <= 0 & bm == mu) = 0;
T(bs <= 0 & bm > mu) = Inf;
T(bs <= 0 & bm < mu) = -Inf;
T = sort(T);
btLower = mu - T(kPlus) * se;
btUpper = mu - T(kMinus) * se;
ci(3) = localPackCi("Bootstrap-t", btLower, btUpper, truth, false);

ratio = ci(3).length / max(ci(2).length, realmin);
overlong = ~ci(3).infinite && ratio > lambda;
trigger = ci(3).infinite || overlong;
if trigger
    ci(4) = localPackCi("Bootstrap-t fallback-rule", ci(2).lower, ci(2).upper, truth, true);
else
    ci(4) = localPackCi("Bootstrap-t fallback-rule", ci(3).lower, ci(3).upper, truth, false);
end

positive = bs(bs > 0 & isfinite(bs));
if isempty(positive)
    minPositive = NaN;
else
    minPositive = min(positive);
end
diag = struct("fallback_trigger", trigger, ...
    "bootstrap_t_infinite", ci(3).infinite, ...
    "overlong_trigger", overlong, "zero_sd_count", sum(bs <= 0), ...
    "min_positive_bootstrap_sd", minPositive, ...
    "bt_to_percentile_length_ratio", ratio, "sample_sd", s);
end

function ci = localPackCi(name, lower, upper, truth, fallbackTrigger)
ci = struct("name", string(name), "lower", lower, "upper", upper, ...
    "length", upper - lower, "infinite", ~(isfinite(lower) && isfinite(upper)), ...
    "contains", false, "left_miss", false, "right_miss", false, ...
    "fallback_trigger", fallbackTrigger);
ci.contains = ~ci.infinite && lower <= truth && truth <= upper;
ci.left_miss = ~ci.infinite && truth < lower;
ci.right_miss = ~ci.infinite && truth > upper;
end

function [baseline, baselineType] = localFormulaBaseline(methodName, alpha, B)
if string(methodName) == "Student-t"
    baseline = 1 - alpha;
    baselineType = "nominal_1_minus_alpha";
else
    [kMinus, kPlus] = localFiniteBRanks(alpha, B);
    baseline = (kPlus - kMinus) / (B + 1);
    baselineType = "finite_B_C0B";
end
end

function [kMinus, kPlus] = localFiniteBRanks(alpha, B)
kMinus = max(1, floor((alpha / 2) * (B + 1)));
kPlus = min(B, ceil((1 - alpha / 2) * (B + 1)));
if kPlus <= kMinus
    error("Invalid finite-B ranks: kMinus=%d, kPlus=%d, B=%d.", kMinus, kPlus, B);
end
end

function candidates = localSelectCandidates(tuning, cfg, expSpec, pointMethod)
raw = tuning(tuning.method == "Bootstrap-t", :);
fallbackRows = tuning(tuning.method == "Bootstrap-t fallback-rule", :);
rows = cell(0, 1);
yValues = double(expSpec.y_points(:)).';
for ir = 1:numel(cfg.R_list)
    R = cfg.R_list(ir);
    for iy = 1:numel(yValues)
        y0 = yValues(iy);
        scores = nan(numel(cfg.h_list), 1);
        rankScores = nan(numel(cfg.h_list), 1);
        stable = false(numel(cfg.h_list), 1);
        for ih = 1:numel(cfg.h_list)
            q = raw(raw.R == R & raw.h_index == ih & raw.y0 == y0, :);
            f = fallbackRows(fallbackRows.R == R & fallbackRows.h_index == ih & fallbackRows.y0 == y0, :);
            if isempty(q) || isempty(f)
                continue;
            end
            scores(ih) = mean(q.abs_coverage_error);
            rankScores(ih) = scores(ih) + cfg.selection_h_anchor_penalty * ...
                abs(log(max(cfg.h_list(ih), realmin) / expSpec.anchor_h));
            stable(ih) = max(q.bootstrap_t_inf_rate) <= cfg.selection_max_inf ...
                && max(f.fallback_rate) <= cfg.selection_max_fallback ...
                && min(q.coverage) >= cfg.selection_min_location_coverage ...
                && mean(q.abs_coverage_error) <= cfg.selection_max_mean_abs_error;
        end
        chosen = [];
        reasons = strings(0, 1);
        eligible = find(stable);
        if ~isempty(eligible)
            [~, loc] = min(rankScores(eligible));
            chosen(end+1) = eligible(loc); %#ok<AGROW>
            reasons(end+1,1) = "tuning_selected_stable"; %#ok<AGROW>
        end
        [~, best] = min(rankScores);
        if isfinite(rankScores(best))
            chosen(end+1) = best; %#ok<AGROW>
            reasons(end+1,1) = "tuning_min_abs_error"; %#ok<AGROW>
        end
        [~, anchor] = min(abs(cfg.h_list - expSpec.anchor_h));
        chosen(end+1) = anchor; %#ok<AGROW>
        reasons(end+1,1) = "prespecified_anchor"; %#ok<AGROW>
        [chosen, ia] = unique(chosen, "stable");
        reasons = reasons(ia);
        keep = 1:min(numel(chosen), cfg.confirmation_h_per_experiment);
        for k = keep
            ih = chosen(k);
            q = raw(raw.R == R & raw.h_index == ih & raw.y0 == y0, :);
            f = fallbackRows(fallbackRows.R == R & fallbackRows.h_index == ih & fallbackRows.y0 == y0, :);
            rows{end+1,1} = table(string(expSpec.name), string(pointMethod), R, y0, ih, ...
                cfg.h_list(ih), reasons(k), stable(ih), mean(q.coverage), ...
                min(q.coverage), mean(q.abs_coverage_error), rankScores(ih), max(f.fallback_rate), ...
                max(q.bootstrap_t_inf_rate), ...
                'VariableNames', {'experiment','point_method','R','y0','h_index','h', ...
                'selection_reason','tuning_stable','tuning_mean_coverage', ...
                'tuning_min_coverage','tuning_mean_abs_error','tuning_rank_score', ...
                'tuning_max_fallback_rate','tuning_max_inf_rate'});
        end
    end
end
candidates = vertcat(rows{:});
end

function output = localRunE8BandDiagnostic(cfg)
outDir = fullfile(cfg.results_root, "E8_plate_RFEM_simultaneous_band_weighted_diagnostic");
dpimnumeric.ensureDir(outDir);
dpimnumeric.ensureDir(fullfile(outDir, "figures"));
expSpec = struct("name", "E8_plate_RFEM_simultaneous_band_weighted_diagnostic", ...
    "d", 5, "y_points", [0.8, 1.0, 1.25, 1.5, 1.8], ...
    "anchor_h", 0.00183298071083244, "paper_role", "diagnostic_only");
problem = localPlateProblem(expSpec);
expCfg = localExperimentCfg(cfg, outDir);
dpimnumeric.writeJson(fullfile(outDir, "config.json"), localJsonSafe(expCfg, expSpec));

hCandidates = unique([1, numel(cfg.h_list), ...
    localNearestIndex(cfg.h_list, expSpec.anchor_h), ...
    localNearestIndex(cfg.h_list, 4e-4)]);
hCandidates = hCandidates(hCandidates >= 1 & hCandidates <= numel(cfg.h_list));
summary = table();
pointwise = table();
weightRows = table();
rawBlocks = cell(0, 1);
for im = 1:numel(cfg.methods)
    method = cfg.methods(im);
    fprintf("Building E8 weighted curve pool: method=%s, curves=%d, n=%d\n", ...
        method, cfg.curve_pool_size, cfg.n);
    curves = dpim_build_weighted_curve_pool(problem, method, cfg.n, expCfg, cfg.curve_pool_size);
    weightRows = [weightRows; localWeightTable(curves, expSpec.name, method)]; %#ok<AGROW>
    yPoints = double(expSpec.y_points(:)).';
    for ih = hCandidates
        h = cfg.h_list(ih);
        truth = zeros(1, numel(yPoints));
        poolMatrix = zeros(numel(curves), numel(yPoints));
        for iy = 1:numel(yPoints)
            truth(iy) = localTruth(problem, yPoints(iy), h, expCfg, 8000 * iy + ih);
            poolMatrix(:, iy) = dpim_curve_point_estimates(curves, yPoints(iy), h);
        end
        seed = localBlockSeed(cfg.seed, 8, im, 1, ih, 1);
        [bandRow, pointRows, raw] = localBandBlock(poolMatrix, truth, expCfg, ...
            expSpec.name, method, yPoints, h, ih, cfg.R_list(1), seed);
        summary = [summary; bandRow]; %#ok<AGROW>
        pointwise = [pointwise; pointRows]; %#ok<AGROW>
        rawBlocks{end+1,1} = raw; %#ok<AGROW>
    end
end
writetable(summary, fullfile(outDir, "summary.csv"));
writetable(pointwise, fullfile(outDir, "pointwise_band_by_h.csv"));
writetable(weightRows, fullfile(outDir, "weight_diagnostics.csv"));
save(fullfile(outDir, "raw_results.mat"), "summary", "pointwise", ...
    "weightRows", "rawBlocks", "cfg", "expSpec", "-v7.3");
output = struct("output_dir", outDir, "summary_csv", fullfile(outDir, "summary.csv"));
end

function problem = localPlateProblem(expSpec)
problem = struct();
problem.name = char(expSpec.name);
problem.short_name = "E8";
problem.d = expSpec.d;
problem.target_distribution = "standard_normal";
problem.center_transform = "normal_icdf";
problem.provider = "voronoi_ci_probability_weights_provider";
coeff = 0.20 * exp(-0.30 * ((1:problem.d)' - 1));
problem.response_fun = @(theta) 1.0 + theta * coeff(:) ...
    + 0.08 * theta(:,1).^2 - 0.05 * theta(:,2) .* theta(:,3);
end

function [summary, pointwise, raw] = localBandBlock(poolMatrix, truth, cfg, ...
    experiment, pointMethod, yPoints, h, hIndex, R, blockSeed)
rng(blockSeed, "twister");
nGrid = numel(yPoints);
hit = false(cfg.M, 1);
infFlag = false(cfg.M, 1);
width = nan(cfg.M, 1);
pointHit = false(cfg.M, nGrid);
for m = 1:cfg.M
    idx = randi(size(poolMatrix, 1), R, 1);
    X = poolMatrix(idx, :);
    mu = mean(X, 1);
    se = std(X, 0, 1) / sqrt(R);
    bootIdx = randi(R, cfg.B, R);
    maxStats = nan(cfg.B, 1);
    for b = 1:cfg.B
        xb = X(bootIdx(b,:), :);
        s = std(xb, 0, 1);
        t = sqrt(R) * (mean(xb, 1) - mu) ./ s;
        t(~isfinite(t)) = Inf;
        maxStats(b) = max(abs(t));
    end
    crit = quantile(maxStats, 1 - cfg.alpha);
    lower = mu - crit * se;
    upper = mu + crit * se;
    infFlag(m) = ~isfinite(crit) || any(~isfinite(lower) | ~isfinite(upper));
    thisHit = lower <= truth & truth <= upper;
    pointHit(m,:) = ~infFlag(m) & thisHit;
    hit(m) = ~infFlag(m) && all(thisHit);
    if ~infFlag(m)
        width(m) = mean(upper - lower, "omitnan");
    end
end
coverage = mean(hit);
summary = table(string(experiment), string(pointMethod), ...
    "Bootstrap-t simultaneous-band", hIndex, h, R, cfg.M, cfg.B, coverage, ...
    cfg.nominal, coverage - cfg.nominal, abs(coverage - cfg.nominal), ...
    sum(infFlag), mean(infFlag), mean(width, "omitnan"), ...
    median(width, "omitnan"), nGrid, blockSeed, ...
    'VariableNames', {'experiment','point_method','method','h_index','h', ...
    'R','M','B','coverage','nominal_coverage','coverage_error', ...
    'abs_coverage_error','bootstrap_t_inf_count','bootstrap_t_inf_rate', ...
    'mean_band_length','median_band_length','grid_count','block_seed'});
pointwise = table(repmat(string(experiment), nGrid, 1), ...
    repmat(string(pointMethod), nGrid, 1), yPoints(:), repmat(hIndex, nGrid, 1), ...
    repmat(h, nGrid, 1), truth(:), mean(pointHit, 1)', ...
    repmat(sum(infFlag), nGrid, 1), repmat(mean(infFlag), nGrid, 1), ...
    'VariableNames', {'experiment','point_method','y0','h_index','h','truth', ...
    'pointwise_band_coverage','invalid_inf_count','invalid_inf_rate'});
raw = struct("hit", hit, "infFlag", infFlag, "pointHit", pointHit, ...
    "h", h, "h_index", hIndex, "R", R, "yPoints", yPoints, "truth", truth);
end

function output = localIndexOriginalPlate(cfg)
outPath = fullfile(cfg.results_root, "original_plate_reference.csv");
base = fullfile(cfg.project_root, "results");
runs = dir(fullfile(base, cfg.original_plate_glob));
if isempty(runs)
    T = table("original_plate", "missing", "", "", ...
        'VariableNames', {'source','status','results_root','summary_root'});
    writetable(T, outPath);
    output = struct("csv", outPath);
    return;
end
[~, idx] = max([runs.datenum]);
root = fullfile(runs(idx).folder, runs(idx).name);
summaryRoot = fullfile(root, "_summary");
status = "indexed_existing_result";
T = table("original_plate", status, string(root), string(summaryRoot), ...
    isfile(fullfile(summaryRoot, "overall_method_summary.csv")), ...
    isfile(fullfile(summaryRoot, "best_h_active_grid.csv")), ...
    'VariableNames', {'source','status','results_root','summary_root', ...
    'has_overall_method_summary','has_best_h_active_grid'});
writetable(T, outPath);
output = struct("csv", outPath, "results_root", root);
end

function W = localWeightTable(curves, experimentName, method)
rows = cell(numel(curves), 1);
for i = 1:numel(curves)
    rows{i} = dpim_weight_summary_row(curves(i), experimentName, method);
    if isfield(curves(i).weightData, "assignment_backend")
        rows{i}.assignment_backend = string(curves(i).weightData.assignment_backend);
    else
        rows{i}.assignment_backend = "unknown";
    end
end
W = struct2table([rows{:}]);
end

function stats = localAggregateWeightStats(W)
stats = struct( ...
    "mean_rho3_w", mean(double(W.rho3_w), "omitnan"), ...
    "mean_rho4_w", mean(double(W.rho4_w), "omitnan"), ...
    "mean_n2_eff_w", mean(double(W.n2_eff_w), "omitnan"), ...
    "mean_n3_eff_w", mean(double(W.n3_eff_w), "omitnan"), ...
    "mean_n4_eff_w", mean(double(W.n4_eff_w), "omitnan"));
end

function stats = localVectorStats(x)
x = double(x(:));
mu = mean(x);
c = x - mu;
variance = mean(c.^2);
sd = sqrt(max(variance, realmin));
stats = struct("mean", mu, "variance", variance, "std", sd, ...
    "skewness", mean(c.^3) / sd^3, ...
    "excess_kurtosis", mean(c.^4) / sd^4 - 3);
end

function BList = localProtocolBList(cfg, experimentName)
if string(experimentName) == "E2_scalar_normal_R_order_weighted"
    BList = unique([cfg.protocol_B_list(cfg.protocol_B_list <= cfg.B), cfg.B]);
else
    BList = cfg.B;
end
BList = sort(double(BList(:))).';
end

function T = localGaussianMomentValidation(yPoints, h, cfg, experimentName)
N = max(10000, double(cfg.truth_N));
rng(cfg.seed + 220001, "twister");
z = randn(N, 1);
rows = cell(numel(yPoints), 1);
for iy = 1:numel(yPoints)
    y0 = yPoints(iy);
    kernelValues = dpim_gaussian_kernel(y0 - z, h);
    analytic = zeros(1, 4);
    numeric = zeros(1, 4);
    simulated = zeros(1, 4);
    simulatedSe = zeros(1, 4);
    for j = 1:4
        analytic(j) = localGaussianKernelRawMoment(y0, h, 1, j);
        integrand = @(x) exp(-0.5 * x.^2) ./ sqrt(2*pi) ...
            .* dpim_gaussian_kernel(y0 - x, h).^j;
        numeric(j) = integral(integrand, -Inf, Inf, ...
            "AbsTol", 1e-12, "RelTol", 1e-10);
        values = kernelValues.^j;
        simulated(j) = mean(values);
        simulatedSe(j) = std(values, 0) / sqrt(N);
    end
    rows{iy} = table(string(experimentName), y0, h, N, ...
        analytic(1), numeric(1), simulated(1), simulatedSe(1), ...
        analytic(2), numeric(2), simulated(2), simulatedSe(2), ...
        analytic(3), numeric(3), simulated(3), simulatedSe(3), ...
        analytic(4), numeric(4), simulated(4), simulatedSe(4), ...
        max(abs(analytic - numeric)), ...
        max(abs(analytic - simulated) ./ max(simulatedSe, realmin)), ...
        "Gaussian response; analytic formula, adaptive quadrature, and independent Monte Carlo.", ...
        'VariableNames', {'experiment','y0','h','simulation_count', ...
        'A1_analytic','A1_numeric','A1_simulated','A1_simulation_se', ...
        'A2_analytic','A2_numeric','A2_simulated','A2_simulation_se', ...
        'A3_analytic','A3_numeric','A3_simulated','A3_simulation_se', ...
        'A4_analytic','A4_numeric','A4_simulated','A4_simulation_se', ...
        'max_abs_analytic_numeric_error','max_simulation_z_error','validation_note'});
end
T = vertcat(rows{:});
end

function value = localGaussianKernelRawMoment(y0, h, sigma2, order)
value = (2*pi)^(-order/2) * h^(-order) ...
    / sqrt(1 + order * sigma2 / h^2) ...
    * exp(-order * y0^2 / (2 * (h^2 + order * sigma2)));
end

function localPlotPointwise(outDir, tuning, confirmation, name)
if isempty(tuning)
    return;
end
fig = figure("Visible", "off", "Color", "w");
tiledlayout(1, 2, "TileSpacing", "compact");
nexttile;
localPlotPhase(tuning(tuning.method == "Bootstrap-t", :), "tuning");
nexttile;
localPlotPhase(confirmation(confirmation.method == "Bootstrap-t", :), "confirmation");
sgtitle(name, "Interpreter", "none");
try
    exportgraphics(fig, fullfile(outDir, "figures", "bootstrap_t_coverage_by_h.png"), "Resolution", 180);
catch
    saveas(fig, fullfile(outDir, "figures", "bootstrap_t_coverage_by_h.png"));
end
close(fig);
end

function localPlotPhase(T, ttl)
if isempty(T)
    title(ttl);
    return;
end
G = groupsummary(T, "h", "mean", "coverage");
semilogx(G.h, G.mean_coverage, "-o", "LineWidth", 1.2);
yline(0.95, "k--");
ylim([0 1]);
grid on;
xlabel("h");
ylabel("coverage");
title(ttl);
end

function seed = localBlockSeed(baseSeed, phaseCode, im, iy, ih, ir)
if im >= 10 || iy >= 10 || ih >= 1000 || ir >= 100
    error("Seed index exceeds the collision-free packing range.");
end
seed = baseSeed + 10000000 * phaseCode + 1000000 * im ...
    + 100000 * iy + 100 * ih + ir;
if seed > double(intmax("uint32"))
    error("Packed seed exceeds the uint32 range accepted by MATLAB RNG.");
end
end

function idx = localNearestIndex(values, target)
[~, idx] = min(abs(values - target));
end

function q = localPercentile(x, p)
x = x(isfinite(x));
if isempty(x)
    q = NaN;
else
    q = quantile(x, p);
end
end
