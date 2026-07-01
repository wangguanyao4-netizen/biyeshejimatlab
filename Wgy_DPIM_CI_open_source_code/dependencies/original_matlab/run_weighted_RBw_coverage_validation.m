function outputs = run_weighted_RBw_coverage_validation(runMode, resultsRoot)
%RUN_WEIGHTED_RBW_COVERAGE_VALIDATION Validate probability-weighted coverage terms.
%
% This script targets the paper formula labelled
% thm:weighted-R-B-w-coverage. It is not a replacement for the main paper
% evidence suite. It creates controlled probability-weighted curve pools,
% reports s2/s3/s4/rho3/rho4, and checks empirical coverage against the
% formula baselines C0B and the explicit Student-t second-order term.

if nargin < 1 || strlength(string(runMode)) == 0
    runMode = "diagnostic";
end

projectRoot = fileparts(mfilename("fullpath"));
if nargin < 2 || strlength(string(resultsRoot)) == 0
    stamp = string(datetime("now", "Format", "yyyyMMdd_HHmmss"));
    resultsRoot = fullfile(projectRoot, "results", ...
        "weighted_RBw_coverage_validation_" + lower(string(runMode)) + "_" + stamp);
end

cfg = weighted_paper_config(runMode, projectRoot, resultsRoot);
localAddPaths(cfg);
vcfg = localValidationConfig(cfg, runMode);

dpimnumeric.ensureDir(resultsRoot);
dpimnumeric.ensureDir(fullfile(resultsRoot, "scheme_details"));
dpimnumeric.ensureDir(fullfile(resultsRoot, "figures"));
dpimnumeric.writeJson(fullfile(resultsRoot, "validation_config.json"), vcfg);

problem = localScalarNormalProblem();
schemes = localWeightSchemes(vcfg.n, vcfg.seed, vcfg.power_list);
[kMinus, kPlus, C0B] = localFiniteBBaseline(vcfg.alpha, vcfg.B);

coverageRows = table();
weightRows = table();
schemeRows = table();

for rep = 1:vcfg.replicate_count
    for i = 1:numel(schemes)
        scheme = schemes(i);
        fprintf("Building RBw validation pool: replicate=%d/%d, scheme=%s, curves=%d, n=%d\n", ...
            rep, vcfg.replicate_count, scheme.name, vcfg.curve_pool_size, vcfg.n);
        curves = localBuildCurves(problem, scheme, cfg, vcfg, rep);
        W = localWeightTable(curves, scheme.name, vcfg.point_method, rep);
        detailName = sprintf("%s_rep%02d_weight_diagnostics.csv", scheme.name, rep);
        writetable(W, fullfile(resultsRoot, "scheme_details", detailName));
        weightRows = [weightRows; W]; %#ok<AGROW>
        schemeRows = [schemeRows; localSchemeSummary(W, scheme.name, scheme.kind, rep)]; %#ok<AGROW>

        for ih = 1:numel(vcfg.h_list)
            h = vcfg.h_list(ih);
            for iy = 1:numel(vcfg.y_points)
                y0 = vcfg.y_points(iy);
                poolVals = dpim_curve_point_estimates(curves, y0, h);
                truth = dpim_gaussian_kernel(y0, sqrt(1 + h^2));
                poolStats = localPoolStats(poolVals);
                for ir = 1:numel(vcfg.R_list)
                    R = vcfg.R_list(ir);
                    blockSeed = vcfg.seed + 100000000 * rep + 1000000 * i + 10000 * ih + 100 * iy + ir;
                    block = localCoverageBlock(poolVals, truth, poolStats, W, vcfg, ...
                        scheme, rep, y0, h, ih, R, blockSeed, C0B, kMinus, kPlus);
                    coverageRows = [coverageRows; block]; %#ok<AGROW>
                end
            end
        end
    end
end

writetable(weightRows, fullfile(resultsRoot, "weight_diagnostics.csv"));
writetable(schemeRows, fullfile(resultsRoot, "weight_moment_summary.csv"));
writetable(coverageRows, fullfile(resultsRoot, "coverage_formula_validation.csv"));
replicateStability = localReplicateStabilityTable(coverageRows);
writetable(replicateStability, fullfile(resultsRoot, "weighted_RBw_replicate_stability.csv"));

[fitSummary, fitCoefficients] = localRegressionTables(coverageRows);
fitPredictions = localPredictionTables(coverageRows);
writetable(fitSummary, fullfile(resultsRoot, "weighted_RBw_fit_summary.csv"));
writetable(fitCoefficients, fullfile(resultsRoot, "weighted_RBw_fit_coefficients.csv"));
writetable(fitPredictions, fullfile(resultsRoot, "weighted_RBw_fit_predictions.csv"));
localPlotValidationFigures(resultsRoot, schemeRows, fitPredictions);

reportPath = fullfile(resultsRoot, "weighted_RBw_coverage_validation_report.md");
localWriteReport(reportPath, vcfg, C0B, kMinus, kPlus, schemeRows, coverageRows, fitSummary);

outputs = struct();
outputs.results_root = resultsRoot;
outputs.coverage_csv = fullfile(resultsRoot, "coverage_formula_validation.csv");
outputs.weight_summary_csv = fullfile(resultsRoot, "weight_moment_summary.csv");
outputs.replicate_stability_csv = fullfile(resultsRoot, "weighted_RBw_replicate_stability.csv");
outputs.fit_summary_csv = fullfile(resultsRoot, "weighted_RBw_fit_summary.csv");
outputs.fit_predictions_csv = fullfile(resultsRoot, "weighted_RBw_fit_predictions.csv");
outputs.report_md = reportPath;
fprintf("Weighted R-B-w coverage validation completed: %s\n", resultsRoot);
end

function localAddPaths(cfg)
addpath(cfg.project_root, "-begin");
addpath(fullfile(cfg.weighted_root, "external_weight_providers"), "-begin");
addpath(fullfile(cfg.weighted_root, "common"), "-begin");
end

function vcfg = localValidationConfig(cfg, runMode)
runMode = lower(string(runMode));
vcfg = struct();
vcfg.schema_version = "weighted_RBw_coverage_validation_v1";
vcfg.run_mode = char(runMode);
vcfg.project_root = cfg.project_root;
vcfg.results_root = cfg.results_root;
vcfg.seed = cfg.seed + 770001;
vcfg.alpha = cfg.alpha;
vcfg.lambda = cfg.lambda;
vcfg.point_method = cfg.main_methods;
vcfg.assignment_backend = cfg.assignment_backend;
switch runMode
    case "formal"
        vcfg.n = cfg.n;
        vcfg.curve_pool_size = max(cfg.curve_pool_size, 1200);
        vcfg.R_list = [16, 32, 64, 128];
        vcfg.M = 1000;
        vcfg.B = 999;
        vcfg.y_points = [-1, 0, 1];
        vcfg.h_list = [0.10, 0.25, 0.50];
        vcfg.replicate_count = 8;
        vcfg.power_list = [0.25, 0.50, 0.75, 1.00, 1.25, 1.50, 2.00];
    case "full"
        vcfg.n = cfg.n;
        vcfg.curve_pool_size = max(cfg.curve_pool_size, 600);
        vcfg.R_list = cfg.R_list;
        vcfg.M = cfg.M;
        vcfg.B = cfg.B;
        vcfg.y_points = [-1, 0, 1];
        vcfg.h_list = [0.10, 0.25, 0.50];
        vcfg.replicate_count = 5;
        vcfg.power_list = [0.25, 0.50, 0.75, 1.00, 1.25, 1.50];
    case "medium"
        vcfg.n = 384;
        vcfg.curve_pool_size = 360;
        vcfg.R_list = [20, 40, 80];
        vcfg.M = 180;
        vcfg.B = 199;
        vcfg.y_points = [-1, 0, 1];
        vcfg.h_list = [0.15, 0.35];
        vcfg.replicate_count = 3;
        vcfg.power_list = [0.25, 0.50, 0.75, 1.00, 1.25, 1.50];
    case "pilot"
        vcfg.n = 256;
        vcfg.curve_pool_size = 300;
        vcfg.R_list = [12, 24, 48];
        vcfg.M = 160;
        vcfg.B = 199;
        vcfg.y_points = [-1, 0, 1];
        vcfg.h_list = [0.15, 0.35];
        vcfg.replicate_count = 2;
        vcfg.power_list = [0.50, 0.75, 1.00, 1.50];
    case "small"
        vcfg.n = 192;
        vcfg.curve_pool_size = 240;
        vcfg.R_list = [12, 24, 48];
        vcfg.M = 120;
        vcfg.B = 199;
        vcfg.y_points = [-1, 0, 1];
        vcfg.h_list = 0.25;
        vcfg.replicate_count = 2;
        vcfg.power_list = [0.50, 0.75, 1.00, 1.50];
    otherwise
        vcfg.n = 96;
        vcfg.curve_pool_size = 120;
        vcfg.R_list = [8, 16, 32];
        vcfg.M = 60;
        vcfg.B = 99;
        vcfg.y_points = [-1, 0, 1];
        vcfg.h_list = 0.25;
        vcfg.replicate_count = 1;
        vcfg.power_list = [0.75, 1.50];
end
vcfg.voronoi_aux_sample_count = min(cfg.voronoi_aux_sample_count, max(5000, 40 * vcfg.n));
vcfg.voronoi_block_size = min(cfg.voronoi_block_size, 1000);
vcfg.note = "Uses explicit finite-B order statistic endpoints; bootstrap methods are compared to C0B, Student-t to its nominal and explicit second-order diagnostic.";
end

function problem = localScalarNormalProblem()
problem = struct();
problem.name = "RBw_scalar_normal_weighted_validation";
problem.short_name = "RBwScalarNormal";
problem.d = 1;
problem.distribution = "standard_normal";
problem.target_distribution = "standard_normal";
problem.center_transform = "normal_icdf";
problem.provider = "voronoi_ci_probability_weights_provider";
problem.response_fun = @(theta) theta(:,1);
end

function schemes = localWeightSchemes(n, seed, powers)
schemeCount = numel(powers) + 2;
schemes = repmat(struct("name", "", "kind", "", "weights", [], "power", NaN), schemeCount, 1);
schemes(1) = struct("name", "equal_fixed", "kind", "fixed", ...
    "weights", ones(n, 1) / n, "power", 0);
for i = 1:numel(powers)
    p = powers(i);
    name = "power_fixed_" + replace(string(sprintf("%.2f", p)), ".", "p");
    schemes(i + 1) = struct("name", name, "kind", "fixed", ...
        "weights", localPowerWeights(n, p, seed + 11 + 37 * i), "power", p);
end
schemes(end) = struct("name", "voronoi_rqmc", "kind", "voronoi", ...
    "weights", [], "power", NaN);
end

function w = localPowerWeights(n, p, seed)
base = ((1:n)').^(-p);
base = base / sum(base);
oldState = rng;
rng(seed, "twister");
cleanupObj = onCleanup(@() rng(oldState)); %#ok<NASGU>
perm = randperm(n);
w = base(perm);
end

function curves = localBuildCurves(problem, scheme, cfg, vcfg, replicateIndex)
if scheme.kind == "voronoi"
    expCfg = cfg;
    expCfg.seed = vcfg.seed + 50000 + 10000000 * replicateIndex;
    expCfg.weighting_cfg.voronoi_aux_sample_count = vcfg.voronoi_aux_sample_count;
    expCfg.weighting_cfg.voronoi_block_size = vcfg.voronoi_block_size;
    expCfg.weighting_cfg.voronoi_assignment_backend = vcfg.assignment_backend;
    expCfg.weighting_cfg.voronoi_output_dir = fullfile(vcfg.results_root, ...
        "scheme_details", sprintf("voronoi_probability_weights_rep%02d", replicateIndex));
    curves = dpim_build_weighted_curve_pool(problem, vcfg.point_method, ...
        vcfg.n, expCfg, vcfg.curve_pool_size);
    return;
end

curves(1,1) = localFixedCurve(problem, scheme, vcfg, replicateIndex, 1);
for r = 2:vcfg.curve_pool_size
    curves(r,1) = localFixedCurve(problem, scheme, vcfg, replicateIndex, r);
end
end

function curve = localFixedCurve(problem, scheme, vcfg, replicateIndex, curveIndex)
seed = vcfg.seed + 10000000 * replicateIndex + 100000 * curveIndex ...
    + 7919 * double(sum(char(scheme.name)));
[U, pointInfo] = dpim_generate_point_set(vcfg.n, problem.d, vcfg.point_method, seed);
theta = localNormInv(U);
y = problem.response_fun(theta);
w = double(scheme.weights(:));
curve = struct();
curve.problem = problem.name;
curve.method = char(vcfg.point_method);
curve.curve_index = curveIndex;
curve.pointInfo = pointInfo;
curve.point_requested_method = char(pointInfo.requested_method);
curve.point_actual_method = char(pointInfo.actual_method);
curve.point_fallback_used = logical(pointInfo.fallback_used);
curve.point_message = char(pointInfo.message);
curve.U = U;
curve.theta = theta;
curve.y = y(:);
curve.weights = w;
curve.weightData = struct("weight_source", "fixed:" + string(scheme.name), ...
    "empty_cell_count", 0, "auxiliary_sample_count", 0, ...
    "assignment_backend", "not_applicable");
curve.weight_source = char("fixed:" + string(scheme.name));
curve.sum_weights = sum(w);
curve.weight_ess = 1 / sum(w.^2);
curve.weight_cv = std(w) / mean(w);
curve.l1_from_equal = sum(abs(w - 1 / numel(w)));
end

function x = localNormInv(u)
u = min(max(double(u), 1e-15), 1 - 1e-15);
x = sqrt(2) * erfinv(2 * u - 1);
end

function W = localWeightTable(curves, schemeName, pointMethod, replicateIndex)
rows = cell(numel(curves), 1);
for i = 1:numel(curves)
    rows{i} = dpim_weight_summary_row(curves(i), schemeName, pointMethod);
    rows{i}.replicate_index = replicateIndex;
    if isfield(curves(i).weightData, "assignment_backend")
        rows{i}.assignment_backend = string(curves(i).weightData.assignment_backend);
    else
        rows{i}.assignment_backend = "unknown";
    end
end
W = struct2table([rows{:}]);
end

function row = localSchemeSummary(W, schemeName, schemeKind, replicateIndex)
row = table(string(schemeName), string(schemeKind), replicateIndex, height(W), W.n(1), ...
    mean(W.s2_w), std(W.s2_w), mean(W.s3_w), std(W.s3_w), ...
    mean(W.s4_w), std(W.s4_w), mean(W.rho3_w), std(W.rho3_w), ...
    mean(W.rho4_w), std(W.rho4_w), mean(W.n2_eff_w), ...
    mean(W.n3_eff_w), mean(W.n4_eff_w), mean(W.ess_ratio), ...
    mean(W.max_over_equal), mean(W.l1_from_equal), ...
    'VariableNames', {'scheme','scheme_kind','replicate_index','curve_count','n', ...
    'mean_s2_w','sd_s2_w','mean_s3_w','sd_s3_w','mean_s4_w','sd_s4_w', ...
    'mean_rho3_w','sd_rho3_w','mean_rho4_w','sd_rho4_w', ...
    'mean_n2_eff_w','mean_n3_eff_w','mean_n4_eff_w', ...
    'mean_ess_ratio','mean_max_over_equal','mean_l1_from_equal'});
end

function stats = localPoolStats(x)
x = double(x(:));
mu = mean(x);
centered = x - mu;
v = mean(centered.^2);
sd = sqrt(max(v, realmin));
stats = struct();
stats.mean = mu;
stats.variance = v;
stats.sd = sd;
stats.skewness = mean(centered.^3) / sd^3;
stats.excess_kurtosis = mean(centered.^4) / sd^4 - 3;
end

function rows = localCoverageBlock(poolVals, truth, poolStats, W, vcfg, scheme, ...
    replicateIndex, y0, h, hIndex, R, blockSeed, C0B, kMinus, kPlus)
poolVals = double(poolVals(:));
rng(blockSeed, "twister");
methodNames = ["Student-t"; "Percentile bootstrap"; "Bootstrap-t"; "Bootstrap-t fallback-rule"];
nMethods = numel(methodNames);
cover = false(vcfg.M, nMethods);
leftMiss = false(vcfg.M, nMethods);
rightMiss = false(vcfg.M, nMethods);
isInf = false(vcfg.M, nMethods);
lengths = nan(vcfg.M, nMethods);
fallback = false(vcfg.M, 1);
btInf = false(vcfg.M, 1);
zeroSdCount = zeros(vcfg.M, 1);
poolSize = numel(poolVals);

for m = 1:vcfg.M
    idx = randi(poolSize, R, 1);
    x = poolVals(idx);
    [ci, diag] = localCiMethodsOrder(x, truth, vcfg.alpha, vcfg.B, ...
        vcfg.lambda, kMinus, kPlus);
    for j = 1:nMethods
        cover(m,j) = ci(j).contains;
        leftMiss(m,j) = ci(j).left_miss;
        rightMiss(m,j) = ci(j).right_miss;
        isInf(m,j) = ci(j).infinite;
        lengths(m,j) = ci(j).length;
    end
    fallback(m) = diag.fallback_trigger;
    btInf(m) = diag.bootstrap_t_infinite;
    zeroSdCount(m) = diag.zero_sd_count;
end

z = localNormInv(1 - vcfg.alpha / 2);
phi = exp(-0.5 * z^2) / sqrt(2 * pi);
tEdgeworthPrediction = 1 - vcfg.alpha + (2 * z * phi / R) * ...
    (((z^2 - 3) / 12) * poolStats.excess_kurtosis ...
    - ((z^4 + 2 * z^2 - 3) / 18) * poolStats.skewness^2);

rho3 = mean(W.rho3_w);
rho4 = mean(W.rho4_w);
finiteLength = isfinite(lengths);
rowCells = cell(nMethods, 1);
for j = 1:nMethods
    coverage = mean(cover(:,j));
    if methodNames(j) == "Student-t"
        baseline = 1 - vcfg.alpha;
        baselineType = "nominal_1_minus_alpha";
        explicitPrediction = tEdgeworthPrediction;
    else
        baseline = C0B;
        baselineType = "finite_B_C0B";
        explicitPrediction = NaN;
    end
    rowCells{j} = table(string(scheme.name), string(scheme.kind), ...
        replicateIndex, string(vcfg.point_method), methodNames(j), y0, hIndex, h, R, ...
        poolSize, vcfg.M, vcfg.B, kMinus, kPlus, C0B, truth, ...
        poolStats.mean, poolStats.mean - truth, poolStats.skewness, ...
        poolStats.excess_kurtosis, rho3, rho4, rho3 / sqrt(R), ...
        rho4 / R, rho3^2 / R, poolStats.skewness / sqrt(R), ...
        poolStats.excess_kurtosis / R, poolStats.skewness^2 / R, ...
        coverage, baseline, string(baselineType), coverage - baseline, ...
        abs(coverage - baseline), coverage - C0B, abs(coverage - C0B), ...
        explicitPrediction, coverage - explicitPrediction, ...
        sqrt(max(coverage * (1 - coverage), realmin) / vcfg.M), ...
        mean(leftMiss(:,j)), mean(rightMiss(:,j)), mean(fallback), ...
        mean(btInf), mean(isInf(:,j)), mean(lengths(finiteLength(:,j),j), "omitnan"), ...
        median(lengths(finiteLength(:,j),j), "omitnan"), sum(finiteLength(:,j)), ...
        sum(finiteLength(:,j)) / vcfg.M, sum(zeroSdCount), ...
        sum(zeroSdCount) / max(vcfg.M * vcfg.B, 1), blockSeed, ...
        'VariableNames', {'scheme','scheme_kind','replicate_index','point_method','method', ...
        'y0','h_index','h','R','pool_size','M','B','k_minus','k_plus', ...
        'C0B','truth','pool_mean','estimator_bias','pool_skewness', ...
        'pool_excess_kurtosis','mean_rho3_w','mean_rho4_w', ...
        'rho3_over_sqrtR','rho4_over_R','rho3sq_over_R', ...
        'pool_skewness_over_sqrtR','pool_kurtosis_over_R', ...
        'pool_skewness_sq_over_R','coverage','formula_baseline', ...
        'formula_baseline_type','coverage_error_to_formula_baseline', ...
        'abs_error_to_formula_baseline','coverage_error_to_C0B', ...
        'abs_error_to_C0B','student_t_edgeworth_prediction', ...
        'student_t_edgeworth_residual','coverage_mc_se','left_miss', ...
        'right_miss','fallback_rate','bootstrap_t_inf_rate', ...
        'interval_inf_rate','mean_interval_length','median_interval_length', ...
        'finite_length_count','finite_length_rate','zero_bootstrap_sd_count', ...
        'zero_bootstrap_sd_rate','block_seed'});
end
rows = vertcat(rowCells{:});
end

function [ci, diag] = localCiMethodsOrder(x, truth, alpha, B, lambda, kMinus, kPlus)
x = double(x(:));
R = numel(x);
mu = mean(x);
s = std(x, 0);
se = s / sqrt(R);
tCrit = tinv(1 - alpha / 2, max(R - 1, 1));
ci(1) = localPackCi("Student-t", mu - tCrit * se, mu + tCrit * se, truth, false);

bootIdx = randi(R, B, R);
boot = x(bootIdx);
bm = mean(boot, 2);
bmSorted = sort(bm);
ci(2) = localPackCi("Percentile bootstrap", bmSorted(kMinus), bmSorted(kPlus), truth, false);

bs = std(boot, 0, 2);
T = sqrt(R) * (bm - mu) ./ bs;
T(bs <= 0 & bm == mu) = 0;
T(bs <= 0 & bm > mu) = Inf;
T(bs <= 0 & bm < mu) = -Inf;
TSorted = sort(T);
btLower = mu - TSorted(kPlus) * se;
btUpper = mu - TSorted(kMinus) * se;
ci(3) = localPackCi("Bootstrap-t", btLower, btUpper, truth, false);

ratio = ci(3).length / max(ci(2).length, realmin);
overlong = ~ci(3).infinite && ratio > lambda;
trigger = ci(3).infinite || overlong;
if trigger
    ci(4) = localPackCi("Bootstrap-t fallback-rule", ci(2).lower, ci(2).upper, truth, true);
else
    ci(4) = localPackCi("Bootstrap-t fallback-rule", ci(3).lower, ci(3).upper, truth, false);
end

diag = struct("fallback_trigger", trigger, ...
    "bootstrap_t_infinite", ci(3).infinite, ...
    "zero_sd_count", sum(bs <= 0));
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

function [kMinus, kPlus, C0B] = localFiniteBBaseline(alpha, B)
kMinus = max(1, floor((alpha / 2) * (B + 1)));
kPlus = min(B, ceil((1 - alpha / 2) * (B + 1)));
if kPlus <= kMinus
    error("Invalid finite-B endpoint ranks: kMinus=%d, kPlus=%d, B=%d.", kMinus, kPlus, B);
end
C0B = (kPlus - kMinus) / (B + 1);
end

function [fitSummary, fitCoefficients] = localRegressionTables(T)
summaryRows = cell(0, 1);
coefRows = cell(0, 1);

[summaryRows, coefRows] = localAddFit(T, summaryRows, coefRows, ...
    "Percentile bootstrap", "Cp_minus_C0B_weight_terms", ...
    ["intercept", "rho3_over_sqrtR", "rho4_over_R", "rho3sq_over_R"], ...
    ["rho3_over_sqrtR", "rho4_over_R", "rho3sq_over_R"], "coverage_error_to_C0B");

[summaryRows, coefRows] = localAddFit(T, summaryRows, coefRows, ...
    "Percentile bootstrap", "Cp_minus_C0B_empirical_cumulants", ...
    ["intercept", "pool_skewness_over_sqrtR", "pool_kurtosis_over_R", ...
    "pool_skewness_sq_over_R"], ...
    ["pool_skewness_over_sqrtR", "pool_kurtosis_over_R", ...
    "pool_skewness_sq_over_R"], "coverage_error_to_C0B");

[summaryRows, coefRows] = localAddStratifiedFit(T, summaryRows, coefRows, ...
    "Percentile bootstrap", "Cp_minus_C0B_empirical_cumulants_yh_fixed", ...
    ["pool_skewness_over_sqrtR", "pool_kurtosis_over_R", ...
    "pool_skewness_sq_over_R"], "coverage_error_to_C0B");

[summaryRows, coefRows] = localAddFit(T, summaryRows, coefRows, ...
    "Bootstrap-t", "Cbt_minus_C0B_weight_terms", ...
    ["intercept", "rho4_over_R", "rho3sq_over_R", "fallback_rate"], ...
    ["rho4_over_R", "rho3sq_over_R", "fallback_rate"], "coverage_error_to_C0B");

[summaryRows, coefRows] = localAddFit(T, summaryRows, coefRows, ...
    "Bootstrap-t", "Cbt_minus_C0B_empirical_cumulants", ...
    ["intercept", "pool_kurtosis_over_R", "pool_skewness_sq_over_R", ...
    "fallback_rate"], ...
    ["pool_kurtosis_over_R", "pool_skewness_sq_over_R", "fallback_rate"], ...
    "coverage_error_to_C0B");

[summaryRows, coefRows] = localAddStratifiedFit(T, summaryRows, coefRows, ...
    "Bootstrap-t", "Cbt_minus_C0B_empirical_cumulants_yh_fixed", ...
    ["pool_kurtosis_over_R", "pool_skewness_sq_over_R", "fallback_rate"], ...
    "coverage_error_to_C0B");

S = T(T.method == "Student-t" & isfinite(T.student_t_edgeworth_residual), :);
if isempty(S)
    summaryRows{end+1,1} = localFitSummaryRow("Student-t", ...
        "explicit_t_edgeworth_prediction", 0, 0, NaN, NaN, NaN, ...
        "No finite Student-t prediction rows."); %#ok<AGROW>
else
    residual = S.student_t_edgeworth_residual;
    summaryRows{end+1,1} = localFitSummaryRow("Student-t", ...
        "explicit_t_edgeworth_prediction", height(S), NaN, ...
        sqrt(mean(residual.^2, "omitnan")), mean(abs(residual), "omitnan"), ...
        NaN, "Uses paper Ct second-order term with pool skewness/kurtosis."); %#ok<AGROW>
end

fitSummary = vertcat(summaryRows{:});
if isempty(coefRows)
    fitCoefficients = table();
else
    fitCoefficients = vertcat(coefRows{:});
end
end

function [summaryRows, coefRows] = localAddFit(T, summaryRows, coefRows, methodName, ...
    modelName, allTerms, predictorTerms, responseName)
S = T(T.method == methodName, :);
if isempty(S)
    summaryRows{end+1,1} = localFitSummaryRow(methodName, modelName, 0, 0, ...
        NaN, NaN, NaN, "No rows for method."); %#ok<AGROW>
    return;
end
y = S.(char(responseName));
X = ones(height(S), numel(allTerms));
for j = 1:numel(predictorTerms)
    X(:, j + 1) = S.(char(predictorTerms(j)));
end
finiteRows = isfinite(y) & all(isfinite(X), 2);
y = y(finiteRows);
X = X(finiteRows, :);
if numel(y) < 2
    summaryRows{end+1,1} = localFitSummaryRow(methodName, modelName, numel(y), ...
        rank(X), NaN, NaN, NaN, "Too few finite rows."); %#ok<AGROW>
    return;
end
beta = pinv(X) * y;
yhat = X * beta;
resid = y - yhat;
ssTot = sum((y - mean(y)).^2);
if ssTot <= 0
    r2 = NaN;
else
    r2 = 1 - sum(resid.^2) / ssTot;
end
summaryRows{end+1,1} = localFitSummaryRow(methodName, modelName, numel(y), ...
    rank(X), sqrt(mean(resid.^2)), mean(abs(resid)), r2, ...
    "Empirical linear diagnostic; coefficients are not theoretical proof."); %#ok<AGROW>
for j = 1:numel(allTerms)
    coefRows{end+1,1} = table(string(methodName), string(modelName), ...
        string(allTerms(j)), beta(j), numel(y), rank(X), ...
        'VariableNames', {'method','model','term','coefficient','row_count','design_rank'}); %#ok<AGROW>
end
end

function S = localReplicateStabilityTable(T)
if isempty(T) || ~ismember("replicate_index", string(T.Properties.VariableNames))
    S = table();
    return;
end
[G, scheme, method, y0, hIndex, h, R] = findgroups(T.scheme, T.method, ...
    T.y0, T.h_index, T.h, T.R);
replicateCount = splitapply(@(x) numel(unique(x)), T.replicate_index, G);
meanCoverage = splitapply(@mean, T.coverage, G);
sdCoverage = splitapply(@std, T.coverage, G);
meanError = splitapply(@mean, T.coverage_error_to_C0B, G);
sdError = splitapply(@std, T.coverage_error_to_C0B, G);
maxAbsError = splitapply(@(x) max(abs(x)), T.coverage_error_to_C0B, G);
meanMcse = splitapply(@mean, T.coverage_mc_se, G);
meanFallback = splitapply(@mean, T.fallback_rate, G);
maxInf = splitapply(@max, T.interval_inf_rate, G);
S = table(scheme, method, y0, hIndex, h, R, replicateCount, meanCoverage, ...
    sdCoverage, meanError, sdError, maxAbsError, meanMcse, meanFallback, maxInf, ...
    'VariableNames', {'scheme','method','y0','h_index','h','R','replicate_count', ...
    'mean_coverage','sd_coverage','mean_error_to_C0B','sd_error_to_C0B', ...
    'max_abs_error_to_C0B','mean_coverage_mc_se','mean_fallback_rate', ...
    'max_interval_inf_rate'});
end

function P = localPredictionTables(T)
rows = cell(0, 1);
rows = localAddPredictionRows(T, rows, "Percentile bootstrap", ...
    "Cp_empirical_cumulants_yh_fixed_prediction", ...
    ["pool_skewness_over_sqrtR", "pool_kurtosis_over_R", ...
    "pool_skewness_sq_over_R"], "coverage_error_to_C0B");
rows = localAddPredictionRows(T, rows, "Bootstrap-t", ...
    "Cbt_empirical_cumulants_yh_fixed_prediction", ...
    ["pool_kurtosis_over_R", "pool_skewness_sq_over_R", "fallback_rate"], ...
    "coverage_error_to_C0B");
if isempty(rows)
    P = table();
else
    P = vertcat(rows{:});
end
end

function rows = localAddPredictionRows(T, rows, methodName, modelName, predictorTerms, responseName)
S = T(T.method == methodName, :);
if isempty(S)
    return;
end
[X, termNames, finiteRows] = localStratifiedDesign(S, predictorTerms);
yAll = S.(char(responseName));
y = yAll(finiteRows);
Xfit = X(finiteRows, :);
if numel(y) < 2
    return;
end
beta = pinv(Xfit) * y;
pred = nan(height(S), 1);
pred(finiteRows) = Xfit * beta;
resid = yAll - pred;
S.model = repmat(string(modelName), height(S), 1);
S.predicted_error_to_C0B = pred;
S.prediction_residual = resid;
S.prediction_abs_residual = abs(resid);
S.model_design_rank = repmat(rank(Xfit), height(S), 1);
S.model_term_count = repmat(numel(termNames), height(S), 1);
keep = ["model","scheme","scheme_kind","replicate_index","method","y0","h_index","h","R", ...
    "coverage","C0B","coverage_error_to_C0B","predicted_error_to_C0B", ...
    "prediction_residual","prediction_abs_residual","coverage_mc_se", ...
    "mean_rho3_w","mean_rho4_w","pool_skewness","pool_excess_kurtosis", ...
    "fallback_rate","interval_inf_rate","model_design_rank","model_term_count"];
rows{end+1,1} = S(:, keep); %#ok<AGROW>
end

function [X, termNames, finiteRows] = localStratifiedDesign(S, predictorTerms)
strata = "h" + string(S.h_index) + "_y" + string(S.y0);
uniqueStrata = unique(strata, "stable");
X = zeros(height(S), numel(uniqueStrata) + numel(predictorTerms));
termNames = strings(1, size(X, 2));
for i = 1:numel(uniqueStrata)
    X(:, i) = double(strata == uniqueStrata(i));
    termNames(i) = "stratum_" + uniqueStrata(i);
end
for j = 1:numel(predictorTerms)
    X(:, numel(uniqueStrata) + j) = S.(char(predictorTerms(j)));
    termNames(numel(uniqueStrata) + j) = predictorTerms(j);
end
finiteRows = all(isfinite(X), 2);
end

function [summaryRows, coefRows] = localAddStratifiedFit(T, summaryRows, coefRows, ...
    methodName, modelName, predictorTerms, responseName)
S = T(T.method == methodName, :);
if isempty(S)
    summaryRows{end+1,1} = localFitSummaryRow(methodName, modelName, 0, 0, ...
        NaN, NaN, NaN, "No rows for method."); %#ok<AGROW>
    return;
end
y = S.(char(responseName));
[X, termNames, finiteRows] = localStratifiedDesign(S, predictorTerms);
finiteRows = finiteRows & isfinite(y);
y = y(finiteRows);
X = X(finiteRows, :);
if numel(y) < 2
    summaryRows{end+1,1} = localFitSummaryRow(methodName, modelName, numel(y), ...
        rank(X), NaN, NaN, NaN, "Too few finite rows."); %#ok<AGROW>
    return;
end
beta = pinv(X) * y;
yhat = X * beta;
resid = y - yhat;
ssTot = sum((y - mean(y)).^2);
if ssTot <= 0
    r2 = NaN;
else
    r2 = 1 - sum(resid.^2) / ssTot;
end
summaryRows{end+1,1} = localFitSummaryRow(methodName, modelName, numel(y), ...
    rank(X), sqrt(mean(resid.^2)), mean(abs(resid)), r2, ...
    "Stratified diagnostic with y/h intercepts; coefficients are not theoretical proof."); %#ok<AGROW>
for j = 1:numel(termNames)
    coefRows{end+1,1} = table(string(methodName), string(modelName), ...
        string(termNames(j)), beta(j), numel(y), rank(X), ...
        'VariableNames', {'method','model','term','coefficient','row_count','design_rank'}); %#ok<AGROW>
end
end

function row = localFitSummaryRow(methodName, modelName, rowCount, designRank, rmse, mae, r2, note)
row = table(string(methodName), string(modelName), rowCount, designRank, ...
    rmse, mae, r2, string(note), ...
    'VariableNames', {'method','model','row_count','design_rank','rmse','mae','r_squared','note'});
end

function localPlotValidationFigures(resultsRoot, schemeRows, fitPredictions)
figDir = fullfile(resultsRoot, "figures");
dpimnumeric.ensureDir(figDir);
if ~isempty(schemeRows)
    plotRows = localAggregateSchemeRows(schemeRows);
    fig = figure("Visible", "off");
    tiledlayout(fig, 1, 2);
    nexttile;
    bar(categorical(plotRows.scheme), [plotRows.mean_rho3_w, plotRows.mean_rho4_w]);
    ylabel("rho value");
    legend(["rho3", "rho4"], "Location", "best");
    title("Weight high-order factors", "Interpreter", "none");
    grid on;
    nexttile;
    bar(categorical(plotRows.scheme), [plotRows.mean_n2_eff_w, ...
        plotRows.mean_n3_eff_w, plotRows.mean_n4_eff_w]);
    ylabel("effective order");
    legend(["n2_eff", "n3_eff", "n4_eff"], "Location", "best");
    title("Weight effective orders", "Interpreter", "none");
    grid on;
    exportgraphics(fig, fullfile(figDir, "weight_moment_ranges.png"), "Resolution", 180);
    close(fig);
end

if ~isempty(fitPredictions)
    fig = figure("Visible", "off");
    methods = unique(fitPredictions.method, "stable");
    tiledlayout(fig, 1, numel(methods));
    for i = 1:numel(methods)
        ax = nexttile;
        S = fitPredictions(fitPredictions.method == methods(i), :);
        scatter(ax, S.coverage_error_to_C0B, S.predicted_error_to_C0B, 28, "filled");
        hold(ax, "on");
        lim = max(abs([S.coverage_error_to_C0B; S.predicted_error_to_C0B]), [], "omitnan");
        if isempty(lim) || ~isfinite(lim) || lim <= 0
            lim = 0.1;
        end
        plot(ax, [-lim, lim], [-lim, lim], "k--");
        xlim(ax, [-lim, lim]);
        ylim(ax, [-lim, lim]);
        xlabel(ax, "observed C-C0B");
        ylabel(ax, "predicted C-C0B");
        title(ax, methods(i), "Interpreter", "none");
        grid(ax, "on");
    end
    exportgraphics(fig, fullfile(figDir, "observed_vs_predicted_coverage_error.png"), "Resolution", 180);
    close(fig);
end
end

function A = localAggregateSchemeRows(S)
schemes = unique(S.scheme, "stable");
rows = cell(numel(schemes), 1);
for i = 1:numel(schemes)
    q = S(S.scheme == schemes(i), :);
    rows{i} = table(q.scheme(1), q.scheme_kind(1), height(q), ...
        mean(q.mean_rho3_w), mean(q.mean_rho4_w), mean(q.mean_n2_eff_w), ...
        mean(q.mean_n3_eff_w), mean(q.mean_n4_eff_w), ...
        'VariableNames', {'scheme','scheme_kind','replicate_count', ...
        'mean_rho3_w','mean_rho4_w','mean_n2_eff_w','mean_n3_eff_w','mean_n4_eff_w'});
end
A = vertcat(rows{:});
end

function localWriteReport(path, vcfg, C0B, kMinus, kPlus, schemeRows, coverageRows, fitSummary)
lines = strings(0, 1);
lines(end+1) = "# Weighted R-B-w Coverage Validation";
lines(end+1) = "";
lines(end+1) = "This run validates the probability-weighted coverage-expression diagnostics for the paper theorem `weighted-R-B-w-coverage`.";
lines(end+1) = "It uses explicit finite-B order statistic endpoints instead of MATLAB `quantile` interpolation.";
lines(end+1) = "";
lines(end+1) = sprintf("- run_mode: `%s`", string(vcfg.run_mode));
lines(end+1) = sprintf("- n: `%d`, curve_pool_size: `%d`, M: `%d`, B: `%d`", ...
    vcfg.n, vcfg.curve_pool_size, vcfg.M, vcfg.B);
lines(end+1) = sprintf("- replicate_count: `%d`; power_list: `%s`", ...
    vcfg.replicate_count, mat2str(vcfg.power_list));
lines(end+1) = sprintf("- endpoint ranks: `k_minus=%d`, `k_plus=%d`, `C0B=%.12g`", ...
    kMinus, kPlus, C0B);
lines(end+1) = sprintf("- R list: `%s`; y points: `%s`; h list: `%s`", ...
    mat2str(vcfg.R_list), mat2str(vcfg.y_points), mat2str(vcfg.h_list));
lines(end+1) = "";
lines(end+1) = "## Weight Moment Range";
plotRows = localAggregateSchemeRows(schemeRows);
for i = 1:height(plotRows)
    lines(end+1) = sprintf("- `%s`: mean rho3=%.6g, mean rho4=%.6g, mean n2_eff=%.6g, mean n3_eff=%.6g, mean n4_eff=%.6g", ...
        plotRows.scheme(i), plotRows.mean_rho3_w(i), plotRows.mean_rho4_w(i), ...
        plotRows.mean_n2_eff_w(i), plotRows.mean_n3_eff_w(i), plotRows.mean_n4_eff_w(i));
end
lines(end+1) = "";
lines(end+1) = "## Coverage Diagnostics";
for i = 1:height(fitSummary)
    lines(end+1) = sprintf("- `%s` / `%s`: rows=%d, rank=%g, RMSE=%.6g, MAE=%.6g, R2=%.6g", ...
        fitSummary.method(i), fitSummary.model(i), fitSummary.row_count(i), ...
        fitSummary.design_rank(i), fitSummary.rmse(i), fitSummary.mae(i), fitSummary.r_squared(i));
end
bt = coverageRows(coverageRows.method == "Bootstrap-t", :);
if ~isempty(bt)
    lines(end+1) = "";
    lines(end+1) = sprintf("- Bootstrap-t max fallback rate: %.6g", max(bt.fallback_rate));
    lines(end+1) = sprintf("- Bootstrap-t max interval Inf rate: %.6g", max(bt.interval_inf_rate));
end
lines(end+1) = "";
lines(end+1) = "## Interpretation Boundary";
lines(end+1) = "This is a numerical diagnostic for the formula structure: C0B baseline, rho3/rho4 weight terms, Student-t explicit second-order term, and A_lambda instability frequency.";
lines(end+1) = "It is not a proof of the Edgeworth expansion and does not estimate the unknown mathfrak A coefficients from first principles.";
dpimnumeric.writeText(path, strjoin(lines, newline));
end
