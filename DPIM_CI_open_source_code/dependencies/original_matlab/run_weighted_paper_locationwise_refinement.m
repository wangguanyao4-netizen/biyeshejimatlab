function outputs = run_weighted_paper_locationwise_refinement(runMode, resultsRoot, candidateSourceRoot)
%RUN_WEIGHTED_PAPER_LOCATIONWISE_REFINEMENT Location-wise h refinement.
%
% This is a follow-up experiment for E2/E4 after the global-h full suite.
% It keeps probability weights, GPU Voronoi assignment, and independent
% tuning/confirmation splits, but selects h for each point_method/R/y0.
% If candidateSourceRoot is provided, the script runs locked validation with
% fixed h values selected from that source result directory.

if nargin < 1 || strlength(string(runMode)) == 0
    runMode = "full";
end

projectRoot = fileparts(mfilename("fullpath"));
if nargin < 2 || strlength(string(resultsRoot)) == 0
    stamp = string(datetime("now", "Format", "yyyyMMdd_HHmmss"));
    if nargin >= 3 && strlength(string(candidateSourceRoot)) > 0
        prefix = "weighted_paper_locationwise_locked_validation_";
    else
        prefix = "weighted_paper_locationwise_refinement_";
    end
    resultsRoot = fullfile(projectRoot, "results", prefix + lower(string(runMode)) + "_" + stamp);
end
if nargin < 3
    candidateSourceRoot = "";
end

cfg = weighted_paper_config(runMode, projectRoot, resultsRoot);
cfg.refinement_scope = "E2_E4_locationwise_h";
cfg.refinement_max_candidates_per_location = 2;
cfg.refinement_note = "Exploratory location-wise h refinement; locked validation should be rerun before main-text claims.";
cfg.candidate_source_root = char(string(candidateSourceRoot));
cfg.locked_validation = strlength(string(candidateSourceRoot)) > 0;
if cfg.locked_validation
    cfg.methods = cfg.main_methods;
end

dpimnumeric.ensureDir(cfg.results_root);
localAddPaths(cfg);
preflight = localPreflight(cfg);
writetable(preflight, fullfile(cfg.results_root, "preflight.csv"));
dpimnumeric.writeJson(fullfile(cfg.results_root, "weighted_paper_locationwise_config.json"), cfg);
if any(preflight.status == "FAIL")
    error("Preflight failed. Inspect %s.", fullfile(cfg.results_root, "preflight.csv"));
end

if cfg.locked_validation
    fprintf("Weighted paper location-wise locked validation: %s\n", cfg.run_mode);
    fprintf("Candidate source: %s\n", cfg.candidate_source_root);
else
    fprintf("Weighted paper location-wise refinement: %s\n", cfg.run_mode);
end
fprintf("Output root: %s\n", cfg.results_root);

specs = localPointExperiments();
outputs = struct();
outputs.results_root = cfg.results_root;
auditRows = cell(0, 1);
for i = 1:numel(specs)
    expSpec = specs(i);
    fprintf("\n=== %s ===\n", expSpec.name);
    result = localRunExperiment(cfg, expSpec);
    outputs.(matlab.lang.makeValidName(expSpec.name)) = result;
    auditRows{end+1,1} = localAuditExperiment(result.output_dir, expSpec.name); %#ok<AGROW>
end

audit = vertcat(auditRows{:});
writetable(audit, fullfile(cfg.results_root, "locationwise_audit.csv"));
summary = groupsummary(audit, "claim_status");
writetable(summary, fullfile(cfg.results_root, "locationwise_status_summary.csv"));
save(fullfile(cfg.results_root, "locationwise_outputs.mat"), "outputs", "cfg", "preflight", "audit", "-v7.3");

fprintf("\nLocation-wise refinement completed: %s\n", cfg.results_root);
end

function localAddPaths(cfg)
addpath(cfg.project_root, "-begin");
addpath(fullfile(cfg.weighted_root, "external_weight_providers"), "-begin");
addpath(fullfile(cfg.weighted_root, "common"), "-begin");
end

function T = localPreflight(cfg)
rows = cell(0, 1);
rows{end+1,1} = localCheckRow("paper_tex", isfile(cfg.paper_tex), cfg.paper_tex);
rows{end+1,1} = localCheckRow("thesis_pdf", isfile(cfg.thesis_pdf), cfg.thesis_pdf);
rows{end+1,1} = localCheckRow("weighted_root", isfolder(cfg.weighted_root), cfg.weighted_root);
for i = 1:numel(cfg.required_functions)
    name = cfg.required_functions(i);
    rows{end+1,1} = localCheckRow("function:" + name, exist(name, "file") == 2, which(name)); %#ok<AGROW>
end
try
    dev = gpuDevice;
    gpuMsg = sprintf("%s, available %.2f GiB", dev.Name, dev.AvailableMemory / 2^30);
    gpuOk = true;
catch ME
    gpuMsg = ME.message;
    gpuOk = false;
end
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

function specs = localPointExperiments()
specs = [ ...
    struct("name", "E2_scalar_normal_R_order_weighted", ...
    "d", 1, "y_points", [-1, 0, 1], "anchor_h", 0.233572146909012), ...
    struct("name", "E4_standard_normal_nonlinear_tail_weighted", ...
    "d", 5, "y_points", [1.6, 2.0, 2.4], "anchor_h", 4e-4) ...
    ];
end

function output = localRunExperiment(cfg, expSpec)
outDir = fullfile(cfg.results_root, expSpec.name);
dpimnumeric.ensureDir(outDir);
problem = localProblem(expSpec);
expCfg = cfg;
expCfg.results_root = outDir;
expCfg.weighting_cfg.voronoi_output_dir = fullfile(outDir, "ci_probability_weights");
dpimnumeric.writeJson(fullfile(outDir, "config.json"), localJsonSafe(expCfg, expSpec));

tuning = table();
confirmation = table();
candidateRows = table();
weightRows = table();

for im = 1:numel(cfg.methods)
    method = cfg.methods(im);
    fprintf("Building weighted curve pool: %s, method=%s, curves=%d, n=%d\n", ...
        expSpec.name, method, cfg.curve_pool_size, cfg.n);
    curves = dpim_build_weighted_curve_pool(problem, method, cfg.n, expCfg, cfg.curve_pool_size);
    weightRows = [weightRows; localWeightTable(curves, expSpec.name, method)]; %#ok<AGROW>
    yPoints = double(expSpec.y_points(:)).';
    tuningIdx = 1:cfg.tuning_pool_size;
    confirmationIdx = (cfg.tuning_pool_size + 1):cfg.curve_pool_size;

    if cfg.locked_validation
        candidates = localLoadLockedCandidates(cfg, expSpec, method);
        hEval = unique(candidates.h_index);
    else
        candidates = table();
        hEval = 1:numel(cfg.h_list);
    end
    poolEst = localEvaluatePool(curves, yPoints, cfg.h_list, hEval);
    truth = localTruthGrid(problem, yPoints, cfg.h_list, expCfg, hEval);

    if cfg.locked_validation
        % Locked validation reuses fixed h values from an earlier exploratory run.
    else
        for iy = 1:numel(yPoints)
            for ih = 1:numel(cfg.h_list)
                for ir = 1:numel(cfg.R_list)
                    R = cfg.R_list(ir);
                    seed = localBlockSeed(cfg.seed, 11, im, iy, ih, ir);
                    rows = localCoverageBlock(squeeze(poolEst(tuningIdx, iy, ih)), ...
                        truth(iy, ih), expCfg, expSpec.name, method, yPoints(iy), ...
                        cfg.h_list(ih), ih, R, "tuning", seed);
                    tuning = [tuning; rows]; %#ok<AGROW>
                end
            end
        end
        candidates = localSelectLocationCandidates(tuning(tuning.point_method == method, :), cfg, expSpec, method);
    end
    candidateRows = [candidateRows; candidates]; %#ok<AGROW>
    for ic = 1:height(candidates)
        iy = find(yPoints == candidates.y0(ic), 1, "first");
        ih = candidates.h_index(ic);
        R = candidates.R(ic);
        ir = find(cfg.R_list == R, 1, "first");
        seed = localBlockSeed(cfg.seed, 12 + ic, im, iy, ih, ir);
        rows = localCoverageBlock(squeeze(poolEst(confirmationIdx, iy, ih)), ...
            truth(iy, ih), expCfg, expSpec.name, method, yPoints(iy), ...
            cfg.h_list(ih), ih, R, "confirmation", seed);
        rows.selection_reason = repmat(candidates.selection_reason(ic), height(rows), 1);
        confirmation = [confirmation; rows]; %#ok<AGROW>
    end
end

writetable(tuning, fullfile(outDir, "locationwise_tuning_summary.csv"));
writetable(candidateRows, fullfile(outDir, "locationwise_candidate_h.csv"));
writetable(confirmation, fullfile(outDir, "locationwise_confirmation_summary.csv"));
writetable(weightRows, fullfile(outDir, "weight_diagnostics.csv"));
save(fullfile(outDir, "locationwise_raw_results.mat"), "cfg", "expSpec", ...
    "tuning", "confirmation", "candidateRows", "weightRows", "-v7.3");
output = struct("output_dir", outDir, ...
    "candidate_csv", fullfile(outDir, "locationwise_candidate_h.csv"), ...
    "confirmation_csv", fullfile(outDir, "locationwise_confirmation_summary.csv"));
end

function S = localJsonSafe(expCfg, expSpec)
S = expCfg;
S.experiment = expSpec.name;
S.y_points = expSpec.y_points;
S.anchor_h = expSpec.anchor_h;
end

function problem = localProblem(expSpec)
problem = struct();
problem.name = char(expSpec.name);
problem.short_name = char(extractBefore(string(expSpec.name) + "_", "_"));
problem.d = expSpec.d;
problem.target_distribution = "standard_normal";
problem.center_transform = "normal_icdf";
problem.provider = "voronoi_ci_probability_weights_provider";
if string(expSpec.name) == "E2_scalar_normal_R_order_weighted"
    problem.response_fun = @(theta) theta(:,1);
else
    problem.response_fun = @(theta) localNormalNonlinearResponse(theta);
end
end

function y = localNormalNonlinearResponse(theta)
y = 2.0 + 0.12 * sum(theta, 2) ...
    + 0.25 * tanh(theta(:,1) .* theta(:,2)) ...
    + 0.08 * (theta(:,3).^2 - 1) ...
    + 0.06 * sin(theta(:,4) .* theta(:,5));
end

function poolEst = localEvaluatePool(curves, yPoints, hList, hEval)
if nargin < 4 || isempty(hEval)
    hEval = 1:numel(hList);
end
poolEst = nan(numel(curves), numel(yPoints), numel(hList));
for iy = 1:numel(yPoints)
    for ih = hEval(:).'
        poolEst(:, iy, ih) = dpim_curve_point_estimates(curves, yPoints(iy), hList(ih));
    end
end
end

function truth = localTruthGrid(problem, yPoints, hList, cfg, hEval)
if nargin < 5 || isempty(hEval)
    hEval = 1:numel(hList);
end
truth = nan(numel(yPoints), numel(hList));
for iy = 1:numel(yPoints)
    for ih = hEval(:).'
        truth(iy, ih) = localTruth(problem, yPoints(iy), hList(ih), cfg, 9000 * iy + ih);
    end
end
end

function truth = localTruth(problem, y0, h, cfg, seedOffset)
if string(problem.name) == "E2_scalar_normal_R_order_weighted"
    truth = dpim_gaussian_kernel(y0, sqrt(1 + h^2));
else
    truth = dpim_truth_smoothed_density(problem, y0, h, cfg, seedOffset);
end
end

function rows = localCoverageBlock(poolVals, truth, cfg, experiment, pointMethod, ...
    y0, h, hIndex, R, phase, blockSeed)
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
poolSize = numel(poolVals);

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
end

rowCells = cell(nMethods, 1);
for j = 1:nMethods
    coverage = mean(cover(:,j));
    finiteLength = isfinite(lengths(:,j));
    [formulaBaseline, formulaBaselineType] = localFormulaBaseline(methodNames(j), cfg.alpha, cfg.B);
    rowCells{j} = table(string(phase), string(experiment), string(pointMethod), ...
        methodNames(j), y0, hIndex, h, R, poolSize, cfg.M, cfg.B, truth, ...
        mean(poolVals), mean(poolVals) - truth, coverage, cfg.nominal, ...
        formulaBaseline, string(formulaBaselineType), ...
        coverage - formulaBaseline, abs(coverage - formulaBaseline), ...
        sqrt(max(coverage * (1 - coverage), realmin) / cfg.M), ...
        mean(leftMiss(:,j)), mean(rightMiss(:,j)), sum(fallback), ...
        mean(fallback), sum(btInf), mean(btInf), sum(isInf(:,j)), ...
        mean(isInf(:,j)), mean(lengths(finiteLength,j), "omitnan"), ...
        median(lengths(finiteLength,j), "omitnan"), sum(finiteLength), ...
        sum(finiteLength) / cfg.M, blockSeed, ...
        'VariableNames', {'phase','experiment','point_method','method','y0', ...
        'h_index','h','R','pool_size','M','B','truth','pool_mean', ...
        'estimator_bias','coverage','nominal_coverage','formula_baseline', ...
        'formula_baseline_type','coverage_error', ...
        'abs_coverage_error','coverage_mc_se','left_miss','right_miss', ...
        'fallback_count','fallback_rate','bootstrap_t_inf_count', ...
        'bootstrap_t_inf_rate','interval_inf_count','interval_inf_rate', ...
        'mean_interval_length','median_interval_length','finite_length_count', ...
        'finite_length_rate','block_seed'});
end
rows = vertcat(rowCells{:});
end

function [ci, diag] = localCiMethods(x, truth, alpha, B, lambda)
x = x(:);
R = numel(x);
mu = mean(x);
s = std(x, 0);
se = s / sqrt(R);
tCrit = tinv(1 - alpha / 2, max(R - 1, 1));
ci(1) = localPackCi("Student-t", mu - tCrit * se, mu + tCrit * se, truth);

bootIdx = randi(R, B, R);
boot = x(bootIdx);
bm = mean(boot, 2);
[kMinus, kPlus] = localFiniteBRanks(alpha, B);
bmSorted = sort(bm);
ci(2) = localPackCi("Percentile bootstrap", bmSorted(kMinus), bmSorted(kPlus), truth);

bs = std(boot, 0, 2);
T = sqrt(R) * (bm - mu) ./ bs;
T(bs <= 0 & bm == mu) = 0;
T(bs <= 0 & bm > mu) = Inf;
T(bs <= 0 & bm < mu) = -Inf;
T = sort(T);
btLower = mu - T(kPlus) * se;
btUpper = mu - T(kMinus) * se;
ci(3) = localPackCi("Bootstrap-t", btLower, btUpper, truth);

ratio = ci(3).length / max(ci(2).length, realmin);
trigger = ci(3).infinite || (~ci(3).infinite && ratio > lambda);
if trigger
    ci(4) = localPackCi("Bootstrap-t fallback-rule", ci(2).lower, ci(2).upper, truth);
else
    ci(4) = localPackCi("Bootstrap-t fallback-rule", ci(3).lower, ci(3).upper, truth);
end
diag = struct("fallback_trigger", trigger, "bootstrap_t_infinite", ci(3).infinite);
end

function ci = localPackCi(name, lower, upper, truth)
ci = struct("name", string(name), "lower", lower, "upper", upper, ...
    "length", upper - lower, "infinite", ~(isfinite(lower) && isfinite(upper)), ...
    "contains", false, "left_miss", false, "right_miss", false);
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

function candidates = localSelectLocationCandidates(tuning, cfg, expSpec, pointMethod)
raw = tuning(tuning.method == "Bootstrap-t", :);
fallbackRows = tuning(tuning.method == "Bootstrap-t fallback-rule", :);
rows = cell(0, 1);
for ir = 1:numel(cfg.R_list)
    R = cfg.R_list(ir);
    for iy = 1:numel(expSpec.y_points)
        y0 = expSpec.y_points(iy);
        scores = nan(numel(cfg.h_list), 1);
        stable = false(numel(cfg.h_list), 1);
        for ih = 1:numel(cfg.h_list)
            q = raw(raw.R == R & raw.y0 == y0 & raw.h_index == ih, :);
            f = fallbackRows(fallbackRows.R == R & fallbackRows.y0 == y0 & fallbackRows.h_index == ih, :);
            if isempty(q) || isempty(f)
                continue;
            end
            scores(ih) = q.abs_coverage_error + 0.25 * max(0, 0.93 - q.coverage) ...
                + 0.25 * f.fallback_rate + 0.25 * q.bootstrap_t_inf_rate;
            stable(ih) = q.bootstrap_t_inf_rate <= cfg.selection_max_inf ...
                && f.fallback_rate <= cfg.selection_max_fallback ...
                && q.coverage >= 0.93 ...
                && q.abs_coverage_error <= cfg.selection_max_mean_abs_error;
        end
        chosen = [];
        reasons = strings(0, 1);
        eligible = find(stable);
        if ~isempty(eligible)
            [~, loc] = min(scores(eligible));
            chosen(end+1) = eligible(loc); %#ok<AGROW>
            reasons(end+1,1) = "location_tuning_selected_stable"; %#ok<AGROW>
            [~, loc] = max(cfg.h_list(eligible));
            chosen(end+1) = eligible(loc); %#ok<AGROW>
            reasons(end+1,1) = "location_largest_stable_h"; %#ok<AGROW>
        else
            [~, best] = min(scores);
            if isfinite(scores(best))
                chosen(end+1) = best; %#ok<AGROW>
                reasons(end+1,1) = "location_min_penalty_no_stable"; %#ok<AGROW>
            end
        end
        [chosen, ia] = unique(chosen, "stable");
        reasons = reasons(ia);
        keep = 1:min(numel(chosen), cfg.refinement_max_candidates_per_location);
        for k = keep
            ih = chosen(k);
            q = raw(raw.R == R & raw.y0 == y0 & raw.h_index == ih, :);
            f = fallbackRows(fallbackRows.R == R & fallbackRows.y0 == y0 & fallbackRows.h_index == ih, :);
            rows{end+1,1} = table(string(expSpec.name), string(pointMethod), ...
                y0, R, ih, cfg.h_list(ih), reasons(k), stable(ih), q.coverage, ...
                q.abs_coverage_error, f.fallback_rate, q.bootstrap_t_inf_rate, ...
                'VariableNames', {'experiment','point_method','y0','R','h_index','h', ...
                'selection_reason','tuning_stable','tuning_coverage', ...
                'tuning_abs_error','tuning_fallback_rate','tuning_inf_rate'});
        end
    end
end
candidates = vertcat(rows{:});
end

function candidates = localLoadLockedCandidates(cfg, expSpec, pointMethod)
sourcePath = fullfile(cfg.candidate_source_root, expSpec.name, "locationwise_confirmation_summary.csv");
if ~isfile(sourcePath)
    error("Locked candidate source is missing: %s", sourcePath);
end
T = readtable(sourcePath, TextType="string", Delimiter=",");
T = T(T.method == "Bootstrap-t" & T.point_method == string(pointMethod), :);
if isempty(T)
    error("No Bootstrap-t source rows for %s/%s.", expSpec.name, string(pointMethod));
end
rows = cell(0, 1);
for ir = 1:numel(cfg.R_list)
    R = cfg.R_list(ir);
    for iy = 1:numel(expSpec.y_points)
        y0 = expSpec.y_points(iy);
        q = T(T.R == R & T.y0 == y0, :);
        if isempty(q)
            continue;
        end
        penalty = q.abs_coverage_error ...
            + 0.25 * max(0, 0.90 - q.coverage) ...
            + 0.25 * q.fallback_rate ...
            + 0.25 * q.bootstrap_t_inf_rate;
        [~, loc] = min(penalty);
        best = q(loc, :);
        ih = localNearestIndex(cfg.h_list, best.h);
        rows{end+1,1} = table(string(expSpec.name), string(pointMethod), ...
            y0, R, ih, cfg.h_list(ih), "locked_best_confirmed_locationwise", true, ...
            best.coverage, best.abs_coverage_error, best.fallback_rate, ...
            best.bootstrap_t_inf_rate, ...
            'VariableNames', {'experiment','point_method','y0','R','h_index','h', ...
            'selection_reason','tuning_stable','tuning_coverage', ...
            'tuning_abs_error','tuning_fallback_rate','tuning_inf_rate'}); %#ok<AGROW>
    end
end
if isempty(rows)
    candidates = table();
else
    candidates = vertcat(rows{:});
end
end

function audit = localAuditExperiment(outDir, expName)
path = fullfile(outDir, "locationwise_confirmation_summary.csv");
if ~isfile(path)
    audit = table(string(expName), "missing", NaN, NaN, NaN, NaN, NaN, ...
        "locationwise_confirmation_summary.csv missing.", ...
        'VariableNames', {'experiment','claim_status','mean_coverage','min_coverage', ...
        'mean_abs_coverage_error','max_fallback_rate','max_bootstrap_t_inf_rate','note'});
    return;
end
T = readtable(path, TextType="string", Delimiter=",");
if ismember("point_method", string(T.Properties.VariableNames))
    T = T(T.point_method == "sobol_scrambled", :);
end
mainReasons = ["location_tuning_selected_stable", "locked_best_confirmed_locationwise"];
bt = T(T.method == "Bootstrap-t" & ismember(T.selection_reason, mainReasons), :);
fb = T(T.method == "Bootstrap-t fallback-rule" & ismember(T.selection_reason, mainReasons), :);
if isempty(bt)
    status = "supplement_only";
    note = "No location_tuning_selected_stable confirmation rows.";
    audit = table(string(expName), status, NaN, NaN, NaN, NaN, NaN, note, ...
        'VariableNames', {'experiment','claim_status','mean_coverage','min_coverage', ...
        'mean_abs_coverage_error','max_fallback_rate','max_bootstrap_t_inf_rate','note'});
    return;
end
meanCov = mean(bt.coverage);
minCov = min(bt.coverage);
meanAbs = mean(bt.abs_coverage_error);
maxInf = max(bt.bootstrap_t_inf_rate);
maxFallback = max(fb.fallback_rate);
stable = minCov >= 0.90 && meanAbs <= 0.05 && maxInf <= 0.01 && maxFallback <= 0.05;
isLocked = any(T.selection_reason == "locked_best_confirmed_locationwise");
if stable
    if isLocked
        status = "main_text_ready";
        note = "Locked location-wise validation satisfies the engineering coverage gate.";
    else
        status = "locationwise_validation_candidate";
        note = "Location-wise h passed exploratory confirmation; rerun locked validation before main-text use.";
    end
else
    status = "supplement_only";
    if isLocked
        note = "Locked location-wise validation still fails at least one gate.";
    else
        note = "Location-wise h still fails at least one exploratory gate.";
    end
end
audit = table(string(expName), status, meanCov, minCov, meanAbs, maxFallback, maxInf, string(note), ...
    'VariableNames', {'experiment','claim_status','mean_coverage','min_coverage', ...
    'mean_abs_coverage_error','max_fallback_rate','max_bootstrap_t_inf_rate','note'});
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

function seed = localBlockSeed(baseSeed, phaseCode, im, iy, ih, ir)
seed = baseSeed + 10000000 * phaseCode + 1000000 * im ...
    + 10000 * iy + 100 * ih + ir;
end

function idx = localNearestIndex(values, target)
[~, idx] = min(abs(values - target));
end
