function outputs = run_weighted_paper_formal_campaign(runMode, resultsRoot)
%RUN_WEIGHTED_PAPER_FORMAL_CAMPAIGN Orchestrate paper-aligned E1--E8 evidence.
%
% Usage:
%   run_weighted_paper_formal_campaign("preflight")
%   run_weighted_paper_formal_campaign("diagnostic")
%   run_weighted_paper_formal_campaign("pilot")
%   run_weighted_paper_formal_campaign("medium")
%   run_weighted_paper_formal_campaign("formal")
%
% This wrapper does not change the paper TeX. It creates a locked campaign
% directory, runs the probability-weighted evidence suite plus the R-B-w
% formula validation, and writes paper-numbered audit artifacts.

if nargin < 1 || strlength(string(runMode)) == 0
    runMode = "diagnostic";
end

runMode = lower(string(runMode));
validModes = ["preflight", "diagnostic", "small", "pilot", "medium", "full", "formal"];
if ~any(runMode == validModes)
    error("runMode must be one of: %s.", strjoin(validModes, ", "));
end

projectRoot = fileparts(mfilename("fullpath"));
if nargin < 2 || strlength(string(resultsRoot)) == 0
    stamp = string(datetime("now", "Format", "yyyyMMdd_HHmmss"));
    resultsRoot = fullfile(projectRoot, "results", ...
        "weighted_paper_formal_" + runMode + "_" + stamp);
end

cfgMode = localUnderlyingMode(runMode);
evidenceRoot = fullfile(resultsRoot, "evidence_suite");
rbwRoot = fullfile(resultsRoot, "coverage_formula_validation");
reportRoot = fullfile(resultsRoot, "_formal_report");
tablesRoot = fullfile(reportRoot, "tables");

cfg = weighted_paper_config(cfgMode, projectRoot, evidenceRoot);
localAddPaths(cfg);

formalCfg = localCampaignConfig(runMode, cfgMode, resultsRoot, evidenceRoot, ...
    rbwRoot, reportRoot, tablesRoot, cfg);

dpimnumeric.ensureDir(resultsRoot);
dpimnumeric.ensureDir(reportRoot);
dpimnumeric.ensureDir(tablesRoot);

dpimnumeric.writeJson(fullfile(resultsRoot, "config_lock.json"), formalCfg);
dpimnumeric.writeJson(fullfile(resultsRoot, "seed_lock.json"), localSeedLock(cfg, cfgMode));

preflight = localPreflight(cfg, formalCfg);
writetable(preflight, fullfile(resultsRoot, "preflight.csv"));

manifest = localManifest(cfg, formalCfg);
dpimnumeric.writeJson(fullfile(resultsRoot, "manifest.json"), manifest);

fprintf("Weighted paper formal campaign: %s\n", runMode);
fprintf("Campaign root: %s\n", resultsRoot);

if any(preflight.status == "FAIL")
    error("Formal campaign preflight failed. Inspect %s.", fullfile(resultsRoot, "preflight.csv"));
end

outputs = struct();
outputs.results_root = resultsRoot;
outputs.config_lock = fullfile(resultsRoot, "config_lock.json");
outputs.seed_lock = fullfile(resultsRoot, "seed_lock.json");
outputs.preflight_csv = fullfile(resultsRoot, "preflight.csv");
outputs.manifest_json = fullfile(resultsRoot, "manifest.json");

if runMode ~= "preflight"
    outputs.evidence_suite = run_weighted_paper_evidence_suite(cfgMode, evidenceRoot);
    outputs.coverage_formula_validation = run_weighted_RBw_coverage_validation(cfgMode, rbwRoot);
end

outputs.assembly = localAssembleFormalOutputs(formalCfg, cfg);
save(fullfile(resultsRoot, "formal_campaign_outputs.mat"), ...
    "outputs", "cfg", "formalCfg", "preflight", "-v7.3");

fprintf("Formal campaign completed: %s\n", resultsRoot);
end

function mode = localUnderlyingMode(runMode)
runMode = lower(string(runMode));
if runMode == "preflight"
    mode = "diagnostic";
else
    mode = runMode;
end
end

function localAddPaths(cfg)
addpath(cfg.project_root, "-begin");
addpath(fullfile(cfg.weighted_root, "external_weight_providers"), "-begin");
addpath(fullfile(cfg.weighted_root, "common"), "-begin");
end

function formalCfg = localCampaignConfig(runMode, cfgMode, resultsRoot, evidenceRoot, ...
    rbwRoot, reportRoot, tablesRoot, cfg)
formalCfg = struct();
formalCfg.schema_version = "weighted_paper_formal_campaign_v1";
formalCfg.created_at = char(datetime("now"));
formalCfg.run_mode = char(runMode);
formalCfg.underlying_mode = char(cfgMode);
formalCfg.project_root = cfg.project_root;
formalCfg.results_root = char(resultsRoot);
formalCfg.evidence_root = char(evidenceRoot);
formalCfg.coverage_formula_root = char(rbwRoot);
formalCfg.report_root = char(reportRoot);
formalCfg.tables_root = char(tablesRoot);
formalCfg.paper_tex = cfg.paper_tex;
formalCfg.paper_figure_dir = localPaperFigureDir(cfg);
formalCfg.core_policy = "Main evidence must be probability-weighted; equal weights are controls or fixed-weight special cases only.";
formalCfg.sample_pool_policy = "Point pools use randomized QMC by default; large MC is reserved for truth/reference integration.";
formalCfg.coverage_policy = "Bootstrap coverage is compared against finite-B C0B or formula baselines, not blindly against 0.95.";
formalCfg.h_protocol = "Tune h on tuning curves, lock selected h, then validate on independent confirmation curves.";
formalCfg.gpu_policy = "GPU may accelerate Voronoi/distance/bootstrap work only after CPU/GPU parity smoke test passes.";
formalCfg.stage_sequence = ["preflight", "diagnostic", "pilot", "medium", "formal"];
formalCfg.acceptance = struct();
formalCfg.acceptance.formal_min_M = 1000;
formalCfg.acceptance.formal_min_B = 999;
formalCfg.acceptance.positive_max_fallback_rate = cfg.selection_max_fallback;
formalCfg.acceptance.positive_max_inf_rate = cfg.selection_max_inf;
formalCfg.acceptance.positive_max_mean_abs_error = cfg.selection_max_mean_abs_error;
formalCfg.acceptance.positive_min_location_coverage = cfg.selection_min_location_coverage;
formalCfg.acceptance.nonformal_status = "pilot_evidence_only";
formalCfg.warning = "Non-formal modes validate the pipeline only; do not copy their positive statuses into the paper as final evidence.";
end

function figureDir = localPaperFigureDir(cfg)
if isfile(cfg.paper_tex)
    figureDir = fullfile(fileparts(cfg.paper_tex), ...
        "DPIM_CI_full_integrated_figures", "formal_e1e8");
else
    figureDir = fullfile(cfg.project_root, "unresolved_paper_figures", "formal_e1e8");
end
figureDir = char(figureDir);
end

function lock = localSeedLock(cfg, cfgMode)
lock = struct();
lock.schema_version = "weighted_paper_seed_lock_v1";
lock.underlying_mode = char(cfgMode);
lock.base_seed = cfg.seed;
lock.evidence_suite_seed = cfg.seed;
lock.rbw_validation_seed = cfg.seed + 770001;
lock.gpu_parity_seed = cfg.seed + 990001;
lock.note = "Validation stages must use seeds derived from this lock; h selected during tuning must not be changed using confirmation results.";
end

function manifest = localManifest(cfg, formalCfg)
manifest = struct();
manifest.schema_version = formalCfg.schema_version;
manifest.created_at = formalCfg.created_at;
manifest.run_mode = formalCfg.run_mode;
manifest.underlying_mode = formalCfg.underlying_mode;
manifest.project_root = formalCfg.project_root;
manifest.results_root = formalCfg.results_root;
manifest.evidence_root = formalCfg.evidence_root;
manifest.coverage_formula_root = formalCfg.coverage_formula_root;
manifest.report_root = formalCfg.report_root;
manifest.paper_tex = cfg.paper_tex;
manifest.paper_figure_dir = formalCfg.paper_figure_dir;
manifest.assumption = "The core paper is the TeX file under the desktop journal-paper folder; this campaign does not edit theory text.";
manifest.main_methods = cfg.main_methods;
manifest.methods = cfg.methods;
manifest.assignment_backend = cfg.assignment_backend;
manifest.n = cfg.n;
manifest.curve_pool_size = cfg.curve_pool_size;
manifest.tuning_pool_size = cfg.tuning_pool_size;
manifest.R_list = cfg.R_list;
manifest.M = cfg.M;
manifest.B = cfg.B;
manifest.h_count = numel(cfg.h_list);
manifest.h_min = min(cfg.h_list);
manifest.h_max = max(cfg.h_list);
end

function T = localPreflight(cfg, formalCfg)
rows = cell(0, 1);
rows{end+1,1} = localCheckRow("paper_tex", isfile(cfg.paper_tex), cfg.paper_tex); %#ok<AGROW>
rows{end+1,1} = localCheckRow("paper_figure_parent", isfolder(fileparts(formalCfg.paper_figure_dir)), fileparts(formalCfg.paper_figure_dir)); %#ok<AGROW>
rows{end+1,1} = localCheckRow("weighted_root", isfolder(cfg.weighted_root), cfg.weighted_root); %#ok<AGROW>
rows{end+1,1} = localCheckRow("main_method_is_rqmc", cfg.main_methods == "sobol_scrambled" && all(string(cfg.methods) == cfg.main_methods), string(cfg.methods)); %#ok<AGROW>
rows{end+1,1} = localCheckRow("no_equal_weight_main_controls", ~cfg.include_control_methods, "include_control_methods=" + string(cfg.include_control_methods)); %#ok<AGROW>
rows{end+1,1} = localCheckRow("h_scan_starts_small", min(cfg.h_list) <= 1e-5 && numel(cfg.h_list) >= 8, sprintf("h_min=%.4g, h_count=%d", min(cfg.h_list), numel(cfg.h_list))); %#ok<AGROW>
rows{end+1,1} = localCheckRow("finite_B_baseline_present", cfg.B >= 39, "B=" + string(cfg.B)); %#ok<AGROW>
if string(formalCfg.run_mode) == "formal"
    rows{end+1,1} = localCheckRow("formal_scale_gate", cfg.M >= 1000 && cfg.B >= 999, sprintf("M=%d, B=%d", cfg.M, cfg.B)); %#ok<AGROW>
else
    rows{end+1,1} = localStatusRow("formal_scale_gate", "WARN", sprintf("run_mode=%s is not final-scale evidence.", formalCfg.run_mode)); %#ok<AGROW>
end

for i = 1:numel(cfg.required_functions)
    name = cfg.required_functions(i);
    rows{end+1,1} = localCheckRow("function:" + name, exist(name, "file") == 2, which(name)); %#ok<AGROW>
end

rows{end+1,1} = localRqmcSmoke(cfg); %#ok<AGROW>
rows{end+1,1} = localGpuParitySmoke(cfg, formalCfg); %#ok<AGROW>
T = vertcat(rows{:});
end

function row = localCheckRow(name, ok, detail)
if ok
    row = localStatusRow(name, "PASS", detail);
else
    row = localStatusRow(name, "FAIL", detail);
end
end

function row = localStatusRow(name, status, detail)
row = table(string(name), string(status), string(detail), ...
    'VariableNames', {'check','status','detail'});
end

function row = localRqmcSmoke(cfg)
try
    [~, info] = dpim_generate_point_set(16, 2, cfg.main_methods, cfg.seed + 990002);
    if string(info.actual_method) == "sobol_scrambled" && ~logical(info.fallback_used)
        row = localStatusRow("rqmc_point_pool_smoke", "PASS", info.message);
    else
        row = localStatusRow("rqmc_point_pool_smoke", "WARN", ...
            "Requested sobol_scrambled but got " + string(info.actual_method) + ": " + string(info.message));
    end
catch ME
    row = localStatusRow("rqmc_point_pool_smoke", "FAIL", ME.message);
end
end

function row = localGpuParitySmoke(cfg, formalCfg)
backend = lower(string(cfg.assignment_backend));
if backend == "cpu"
    row = localStatusRow("gpu_parity_smoke", "PASS", "CPU backend requested; GPU parity is not required.");
    return;
end

try
    dev = gpuDevice;
catch ME
    if backend == "gpu"
        row = localStatusRow("gpu_parity_smoke", "FAIL", "GPU backend requested but gpuDevice failed: " + string(ME.message));
    else
        row = localStatusRow("gpu_parity_smoke", "WARN", "GPU unavailable; auto backend will fall back to CPU: " + string(ME.message));
    end
    return;
end

try
    n = 16;
    d = 2;
    seed = cfg.seed + 990001;
    [U, ~] = dpim_generate_point_set(n, d, "sobol_scrambled", seed);
    ctx = cfg;
    ctx.project_dir = fullfile(formalCfg.results_root, "_preflight", "gpu_parity");
    ctx.point_seed = seed;
    ctx.n = n;
    ctx.d = d;
    ctx.weighting_cfg.voronoi_aux_sample_count = 512;
    ctx.weighting_cfg.voronoi_block_size = 128;
    ctx.weighting_cfg.voronoi_enable_cache = false;
    ctx.weighting_cfg.voronoi_save_outputs = false;

    ctx.weighting_cfg.voronoi_assignment_backend = "cpu";
    cpuWeights = voronoi_ci_probability_weights_provider(n, d, "sobol_scrambled", "gpu_parity_smoke", U, ctx);

    ctx.weighting_cfg.voronoi_assignment_backend = "gpu";
    gpuWeights = voronoi_ci_probability_weights_provider(n, d, "sobol_scrambled", "gpu_parity_smoke", U, ctx);

    diff = abs(double(cpuWeights.weights(:)) - double(gpuWeights.weights(:)));
    l1 = sum(diff);
    linf = max(diff);
    detail = sprintf("%s; L1=%.3g, Linf=%.3g", dev.Name, l1, linf);
    if l1 <= 1e-10 && linf <= 1e-10
        row = localStatusRow("gpu_parity_smoke", "PASS", detail);
    else
        row = localStatusRow("gpu_parity_smoke", "FAIL", "CPU/GPU Voronoi weights differ: " + string(detail));
    end
catch ME
    row = localStatusRow("gpu_parity_smoke", "FAIL", ME.message);
end
end

function outputs = localAssembleFormalOutputs(formalCfg, cfg)
dpimnumeric.ensureDir(formalCfg.report_root);
dpimnumeric.ensureDir(formalCfg.tables_root);

supplemental = localRunSupplementalFormalExperiments(formalCfg, cfg);
map = localPaperExperimentMap(formalCfg);
writetable(map, fullfile(formalCfg.results_root, "paper_experiment_map.csv"));

methodAudit = localBuildFormalMethodAudit(formalCfg);
writetable(methodAudit, fullfile(formalCfg.results_root, "method_audit.csv"));

audit = localBuildFormalClaimAudit(formalCfg, cfg, map, methodAudit);
writetable(audit, fullfile(formalCfg.results_root, "claim_audit.csv"));

localCopyCoreArtifacts(formalCfg);
latexTables = localWriteLatexTables(formalCfg, audit);
figureInventory = localCopyPaperFigures(formalCfg);
notesPath = localWritePaperUpdateNotes(formalCfg, audit, figureInventory);
report = localWriteFormalReport(formalCfg, audit, map, figureInventory, methodAudit);

outputs = struct();
outputs.paper_experiment_map_csv = fullfile(formalCfg.results_root, "paper_experiment_map.csv");
outputs.claim_audit_csv = fullfile(formalCfg.results_root, "claim_audit.csv");
outputs.method_audit_csv = fullfile(formalCfg.results_root, "method_audit.csv");
outputs.paper_update_notes_md = notesPath;
outputs.report_md = report.markdown;
outputs.report_html = report.html;
outputs.latex_tables = latexTables;
outputs.figure_inventory_csv = fullfile(formalCfg.report_root, "paper_figure_inventory.csv");
outputs.paper_figure_dir = formalCfg.paper_figure_dir;
outputs.supplemental = supplemental;
end

function outputs = localRunSupplementalFormalExperiments(formalCfg, cfg)
root = fullfile(formalCfg.results_root, "formal_experiments");
dpimnumeric.ensureDir(root);
outputs = struct("root", root);
if string(formalCfg.run_mode) == "preflight"
    dpimnumeric.writeText(fullfile(root, "README.md"), ...
        "Supplemental formal experiments are skipped in preflight mode." + newline);
    return;
end

[outputs.E5, densityCurves] = localRunFormalE5Density(formalCfg, cfg, root);
outputs.E6 = localRunFormalE6RqmcDiagnosticV2(formalCfg, cfg, root);
outputs.E7 = localRunFormalE7FiniteGridBand(formalCfg, cfg, root, densityCurves);
outputs.E8 = localRunFormalE8HoldoutFormula(formalCfg, root);
end

function [output, curves] = localRunFormalE5Density(formalCfg, cfg, root)
outDir = fullfile(root, "E5_probability_weighted_density");
dpimnumeric.ensureDir(outDir);
dpimnumeric.ensureDir(fullfile(outDir, "figures"));
problem = localSupplementalNormalProblem("E5_probability_weighted_density");
expCfg = localSupplementalExperimentCfg(cfg, outDir);
pointMethod = cfg.main_methods;
fprintf("Building formal E5 weighted density pool: method=%s, curves=%d, n=%d\n", ...
    pointMethod, cfg.curve_pool_size, cfg.n);
curves = dpim_build_weighted_curve_pool(problem, pointMethod, cfg.n, expCfg, cfg.curve_pool_size);
weights = localCurveWeightTable(curves, "E5_probability_weighted_density", pointMethod);
weightStats = localAggregateWeightStats(weights);

yPoints = [-1, 0, 1];
hList = localSupplementalHList(cfg, [0.10, 0.25, 0.50]);
summary = table();
rawBlocks = cell(0, 1);
for ih = 1:numel(hList)
    h = hList(ih);
    for iy = 1:numel(yPoints)
        y0 = yPoints(iy);
        poolVals = dpim_curve_point_estimates(curves, y0, h);
        truth = dpim_gaussian_kernel(y0, sqrt(1 + h^2));
        for ir = 1:numel(cfg.R_list)
            R = cfg.R_list(ir);
            seed = cfg.seed + 51000000 + 10000 * ih + 100 * iy + ir;
            [rows, raw] = localCoverageRowsFromPool(poolVals, truth, cfg, ...
                "E5", "E5_probability_weighted_density", pointMethod, y0, h, R, seed, weightStats);
            summary = [summary; rows]; %#ok<AGROW>
            rawBlocks{end+1,1} = raw; %#ok<AGROW>
        end
    end
end

writetable(summary, fullfile(outDir, "summary.csv"));
writetable(weights, fullfile(outDir, "weight_diagnostics.csv"));
localPlotCoverageBars(outDir, summary, "E5 probability-weighted density");
save(fullfile(outDir, "raw_results.mat"), "summary", "weights", "rawBlocks", "cfg", "formalCfg", "-v7.3");
output = struct("output_dir", outDir, "summary_csv", fullfile(outDir, "summary.csv"), ...
    "weight_csv", fullfile(outDir, "weight_diagnostics.csv"));
end

function output = localRunFormalE6RqmcDiagnostic(~, cfg, root)
outDir = fullfile(root, "E6_rqmc_effective_order");
dpimnumeric.ensureDir(outDir);
problem = localSupplementalNormalProblem("E6_rqmc_effective_order");
methods = ["mc", "sobol_qmc", "sobol_scrambled"];
summaryRows = cell(0, 1);
coverageRows = table();
weightRows = table();
poolVars = nan(numel(methods), 1);
for im = 1:numel(methods)
    method = methods(im);
    expCfg = localSupplementalExperimentCfg(cfg, fullfile(outDir, method));
    fprintf("Building formal E6 point-pool diagnostic: method=%s, curves=%d, n=%d\n", ...
        method, cfg.curve_pool_size, cfg.n);
    curves = dpim_build_weighted_curve_pool(problem, method, cfg.n, expCfg, cfg.curve_pool_size);
    W = localCurveWeightTable(curves, "E6_rqmc_effective_order", method);
    weightRows = [weightRows; W]; %#ok<AGROW>
    weightStats = localAggregateWeightStats(W);

    h = localSupplementalHList(cfg, 0.25);
    h = h(1);
    y0 = 0;
    truth = dpim_gaussian_kernel(y0, sqrt(1 + h^2));
    poolVals = dpim_curve_point_estimates(curves, y0, h);
    stats = localVectorStats(poolVals);
    poolVars(im) = stats.variance;
    seed = cfg.seed + 62000000 + 1000 * im;
    [covRows, ~] = localCoverageRowsFromPool(poolVals, truth, cfg, ...
        "E6", "E6_rqmc_effective_order", method, y0, h, cfg.R_list(1), seed, weightStats);
    coverageRows = [coverageRows; covRows]; %#ok<AGROW>

    summaryRows{end+1,1} = table(method, strjoin(unique(string(W.point_actual_method), "stable"), "/"), ...
        mean(localTruthy(W.point_fallback_used)), cfg.curve_pool_size, cfg.n, ...
        mean(poolVals), stats.variance, stats.skewness, stats.excess_kurtosis, ...
        mean(W.rho3_w, "omitnan"), mean(W.rho4_w, "omitnan"), mean(W.ess_ratio, "omitnan"), NaN, ...
        'VariableNames', {'point_method','actual_methods','point_fallback_rate', ...
        'curve_count','n','mean_estimate','variance_estimate','skewness_estimate', ...
        'excess_kurtosis_estimate','mean_rho3_w','mean_rho4_w','mean_ess_ratio', ...
        'variance_ratio_vs_mc'}); %#ok<AGROW>
end
summary = vertcat(summaryRows{:});
mcVar = poolVars(find(methods == "mc", 1, "first"));
if isfinite(mcVar) && mcVar > 0
    summary.variance_ratio_vs_mc = summary.variance_estimate ./ mcVar;
end
writetable(summary, fullfile(outDir, "summary.csv"));
writetable(coverageRows, fullfile(outDir, "coverage_by_method.csv"));
writetable(weightRows, fullfile(outDir, "weight_diagnostics.csv"));
output = struct("output_dir", outDir, "summary_csv", fullfile(outDir, "summary.csv"), ...
    "coverage_csv", fullfile(outDir, "coverage_by_method.csv"));
end

function output = localRunFormalE6RqmcDiagnosticV2(~, cfg, root)
outDir = fullfile(root, "E6_rqmc_effective_order");
dpimnumeric.ensureDir(outDir);
problem = localSupplementalNormalProblem("E6_rqmc_effective_order");
methods = ["mc", "sobol_qmc", "sobol_scrambled"];
nList = localE6NList(cfg);
summaryRows = cell(0, 1);
coverageRows = table();
weightRows = table();
for in = 1:numel(nList)
    nVal = nList(in);
    for im = 1:numel(methods)
        method = methods(im);
        expCfg = localSupplementalExperimentCfg(cfg, fullfile(outDir, method + "_n" + string(nVal)));
        fprintf("Building formal E6 point-pool diagnostic V2: method=%s, curves=%d, n=%d\n", ...
            method, cfg.curve_pool_size, nVal);
        curves = dpim_build_weighted_curve_pool(problem, method, nVal, expCfg, cfg.curve_pool_size);
        W = localCurveWeightTable(curves, "E6_rqmc_effective_order", method);
        W.point_pool_n = repmat(nVal, height(W), 1);
        weightRows = [weightRows; W]; %#ok<AGROW>
        weightStats = localAggregateWeightStats(W);

        h = localSupplementalHList(cfg, 0.25);
        h = h(1);
        y0 = 0;
        truth = dpim_gaussian_kernel(y0, sqrt(1 + h^2));
        poolVals = dpim_curve_point_estimates(curves, y0, h);
        stats = localVectorStats(poolVals);
        coverageRole = "variance_and_coverage_probe";
        if method == "sobol_qmc"
            coverageRole = "deterministic_error_only";
        else
            seed = cfg.seed + 62000000 + 10000 * in + 1000 * im;
            [covRows, ~] = localCoverageRowsFromPool(poolVals, truth, cfg, ...
                "E6", "E6_rqmc_effective_order", method, y0, h, cfg.R_list(1), seed, weightStats);
            coverageRows = [coverageRows; covRows]; %#ok<AGROW>
        end

        summaryRows{end+1,1} = table(method, strjoin(unique(string(W.point_actual_method), "stable"), "/"), ...
            string(coverageRole), mean(localTruthy(W.point_fallback_used)), cfg.curve_pool_size, nVal, h, truth, ...
            mean(poolVals), abs(mean(poolVals) - truth), stats.variance, stats.skewness, stats.excess_kurtosis, ...
            mean(W.rho3_w, "omitnan"), mean(W.rho4_w, "omitnan"), mean(W.ess_ratio, "omitnan"), NaN, ...
            'VariableNames', {'point_method','actual_methods','coverage_role','point_fallback_rate', ...
            'curve_count','n','h','truth','mean_estimate','abs_bias','variance_estimate', ...
            'skewness_estimate','excess_kurtosis_estimate','mean_rho3_w','mean_rho4_w', ...
            'mean_ess_ratio','variance_ratio_vs_mc_same_n'}); %#ok<AGROW>
    end
end
summary = vertcat(summaryRows{:});
for in = 1:numel(nList)
    nVal = nList(in);
    mcRow = summary(summary.point_method == "mc" & summary.n == nVal, :);
    if ~isempty(mcRow) && isfinite(mcRow.variance_estimate(1)) && mcRow.variance_estimate(1) > 0
        mask = summary.n == nVal;
        summary.variance_ratio_vs_mc_same_n(mask) = summary.variance_estimate(mask) ./ mcRow.variance_estimate(1);
    end
end
effectiveOrder = localE6EffectiveOrderSummary(summary);
writetable(summary, fullfile(outDir, "summary.csv"));
writetable(coverageRows, fullfile(outDir, "coverage_by_method.csv"));
writetable(weightRows, fullfile(outDir, "weight_diagnostics.csv"));
writetable(effectiveOrder, fullfile(outDir, "effective_order_summary.csv"));
output = struct("output_dir", outDir, "summary_csv", fullfile(outDir, "summary.csv"), ...
    "coverage_csv", fullfile(outDir, "coverage_by_method.csv"), ...
    "effective_order_csv", fullfile(outDir, "effective_order_summary.csv"));
end

function output = localRunFormalE7FiniteGridBand(formalCfg, cfg, root, curves)
outDir = fullfile(root, "E7_finite_grid_band");
dpimnumeric.ensureDir(outDir);
dpimnumeric.ensureDir(fullfile(outDir, "figures"));
if nargin < 4 || isempty(curves)
    problem = localSupplementalNormalProblem("E7_finite_grid_band");
    expCfg = localSupplementalExperimentCfg(cfg, outDir);
    curves = dpim_build_weighted_curve_pool(problem, cfg.main_methods, cfg.n, expCfg, cfg.curve_pool_size);
end
h = localSupplementalHList(cfg, 0.25);
h = h(1);

summary = table();
pointwise = table();
rawBlocks = cell(0, 1);
gridCounts = unique(double(cfg.e7_grid_counts(:))).';
for ig = 1:numel(gridCounts)
    gridCount = gridCounts(ig);
    yGrid = linspace(-1, 1, gridCount);
    truth = dpim_gaussian_kernel(yGrid(:), sqrt(1 + h^2)).';
    poolMatrix = zeros(numel(curves), numel(yGrid));
    for iy = 1:numel(yGrid)
        poolMatrix(:, iy) = dpim_curve_point_estimates(curves, yGrid(iy), h);
    end
    for ir = 1:numel(cfg.R_list)
        R = cfg.R_list(ir);
        seed = cfg.seed + 73000000 + 100 * ig + ir;
        [rows, pointRows, raw] = localFiniteGridBandRows(poolMatrix, truth, cfg, ...
            "E7", cfg.main_methods, yGrid, h, R, seed);
        summary = [summary; rows]; %#ok<AGROW>
        pointwise = [pointwise; pointRows]; %#ok<AGROW>
        rawBlocks{end+1,1} = raw; %#ok<AGROW>
    end
end
writetable(summary, fullfile(outDir, "summary.csv"));
writetable(pointwise, fullfile(outDir, "pointwise_summary.csv"));
localPlotCoverageBars(outDir, summary, "E7 finite-grid simultaneous band");
save(fullfile(outDir, "raw_results.mat"), "summary", "pointwise", "rawBlocks", "cfg", "formalCfg", "-v7.3");
output = struct("output_dir", outDir, "summary_csv", fullfile(outDir, "summary.csv"), ...
    "pointwise_csv", fullfile(outDir, "pointwise_summary.csv"));
end

function output = localRunFormalE8HoldoutFormula(formalCfg, root)
outDir = fullfile(root, "E8_weighted_formula_holdout");
dpimnumeric.ensureDir(outDir);
coveragePath = fullfile(formalCfg.coverage_formula_root, "coverage_formula_validation.csv");
T = localReadCsv(coveragePath);
if isempty(T)
    summary = table("missing", "missing", 0, 0, NaN, NaN, NaN, NaN, NaN, ...
        "RBw coverage_formula_validation.csv is missing.", ...
        'VariableNames', {'method','model','train_count','test_count','rmse','mae', ...
        'r_squared','max_fallback_rate','max_inf_rate','note'});
    predictions = table();
else
    [summary, predictions] = localHoldoutFormulaTables(T);
end
writetable(summary, fullfile(outDir, "holdout_summary.csv"));
writetable(predictions, fullfile(outDir, "holdout_predictions.csv"));
output = struct("output_dir", outDir, "summary_csv", fullfile(outDir, "holdout_summary.csv"), ...
    "predictions_csv", fullfile(outDir, "holdout_predictions.csv"));
end

function problem = localSupplementalNormalProblem(name)
problem = struct();
problem.name = char(string(name));
problem.short_name = "FormalNormal";
problem.d = 1;
problem.target_distribution = "standard_normal";
problem.center_transform = "normal_icdf";
problem.provider = "voronoi_ci_probability_weights_provider";
problem.response_fun = @(theta) theta(:,1);
end

function expCfg = localSupplementalExperimentCfg(cfg, outDir)
expCfg = cfg;
expCfg.results_root = outDir;
expCfg.weighting_cfg = cfg.weighting_cfg;
expCfg.weighting_cfg.voronoi_output_dir = fullfile(outDir, "ci_probability_weights");
end

function nList = localE6NList(cfg)
raw = round([cfg.n / 4, cfg.n / 2, cfg.n]);
raw(raw < 16) = 16;
nList = unique(raw, "stable");
if numel(nList) < 3
    nList = unique([max(16, round(cfg.n / 3)), max(16, round(2 * cfg.n / 3)), cfg.n], "stable");
end
end

function T = localE6EffectiveOrderSummary(summary)
methods = unique(string(summary.point_method), "stable");
rows = cell(numel(methods), 1);
for i = 1:numel(methods)
    method = methods(i);
    S = summary(summary.point_method == method, :);
    valid = isfinite(S.variance_estimate) & S.variance_estimate > 0 & S.n > 0;
    slope = NaN;
    order = NaN;
    if sum(valid) >= 2
        coeff = polyfit(log(double(S.n(valid))), log(S.variance_estimate(valid)), 1);
        slope = coeff(1);
        order = -slope;
    end
    maxN = max(S.n);
    atMax = S(S.n == maxN, :);
    ratioAtMax = NaN;
    if ~isempty(atMax)
        ratioAtMax = atMax.variance_ratio_vs_mc_same_n(1);
    end
    note = "Variance slope is estimated from point-pool replicates; it is numerical evidence, not a proof.";
    if method == "sobol_qmc"
        note = "Plain QMC is deterministic here; use abs bias and weight moments, not coverage frequency.";
    end
    rows{i} = table(method, height(S), min(S.n), maxN, slope, order, ...
        mean(S.abs_bias, "omitnan"), max(S.abs_bias, [], "omitnan"), ...
        max(S.point_fallback_rate, [], "omitnan"), mean(S.mean_rho3_w, "omitnan"), ...
        mean(S.mean_rho4_w, "omitnan"), mean(S.mean_ess_ratio, "omitnan"), ratioAtMax, string(note), ...
        'VariableNames', {'method','n_count','n_min','n_max','slope_log_variance_vs_log_n', ...
        'effective_order_variance','mean_abs_bias','max_abs_bias','max_fallback_rate', ...
        'mean_rho3_w','mean_rho4_w','mean_ess_ratio','variance_ratio_at_max_n_vs_mc','note'});
end
T = vertcat(rows{:});
end

function h = localSupplementalHList(cfg, targets)
targets = double(targets(:)).';
h = zeros(size(targets));
for i = 1:numel(targets)
    [~, idx] = min(abs(cfg.h_list - targets(i)));
    h(i) = cfg.h_list(idx);
end
h = unique(h, "stable");
end

function W = localCurveWeightTable(curves, experimentName, method)
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

function [rows, raw] = localCoverageRowsFromPool(poolVals, truth, cfg, ...
    paperExperiment, sourceExperiment, pointMethod, y0, h, R, blockSeed, weightStats)
poolVals = double(poolVals(:));
rng(blockSeed, "twister");
methodNames = ["Student-t"; "Percentile bootstrap"; "Bootstrap-t"; "Hybrid"];
nMethods = numel(methodNames);
cover = false(cfg.M, nMethods);
leftMiss = false(cfg.M, nMethods);
rightMiss = false(cfg.M, nMethods);
isInf = false(cfg.M, nMethods);
lengths = nan(cfg.M, nMethods);
fallback = false(cfg.M, 1);
btInf = false(cfg.M, 1);
poolSize = numel(poolVals);
poolStats = localVectorStats(poolVals);
poolSe = sqrt(max(poolStats.variance, realmin)) / sqrt(R);
biasSeRatio = (poolStats.mean - truth) / max(poolSe, realmin);
for m = 1:cfg.M
    idx = randi(poolSize, R, 1);
    x = poolVals(idx);
    [ci, diag] = localFormalCiMethods(x, truth, cfg.alpha, cfg.B, cfg.lambda);
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
[kMinus, kPlus, C0B] = localFormalFiniteBBaseline(cfg.alpha, cfg.B);
rowCells = cell(nMethods, 1);
for j = 1:nMethods
    coverage = mean(cover(:,j));
    if methodNames(j) == "Student-t"
        baseline = 1 - cfg.alpha;
        baselineType = "nominal_1_minus_alpha";
    else
        baseline = C0B;
        baselineType = "finite_B_C0B";
    end
    finiteLength = isfinite(lengths(:,j));
    rowCells{j} = table(string(paperExperiment), string(sourceExperiment), ...
        string(pointMethod), methodNames(j), y0, h, R, cfg.M, cfg.B, ...
        kMinus, kPlus, C0B, truth, poolStats.mean, poolStats.mean - truth, ...
        sqrt(max(poolStats.variance, realmin)), poolStats.skewness, ...
        poolStats.excess_kurtosis, biasSeRatio, ...
        weightStats.mean_rho3_w, weightStats.mean_rho4_w, ...
        weightStats.mean_n2_eff_w, weightStats.mean_n3_eff_w, ...
        weightStats.mean_n4_eff_w, ...
        coverage, baseline, string(baselineType), coverage - baseline, ...
        abs(coverage - baseline), sqrt(max(coverage * (1 - coverage), realmin) / cfg.M), ...
        mean(leftMiss(:,j)), mean(rightMiss(:,j)), mean(fallback), mean(btInf), ...
        mean(isInf(:,j)), mean(lengths(finiteLength,j), "omitnan"), ...
        median(lengths(finiteLength,j), "omitnan"), sum(finiteLength) / cfg.M, blockSeed, ...
        'VariableNames', {'paper_experiment','source_experiment','point_method', ...
        'method','y0','h','R','M','B','k_minus','k_plus','C0B','truth', ...
        'pool_mean','estimator_bias','pool_sd','gamma','kappa','bias_se_ratio', ...
        'mean_rho3_w','mean_rho4_w','mean_n2_eff_w','mean_n3_eff_w', ...
        'mean_n4_eff_w','coverage','formula_baseline', ...
        'formula_baseline_type','coverage_error','abs_coverage_error', ...
        'coverage_mc_se','left_miss_rate','right_miss_rate','fallback_rate', ...
        'bootstrap_t_inf_rate','interval_inf_rate','mean_interval_length', ...
        'median_interval_length','finite_length_rate','block_seed'});
end
rows = vertcat(rowCells{:});
raw = struct("paper_experiment", string(paperExperiment), "source_experiment", string(sourceExperiment), ...
    "point_method", string(pointMethod), "y0", y0, "h", h, "R", R, ...
    "block_seed", blockSeed, "cover", cover, "interval_inf", isInf, ...
    "interval_lengths", single(lengths), "fallback", fallback, ...
    "bootstrap_t_inf", btInf, "pool_stats", poolStats, ...
    "weight_stats", weightStats);
end

function [ci, diag] = localFormalCiMethods(x, truth, alpha, B, lambda)
x = double(x(:));
R = numel(x);
mu = mean(x);
s = std(x, 0);
se = s / sqrt(R);
tCrit = tinv(1 - alpha / 2, max(R - 1, 1));
ci(1) = localFormalPackCi("Student-t", mu - tCrit * se, mu + tCrit * se, truth);

bootIdx = randi(R, B, R);
boot = x(bootIdx);
bm = mean(boot, 2);
[kMinus, kPlus] = localFormalFiniteBRanks(alpha, B);
bmSorted = sort(bm);
ci(2) = localFormalPackCi("Percentile bootstrap", bmSorted(kMinus), bmSorted(kPlus), truth);

bs = std(boot, 0, 2);
T = sqrt(R) * (bm - mu) ./ bs;
T(bs <= 0 & bm == mu) = 0;
T(bs <= 0 & bm > mu) = Inf;
T(bs <= 0 & bm < mu) = -Inf;
TSorted = sort(T);
btLower = mu - TSorted(kPlus) * se;
btUpper = mu - TSorted(kMinus) * se;
ci(3) = localFormalPackCi("Bootstrap-t", btLower, btUpper, truth);

ratio = ci(3).length / max(ci(2).length, realmin);
trigger = ci(3).infinite || (~ci(3).infinite && ratio > lambda);
if trigger
    ci(4) = localFormalPackCi("Hybrid", ci(2).lower, ci(2).upper, truth);
else
    ci(4) = localFormalPackCi("Hybrid", ci(3).lower, ci(3).upper, truth);
end
diag = struct("fallback_trigger", trigger, "bootstrap_t_infinite", ci(3).infinite);
end

function ci = localFormalPackCi(name, lower, upper, truth)
ci = struct("name", string(name), "lower", lower, "upper", upper, ...
    "length", upper - lower, "infinite", ~(isfinite(lower) && isfinite(upper)), ...
    "contains", false, "left_miss", false, "right_miss", false);
ci.contains = ~ci.infinite && lower <= truth && truth <= upper;
ci.left_miss = ~ci.infinite && truth < lower;
ci.right_miss = ~ci.infinite && truth > upper;
end

function [kMinus, kPlus, C0B] = localFormalFiniteBBaseline(alpha, B)
[kMinus, kPlus] = localFormalFiniteBRanks(alpha, B);
C0B = (kPlus - kMinus) / (B + 1);
end

function [kMinus, kPlus] = localFormalFiniteBRanks(alpha, B)
kMinus = max(1, floor((alpha / 2) * (B + 1)));
kPlus = min(B, ceil((1 - alpha / 2) * (B + 1)));
if kPlus <= kMinus
    error("Invalid finite-B ranks: kMinus=%d, kPlus=%d, B=%d.", kMinus, kPlus, B);
end
end

function stats = localVectorStats(x)
x = double(x(:));
mu = mean(x);
c = x - mu;
v = mean(c.^2);
sd = sqrt(max(v, realmin));
stats = struct("mean", mu, "variance", v, ...
    "skewness", mean(c.^3) / sd^3, ...
    "excess_kurtosis", mean(c.^4) / sd^4 - 3);
end

function [summary, pointwise, raw] = localFiniteGridBandRows(poolMatrix, truth, cfg, ...
    paperExperiment, pointMethod, yGrid, h, R, blockSeed)
rng(blockSeed, "twister");
nGrid = size(poolMatrix, 2);
if nGrid > 1
    gridSpacing = max(diff(double(yGrid)));
else
    gridSpacing = NaN;
end
methodNames = ["Percentile max-deviation band"; "Bootstrap-t max-stat band"; "Hybrid max-stat band"];
nMethods = numel(methodNames);
hit = false(cfg.M, nMethods);
infFlag = false(cfg.M, nMethods);
width = nan(cfg.M, nMethods);
pointHit = false(cfg.M, nGrid, nMethods);
k = min(cfg.B, max(1, ceil((1 - cfg.alpha) * (cfg.B + 1))));
C0BMax = k / (cfg.B + 1);
criticalValue = nan(cfg.M, nMethods);
for m = 1:cfg.M
    idx = randi(size(poolMatrix, 1), R, 1);
    X = poolMatrix(idx, :);
    mu = mean(X, 1);
    se = std(X, 0, 1) / sqrt(R);
    bootIdx = randi(R, cfg.B, R);
    maxDev = nan(cfg.B, 1);
    maxT = nan(cfg.B, 1);
    for b = 1:cfg.B
        xb = X(bootIdx(b,:), :);
        bm = mean(xb, 1);
        bs = std(xb, 0, 1);
        maxDev(b) = max(abs(bm - mu));
        t = sqrt(R) * (bm - mu) ./ bs;
        t(~isfinite(t)) = Inf;
        maxT(b) = max(abs(t));
    end
    devCrit = sort(maxDev);
    devCrit = devCrit(k);
    tCrits = sort(maxT);
    tCrit = tCrits(k);
    criticalValue(m,1) = devCrit;
    criticalValue(m,2) = tCrit;
    [hit(m,1), infFlag(m,1), width(m,1), pointHit(m,:,1)] = ...
        localBandHit(mu - devCrit, mu + devCrit, truth);
    [hit(m,2), infFlag(m,2), width(m,2), pointHit(m,:,2)] = ...
        localBandHit(mu - tCrit * se, mu + tCrit * se, truth);
    if infFlag(m,2) || width(m,2) > cfg.lambda * max(width(m,1), realmin)
        hit(m,3) = hit(m,1);
        infFlag(m,3) = infFlag(m,1);
        width(m,3) = width(m,1);
        pointHit(m,:,3) = pointHit(m,:,1);
        criticalValue(m,3) = devCrit;
    else
        hit(m,3) = hit(m,2);
        infFlag(m,3) = infFlag(m,2);
        width(m,3) = width(m,2);
        pointHit(m,:,3) = pointHit(m,:,2);
        criticalValue(m,3) = tCrit;
    end
end

rows = cell(nMethods, 1);
pointRows = cell(nMethods, 1);
for j = 1:nMethods
    coverage = mean(hit(:,j));
    rows{j} = table(string(paperExperiment), string(pointMethod), methodNames(j), ...
        h, R, cfg.M, cfg.B, nGrid, gridSpacing, k, C0BMax, coverage, ...
        coverage - C0BMax, abs(coverage - C0BMax), ...
        sqrt(max(coverage * (1 - coverage), realmin) / cfg.M), ...
        mean(infFlag(:,j)), mean(width(:,j), "omitnan"), median(width(:,j), "omitnan"), ...
        mean(criticalValue(:,j), "omitnan"), median(criticalValue(:,j), "omitnan"), ...
        NaN, "not_estimated_finite_grid_only", "finite_grid_simultaneous", false, blockSeed, ...
        'VariableNames', {'paper_experiment','point_method','method','h','R','M','B', ...
        'grid_count','grid_spacing','max_stat_order_index','formula_baseline', ...
        'coverage','coverage_error','abs_coverage_error','coverage_mc_se', ...
        'interval_inf_rate','mean_band_width','median_band_width', ...
        'mean_max_stat_quantile','max_stat_quantile','interpolation_remainder', ...
        'interpolation_statement','coverage_scope','continuous_coverage_claimed','block_seed'});
    pointRows{j} = table(repmat(string(paperExperiment), nGrid, 1), ...
        repmat(string(pointMethod), nGrid, 1), repmat(methodNames(j), nGrid, 1), ...
        yGrid(:), repmat(h, nGrid, 1), repmat(R, nGrid, 1), ...
        repmat(nGrid, nGrid, 1), repmat(gridSpacing, nGrid, 1), truth(:), ...
        squeeze(mean(pointHit(:,:,j), 1))', ...
        'VariableNames', {'paper_experiment','point_method','method','y0','h','R', ...
        'grid_count','grid_spacing','truth','pointwise_coverage'});
end
summary = vertcat(rows{:});
pointwise = vertcat(pointRows{:});
raw = struct("hit", hit, "infFlag", infFlag, "width", width, "pointHit", pointHit, ...
    "criticalValue", criticalValue, "truth", truth, "yGrid", yGrid, ...
    "gridSpacing", gridSpacing, "maxStatOrderIndex", k, "C0BMax", C0BMax, ...
    "h", h, "R", R, "blockSeed", blockSeed);
end

function [hit, infFlag, width, pointHit] = localBandHit(lower, upper, truth)
infFlag = any(~isfinite(lower) | ~isfinite(upper));
pointHit = ~infFlag & lower <= truth & truth <= upper;
hit = ~infFlag && all(pointHit);
if infFlag
    width = NaN;
else
    width = mean(upper - lower, "omitnan");
end
end

function [summary, predictions] = localHoldoutFormulaTables(T)
methods = ["Percentile bootstrap", "Bootstrap-t"];
summaryRows = cell(0, 1);
predictionRows = cell(0, 1);
for i = 1:numel(methods)
    methodName = methods(i);
    S = T(T.method == methodName, :);
    if isempty(S)
        summaryRows{end+1,1} = localHoldoutSummaryRow(methodName, "heldout_empirical_terms", 0, 0, ...
            NaN, NaN, NaN, NaN, NaN, "No rows for method."); %#ok<AGROW>
        continue;
    end
    if methodName == "Percentile bootstrap"
        predictors = ["pool_skewness_over_sqrtR", "pool_kurtosis_over_R", "pool_skewness_sq_over_R"];
    else
        predictors = ["pool_kurtosis_over_R", "pool_skewness_sq_over_R", "fallback_rate"];
    end
    [X, y, finiteRows] = localHoldoutDesign(S, predictors, "coverage_error_to_C0B");
    S = S(finiteRows, :);
    if numel(y) < numel(predictors) + 3
        summaryRows{end+1,1} = localHoldoutSummaryRow(methodName, "heldout_empirical_terms", ...
            0, numel(y), NaN, NaN, NaN, max(S.fallback_rate), max(S.interval_inf_rate), ...
            "Too few finite rows for held-out fit."); %#ok<AGROW>
        continue;
    end
    splitKey = double(S.block_seed);
    trainMask = mod(splitKey, 2) == 0;
    if sum(trainMask) < numel(predictors) + 1 || sum(~trainMask) < 2
        trainMask = false(height(S), 1);
        trainMask(1:2:height(S)) = true;
    end
    beta = pinv(X(trainMask, :)) * y(trainMask);
    pred = X * beta;
    resid = y - pred;
    testResid = resid(~trainMask);
    testY = y(~trainMask);
    ssTot = sum((testY - mean(testY)).^2);
    if ssTot <= 0
        r2 = NaN;
    else
        r2 = 1 - sum(testResid.^2) / ssTot;
    end
    summaryRows{end+1,1} = localHoldoutSummaryRow(methodName, "heldout_empirical_terms", ...
        sum(trainMask), sum(~trainMask), sqrt(mean(testResid.^2)), mean(abs(testResid)), ...
        r2, max(S.fallback_rate), max(S.interval_inf_rate), ...
        "Held-out split by block_seed parity; empirical diagnostic, not a theorem proof."); %#ok<AGROW>
    P = S(:, intersect(["scheme","scheme_kind","replicate_index","method","y0","h","R", ...
        "coverage_error_to_C0B","coverage_mc_se","fallback_rate","interval_inf_rate"], ...
        string(S.Properties.VariableNames), "stable"));
    P.model = repmat("heldout_empirical_terms", height(S), 1);
    P.is_train = trainMask;
    P.predicted_error_to_C0B = pred;
    P.prediction_residual = resid;
    P.prediction_abs_residual = abs(resid);
    predictionRows{end+1,1} = P; %#ok<AGROW>
end
summary = vertcat(summaryRows{:});
if isempty(predictionRows)
    predictions = table();
else
    predictions = vertcat(predictionRows{:});
end
end

function [X, y, finiteRows] = localHoldoutDesign(S, predictors, responseName)
X = ones(height(S), numel(predictors) + 1);
for j = 1:numel(predictors)
    X(:, j + 1) = S.(char(predictors(j)));
end
y = S.(char(responseName));
finiteRows = isfinite(y) & all(isfinite(X), 2);
X = X(finiteRows, :);
y = y(finiteRows);
end

function row = localHoldoutSummaryRow(methodName, modelName, trainCount, testCount, ...
    rmse, mae, r2, maxFallback, maxInf, note)
row = table(string(methodName), string(modelName), trainCount, testCount, ...
    rmse, mae, r2, maxFallback, maxInf, string(note), ...
    'VariableNames', {'method','model','train_count','test_count','rmse','mae', ...
    'r_squared','max_fallback_rate','max_inf_rate','note'});
end

function localPlotCoverageBars(outDir, T, ttl)
if isempty(T) || ~ismember("method", string(T.Properties.VariableNames)) || ~ismember("coverage", string(T.Properties.VariableNames))
    return;
end
fig = figure("Visible", "off", "Color", "w");
G = groupsummary(T, "method", "mean", "coverage");
bar(categorical(G.method), G.mean_coverage);
yline(0.95, "k--");
ylim([0 1]);
grid on;
ylabel("mean coverage");
title(ttl, "Interpreter", "none");
try
    exportgraphics(fig, fullfile(outDir, "figures", "method_mean_coverage.png"), "Resolution", 180);
catch
    saveas(fig, fullfile(outDir, "figures", "method_mean_coverage.png"));
end
close(fig);
end

function methodAudit = localBuildFormalMethodAudit(formalCfg)
rows = cell(0, 1);
rows = localAppendMethodAuditRows(rows, "E2", "evidence_confirmation", ...
    fullfile(formalCfg.evidence_root, "E2_scalar_normal_R_order_weighted", "confirmation_summary.csv"));
rows = localAppendMethodAuditRows(rows, "E3", "evidence_confirmation", ...
    fullfile(formalCfg.evidence_root, "E4_standard_normal_nonlinear_tail_weighted", "confirmation_summary.csv"));
rows = localAppendMethodAuditRows(rows, "E4", "instability_confirmation", ...
    fullfile(formalCfg.evidence_root, "E5_standard_normal_bootstrap_t_instability_weighted", "confirmation_summary.csv"));
rows = localAppendMethodAuditRows(rows, "E5", "formal_density", ...
    fullfile(formalCfg.results_root, "formal_experiments", "E5_probability_weighted_density", "summary.csv"));
rows = localAppendE6EffectiveOrderAuditRows(rows, ...
    fullfile(formalCfg.results_root, "formal_experiments", "E6_rqmc_effective_order", "effective_order_summary.csv"));
rows = localAppendMethodAuditRows(rows, "E7", "finite_grid_band", ...
    fullfile(formalCfg.results_root, "formal_experiments", "E7_finite_grid_band", "summary.csv"));
rows = localAppendHoldoutAuditRows(rows, fullfile(formalCfg.results_root, ...
    "formal_experiments", "E8_weighted_formula_holdout", "holdout_summary.csv"));
if isempty(rows)
    methodAudit = table();
else
    methodAudit = vertcat(rows{:});
end
end

function rows = localAppendMethodAuditRows(rows, paperExperiment, sourceName, path)
T = localReadCsv(path);
if isempty(T) || ~ismember("method", string(T.Properties.VariableNames))
    rows{end+1,1} = localMethodAuditRow(paperExperiment, sourceName, "missing", 0, ...
        NaN, NaN, NaN, NaN, NaN, NaN, "missing", "Source table missing or has no method column."); %#ok<AGROW>
    return;
end
if string(paperExperiment) ~= "E4" && ismember("selection_reason", string(T.Properties.VariableNames))
    stableMask = T.selection_reason == "tuning_selected_stable";
    if any(stableMask)
        T = T(stableMask, :);
    end
end
methods = unique(string(T.method), "stable");
for i = 1:numel(methods)
    S = T(string(T.method) == methods(i), :);
    meanCoverage = localMeanIfPresent(S, "coverage");
    minCoverage = localMinIfPresent(S, "coverage");
    meanAbs = localMeanFirstPresent(S, ["abs_coverage_error", "abs_error_to_formula_baseline", "abs_error_to_C0B"]);
    maxFallback = localMaxFirstPresent(S, ["fallback_rate", "max_fallback_rate"]);
    maxInf = localMaxFirstPresent(S, ["interval_inf_rate", "bootstrap_t_inf_rate", "max_inf_rate"]);
    maxMcse = localMaxFirstPresent(S, ["coverage_mc_se", "mean_coverage_mc_se"]);
    methodLabel = methods(i);
    if ~any(contains(lower(methodLabel), ["bootstrap-t", "hybrid", "fallback"]))
        maxFallback = 0;
    end
    status = localMethodGateStatus(paperExperiment, methodLabel, height(S), meanCoverage, minCoverage, meanAbs, maxFallback, maxInf, maxMcse);
    rows{end+1,1} = localMethodAuditRow(paperExperiment, sourceName, methods(i), height(S), ...
        meanCoverage, minCoverage, meanAbs, maxFallback, maxInf, maxMcse, status, ...
        "Method-level gate uses formula/C0B baseline columns when available."); %#ok<AGROW>
end
end

function rows = localAppendHoldoutAuditRows(rows, path)
T = localReadCsv(path);
if isempty(T)
    rows{end+1,1} = localMethodAuditRow("E8", "formula_holdout", "missing", 0, ...
        NaN, NaN, NaN, NaN, NaN, NaN, "missing", "Holdout summary missing."); %#ok<AGROW>
    return;
end
for i = 1:height(T)
    status = "supplement_only";
    if isfinite(T.mae(i)) && T.mae(i) <= 0.025 && T.max_inf_rate(i) <= 0.01
        status = "formula_diagnostic_ready";
    end
    rows{end+1,1} = localMethodAuditRow("E8", "formula_holdout", T.method(i), T.test_count(i), ...
        NaN, NaN, T.mae(i), T.max_fallback_rate(i), T.max_inf_rate(i), NaN, status, T.note(i)); %#ok<AGROW>
end
end

function rows = localAppendE6EffectiveOrderAuditRows(rows, path)
T = localReadCsv(path);
if isempty(T)
    rows{end+1,1} = localMethodAuditRow("E6", "rqmc_effective_order", "missing", 0, ...
        NaN, NaN, NaN, NaN, NaN, NaN, "missing", "E6 effective-order summary missing."); %#ok<AGROW>
    return;
end
for i = 1:height(T)
    method = string(T.method(i));
    status = "supplement_only";
    strongRqmcOrder = T.effective_order_variance(i) >= 0.1 ...
        && T.variance_ratio_at_max_n_vs_mc(i) <= 0.9;
    if method == "sobol_scrambled" && T.n_count(i) >= 3 && T.max_fallback_rate(i) == 0 ...
            && isfinite(T.effective_order_variance(i)) && strongRqmcOrder
        status = "effective_order_ready";
    elseif method == "sobol_scrambled" && T.max_fallback_rate(i) == 0
        status = "rqmc_constructed_weak_order";
    elseif method == "sobol_qmc"
        status = "deterministic_qmc_control";
    end
    note = sprintf("n=%g..%g; variance order=%.4g; variance ratio at max n=%.4g; %s", ...
        T.n_min(i), T.n_max(i), T.effective_order_variance(i), ...
        T.variance_ratio_at_max_n_vs_mc(i), string(T.note(i)));
    rows{end+1,1} = localMethodAuditRow("E6", "rqmc_effective_order", method, T.n_count(i), ...
        NaN, NaN, T.mean_abs_bias(i), T.max_fallback_rate(i), 0, NaN, status, note); %#ok<AGROW>
end
end

function status = localMethodGateStatus(paperExperiment, methodName, rowCount, meanCoverage, minCoverage, meanAbs, maxFallback, maxInf, maxMcse)
if rowCount == 0
    status = "missing";
elseif string(paperExperiment) == "E4"
    status = "failure_boundary";
elseif isfinite(maxInf) && maxInf > 0.01
    status = "failed_boundary";
elseif isfinite(maxFallback) && maxFallback > 0.05 && any(contains(string(methodName), ["Bootstrap", "Hybrid"]))
    status = "failed_boundary";
elseif isfinite(meanAbs) && isfinite(minCoverage) ...
        && meanAbs <= min(0.05, 0.025 + localFiniteOrZero(maxMcse)) ...
        && minCoverage >= 0.90 - max(0.02, 2 * localFiniteOrZero(maxMcse))
    status = "main_candidate";
elseif isfinite(meanAbs) && meanAbs <= 0.05
    status = "supplement_only";
else
    status = "supplement_or_boundary";
end
end

function x = localFiniteOrZero(x)
if ~isfinite(x)
    x = 0;
end
end

function row = localMethodAuditRow(paperExperiment, sourceName, methodName, rowCount, ...
    meanCoverage, minCoverage, meanAbsError, maxFallback, maxInf, maxMcse, status, note)
row = table(string(paperExperiment), string(sourceName), string(methodName), rowCount, ...
    meanCoverage, minCoverage, meanAbsError, maxFallback, maxInf, maxMcse, ...
    string(status), string(note), ...
    'VariableNames', {'paper_experiment','source','method','row_count','mean_coverage', ...
    'min_coverage','mean_abs_error','max_fallback_rate','max_inf_rate','max_mc_se', ...
    'method_status','note'});
end

function value = localMeanIfPresent(T, name)
if ismember(name, string(T.Properties.VariableNames))
    value = mean(T.(char(name)), "omitnan");
else
    value = NaN;
end
end

function value = localMinIfPresent(T, name)
if ismember(name, string(T.Properties.VariableNames))
    value = min(T.(char(name)), [], "omitnan");
else
    value = NaN;
end
end

function value = localMeanFirstPresent(T, names)
value = NaN;
for i = 1:numel(names)
    if ismember(names(i), string(T.Properties.VariableNames))
        value = mean(T.(char(names(i))), "omitnan");
        return;
    end
end
end

function value = localMaxFirstPresent(T, names)
value = NaN;
for i = 1:numel(names)
    if ismember(names(i), string(T.Properties.VariableNames))
        value = max(T.(char(names(i))), [], "omitnan");
        return;
    end
end
end

function map = localPaperExperimentMap(formalCfg)
rows = cell(8, 1);
rows{1} = localMapRow("E1", "Finite-B order-statistic grid error", ...
    "evidence_suite", "E1_finite_B_rank_audit", ...
    fullfile(formalCfg.evidence_root, "E1_finite_B_rank_audit", "summary.csv"), ...
    "main_or_baseline", "Exact C0B integerization; not a 0.95-only target.");
rows{2} = localMapRow("E2", "Gaussian / linear DPIM closed-form benchmark", ...
    "evidence_suite", "E2_scalar_normal_R_order_weighted", ...
    fullfile(formalCfg.evidence_root, "E2_scalar_normal_R_order_weighted", "confirmation_summary.csv"), ...
    "main_candidate", "Scalar normal is the primary positive closed-form benchmark; linear beam is supplementary.");
rows{3} = localMapRow("E3", "Skewed or nonlinear response coverage correction", ...
    "evidence_suite", "E4_standard_normal_nonlinear_tail_weighted", ...
    fullfile(formalCfg.evidence_root, "E4_standard_normal_nonlinear_tail_weighted", "confirmation_summary.csv"), ...
    "main_candidate", "Uses the current nonlinear-tail weighted experiment as the paper E3 source.");
rows{4} = localMapRow("E4", "Bootstrap-t small-denominator instability", ...
    "evidence_suite", "E5_standard_normal_bootstrap_t_instability_weighted", ...
    fullfile(formalCfg.evidence_root, "E5_standard_normal_bootstrap_t_instability_weighted", "mechanism_summary.csv"), ...
    "failure_boundary", "Negative mechanism evidence; do not present as positive coverage success.");
rows{5} = localMapRow("E5", "Probability-weighted DPIM density evidence", ...
    "formal_experiments", "E5_probability_weighted_density", ...
    fullfile(formalCfg.results_root, "formal_experiments", "E5_probability_weighted_density", "summary.csv"), ...
    "main_candidate_after_scale", "Requires Voronoi probability weights and randomized QMC point pools.");
rows{6} = localMapRow("E6", "RQMC point-pool effective-order diagnostic", ...
    "formal_experiments", "E6_rqmc_effective_order", ...
    fullfile(formalCfg.results_root, "formal_experiments", "E6_rqmc_effective_order", "summary.csv"), ...
    "diagnostic_or_supplement", "Verifies that random point selection is RQMC/randomized QMC, not plain equal-weight MC.");
rows{7} = localMapRow("E7", "Finite-grid simultaneous confidence band", ...
    "formal_experiments", "E7_finite_grid_band", ...
    fullfile(formalCfg.results_root, "formal_experiments", "E7_finite_grid_band", "summary.csv"), ...
    "main_candidate_finite_grid_only", "Finite-grid band only; no continuous-curve coverage claim.");
rows{8} = localMapRow("E8", "Weighted R-B-w coverage-expression validation", ...
    "formal_experiments", "E8_weighted_formula_holdout", ...
    fullfile(formalCfg.results_root, "formal_experiments", "E8_weighted_formula_holdout", "holdout_summary.csv"), ...
    "formula_diagnostic", "Numerically checks C0B, rho3/rho4 terms, and A_lambda boundary diagnostics.");
map = vertcat(rows{:});
end

function row = localMapRow(paperExperiment, paperTitle, sourceComponent, legacyLabel, sourcePath, role, note)
row = table(string(paperExperiment), string(paperTitle), string(sourceComponent), ...
    string(legacyLabel), string(sourcePath), string(role), string(note), ...
    'VariableNames', {'paper_experiment','paper_title','source_component', ...
    'legacy_source_label','source_path','role_in_paper','note'});
end

function audit = localBuildFormalClaimAudit(formalCfg, cfg, map, methodAudit)
evidenceAudit = localReadCsv(fullfile(formalCfg.evidence_root, "_analysis", "paper_claim_audit.csv"));
rbwFit = localReadCsv(fullfile(formalCfg.coverage_formula_root, "weighted_RBw_fit_summary.csv"));
rbwCoverage = localReadCsv(fullfile(formalCfg.coverage_formula_root, "coverage_formula_validation.csv"));

rows = cell(0, 1);
rows{end+1,1} = localClaimFromEvidence("E1", "finite_B_rank_formula", map, evidenceAudit, ...
    "E1_finite_B_rank_audit", formalCfg, "C0B baseline for all bootstrap order-statistic comparisons."); %#ok<AGROW>
rows{end+1,1} = localPointwiseMethodClaim("E2", "pointwise_closed_form_or_high_precision", map, methodAudit, formalCfg, ...
    "Gaussian/linear pointwise evidence is method-aware; unstable methods must be reported separately."); %#ok<AGROW>
rows{end+1,1} = localPointwiseMethodClaim("E3", "skew_or_nonlinear_coverage_correction", map, methodAudit, formalCfg, ...
    "Nonlinear-tail evidence is method-aware; exact coefficient-level proof is not implied."); %#ok<AGROW>
rows{end+1,1} = localClaimFromEvidence("E4", "bootstrap_t_small_denominator_boundary", map, evidenceAudit, ...
    "E5_standard_normal_bootstrap_t_instability_weighted", formalCfg, "Boundary/failure result; instability is useful negative evidence."); %#ok<AGROW>

rows{end+1,1} = localWeightedDpimClaim(formalCfg, cfg, map, methodAudit); %#ok<AGROW>
rows{end+1,1} = localRqmcClaim(formalCfg, map, methodAudit); %#ok<AGROW>

rows{end+1,1} = localFiniteGridBandClaim(formalCfg, map, methodAudit); %#ok<AGROW>
rows{end+1,1} = localRbwClaim(formalCfg, map, rbwFit, rbwCoverage, methodAudit); %#ok<AGROW>

audit = vertcat(rows{:});
end

function row = localClaimFromEvidence(paperExperiment, claimType, map, evidenceAudit, legacyLabel, formalCfg, extraNote)
M = map(map.paper_experiment == paperExperiment, :);
A = localFindEvidenceRow(evidenceAudit, legacyLabel);
if isempty(A)
    row = localClaimRow(M, claimType, "missing", NaN, NaN, NaN, NaN, NaN, ...
        "not_available", "Source evidence audit row is missing. " + string(extraNote));
    return;
end

status = localScaleStatus(A.claim_status(1), formalCfg);
decision = localDecisionText(status);
note = string(A.note(1)) + " " + string(extraNote);
if string(A.claim_status(1)) == "main_text_ready" && string(status) ~= "main_text_ready"
    note = note + " Non-formal mode is pipeline evidence only; rerun formal mode before manuscript insertion.";
end
row = localClaimRow(M, claimType, status, ...
    localValue(A, "mean_coverage"), localValue(A, "min_coverage"), ...
    localValue(A, "mean_abs_coverage_error"), localValue(A, "max_fallback_rate"), ...
    localValue(A, "max_bootstrap_t_inf_rate"), decision, note);
end

function row = localPointwiseMethodClaim(paperExperiment, claimType, map, methodAudit, formalCfg, extraNote)
M = map(map.paper_experiment == paperExperiment, :);
if isempty(methodAudit)
    row = localClaimRow(M, claimType, "missing", NaN, NaN, NaN, NaN, NaN, ...
        "not_ready", "No method-level audit available. " + string(extraNote));
    return;
end
S = methodAudit(methodAudit.paper_experiment == string(paperExperiment), :);
if isempty(S)
    row = localClaimRow(M, claimType, "missing", NaN, NaN, NaN, NaN, NaN, ...
        "not_ready", "No method-level rows available. " + string(extraNote));
    return;
end
[~, bestIdx] = min(S.mean_abs_error);
best = S(bestIdx, :);
candidateMask = S.method_status == "main_candidate";
if any(candidateMask)
    candidate = S(candidateMask, :);
    [~, loc] = min(candidate.mean_abs_error);
    best = candidate(loc, :);
    if string(formalCfg.run_mode) == "formal"
        status = "main_text_ready";
        decision = "eligible_for_main_text_after_wording_review";
    else
        status = string(formalCfg.acceptance.nonformal_status);
        decision = "pipeline_validated_but_not_final_scale";
    end
elseif any(S.method_status == "supplement_only")
    status = "supplement_only";
    decision = "supplement_or_cautious_support_only";
else
    status = "failed_boundary";
    decision = "use_as_boundary_or_failure_mechanism";
end
note = sprintf("Best method `%s`: mean coverage=%.4g, min coverage=%.4g, mean abs error=%.4g, max Inf=%.4g. %s %s", ...
    best.method(1), best.mean_coverage(1), best.min_coverage(1), best.mean_abs_error(1), ...
    best.max_inf_rate(1), extraNote, localMethodAuditNote(methodAudit, paperExperiment));
row = localClaimRow(M, claimType, status, best.mean_coverage(1), best.min_coverage(1), ...
    best.mean_abs_error(1), best.max_fallback_rate(1), best.max_inf_rate(1), decision, note);
end

function status = localScaleStatus(status, formalCfg)
status = string(status);
if status == "main_text_ready" && string(formalCfg.run_mode) ~= "formal"
    status = string(formalCfg.acceptance.nonformal_status);
end
end

function decision = localDecisionText(status)
switch string(status)
    case "main_text_ready"
        decision = "eligible_for_main_text_after_wording_review";
    case "pilot_evidence_only"
        decision = "pipeline_validated_but_not_final_scale";
    case "supplement_only"
        decision = "supplement_or_cautious_support_only";
    case "failed_boundary"
        decision = "use_as_boundary_or_failure_mechanism";
    case "diagnostic_only"
        decision = "diagnostic_not_main_claim";
    otherwise
        decision = "not_ready";
end
end

function row = localWeightedDpimClaim(formalCfg, cfg, map, methodAudit)
M = map(map.paper_experiment == "E5", :);
summaryPath = fullfile(formalCfg.results_root, "formal_experiments", ...
    "E5_probability_weighted_density", "summary.csv");
weightPath = fullfile(formalCfg.results_root, "formal_experiments", ...
    "E5_probability_weighted_density", "weight_diagnostics.csv");
S = localReadCsv(summaryPath);
W = localReadCsv(weightPath);
if isempty(S) || isempty(W)
    row = localClaimRow(M, "probability_weighted_density_evidence", "missing", ...
        NaN, NaN, NaN, NaN, NaN, "not_ready", ...
        "Dedicated E5 probability-weighted density outputs are missing.");
    return;
end
weightedOk = all(contains(lower(string(W.weight_source)), "voronoi")) ...
    || all(contains(lower(string(W.weight_source)), "external"));
rqmcOk = all(string(W.point_actual_method) == "sobol_scrambled") ...
    && ~any(localTruthy(W.point_fallback_used));
M5 = methodAudit(methodAudit.paper_experiment == "E5", :);
if isempty(M5)
    row = localClaimRow(M, "probability_weighted_density_evidence", "missing", ...
        NaN, NaN, NaN, NaN, NaN, "not_ready", ...
        "Method-level E5 audit rows are missing.");
    return;
end
candidateMask = M5.method_status == "main_candidate";
if any(candidateMask)
    candidates = M5(candidateMask, :);
else
    candidates = M5;
end
[~, bestIdx] = min(candidates.mean_abs_error);
best = candidates(bestIdx, :);
meanCov = best.mean_coverage;
minCov = best.min_coverage;
meanAbs = best.mean_abs_error;
maxFallback = best.max_fallback_rate;
maxInf = best.max_inf_rate;
maxMcse = best.max_mc_se;
goodCoverage = any(candidateMask);
if weightedOk && rqmcOk && goodCoverage
    if string(formalCfg.run_mode) == "formal" && cfg.M >= formalCfg.acceptance.formal_min_M
        status = "main_text_ready";
    else
        status = string(formalCfg.acceptance.nonformal_status);
    end
else
    status = "supplement_only";
end
methodNote = localMethodAuditNote(methodAudit, "E5");
note = sprintf("Best stable method `%s`; dedicated E5 rows=%d; weight rows=%d; probability-weighted=%d; RQMC=%d; mean abs error=%.4g; max fallback=%.4g; max Inf=%.4g; max MCSE=%.4g. %s", ...
    best.method(1), height(S), height(W), weightedOk, rqmcOk, meanAbs, maxFallback, maxInf, maxMcse, methodNote);
row = localClaimRow(M, "probability_weighted_density_evidence", status, ...
    meanCov, minCov, meanAbs, maxFallback, maxInf, localDecisionText(status), note);
end

function row = localRqmcClaim(formalCfg, map, methodAudit)
M = map(map.paper_experiment == "E6", :);
summaryPath = fullfile(formalCfg.results_root, "formal_experiments", ...
    "E6_rqmc_effective_order", "summary.csv");
orderPath = fullfile(formalCfg.results_root, "formal_experiments", ...
    "E6_rqmc_effective_order", "effective_order_summary.csv");
S = localReadCsv(summaryPath);
E = localReadCsv(orderPath);
if isempty(S) || isempty(E)
    row = localClaimRow(M, "rqmc_effective_order_diagnostic", "missing", ...
        NaN, NaN, NaN, NaN, NaN, "not_ready", "Dedicated E6 RQMC diagnostic outputs are missing.");
    return;
end
rqmc = E(E.method == "sobol_scrambled", :);
fallbackRate = max(S.point_fallback_rate, [], "omitnan");
meanRho3 = mean(S.mean_rho3_w, "omitnan");
meanRho4 = mean(S.mean_rho4_w, "omitnan");
meanEss = mean(S.mean_ess_ratio, "omitnan");
rqmcConstructed = ~isempty(rqmc) && rqmc.max_fallback_rate(1) == 0 ...
    && isfinite(rqmc.effective_order_variance(1));
rqmcReady = rqmcConstructed && rqmc.effective_order_variance(1) >= 0.1 ...
    && rqmc.variance_ratio_at_max_n_vs_mc(1) <= 0.9;
if rqmcReady
    if string(formalCfg.run_mode) == "formal"
        status = "main_text_ready";
        decision = "rqmc_effective_order_documented";
    else
        status = string(formalCfg.acceptance.nonformal_status);
        decision = "pipeline_validated_but_not_final_scale";
    end
elseif rqmcConstructed
    status = "supplement_only";
    decision = "rqmc_constructed_but_effective_order_weak";
else
    status = "diagnostic_only";
    decision = "rqmc_fallback_needs_review";
end
methodNote = localMethodAuditNote(methodAudit, "E6");
if isempty(rqmc)
    rqmcOrder = NaN;
    rqmcRatio = NaN;
else
    rqmcOrder = rqmc.effective_order_variance(1);
    rqmcRatio = rqmc.variance_ratio_at_max_n_vs_mc(1);
end
note = sprintf("methods=%s; fallback_rate=%.4g; mean rho3=%.4g; mean rho4=%.4g; mean ESS ratio=%.4g; scrambled RQMC variance order=%.4g; ratio at max n=%.4g. Plain QMC is deterministic control only. %s", ...
    strjoin(unique(string(S.point_method), "stable"), "/"), fallbackRate, meanRho3, meanRho4, meanEss, rqmcOrder, rqmcRatio, methodNote);
row = localClaimRow(M, "rqmc_effective_order_diagnostic", status, ...
    NaN, NaN, NaN, NaN, NaN, decision, note);
end

function row = localFiniteGridBandClaim(formalCfg, map, methodAudit)
M = map(map.paper_experiment == "E7", :);
summaryPath = fullfile(formalCfg.results_root, "formal_experiments", ...
    "E7_finite_grid_band", "summary.csv");
S = localReadCsv(summaryPath);
if isempty(S)
    row = localClaimRow(M, "finite_grid_simultaneous_band", "missing", ...
        NaN, NaN, NaN, NaN, NaN, "not_ready", "Dedicated E7 finite-grid band outputs are missing.");
    return;
end
stable = S(S.method == "Percentile max-deviation band" | S.method == "Hybrid max-stat band", :);
if isempty(stable)
    row = localClaimRow(M, "finite_grid_simultaneous_band", "missing", ...
        NaN, NaN, NaN, NaN, NaN, "not_ready", ...
        "Dedicated E7 summary has no stable finite-grid band rows.");
    return;
end
M7 = methodAudit(methodAudit.paper_experiment == "E7", :);
candidateMask = false(height(M7), 1);
if ~isempty(M7)
    candidateMask = M7.method_status == "main_candidate";
end
if any(candidateMask)
    candidates = M7(candidateMask, :);
    [~, bestIdx] = min(candidates.mean_abs_error);
    best = candidates(bestIdx, :);
    meanCov = best.mean_coverage;
    minCov = best.min_coverage;
    meanAbs = best.mean_abs_error;
    maxInf = best.max_inf_rate;
    if string(formalCfg.run_mode) == "formal"
        status = "main_text_ready";
        decision = "finite_grid_band_ready_with_cautious_wording";
    else
        status = string(formalCfg.acceptance.nonformal_status);
        decision = "pipeline_validated_but_not_final_scale";
    end
else
    meanCov = mean(stable.coverage, "omitnan");
    minCov = min(stable.coverage, [], "omitnan");
    meanAbs = mean(stable.abs_coverage_error, "omitnan");
    maxInf = max(stable.interval_inf_rate, [], "omitnan");
    if meanAbs <= 0.05 && minCov >= 0.85 && maxInf <= 0.01
        status = "supplement_only";
        decision = "finite_grid_band_ready_with_cautious_wording";
    else
        status = "failed_boundary";
        decision = "use_as_boundary_or_failure_mechanism";
    end
end
methodNote = localMethodAuditNote(methodAudit, "E7");
note = sprintf("Finite-grid only, grid count=%d; mean abs error=%.4g; max Inf=%.4g. Does not prove continuous-curve coverage. %s", ...
    max(S.grid_count), meanAbs, maxInf, methodNote);
row = localClaimRow(M, "finite_grid_simultaneous_band", status, ...
    meanCov, minCov, meanAbs, NaN, maxInf, decision, note);
end

function row = localRbwClaim(formalCfg, map, rbwFit, rbwCoverage, methodAudit)
M = map(map.paper_experiment == "E8", :);
if isempty(rbwFit) || isempty(rbwCoverage)
    row = localClaimRow(M, "weighted_RBw_formula_diagnostic", "missing", ...
        NaN, NaN, NaN, NaN, NaN, "not_ready", "RBw validation outputs are missing.");
    return;
end
finiteMae = rbwFit.mae(isfinite(rbwFit.mae));
if isempty(finiteMae)
    bestMae = NaN;
else
    bestMae = min(finiteMae);
end
bt = rbwCoverage(rbwCoverage.method == "Bootstrap-t", :);
if isempty(bt)
    maxFallback = NaN;
    maxInf = NaN;
else
    maxFallback = max(bt.fallback_rate);
    maxInf = max(bt.interval_inf_rate);
end
holdout = localReadCsv(fullfile(formalCfg.results_root, "formal_experiments", ...
    "E8_weighted_formula_holdout", "holdout_summary.csv"));
if isempty(holdout)
    holdoutMae = NaN;
else
    finiteHoldoutMae = holdout.mae(isfinite(holdout.mae));
    if isempty(finiteHoldoutMae)
        holdoutMae = NaN;
    else
        holdoutMae = min(finiteHoldoutMae);
    end
end
status = "supplement_only";
decision = "formula_structure_diagnostic_not_proof";
methodNote = localMethodAuditNote(methodAudit, "E8");
note = sprintf("Best in-sample MAE=%.4g; best held-out MAE=%.4g; Bootstrap-t max fallback=%.4g; max Inf=%.4g. This validates formula structure numerically, not the unknown A coefficients from first principles. %s", ...
    bestMae, holdoutMae, maxFallback, maxInf, methodNote);
if string(formalCfg.run_mode) ~= "formal"
    note = note + " Non-formal mode is not final-scale formula evidence.";
end
row = localClaimRow(M, "weighted_RBw_formula_diagnostic", status, ...
    NaN, NaN, holdoutMae, maxFallback, maxInf, decision, note);
end

function row = localClaimRow(mapRow, claimType, status, meanCoverage, minCoverage, ...
    meanAbsError, maxFallback, maxInf, decision, note)
row = table(mapRow.paper_experiment(1), mapRow.paper_title(1), ...
    string(claimType), string(status), string(decision), ...
    string(mapRow.source_component(1)), string(mapRow.legacy_source_label(1)), ...
    string(mapRow.source_path(1)), meanCoverage, minCoverage, meanAbsError, ...
    maxFallback, maxInf, string(note), ...
    'VariableNames', {'paper_experiment','paper_title','claim_type', ...
    'claim_status','main_text_decision','source_component', ...
    'legacy_source_label','source_path','mean_coverage','min_coverage', ...
    'mean_abs_coverage_error','max_fallback_rate','max_inf_rate','note'});
end

function A = localFindEvidenceRow(evidenceAudit, legacyLabel)
if isempty(evidenceAudit) || ~ismember("experiment", string(evidenceAudit.Properties.VariableNames))
    A = table();
    return;
end
A = evidenceAudit(evidenceAudit.experiment == string(legacyLabel), :);
if ~isempty(A)
    A = A(1, :);
end
end

function value = localValue(T, name)
if isempty(T) || ~ismember(name, string(T.Properties.VariableNames))
    value = NaN;
else
    value = T.(char(name))(1);
end
end

function status = localStatusOrMissing(T)
if isempty(T) || ~ismember("claim_status", string(T.Properties.VariableNames))
    status = "missing";
else
    status = string(T.claim_status(1));
end
end

function note = localMethodAuditNote(methodAudit, paperExperiment)
if isempty(methodAudit) || ~ismember("paper_experiment", string(methodAudit.Properties.VariableNames))
    note = "No method-level audit available.";
    return;
end
S = methodAudit(methodAudit.paper_experiment == string(paperExperiment), :);
if isempty(S)
    note = "No method-level rows for this paper experiment.";
    return;
end
statusBits = strings(height(S), 1);
for i = 1:height(S)
    statusBits(i) = string(S.method(i)) + ":" + string(S.method_status(i));
end
note = "Method audit: " + strjoin(statusBits, "; ") + ".";
end

function T = localReadCsv(path)
if ~isfile(path)
    T = table();
    return;
end
T = readtable(path, TextType="string", Delimiter=",");
end

function W = localReadAllWeightDiagnostics(evidenceRoot)
patterns = ["E2_scalar_normal_R_order_weighted", ...
    "E3_linear_beam_weighted", ...
    "E4_standard_normal_nonlinear_tail_weighted", ...
    "E5_standard_normal_bootstrap_t_instability_weighted", ...
    "E8_plate_RFEM_simultaneous_band_weighted_diagnostic"];
blocks = cell(0, 1);
for i = 1:numel(patterns)
    p = fullfile(evidenceRoot, patterns(i), "weight_diagnostics.csv");
    if isfile(p)
        T = readtable(p, TextType="string", Delimiter=",");
        T.formal_source_experiment = repmat(patterns(i), height(T), 1);
        blocks{end+1,1} = T; %#ok<AGROW>
    end
end
if isempty(blocks)
    W = table();
else
    W = vertcat(blocks{:});
end
end

function mask = localTruthy(x)
if isnumeric(x) || islogical(x)
    mask = x ~= 0;
else
    sx = lower(strtrim(string(x)));
    mask = sx == "1" | sx == "true";
end
mask = mask(:);
end

function localCopyCoreArtifacts(formalCfg)
localCopyOrNotRun(fullfile(formalCfg.coverage_formula_root, "coverage_formula_validation.csv"), ...
    fullfile(formalCfg.results_root, "coverage_formula_validation.csv"));
copyPairs = [
    "weighted_RBw_fit_summary.csv"
    "weighted_RBw_fit_coefficients.csv"
    "weighted_RBw_fit_predictions.csv"
    "weighted_RBw_replicate_stability.csv"
    "weight_moment_summary.csv"
    "weight_diagnostics.csv"];
for i = 1:numel(copyPairs)
    localCopyIfFile(fullfile(formalCfg.coverage_formula_root, copyPairs(i)), ...
        fullfile(formalCfg.report_root, copyPairs(i)));
end
localCopyIfFile(fullfile(formalCfg.evidence_root, "_analysis", "paper_claim_audit.csv"), ...
    fullfile(formalCfg.report_root, "legacy_evidence_claim_audit.csv"));
localCopyIfFile(fullfile(formalCfg.evidence_root, "_analysis", "claim_status_summary.csv"), ...
    fullfile(formalCfg.report_root, "legacy_evidence_claim_status_summary.csv"));
end

function localCopyOrNotRun(src, dst)
if isfile(src)
    localCopyIfFile(src, dst);
else
    T = table("not_run", "RBw validation has not been run in this campaign.", ...
        'VariableNames', {'status','note'});
    writetable(T, dst);
end
end

function localCopyIfFile(src, dst)
if ~isfile(src)
    return;
end
[parentDir, ~, ~] = fileparts(dst);
dpimnumeric.ensureDir(parentDir);
copyfile(src, dst, "f");
end

function latexTables = localWriteLatexTables(formalCfg, audit)
claimPath = fullfile(formalCfg.tables_root, "formal_claim_audit_table.tex");
fitPath = fullfile(formalCfg.tables_root, "formal_rbw_fit_summary_table.tex");
methodPath = fullfile(formalCfg.tables_root, "formal_method_audit_table.tex");

lines = strings(0, 1);
lines(end+1) = "\begin{tabular}{llll}";
lines(end+1) = "\hline";
lines(end+1) = "Exp. & Claim & Status & Decision \\";
lines(end+1) = "\hline";
for i = 1:height(audit)
    lines(end+1) = sprintf("%s & %s & %s & %s \\\\", ...
        localTex(audit.paper_experiment(i)), localTex(audit.claim_type(i)), ...
        localTex(audit.claim_status(i)), localTex(audit.main_text_decision(i)));
end
lines(end+1) = "\hline";
lines(end+1) = "\end{tabular}";
dpimnumeric.writeText(claimPath, strjoin(lines, newline));

fit = localReadCsv(fullfile(formalCfg.coverage_formula_root, "weighted_RBw_fit_summary.csv"));
lines = strings(0, 1);
lines(end+1) = "\begin{tabular}{lrrr}";
lines(end+1) = "\hline";
lines(end+1) = "Model & Rows & MAE & $R^2$ \\";
lines(end+1) = "\hline";
if isempty(fit)
    lines(end+1) = "not run & 0 & -- & -- \\";
else
    maxRows = min(height(fit), 8);
    for i = 1:maxRows
        lines(end+1) = sprintf("%s & %d & %s & %s \\\\", ...
            localTex(fit.model(i)), fit.row_count(i), ...
            localNumTex(fit.mae(i)), localNumTex(fit.r_squared(i)));
    end
end
lines(end+1) = "\hline";
lines(end+1) = "\end{tabular}";
dpimnumeric.writeText(fitPath, strjoin(lines, newline));

methodAudit = localReadCsv(fullfile(formalCfg.results_root, "method_audit.csv"));
lines = strings(0, 1);
lines(end+1) = "\begin{tabular}{lllr}";
lines(end+1) = "\hline";
lines(end+1) = "Exp. & Method & Status & Mean abs. err. \\";
lines(end+1) = "\hline";
if isempty(methodAudit)
    lines(end+1) = "not run & -- & missing & -- \\";
else
    maxRows = min(height(methodAudit), 14);
    for i = 1:maxRows
        lines(end+1) = sprintf("%s & %s & %s & %s \\\\", ...
            localTex(methodAudit.paper_experiment(i)), localTex(methodAudit.method(i)), ...
            localTex(methodAudit.method_status(i)), localNumTex(methodAudit.mean_abs_error(i)));
    end
end
lines(end+1) = "\hline";
lines(end+1) = "\end{tabular}";
dpimnumeric.writeText(methodPath, strjoin(lines, newline));

latexTables = [string(claimPath); string(fitPath); string(methodPath)];
end

function s = localTex(x)
s = char(string(x));
s = strrep(s, "\", "\textbackslash{}");
s = strrep(s, "_", "\_");
s = strrep(s, "%", "\%");
s = strrep(s, "&", "\&");
s = strrep(s, "#", "\#");
end

function s = localNumTex(x)
if isnan(x)
    s = "--";
else
    s = sprintf("%.4g", x);
end
end

function figureInventory = localCopyPaperFigures(formalCfg)
dpimnumeric.ensureDir(formalCfg.paper_figure_dir);
rows = cell(0, 1);
figureSpecs = [
    localFigureSpec("E2_bootstrap_t_h_scan", fullfile(formalCfg.evidence_root, "E2_scalar_normal_R_order_weighted", "figures", "bootstrap_t_coverage_by_h.png"))
    localFigureSpec("E3_nonlinear_tail_h_scan", fullfile(formalCfg.evidence_root, "E4_standard_normal_nonlinear_tail_weighted", "figures", "bootstrap_t_coverage_by_h.png"))
    localFigureSpec("E4_bootstrap_t_instability_h_scan", fullfile(formalCfg.evidence_root, "E5_standard_normal_bootstrap_t_instability_weighted", "figures", "bootstrap_t_coverage_by_h.png"))
    localFigureSpec("E5_probability_weighted_density_methods", fullfile(formalCfg.results_root, "formal_experiments", "E5_probability_weighted_density", "figures", "method_mean_coverage.png"))
    localFigureSpec("E7_finite_grid_band_methods", fullfile(formalCfg.results_root, "formal_experiments", "E7_finite_grid_band", "figures", "method_mean_coverage.png"))
    localFigureSpec("E8_weight_moment_ranges", fullfile(formalCfg.coverage_formula_root, "figures", "weight_moment_ranges.png"))
    localFigureSpec("E8_observed_vs_predicted_error", fullfile(formalCfg.coverage_formula_root, "figures", "observed_vs_predicted_coverage_error.png"))];

for i = 1:numel(figureSpecs)
    spec = figureSpecs(i);
    dst = fullfile(formalCfg.paper_figure_dir, spec.name + ".png");
    copied = false;
    if isfile(spec.source)
        copyfile(spec.source, dst, "f");
        copied = true;
    end
    rows{end+1,1} = table(spec.name, string(spec.source), string(dst), copied, ...
        'VariableNames', {'figure_id','source_path','paper_path','copied'}); %#ok<AGROW>
end
figureInventory = vertcat(rows{:});
writetable(figureInventory, fullfile(formalCfg.report_root, "paper_figure_inventory.csv"));
end

function spec = localFigureSpec(name, source)
spec = struct("name", string(name), "source", string(source));
end

function notesPath = localWritePaperUpdateNotes(formalCfg, audit, figureInventory)
notesPath = fullfile(formalCfg.results_root, "paper_update_notes.md");
lines = strings(0, 1);
lines(end+1) = "# Paper Update Notes";
lines(end+1) = "";
lines(end+1) = sprintf("- Campaign mode: `%s` (underlying mode `%s`).", formalCfg.run_mode, formalCfg.underlying_mode);
lines(end+1) = sprintf("- Core TeX: `%s`.", formalCfg.paper_tex);
lines(end+1) = sprintf("- Paper figure directory: `%s`.", formalCfg.paper_figure_dir);
lines(end+1) = "";
lines(end+1) = "## Claim Boundaries";
for i = 1:height(audit)
    lines(end+1) = sprintf("- `%s`: status `%s`; decision `%s`; note: %s", ...
        audit.paper_experiment(i), audit.claim_status(i), audit.main_text_decision(i), audit.note(i));
end
lines(end+1) = "";
lines(end+1) = "## Figure Copy Status";
for i = 1:height(figureInventory)
    lines(end+1) = sprintf("- `%s`: copied=%d; target `%s`.", ...
        figureInventory.figure_id(i), figureInventory.copied(i), figureInventory.paper_path(i));
end
lines(end+1) = "";
lines(end+1) = "## Do Not Overclaim";
lines(end+1) = "- Non-formal modes are pipeline checks only.";
lines(end+1) = "- RBw formula diagnostics are numerical support, not a proof of the unknown expansion coefficients.";
lines(end+1) = "- Failed bootstrap-t or band rows must be written as boundary mechanisms, not positive coverage successes.";
dpimnumeric.writeText(notesPath, strjoin(lines, newline));
end

function report = localWriteFormalReport(formalCfg, audit, map, figureInventory, methodAudit)
reportMd = fullfile(formalCfg.report_root, "weighted_paper_formal_campaign_report.md");
reportHtml = fullfile(formalCfg.report_root, "weighted_paper_formal_campaign_report.html");
statusSummary = groupsummary(audit, "claim_status");
writetable(statusSummary, fullfile(formalCfg.report_root, "claim_status_summary.csv"));

lines = strings(0, 1);
lines(end+1) = "# Weighted Paper Formal Campaign Report";
lines(end+1) = "";
lines(end+1) = sprintf("- Results root: `%s`", formalCfg.results_root);
lines(end+1) = sprintf("- Run mode: `%s`; underlying mode: `%s`", formalCfg.run_mode, formalCfg.underlying_mode);
lines(end+1) = sprintf("- Core paper: `%s`", formalCfg.paper_tex);
lines(end+1) = "";
lines(end+1) = "## Status Summary";
for i = 1:height(statusSummary)
    lines(end+1) = sprintf("- `%s`: %d", statusSummary.claim_status(i), statusSummary.GroupCount(i));
end
lines(end+1) = "";
lines(end+1) = "## Paper E1--E8 Map";
for i = 1:height(map)
    lines(end+1) = sprintf("- `%s`: %s; source `%s`.", ...
        map.paper_experiment(i), map.paper_title(i), map.legacy_source_label(i));
end
lines(end+1) = "";
lines(end+1) = "## Claim Audit";
for i = 1:height(audit)
    lines(end+1) = sprintf("- `%s`: `%s`, decision `%s`.", ...
        audit.paper_experiment(i), audit.claim_status(i), audit.main_text_decision(i));
end
lines(end+1) = "";
lines(end+1) = "## Method Audit";
if isempty(methodAudit)
    lines(end+1) = "- No method-level audit rows were generated.";
else
    for i = 1:height(methodAudit)
        lines(end+1) = sprintf("- `%s` / `%s`: `%s`, mean_abs_error=%s, max_inf=%s.", ...
            methodAudit.paper_experiment(i), methodAudit.method(i), methodAudit.method_status(i), ...
            localMdNum(methodAudit.mean_abs_error(i)), localMdNum(methodAudit.max_inf_rate(i)));
    end
end
lines(end+1) = "";
lines(end+1) = "## Figures";
for i = 1:height(figureInventory)
    lines(end+1) = sprintf("- `%s`: copied=%d.", figureInventory.figure_id(i), figureInventory.copied(i));
end
dpimnumeric.writeText(reportMd, strjoin(lines, newline));

htmlLines = strings(0, 1);
htmlLines(end+1) = "<!doctype html><html><head><meta charset=""utf-8""><title>Weighted Paper Formal Campaign</title>";
htmlLines(end+1) = "<style>body{font-family:Segoe UI,Microsoft YaHei,sans-serif;max-width:1180px;margin:32px auto;line-height:1.55;color:#1b211b}table{border-collapse:collapse;width:100%;font-size:14px}td,th{border-bottom:1px solid #ddd;padding:7px;text-align:left}th{background:#f4efe5}.main_text_ready{color:#137333;font-weight:700}.pilot_evidence_only{color:#8a5a00;font-weight:700}.supplement_only{color:#8a5a00;font-weight:700}.failed_boundary{color:#9b2f24;font-weight:700}.diagnostic_only{color:#315f8c;font-weight:700}.missing{color:#777;font-weight:700}code{background:#f4efe5;padding:2px 4px}</style></head><body>";
htmlLines(end+1) = "<h1>Weighted Paper Formal Campaign</h1>";
htmlLines(end+1) = "<p>Results root: <code>" + localHtml(formalCfg.results_root) + "</code></p>";
htmlLines(end+1) = "<p>Run mode: <code>" + localHtml(formalCfg.run_mode) + "</code>; underlying mode: <code>" + localHtml(formalCfg.underlying_mode) + "</code></p>";
htmlLines(end+1) = "<h2>Claim Audit</h2><table><tr><th>Paper Exp.</th><th>Claim</th><th>Status</th><th>Decision</th><th>Note</th></tr>";
for i = 1:height(audit)
    cls = audit.claim_status(i);
    htmlLines(end+1) = "<tr><td>" + localHtml(audit.paper_experiment(i)) + "</td><td>" + ...
        localHtml(audit.claim_type(i)) + "</td><td class=""" + cls + """>" + ...
        localHtml(cls) + "</td><td>" + localHtml(audit.main_text_decision(i)) + ...
        "</td><td>" + localHtml(audit.note(i)) + "</td></tr>";
end
htmlLines(end+1) = "</table><h2>Method Audit</h2><table><tr><th>Paper Exp.</th><th>Method</th><th>Status</th><th>Mean Abs Error</th><th>Max Inf</th></tr>";
if ~isempty(methodAudit)
    for i = 1:height(methodAudit)
        htmlLines(end+1) = "<tr><td>" + localHtml(methodAudit.paper_experiment(i)) + ...
            "</td><td>" + localHtml(methodAudit.method(i)) + "</td><td>" + ...
            localHtml(methodAudit.method_status(i)) + "</td><td>" + ...
            localHtml(localMdNum(methodAudit.mean_abs_error(i))) + "</td><td>" + ...
            localHtml(localMdNum(methodAudit.max_inf_rate(i))) + "</td></tr>";
    end
end
htmlLines(end+1) = "</table><h2>Use Boundary</h2><p>Non-formal modes validate the pipeline only. Do not insert positive claims into the manuscript until the formal scale has been run and re-audited.</p>";
htmlLines(end+1) = "</body></html>";
dpimnumeric.writeText(reportHtml, strjoin(htmlLines, newline));

report = struct("markdown", reportMd, "html", reportHtml);
end

function s = localMdNum(x)
if isnan(x)
    s = "--";
else
    s = sprintf("%.4g", x);
end
end

function s = localHtml(x)
s = string(x);
s = replace(s, "&", "&amp;");
s = replace(s, "<", "&lt;");
s = replace(s, ">", "&gt;");
s = replace(s, """", "&quot;");
end
