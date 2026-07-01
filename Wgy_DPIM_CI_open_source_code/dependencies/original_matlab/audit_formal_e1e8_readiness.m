function outputs = audit_formal_e1e8_readiness(resultsRoot, paperTex)
%AUDIT_FORMAL_E1E8_READINESS Audit whether a formal campaign matches the paper protocol.
%
% This function is downstream-only. It does not rerun experiments or edit
% the paper. It checks the current formal outputs against the E1--E8
% requirements stated in the core TeX and writes a reproducible gap matrix.

if nargin < 1 || strlength(string(resultsRoot)) == 0
    resultsRoot = localLatestFormalRoot();
end
if nargin < 2 || strlength(string(paperTex)) == 0
    paperTex = "C:\Users\Wangg\Desktop\期刊论文\DPIM_CI_full_integrated_weighted_RBn_natural.tex";
end

resultsRoot = char(string(resultsRoot));
paperTex = char(string(paperTex));
projectRoot = fileparts(mfilename("fullpath"));
outDir = fullfile(resultsRoot, "_readiness_audit");
dpimnumeric.ensureDir(outDir);

S = localLoadSources(resultsRoot);
config = jsondecode(fileread(fullfile(resultsRoot, "config_lock.json")));
manifest = jsondecode(fileread(fullfile(resultsRoot, "manifest.json")));
paperText = string(fileread(paperTex));

integrity = localIntegrityChecks(S, config, manifest, projectRoot);
seedAudit = localSeedAudit(resultsRoot);
requirements = localRequirementMatrix(S, config, paperText, projectRoot);
conflicts = localClaimConflicts(S, requirements);
inventory = localExperimentInventory(S);

writetable(integrity, fullfile(outDir, "data_integrity_checks.csv"));
writetable(seedAudit, fullfile(outDir, "seed_collision_audit.csv"));
writetable(requirements, fullfile(outDir, "experiment_requirement_matrix.csv"));
writetable(conflicts, fullfile(outDir, "claim_conflict_audit.csv"));
writetable(inventory, fullfile(outDir, "experiment_data_inventory.csv"));

acceptableStatuses = ["PASS", "BOUNDARY_ONLY", "SUPPLEMENT_ONLY"];
overallReady = all(ismember(requirements.paper_status, acceptableStatuses)) ...
    && ~any(integrity.status == "FAIL") ...
    && ~any(seedAudit.status == "FAIL") ...
    && ~any(conflicts.severity == "BLOCKER");
summaryPath = fullfile(outDir, "formal_readiness_audit.md");
localWriteSummary(summaryPath, resultsRoot, paperTex, overallReady, ...
    integrity, seedAudit, requirements, conflicts, inventory);

outputs = struct();
outputs.results_root = resultsRoot;
outputs.paper_tex = paperTex;
outputs.audit_root = outDir;
outputs.overall_ready = overallReady;
outputs.data_integrity_csv = fullfile(outDir, "data_integrity_checks.csv");
outputs.seed_audit_csv = fullfile(outDir, "seed_collision_audit.csv");
outputs.requirement_matrix_csv = fullfile(outDir, "experiment_requirement_matrix.csv");
outputs.claim_conflict_csv = fullfile(outDir, "claim_conflict_audit.csv");
outputs.inventory_csv = fullfile(outDir, "experiment_data_inventory.csv");
outputs.summary_md = summaryPath;
dpimnumeric.writeJson(fullfile(outDir, "audit_manifest.json"), outputs);

fprintf("Formal readiness audit completed. overall_ready=%d\n%s\n", overallReady, outDir);
end

function root = localLatestFormalRoot()
runs = dir(fullfile(pwd, "results", "weighted_paper_formal_formal_*"));
assert(~isempty(runs), "No weighted_paper_formal_formal_* directory found.");
[~, idx] = max([runs.datenum]);
root = fullfile(runs(idx).folder, runs(idx).name);
end

function S = localLoadSources(root)
S = struct();
S.preflight = localRead(fullfile(root, "preflight.csv"));
S.map = localRead(fullfile(root, "paper_experiment_map.csv"));
S.claims = localRead(fullfile(root, "claim_audit.csv"));
S.methods = localRead(fullfile(root, "method_audit.csv"));
S.e1 = localRead(fullfile(root, "evidence_suite", "E1_finite_B_rank_audit", "summary.csv"));
S.e2 = localRead(fullfile(root, "evidence_suite", "E2_scalar_normal_R_order_weighted", "confirmation_summary.csv"));
S.e2Candidates = localRead(fullfile(root, "evidence_suite", "E2_scalar_normal_R_order_weighted", "candidate_h.csv"));
S.e2Protocol = localRead(fullfile(root, "evidence_suite", "E2_scalar_normal_R_order_weighted", "protocol_summary.csv"));
S.e2Moments = localRead(fullfile(root, "evidence_suite", "E2_scalar_normal_R_order_weighted", "kernel_moment_validation.csv"));
S.e3 = localRead(fullfile(root, "evidence_suite", "E4_standard_normal_nonlinear_tail_weighted", "confirmation_summary.csv"));
S.e3Candidates = localRead(fullfile(root, "evidence_suite", "E4_standard_normal_nonlinear_tail_weighted", "candidate_h.csv"));
S.e3Protocol = localRead(fullfile(root, "evidence_suite", "E4_standard_normal_nonlinear_tail_weighted", "protocol_summary.csv"));
S.e4 = localRead(fullfile(root, "evidence_suite", "E5_standard_normal_bootstrap_t_instability_weighted", "confirmation_summary.csv"));
S.e4Mechanism = localRead(fullfile(root, "evidence_suite", "E5_standard_normal_bootstrap_t_instability_weighted", "mechanism_summary.csv"));
S.e5 = localRead(fullfile(root, "formal_experiments", "E5_probability_weighted_density", "summary.csv"));
S.e5Weights = localRead(fullfile(root, "formal_experiments", "E5_probability_weighted_density", "weight_diagnostics.csv"));
S.e6 = localRead(fullfile(root, "formal_experiments", "E6_rqmc_effective_order", "summary.csv"));
S.e6Order = localRead(fullfile(root, "formal_experiments", "E6_rqmc_effective_order", "effective_order_summary.csv"));
S.e7 = localRead(fullfile(root, "formal_experiments", "E7_finite_grid_band", "summary.csv"));
S.e7Pointwise = localRead(fullfile(root, "formal_experiments", "E7_finite_grid_band", "pointwise_summary.csv"));
S.e8 = localRead(fullfile(root, "formal_experiments", "E8_weighted_formula_holdout", "holdout_summary.csv"));
S.rbw = localRead(fullfile(root, "coverage_formula_validation", "coverage_formula_validation.csv"));
S.rbwWeights = localRead(fullfile(root, "coverage_formula_validation", "weight_moment_summary.csv"));
end

function T = localRead(path)
if isfile(path)
    T = readtable(path, TextType="string", Delimiter=",");
else
    T = table();
end
end

function T = localIntegrityChecks(S, config, manifest, projectRoot)
rows = cell(0, 1);
rows{end+1,1} = localCheck("preflight_all_pass", ...
    ~isempty(S.preflight) && all(S.preflight.status == "PASS"), ...
    "Formal preflight must have no WARN/FAIL rows."); %#ok<AGROW>
rows{end+1,1} = localCheck("formal_scale", ...
    string(config.run_mode) == "formal" && config.acceptance.formal_min_M <= 1000 ...
    && config.acceptance.formal_min_B <= 999 && manifest.M >= 1000 && manifest.B >= 999, ...
    sprintf("run_mode=%s, M=%g, B=%g", string(config.run_mode), manifest.M, manifest.B)); %#ok<AGROW>
rows{end+1,1} = localCheck("main_point_pool_is_rqmc", ...
    string(manifest.main_methods) == "sobol_scrambled", ...
    "Main point pools must use independently scrambled Sobol RQMC."); %#ok<AGROW>

rows{end+1,1} = localToleranceCheck("E1_C0B_identity", ...
    max(abs(double(S.e1.C0B) - (double(S.e1.k_plus) - double(S.e1.k_minus)) ./ (double(S.e1.B) + 1))), ...
    1e-14, "C0B=(k_plus-k_minus)/(B+1)."); %#ok<AGROW>
rows{end+1,1} = localToleranceCheck("E1_simulation_within_reported_MC_band", ...
    max(double(S.e1.simulation_error_vs_C0B) - double(S.e1.simulation_ci_half_width_95)), ...
    0, "Absolute simulation error must not exceed the reported 95% half-width."); %#ok<AGROW>

rows{end+1,1} = localCoverageMassCheck("E2_coverage_mass", S.e2, "left_miss", "right_miss", "interval_inf_rate"); %#ok<AGROW>
rows{end+1,1} = localCoverageMassCheck("E3_coverage_mass", S.e3, "left_miss", "right_miss", "interval_inf_rate"); %#ok<AGROW>
rows{end+1,1} = localCoverageMassCheck("E4_coverage_mass", S.e4, "left_miss", "right_miss", "interval_inf_rate"); %#ok<AGROW>
rows{end+1,1} = localCoverageMassCheck("E5_coverage_mass", S.e5, "left_miss_rate", "right_miss_rate", "interval_inf_rate"); %#ok<AGROW>

rows{end+1,1} = localToleranceCheck("E5_weight_sum", ...
    max(abs(double(S.e5Weights.sum_weights) - 1)), 1e-12, ...
    "All probability-weight vectors must sum to one."); %#ok<AGROW>
rows{end+1,1} = localCheck("E5_weights_nonnegative", ...
    min(double(S.e5Weights.min_w)) >= 0, "Probability weights must be nonnegative."); %#ok<AGROW>
rows{end+1,1} = localCheck("E5_no_point_or_weight_fallback", ...
    ~any(localTruthy(S.e5Weights.point_fallback_used)) ...
    && all(string(S.e5Weights.assignment_backend) == "gpu"), ...
    "Formal E5 must use scrambled Sobol points and GPU Voronoi assignment without fallback."); %#ok<AGROW>

suiteCode = string(fileread(fullfile(projectRoot, "run_weighted_paper_evidence_suite.m")));
hasIndependentE1 = contains(suiteCode, "localSimulateFiniteBOrderStatisticCoverage") ...
    && ~contains(suiteCode, "sim = rand(Msim, 1) <= C0B");
rows{end+1,1} = localCheck("E1_independent_order_statistic_simulation", ...
    hasIndependentE1, ...
    "E1 must independently simulate a true Uniform pivot and B bootstrap pivots."); %#ok<AGROW>
T = vertcat(rows{:});
end

function row = localCoverageMassCheck(name, T, leftName, rightName, infName)
vars = string(T.Properties.VariableNames);
if isempty(T) || ~all(ismember(["coverage", leftName, rightName, infName], vars))
    row = localStatusRow(name, "FAIL", "Required coverage decomposition columns are missing.", false);
    return;
end
residual = abs(double(T.coverage) + double(T.(leftName)) + double(T.(rightName)) ...
    + double(T.(infName)) - 1);
row = localToleranceCheck(name, max(residual), 5e-12, ...
    "coverage+left_miss+right_miss+Inf must equal one.");
end

function T = localSeedAudit(root)
experiments = ["E2_scalar_normal_R_order_weighted", ...
    "E4_standard_normal_nonlinear_tail_weighted", ...
    "E5_standard_normal_bootstrap_t_instability_weighted"];
rows = cell(numel(experiments), 1);
for i = 1:numel(experiments)
    base = fullfile(root, "evidence_suite", experiments(i));
    tuning = localRead(fullfile(base, "tuning_summary.csv"));
    confirm = localRead(fullfile(base, "confirmation_summary.csv"));
    tuningSeeds = unique(double(tuning.block_seed));
    confirmSeeds = unique(double(confirm.block_seed));
    expectedTuningBlocks = height(tuning) / numel(unique(tuning.method));
    collisionCount = expectedTuningBlocks - numel(tuningSeeds);
    overlap = intersect(tuningSeeds, confirmSeeds);
    status = "PASS";
    note = "Tuning and confirmation seeds are disjoint and unique.";
    if collisionCount > 0
        status = "FAIL";
        note = "Tuning seed formula collides across y/h blocks; this invalidates a clean parameter-lock audit.";
    elseif ~isempty(overlap)
        status = "FAIL";
        note = "Tuning and confirmation seed sets overlap.";
    end
    rows{i} = table(experiments(i), expectedTuningBlocks, numel(tuningSeeds), ...
        collisionCount, numel(confirmSeeds), numel(overlap), string(status), string(note), ...
        'VariableNames', {'experiment','expected_tuning_blocks','unique_tuning_seeds', ...
        'tuning_seed_collisions','unique_confirmation_seeds','tuning_confirmation_overlap', ...
        'status','note'});
end
T = vertcat(rows{:});
end

function T = localRequirementMatrix(S, config, paperText, projectRoot)
rows = cell(8, 1);

e1B = sort(unique(double(S.e1.B))).';
e1Formula = isequal(e1B, [99 199 399 999 1999]);
suiteCode = string(fileread(fullfile(projectRoot, "run_weighted_paper_evidence_suite.m")));
e1Independent = contains(suiteCode, "localSimulateFiniteBOrderStatisticCoverage") ...
    && ~contains(suiteCode, "sim = rand(Msim, 1) <= C0B");
rows{1} = localRequirement("E1", "TeX 3580--3596", ...
    "B grid, endpoint ranks, C0B, grid error, independent order-statistic simulation", ...
    e1Formula && e1Independent, "PARTIAL", ...
    "B grid and formula are present, but the simulation is Bernoulli(C0B), not an independent order-statistic program check.", ...
    "Replace E1 simulation with independent Uniform/order-statistic coverage and rerun."); %#ok<AGROW>

e2Source = S.e2Protocol;
if isempty(e2Source)
    e2Source = S.e2;
end
e2B = sort(unique(double(e2Source.B))).';
e2FixedH = localSameHInProtocol(e2Source);
momentVars = string(S.e2Moments.Properties.VariableNames);
e2MomentComplete = ~isempty(S.e2Moments) && all(ismember( ...
    ["A1_analytic","A1_numeric","A1_simulated","A4_analytic","A4_numeric","A4_simulated"], momentVars));
e2Complete = numel(e2B) >= 4 && e2FixedH && e2MomentComplete;
rows{2} = localRequirement("E2", "TeX 3598--3629", ...
    "A1--A4 analytic/numerical/simulation comparison; fixed-h R order; B={99,199,399,999}", ...
    e2Complete, "FAIL", ...
    sprintf("Formal E2 protocol B=%s; fixed-h=%d; three-way moment validation=%d.", ...
    mat2str(e2B), e2FixedH, e2MomentComplete), ...
    "Build a dedicated closed-form benchmark with three-way moment validation and factorial fixed-h/fixed-B paths."); %#ok<AGROW>

e3Source = S.e3Protocol;
if isempty(e3Source)
    e3Source = S.e3;
end
e3Fields = ["pool_skewness","pool_excess_kurtosis","A_p1_hat","A_bt1_hat"];
e3Has = all(ismember(e3Fields, string(e3Source.Properties.VariableNames)));
e3FixedH = localSameHInProtocol(e3Source);
rows{3} = localRequirement("E3", "TeX 3631--3655", ...
    "Skewed nonlinear response; gamma/kappa/effective order; Ap1/Abt1 sign relation", ...
    e3Has && e3FixedH, "FAIL", ...
    "The nonlinear response and high-precision smoothed truth exist, but confirmation output omits gamma/kappa and coefficient/sign diagnostics; h changes with R.", ...
    "Rerun with locked h paths, save pool cumulants, and test coverage-error sign against estimated first-order terms."); %#ok<AGROW>

e4Required = ["fallback_rate","minimum_positive_bootstrap_sd", ...
    "median_bt_to_percentile_length_ratio","p95_bt_to_percentile_length_ratio"];
e4Has = all(ismember(e4Required, string(S.e4Mechanism.Properties.VariableNames)));
rows{4} = localRequirement("E4", "TeX 3657--3675", ...
    "A_lambda, length-ratio quantiles, minimum bootstrap SD, raw BT versus hybrid", ...
    e4Has, "BOUNDARY_ONLY", ...
    "Required mechanism arrays exist, but the paper analysis must use the prespecified h subset and derive normalized min-SD diagnostics.", ...
    "Reprocess existing raw MAT blocks; retain E4 only as a failure-boundary result."); %#ok<AGROW>

e5Core = ~isempty(S.e5Weights) && all(string(S.e5Weights.point_actual_method) == "sobol_scrambled") ...
    && all(contains(lower(string(S.e5Weights.weight_source)), "voronoi"));
e5Fields = all(ismember(["gamma","kappa","bias_se_ratio"], string(S.e5.Properties.VariableNames)));
rows{5} = localRequirement("E5", "TeX 3677--3710", ...
    "Probability-weighted RQMC DPIM; n2/n3/n4 effective sizes; gamma/kappa; full diagnostic vector", ...
    e5Core && e5Fields, "PARTIAL", ...
    sprintf("Probability-weighted RQMC construction passes, but the summary omits gamma/kappa and bias/SE; minimum coverage is %.3f.", min(double(S.e5.coverage))), ...
    "Regenerate E5 with saved pool cumulants and bias/SE; keep local undercoverage explicit."); %#ok<AGROW>

e6Methods = unique(string(S.e6.point_method));
e6N = unique(double(S.e6.n));
e6Rqmc = S.e6(S.e6.point_method == "sobol_scrambled", :);
e6Complete = ~isempty(e6Rqmc) && ~isempty(S.e6Order) ...
    && all(double(e6Rqmc.point_fallback_rate) == 0) ...
    && all(ismember(["mean_rho3_w","mean_rho4_w","mean_ess_ratio"], ...
    string(S.e6.Properties.VariableNames)));
rows{6} = localRequirement("E6", "TeX 3713--3725", ...
    "RQMC construction, variance/effective-order and weight-moment diagnostics; no strong acceleration claim", ...
    e6Complete, "SUPPLEMENT_ONLY", ...
    sprintf("Current methods=%s and n=%s; retained only as construction/effective-order diagnostic.", ...
    strjoin(e6Methods, "/"), mat2str(e6N.')), ...
    "Keep E6 in supplement or limitations unless a stronger independently replicated efficiency study is added."); %#ok<AGROW>

e7G = unique(double(S.e7.grid_count));
e7Fields = all(ismember(["grid_spacing","max_stat_quantile","interpolation_remainder"], ...
    string(S.e7.Properties.VariableNames)));
rows{7} = localRequirement("E7", "TeX 3727--3739", ...
    "Pointwise and simultaneous coverage; multiple G; grid spacing; qmax; interpolation statement", ...
    numel(e7G) >= 3 && e7Fields, "FAIL", ...
    sprintf("Only G=%s is present; qmax, grid spacing and G-sensitivity are not saved.", mat2str(e7G.')), ...
    "Rerun finite-grid bands for several G and save critical values; do not claim continuous-curve coverage."); %#ok<AGROW>

e8R2p = localMethodValue(S.e8, "Percentile bootstrap", "r_squared");
e8R2bt = localMethodValue(S.e8, "Bootstrap-t", "r_squared");
e8Complete = ~isempty(S.e8) && ~isempty(S.rbw) ...
    && all(ismember(["mean_rho3_w","mean_rho4_w","coverage_error_to_C0B"], ...
    string(S.rbw.Properties.VariableNames)));
rows{8} = localRequirement("E8", "TeX 3779--3810", ...
    "Observed/predicted holdout error and rho3/rho4 structure diagnostic; no coefficient-proof claim", ...
    e8Complete, "SUPPLEMENT_ONLY", ...
    sprintf("Holdout R2 is %.3f for percentile and %.3f for bootstrap-t; random-weight cross-cumulants remain outside the fitted formula.", e8R2p, e8R2bt), ...
    "Keep E8 diagnostic-only and state that fitted coefficients and random-weight remainder terms are not theoretically proved."); %#ok<AGROW>

T = vertcat(rows{:});
T.paper_status(~T.complete & T.paper_status == "PASS") = "FAIL";
assert(contains(paperText, "\subsection{实验一：有限") ...
    && contains(paperText, "\subsection{实验八：概率权重"), ...
    "Core paper E1--E8 protocol markers were not found.");
end

function tf = localSameHAcrossR(candidates)
tf = false;
if isempty(candidates)
    return;
end
anchor = candidates(candidates.selection_reason == "prespecified_anchor", :);
if isempty(anchor)
    return;
end
yVals = unique(double(anchor.y0));
tf = true;
for i = 1:numel(yVals)
    h = unique(double(anchor.h(anchor.y0 == yVals(i))));
    tf = tf && numel(h) == 1;
end
end

function tf = localSameHInProtocol(T)
tf = false;
if isempty(T) || ~all(ismember(["y0","h","R"], string(T.Properties.VariableNames)))
    return;
end
yVals = unique(double(T.y0));
tf = true;
for i = 1:numel(yVals)
    h = unique(double(T.h(T.y0 == yVals(i))));
    tf = tf && numel(h) == 1;
end
end

function info = localLatestMomentEvidence(projectRoot)
dirs = dir(fullfile(projectRoot, "results", "formula_validation_v2_*"));
info = struct("path", "", "has_triple_comparison", false);
if isempty(dirs)
    return;
end
[~, idx] = max([dirs.datenum]);
path = fullfile(dirs(idx).folder, dirs(idx).name, "linear_gaussian_closed_kernel_moments.csv");
if ~isfile(path)
    return;
end
T = localRead(path);
vars = string(T.Properties.VariableNames);
info.path = path;
info.has_triple_comparison = all(ismember( ...
    ["A1_analytic","A1_numeric","A1_simulated","A4_analytic","A4_numeric","A4_simulated"], vars));
end

function T = localClaimConflicts(S, requirements)
rows = cell(0, 1);
for i = 1:height(requirements)
    expId = requirements.paper_experiment(i);
    claim = S.claims(S.claims.paper_experiment == expId, :);
    map = S.map(S.map.paper_experiment == expId, :);
    if isempty(claim)
        continue;
    end
    if claim.claim_status(1) == "main_text_ready" && requirements.paper_status(i) ~= "PASS"
        rows{end+1,1} = localConflict(expId, "BLOCKER", ...
            "claim_audit marks main_text_ready although the paper protocol is incomplete.", ...
            requirements.gap(i)); %#ok<AGROW>
    end
    if ~isempty(map) && map.role_in_paper(1) == "failure_or_diagnostic" ...
            && claim.claim_status(1) == "main_text_ready"
        rows{end+1,1} = localConflict(expId, "BLOCKER", ...
            "paper_experiment_map and claim_audit assign contradictory roles.", ...
            "Resolve the role before manuscript integration."); %#ok<AGROW>
    end
end

e5Claim = S.claims(S.claims.paper_experiment == "E5", :);
if ~isempty(e5Claim) && double(e5Claim.min_coverage(1)) < 0.90
    rows{end+1,1} = localConflict("E5", "MAJOR", ...
        sprintf("Minimum coverage %.3f is materially below the local acceptance threshold.", double(e5Claim.min_coverage(1))), ...
        "Do not summarize E5 only by mean coverage; report location/h/R cells and MCSE."); %#ok<AGROW>
end
if isempty(rows)
    T = table(strings(0,1), strings(0,1), strings(0,1), strings(0,1), ...
        'VariableNames', {'paper_experiment','severity','conflict','required_action'});
else
    T = vertcat(rows{:});
end
end

function T = localExperimentInventory(S)
ids = ["E1","E2","E3","E4","E5","E6","E7","E8"];
tables = {S.e1,S.e2,S.e3,S.e4,S.e5,S.e6,S.e7,S.rbw};
rows = cell(numel(ids), 1);
for i = 1:numel(ids)
    Q = tables{i};
    vars = string(Q.Properties.VariableNames);
    rows{i} = table(ids(i), height(Q), strjoin(vars, ","), ...
        localUniqueText(Q, "method"), localUniqueNumeric(Q, "R"), ...
        localUniqueNumeric(Q, "B"), localUniqueNumeric(Q, "h"), ...
        localUniqueNumeric(Q, "y0"), ...
        'VariableNames', {'paper_experiment','row_count','columns','methods', ...
        'R_values','B_values','h_values','y_values'});
end
T = vertcat(rows{:});
end

function value = localUniqueText(T, name)
if isempty(T) || ~ismember(name, string(T.Properties.VariableNames))
    value = "";
else
    value = strjoin(unique(string(T.(name)), "stable"), "/");
end
end

function value = localUniqueNumeric(T, name)
if isempty(T) || ~ismember(name, string(T.Properties.VariableNames))
    value = "";
else
    value = string(mat2str(sort(unique(double(T.(name)))).'));
end
end

function row = localRequirement(id, source, requirement, complete, status, gap, action)
row = table(string(id), string(source), string(requirement), logical(complete), ...
    string(status), string(gap), string(action), ...
    'VariableNames', {'paper_experiment','paper_source','requirement','complete', ...
    'paper_status','gap','required_action'});
end

function row = localConflict(id, severity, conflict, action)
row = table(string(id), string(severity), string(conflict), string(action), ...
    'VariableNames', {'paper_experiment','severity','conflict','required_action'});
end

function row = localCheck(name, ok, detail)
if ok
    status = "PASS";
else
    status = "FAIL";
end
row = localStatusRow(name, status, detail, ok);
end

function row = localToleranceCheck(name, value, tolerance, detail)
ok = isfinite(value) && value <= tolerance;
row = localStatusRow(name, localPassFail(ok), ...
    sprintf("%s Observed=%.6g, tolerance=%.6g.", detail, value, tolerance), ok);
end

function value = localMethodValue(T, method, variable)
idx = find(T.method == string(method), 1);
if isempty(idx) || ~ismember(variable, string(T.Properties.VariableNames))
    value = NaN;
else
    value = double(T.(variable)(idx));
end
end

function s = localPassFail(ok)
if ok
    s = "PASS";
else
    s = "FAIL";
end
end

function row = localStatusRow(name, status, detail, ok)
row = table(string(name), string(status), logical(ok), string(detail), ...
    'VariableNames', {'check','status','passed','detail'});
end

function tf = localTruthy(x)
if islogical(x)
    tf = x;
elseif isnumeric(x)
    tf = x ~= 0;
else
    tf = any(lower(string(x)) == ["1","true","yes"], 2);
end
end

function localWriteSummary(path, resultsRoot, paperTex, overallReady, ...
    integrity, seedAudit, requirements, conflicts, inventory)
lines = strings(0, 1);
lines(end+1) = "# Formal E1--E8 Readiness Audit";
lines(end+1) = "";
lines(end+1) = "- Results root: `" + string(resultsRoot) + "`";
lines(end+1) = "- Core paper: `" + string(paperTex) + "`";
lines(end+1) = "- Overall paper-ready: **" + upper(string(overallReady)) + "**";
lines(end+1) = "";
lines(end+1) = "## Explicit Judgment";
lines(end+1) = "";
if overallReady
    lines(end+1) = "The current campaign satisfies the paper protocol and may proceed to manuscript integration.";
else
    lines(end+1) = "The current campaign is not yet suitable for final manuscript integration. Existing outputs are useful pilot/formal-scale diagnostics, but protocol gaps and a tuning-seed collision require correction and selective reruns.";
end
lines(end+1) = "";
lines(end+1) = "## Integrity Checks";
for i = 1:height(integrity)
    lines(end+1) = sprintf("- `%s`: **%s**. %s", integrity.check(i), integrity.status(i), integrity.detail(i));
end
lines(end+1) = "";
lines(end+1) = "## Seed Audit";
for i = 1:height(seedAudit)
    lines(end+1) = sprintf("- `%s`: **%s**; tuning collisions=%d; tuning/confirmation overlap=%d. %s", ...
        seedAudit.experiment(i), seedAudit.status(i), seedAudit.tuning_seed_collisions(i), ...
        seedAudit.tuning_confirmation_overlap(i), seedAudit.note(i));
end
lines(end+1) = "";
lines(end+1) = "## E1--E8 Protocol Matrix";
for i = 1:height(requirements)
    lines(end+1) = sprintf("- `%s`: **%s**. Gap: %s Action: %s", ...
        requirements.paper_experiment(i), requirements.paper_status(i), ...
        requirements.gap(i), requirements.required_action(i));
end
lines(end+1) = "";
lines(end+1) = "## Claim Conflicts";
if isempty(conflicts)
    lines(end+1) = "- No conflicts.";
else
    for i = 1:height(conflicts)
        lines(end+1) = sprintf("- `%s` `%s`: %s Required: %s", ...
            conflicts.paper_experiment(i), conflicts.severity(i), ...
            conflicts.conflict(i), conflicts.required_action(i));
    end
end
lines(end+1) = "";
lines(end+1) = "## Data Inventory";
for i = 1:height(inventory)
    lines(end+1) = sprintf("- `%s`: rows=%d; R=%s; B=%s; h=%s; y=%s.", ...
        inventory.paper_experiment(i), inventory.row_count(i), inventory.R_values(i), ...
        inventory.B_values(i), inventory.h_values(i), inventory.y_values(i));
end
lines(end+1) = "";
lines(end+1) = "## Gate";
lines(end+1) = "";
lines(end+1) = "Do not insert the generated `formal_e1e8_results_section.tex` into the core paper until all BLOCKER conflicts are removed and E1/E2/E3/E5/E7 requirements are rerun or explicitly downgraded.";
dpimnumeric.writeText(path, strjoin(lines, newline));
end
