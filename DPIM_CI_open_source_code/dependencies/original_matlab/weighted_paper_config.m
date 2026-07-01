function cfg = weighted_paper_config(runMode, projectRoot, resultsRoot)
%WEIGHTED_PAPER_CONFIG Central configuration for the paper evidence suite.
%
% The suite is probability-weighted by default. Equal-weight rows are not
% produced by this config and must not be used as main paper evidence.

if nargin < 1 || strlength(string(runMode)) == 0
    runMode = "diagnostic";
end
if nargin < 2 || strlength(string(projectRoot)) == 0
    projectRoot = fileparts(mfilename("fullpath"));
end
if nargin < 3 || strlength(string(resultsRoot)) == 0
    stamp = string(datetime("now", "Format", "yyyyMMdd_HHmmss"));
    resultsRoot = fullfile(projectRoot, "results", ...
        "weighted_paper_evidence_" + lower(string(runMode)) + "_" + stamp);
end

runMode = lower(string(runMode));
if ~any(runMode == ["diagnostic", "small", "pilot", "medium", "full", "formal"])
    error("runMode must be diagnostic, small, pilot, medium, full, or formal.");
end

desktopRoot = fileparts(projectRoot);
paperTex = localFindPaperTex(desktopRoot);
thesisPdf = localFindThesisPdf(desktopRoot);
if strlength(string(thesisPdf)) > 0
    referenceDir = fileparts(fileparts(thesisPdf));
else
    referenceDir = "";
end

cfg = struct();
cfg.schema_version = "weighted_paper_evidence_suite_v1";
cfg.run_mode = char(runMode);
cfg.project_root = char(projectRoot);
cfg.results_root = char(resultsRoot);
cfg.paper_tex = char(paperTex);
cfg.reference_dir = char(referenceDir);
cfg.thesis_pdf = char(thesisPdf);
cfg.source_scope_note = "Core paper is restricted to the desktop journal-paper project; thesis PDF is restricted to the user-provided 2_*20221091006*.pdf outside the DPIM bundle.";
cfg.weighted_root = fullfile(projectRoot, ...
    "dpim_h_weighted_rebuild_package_fixed4", "dpim_h_weighted_rebuild");
cfg.seed = 2026061601;
cfg.alpha = 0.05;
cfg.nominal = 0.95;
cfg.lambda = 5;
cfg.main_methods = "sobol_scrambled";
cfg.control_methods = "mc";
cfg.include_control_methods = false;
cfg.methods = cfg.main_methods;
cfg.assignment_backend = "auto";
cfg.use_parallel = false;
cfg.run_original_plate_summary = true;
cfg.original_plate_glob = "original_plate_full_*";
cfg.selection_max_fallback = 0.05;
cfg.selection_max_inf = 0.01;
cfg.selection_max_mean_abs_error = 0.05;
cfg.selection_min_location_coverage = 0.90;
cfg.selection_h_anchor_penalty = 0.005;
cfg.max_reasonable_mean_length = 1e8;
cfg.protocol_B_list = [99, 199, 399, 999];
cfg.e7_grid_counts = [3, 5, 9];

switch runMode
    case "formal"
        cfg.n = 768;
        cfg.curve_pool_size = 1200;
        cfg.tuning_pool_size = 600;
        cfg.R_list = [16, 32, 64, 128];
        cfg.M = 1000;
        cfg.B = 999;
        cfg.tuning_M = 100;
        cfg.truth_N = 1000000;
        cfg.voronoi_aux_sample_count = 300000;
        cfg.voronoi_block_size = 4000;
        cfg.h_list = localUniqueH(logspace(-7, 0, 161));
    case "full"
        cfg.n = 800;
        cfg.curve_pool_size = 600;
        cfg.tuning_pool_size = 300;
        cfg.R_list = [20, 40, 80];
        cfg.M = 300;
        cfg.B = 399;
        cfg.tuning_M = 120;
        cfg.truth_N = 300000;
        cfg.voronoi_aux_sample_count = 200000;
        cfg.voronoi_block_size = 2000;
        cfg.h_list = localUniqueH(logspace(-6, 0, 121));
    case "medium"
        cfg.n = 384;
        cfg.curve_pool_size = 240;
        cfg.tuning_pool_size = 120;
        cfg.R_list = [20, 40];
        cfg.M = 120;
        cfg.B = 199;
        cfg.truth_N = 100000;
        cfg.voronoi_aux_sample_count = 50000;
        cfg.voronoi_block_size = 1000;
        cfg.h_list = localUniqueH(logspace(-6, -0.5, 61));
    case "pilot"
        cfg.n = 256;
        cfg.curve_pool_size = 180;
        cfg.tuning_pool_size = 90;
        cfg.R_list = [12, 24, 48];
        cfg.M = 100;
        cfg.B = 199;
        cfg.truth_N = 60000;
        cfg.voronoi_aux_sample_count = 25000;
        cfg.voronoi_block_size = 1000;
        cfg.h_list = localUniqueH(logspace(-7, -0.2, 61));
    case "small"
        cfg.n = 192;
        cfg.curve_pool_size = 120;
        cfg.tuning_pool_size = 60;
        cfg.R_list = [12, 20];
        cfg.M = 50;
        cfg.B = 99;
        cfg.truth_N = 30000;
        cfg.voronoi_aux_sample_count = 12000;
        cfg.voronoi_block_size = 1000;
        cfg.h_list = localUniqueH(logspace(-6, -0.3, 41));
    otherwise
        cfg.n = 48;
        cfg.curve_pool_size = 20;
        cfg.tuning_pool_size = 10;
        cfg.R_list = 8;
        cfg.M = 8;
        cfg.B = 39;
        cfg.truth_N = 3000;
        cfg.voronoi_aux_sample_count = 1200;
        cfg.voronoi_block_size = 400;
        cfg.methods = "sobol_scrambled";
        cfg.h_list = localUniqueH([logspace(-5, -0.5, 8), ...
            4e-4, 0.00183298071083244, 0.233572146909012]);
end
if ~isfield(cfg, "tuning_M")
    cfg.tuning_M = cfg.M;
end
cfg.tuning_B = cfg.B;

cfg.confirmation_h_per_experiment = 3;
cfg.required_functions = ["dpim_build_weighted_curve_pool", ...
    "dpim_curve_point_estimates", "dpim_gaussian_kernel", ...
    "dpim_truth_smoothed_density", "voronoi_ci_probability_weights_provider"];
cfg.weighting_cfg = struct();
cfg.weighting_cfg.provider = "voronoi_ci_probability_weights_provider";
cfg.weighting_cfg.voronoi_aux_sample_count = cfg.voronoi_aux_sample_count;
cfg.weighting_cfg.voronoi_block_size = cfg.voronoi_block_size;
cfg.weighting_cfg.voronoi_enable_cache = true;
cfg.weighting_cfg.voronoi_save_outputs = false;
cfg.weighting_cfg.voronoi_assignment_backend = cfg.assignment_backend;
end

function h = localUniqueH(values)
anchor = [4e-4, 0.00183298071083244, 0.233572146909012];
h = unique([values(:); anchor(:)]);
h = h(isfinite(h) & h > 0);
h = sort(h(:)).';
end

function path = localFindPaperTex(rootDir)
hits = dir(fullfile(rootDir, "**", "DPIM_CI_full_integrated_weighted_RBn_natural.tex"));
path = localBestPath(hits, @(h) localPaperScore(h), 30);
end

function score = localPaperScore(hit)
folder = string(hit.folder);
score = 0;
if isfolder(fullfile(hit.folder, "DPIM_CI_full_integrated_figures"))
    score = score + 20;
end
if isfolder(fullfile(hit.folder, ".latex-build"))
    score = score + 10;
end
if contains(folder, "DPIM_CI_final_no_wrong_formula_bundle")
    score = score - 20;
end
end

function path = localFindThesisPdf(rootDir)
hits = dir(fullfile(rootDir, "**", "2 *20221091006*.pdf"));
if ~isempty(hits)
    keep = true(numel(hits), 1);
    for i = 1:numel(hits)
        keep(i) = ~contains(string(hits(i).folder), "DPIM_CI_final_no_wrong_formula_bundle");
    end
    hits = hits(keep);
end
path = localBestPath(hits, @(h) localThesisScore(h), 30);
end

function score = localThesisScore(hit)
folder = string(hit.folder);
score = 0;
if startsWith(string(hit.name), "2 ")
    score = score + 10;
end
if contains(folder, "20221091006")
    score = score + 20;
end
if contains(folder, "DPIM_CI_final_no_wrong_formula_bundle")
    score = score - 20;
end
end

function path = localBestPath(hits, scoreFcn, minScore)
if nargin < 3
    minScore = -Inf;
end
path = "";
if isempty(hits)
    return;
end
hits = hits(~[hits.isdir]);
if isempty(hits)
    return;
end
scores = zeros(numel(hits), 1);
for i = 1:numel(hits)
    scores(i) = scoreFcn(hits(i));
end
[bestScore, idx] = max(scores);
if bestScore < minScore
    return;
end
same = find(scores == bestScore);
if numel(same) > 1
    [~, loc] = max([hits(same).datenum]);
    idx = same(loc);
end
path = string(fullfile(hits(idx).folder, hits(idx).name));
end
