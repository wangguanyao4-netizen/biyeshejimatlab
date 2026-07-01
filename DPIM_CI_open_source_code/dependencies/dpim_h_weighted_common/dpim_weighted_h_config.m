function cfg = dpim_weighted_h_config(runMode, projectRoot)
%DPIM_WEIGHTED_H_CONFIG Configuration for weighted DPIM h studies.
if nargin < 1; runMode = "small"; end
if nargin < 2 || isempty(projectRoot); projectRoot = pwd; end
runMode = string(runMode);

cfg = struct();
cfg.run_mode = runMode;
cfg.project_root = char(projectRoot);
cfg.results_root = fullfile(char(projectRoot), 'results', ['weighted_dpim_h_' char(runMode) '_' char(datetime('now','Format','yyyyMMdd_HHmmss'))]);
cfg.seed = 20260608;
cfg.alpha = 0.05;
cfg.nominal = 0.95;
cfg.coverage_tol = 0.03;
cfg.inf_tol = 0.01;
cfg.min_local_neff = 5;
cfg.draw_with_replacement = true;

switch lower(char(runMode))
    case 'full'
        cfg.n_inner = 256;
        cfg.num_curves = 600;
        cfg.R = 96;
        cfg.M = 300;
        cfg.B = 799;
        cfg.truth_N = 300000;
        cfg.aux_sample_count = 12000;
        cfg.h_list = [0.005 0.01 0.02 0.05 0.1 0.2 0.5 1.0 2.0];
        cfg.methods = ["mc", "sobol_scrambled"];
        cfg.audit_methods = ["mc", "sobol_qmc", "sobol_scrambled"];
    case 'medium'
        cfg.n_inner = 128;
        cfg.num_curves = 250;
        cfg.R = 64;
        cfg.M = 120;
        cfg.B = 399;
        cfg.truth_N = 120000;
        cfg.aux_sample_count = 6000;
        cfg.h_list = [0.01 0.02 0.05 0.1 0.2 0.5 1.0 2.0];
        cfg.methods = ["mc", "sobol_scrambled"];
        cfg.audit_methods = ["mc", "sobol_qmc", "sobol_scrambled"];
    otherwise
        cfg.n_inner = 64;
        cfg.num_curves = 80;
        cfg.R = 32;
        cfg.M = 40;
        cfg.B = 199;
        cfg.truth_N = 30000;
        cfg.aux_sample_count = 2000;
        cfg.h_list = [0.02 0.05 0.1 0.2 0.5 1.0];
        cfg.methods = ["mc", "sobol_scrambled"];
        cfg.audit_methods = ["mc", "sobol_qmc", "sobol_scrambled"];
end

cfg.weighting_cfg = struct();
cfg.weighting_cfg.voronoi_aux_sample_count = cfg.aux_sample_count;
cfg.weighting_cfg.voronoi_block_size = min(1000, cfg.aux_sample_count);
cfg.weighting_cfg.voronoi_enable_cache = true;
cfg.weighting_cfg.voronoi_save_outputs = false;
cfg.weighting_cfg.voronoi_output_dir = fullfile(cfg.results_root, 'ci_probability_weights');

cfg.problems = dpim_weighted_define_problems(runMode);
end
