clc;
clear;
close all;

project_root = fileparts(mfilename('fullpath'));
cd(project_root);

sample_pool_count = 200;
sample_count = 200;
confidence_level = 0.95;
sigma = 1;
base_seed = 1;
ci_sample_size = 10;
coverage_trials = 400;
bootstrap_repetitions = 1600;
bootstrap_seed = 1;
expected_w_grid_count = 800;

reference_files = dir('mc_reference_truth_dual_plate_n*.mat');
if isempty(reference_files)
    error('Missing dual-plate MC reference truth file. Run run_dual_plate_mc_reference_truth_fullflow.m first.');
end
reference_names = {reference_files.name};
[~, newest_idx] = max([reference_files.datenum]);
reference_file = reference_names{newest_idx};

ref_data = load(reference_file);
w_list = ref_data.w_list(:).';
if numel(w_list) ~= expected_w_grid_count
    error('Reference grid has %d points, expected %d. Rerun run_dual_plate_mc_reference_truth_fullflow.m.', numel(w_list), expected_w_grid_count);
end
if ~isfield(ref_data, 'thin') || ~isfield(ref_data.thin, 'p_mc')
    error('Reference file %s is missing thin.p_mc. Regenerate it with run_dual_plate_mc_reference_truth_fullflow.m.', reference_file);
end
if ~isfield(ref_data, 'mindlin') || ~isfield(ref_data.mindlin, 'p_mc')
    error('Reference file %s is missing mindlin.p_mc. Regenerate it with run_dual_plate_mc_reference_truth_fullflow.m.', reference_file);
end
reference_curve_thin = ref_data.thin.p_mc(:).';
reference_curve_mindlin = ref_data.mindlin.p_mc(:).';

rqmc_driver_state = struct();
rqmc_driver_state.config = struct();
rqmc_driver_state.config.sample_count = sample_count;
rqmc_driver_state.config.confidence_level = confidence_level;
rqmc_driver_state.config.sigma = sigma;
rqmc_driver_state.config.w_list = w_list;
rqmc_driver_state.config.base_seed = base_seed;
rqmc_driver_state.config.sample_pool_count = sample_pool_count;
rqmc_driver_state.config.ci_sample_size = ci_sample_size;
rqmc_driver_state.config.coverage_trials = coverage_trials;
rqmc_driver_state.config.bootstrap_repetitions = bootstrap_repetitions;
rqmc_driver_state.config.bootstrap_seed = bootstrap_seed;
rqmc_driver_state.config.voronoi_aux_sample_count = 50000;
rqmc_driver_state.config.voronoi_train_aux_count = 5000;
rqmc_driver_state.config.voronoi_exact_block_size = 1000;
rqmc_driver_state.config.save_figures = false;
rqmc_driver_state.config.sample_fig_on = false;
rqmc_driver_state.config.draw_figures = false;
rqmc_driver_state.config.save_png_figures = false;
rqmc_driver_state.config.plot_gf_figures = false;
rqmc_driver_state.outer_repetitions = sample_pool_count;
rqmc_driver_state.current_rep = 1;
rqmc_driver_state.seed_list = base_seed + (0:(sample_pool_count - 1));
rqmc_driver_state.curve_pool_thin = zeros(sample_pool_count, numel(w_list));
rqmc_driver_state.curve_pool_mindlin = zeros(sample_pool_count, numel(w_list));
rqmc_driver_state.wc_runs_thin = cell(sample_pool_count, 1);
rqmc_driver_state.wc_runs_mindlin = cell(sample_pool_count, 1);
rqmc_driver_state.weight_runs = cell(sample_pool_count, 1);
rqmc_driver_state.cumulative_runs_thin = cell(sample_pool_count, 1);
rqmc_driver_state.cumulative_runs_mindlin = cell(sample_pool_count, 1);
rqmc_driver_state.reference_file = reference_file;
rqmc_driver_state.reference_curve_thin = reference_curve_thin;
rqmc_driver_state.reference_curve_mindlin = reference_curve_mindlin;

ensure_eval_centers_rqmc_file();

while rqmc_driver_state.current_rep <= rqmc_driver_state.outer_repetitions
    fprintf('Running dual-plate RQMC outer sample %d / %d\n', ...
        rqmc_driver_state.current_rep, rqmc_driver_state.outer_repetitions);
    fprintf('  reference file: %s\n', rqmc_driver_state.reference_file);
    run('run_single_dual_plate_rqmc_outer_sample.m');
end

if exist('rqmc_pipeline_config.mat', 'file') == 2
    delete('rqmc_pipeline_config.mat');
end

sample_pool_count = rqmc_driver_state.config.sample_pool_count;
sample_count = rqmc_driver_state.config.sample_count;
confidence_level = rqmc_driver_state.config.confidence_level;
sigma = rqmc_driver_state.config.sigma;
base_seed = rqmc_driver_state.config.base_seed;
ci_sample_size = rqmc_driver_state.config.ci_sample_size;
coverage_trials = rqmc_driver_state.config.coverage_trials;
bootstrap_repetitions = rqmc_driver_state.config.bootstrap_repetitions;
bootstrap_seed = rqmc_driver_state.config.bootstrap_seed;
w_list = rqmc_driver_state.config.w_list;
reference_file = rqmc_driver_state.reference_file;

thin_result = local_build_plate_result( ...
    'thin', ...
    rqmc_driver_state.curve_pool_thin, ...
    rqmc_driver_state.reference_curve_thin, ...
    rqmc_driver_state, ...
    alpha_from_confidence(confidence_level), ...
    'rqmc_dual_plate_thin_curve_summary.csv', ...
    'rqmc_dual_plate_thin_coverage_summary.csv');

mindlin_result = local_build_plate_result( ...
    'mindlin', ...
    rqmc_driver_state.curve_pool_mindlin, ...
    rqmc_driver_state.reference_curve_mindlin, ...
    rqmc_driver_state, ...
    alpha_from_confidence(confidence_level), ...
    'rqmc_dual_plate_mindlin_curve_summary.csv', ...
    'rqmc_dual_plate_mindlin_coverage_summary.csv');

dual_result = struct();
dual_result.sample_pool_count = sample_pool_count;
dual_result.sample_count = sample_count;
dual_result.confidence_level = confidence_level;
dual_result.sigma = sigma;
dual_result.base_seed = base_seed;
dual_result.ci_sample_size = ci_sample_size;
dual_result.coverage_trials = coverage_trials;
dual_result.bootstrap_repetitions = bootstrap_repetitions;
dual_result.bootstrap_seed = bootstrap_seed;
dual_result.w_list = w_list;
dual_result.seed_list = rqmc_driver_state.seed_list;
dual_result.weight_runs = rqmc_driver_state.weight_runs;
dual_result.reference_file = reference_file;
dual_result.thin = thin_result;
dual_result.mindlin = mindlin_result;

save('rqmc_dual_plate_sample_pool_ci_results.mat', '-struct', 'dual_result');

fprintf('Saved dual-plate sample-pool CI result to rqmc_dual_plate_sample_pool_ci_results.mat\n');

function alpha = alpha_from_confidence(confidence_level)
alpha = 1 - confidence_level;
end

function ensure_eval_centers_rqmc_file()
    a = 1.0;
    b = 1.0;
    alpha_deg = 90;
    nx = 20;
    ny = 20;

    alpha = deg2rad(alpha_deg);
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

    writematrix(nodes, 'eval_centers_rqmc.txt', 'Delimiter', 'tab');
end

function plate_result = local_build_plate_result( ...
    plate_tag, curve_pool, reference_curve, state, alpha, curve_summary_file, coverage_summary_file)

w_list = state.config.w_list;
sample_pool_count = state.config.sample_pool_count;
ci_sample_size = state.config.ci_sample_size;
coverage_trials = state.config.coverage_trials;
bootstrap_repetitions = state.config.bootstrap_repetitions;
bootstrap_seed = state.config.bootstrap_seed;
base_seed = state.config.base_seed;

mean_curve = mean(curve_pool, 1);
std_curve = std(curve_pool, 0, 1);
se_curve = std_curve / sqrt(sample_pool_count);

pool_ci_results = local_ci_methods(curve_pool, alpha, bootstrap_repetitions, bootstrap_seed, reference_curve);
coverage_result = local_evaluate_ci_trials( ...
    curve_pool, alpha, bootstrap_repetitions, bootstrap_seed, reference_curve, ...
    ci_sample_size, coverage_trials, base_seed);

plate_result = struct();
plate_result.plate_tag = plate_tag;
plate_result.reference_curve = reference_curve;
plate_result.curve_pool = curve_pool;
plate_result.mean_curve = mean_curve;
plate_result.std_curve = std_curve;
plate_result.se_curve = se_curve;
plate_result.pool_ci_results = pool_ci_results;
plate_result.coverage_result = coverage_result;
if strcmp(plate_tag, 'thin')
    plate_result.wc_runs = state.wc_runs_thin;
    plate_result.cumulative_runs = state.cumulative_runs_thin;
else
    plate_result.wc_runs = state.wc_runs_mindlin;
    plate_result.cumulative_runs = state.cumulative_runs_mindlin;
end

summary_table = table(w_list(:), mean_curve(:), std_curve(:), se_curve(:), ...
    reference_curve(:), ...
    pool_ci_results(1).lower(:), pool_ci_results(1).upper(:), ...
    pool_ci_results(2).lower(:), pool_ci_results(2).upper(:), ...
    pool_ci_results(3).lower(:), pool_ci_results(3).upper(:), ...
    'VariableNames', {'w', 'mean_curve', 'std_curve', 'se_curve', 'reference_curve', ...
    'student_t_lower', 'student_t_upper', ...
    'percentile_bootstrap_lower', 'percentile_bootstrap_upper', ...
    'bootstrap_t_lower', 'bootstrap_t_upper'});
writetable(summary_table, curve_summary_file);

coverage_summary_table = table(w_list(:), reference_curve(:), ...
    coverage_result.coverage(1, :).', coverage_result.mean_length(1, :).', coverage_result.median_length(1, :).', ...
    coverage_result.coverage(2, :).', coverage_result.mean_length(2, :).', coverage_result.median_length(2, :).', ...
    coverage_result.coverage(3, :).', coverage_result.mean_length(3, :).', coverage_result.median_length(3, :).', ...
    'VariableNames', {'w', 'reference_curve', ...
    'student_t_coverage', 'student_t_mean_length', 'student_t_median_length', ...
    'percentile_bootstrap_coverage', 'percentile_bootstrap_mean_length', 'percentile_bootstrap_median_length', ...
    'bootstrap_t_coverage', 'bootstrap_t_mean_length', 'bootstrap_t_median_length'});
writetable(coverage_summary_table, coverage_summary_file);

fprintf('%s average pointwise coverage across %d trials:\n', upper(plate_tag), coverage_trials);
for idx = 1:numel(coverage_result.method_names)
    fprintf('  %s = %.6f\n', ...
        coverage_result.method_names{idx}, mean(coverage_result.coverage(idx, :)));
end
end

function ci_results = local_ci_methods(Y, alpha, B, seed, true_values)
if ~isempty(seed)
    previous_state = rng;
    rng(double(seed), 'twister');
    cleanup_obj = onCleanup(@() rng(previous_state)); %#ok<NASGU>
end

Y = double(Y);
[R, num_grid] = size(Y);
truth_row = double(true_values(:)).';
if numel(truth_row) ~= num_grid
    error('reference_curve length must match the number of grid points.');
end

Ybar = mean(Y, 1);
SR = std(Y, 0, 1);
sqrtR = sqrt(R);

t_crit = local_student_t_quantile(1 - alpha / 2, R - 1);
student_lower = Ybar - t_crit .* SR ./ sqrtR;
student_upper = Ybar + t_crit .* SR ./ sqrtR;

boot_indices = randi(R, R, B);
Ystar_flat = Y(boot_indices(:), :);
Ystar3 = reshape(Ystar_flat, R, B, num_grid);
boot_means = reshape(mean(Ystar3, 1), B, num_grid);
[pct_lower, pct_upper] = local_two_sided_quantiles_matrix(boot_means, alpha);

boot_std = reshape(std(Ystar3, 0, 1), B, num_grid);
Tstar = sqrtR .* (boot_means - Ybar) ./ boot_std;
zero_std_mask = (boot_std == 0);
zero_std_count = sum(zero_std_mask, 1);
if any(zero_std_mask, 'all')
    delta = boot_means - Ybar;
    sign_delta = sign(delta);
    sign_delta(sign_delta == 0) = 1;
    Tstar(zero_std_mask) = sign_delta(zero_std_mask) .* Inf;
end
[t_lower, t_upper] = local_two_sided_quantiles_matrix(Tstar, alpha);
boot_t_lower = Ybar;
boot_t_upper = Ybar;
positive_std_mask = (SR > 0);
boot_t_lower(positive_std_mask) = Ybar(positive_std_mask) - SR(positive_std_mask) .* t_upper(positive_std_mask) ./ sqrtR;
boot_t_upper(positive_std_mask) = Ybar(positive_std_mask) - SR(positive_std_mask) .* t_lower(positive_std_mask) ./ sqrtR;

ci_results = repmat(local_empty_ci(num_grid), 3, 1);
ci_results(1) = local_pack_ci('Student-t', student_lower, student_upper, zeros(1, num_grid), truth_row);
ci_results(2) = local_pack_ci('Percentile bootstrap', pct_lower, pct_upper, zeros(1, num_grid), truth_row);
ci_results(3) = local_pack_ci('Bootstrap-t', boot_t_lower, boot_t_upper, zero_std_count, truth_row);
end

function coverage_result = local_evaluate_ci_trials( ...
    curve_pool, alpha, B, bootstrap_seed, true_values, sample_size, num_trials, base_seed)

pool_size = size(curve_pool, 1);
num_grid = size(curve_pool, 2);
if sample_size > pool_size
    error('ci_sample_size (%d) cannot exceed sample_pool_count (%d).', sample_size, pool_size);
end

method_names = {'Student-t', 'Percentile bootstrap', 'Bootstrap-t'};
num_methods = numel(method_names);
sample_indices = zeros(num_trials, sample_size);
lower = zeros(num_trials, num_methods, num_grid);
upper = zeros(num_trials, num_methods, num_grid);
interval_length = zeros(num_trials, num_methods, num_grid);
contains_true = false(num_trials, num_methods, num_grid);
zero_std_count = zeros(num_trials, num_methods, num_grid);
infinite_flag = false(num_trials, num_methods, num_grid);

progress_stride = max(1, floor(num_trials / 5));
for trial_idx = 1:num_trials
    sample_seed = base_seed + 1000000 + trial_idx;
    ci_seed = bootstrap_seed + 2000000 + trial_idx;
    trial_indices = local_sample_without_replacement(pool_size, sample_size, sample_seed);
    sample_indices(trial_idx, :) = trial_indices;
    trial_ci_results = local_ci_methods(curve_pool(trial_indices, :), alpha, B, ci_seed, true_values);

    for method_idx = 1:num_methods
        lower(trial_idx, method_idx, :) = reshape(trial_ci_results(method_idx).lower, 1, 1, num_grid);
        upper(trial_idx, method_idx, :) = reshape(trial_ci_results(method_idx).upper, 1, 1, num_grid);
        interval_length(trial_idx, method_idx, :) = reshape(trial_ci_results(method_idx).length, 1, 1, num_grid);
        contains_true(trial_idx, method_idx, :) = reshape(trial_ci_results(method_idx).contains_true, 1, 1, num_grid);
        zero_std_count(trial_idx, method_idx, :) = reshape(trial_ci_results(method_idx).zero_std_count, 1, 1, num_grid);
        infinite_flag(trial_idx, method_idx, :) = reshape(trial_ci_results(method_idx).infinite_flag, 1, 1, num_grid);
    end

    if mod(trial_idx, progress_stride) == 0 || trial_idx == num_trials
        fprintf('CI coverage trials: %d / %d\n', trial_idx, num_trials);
    end
end

coverage_result = struct();
coverage_result.method_names = method_names;
coverage_result.sample_indices = sample_indices;
coverage_result.lower = lower;
coverage_result.upper = upper;
coverage_result.length = interval_length;
coverage_result.contains_true = contains_true;
coverage_result.zero_std_count = zero_std_count;
coverage_result.infinite_flag = infinite_flag;
coverage_result.coverage = squeeze(mean(contains_true, 1));
coverage_result.mean_length = squeeze(mean(interval_length, 1));
coverage_result.median_length = squeeze(median(interval_length, 1));
end

function indices = local_sample_without_replacement(pool_size, sample_size, seed)
previous_state = rng;
rng(double(seed), 'twister');
cleanup_obj = onCleanup(@() rng(previous_state)); %#ok<NASGU>
indices = randperm(pool_size, sample_size);
end

function ci = local_empty_ci(num_grid)
ci = struct( ...
    'CI_method', '', ...
    'lower', nan(1, num_grid), ...
    'upper', nan(1, num_grid), ...
    'length', inf(1, num_grid), ...
    'contains_true', false(1, num_grid), ...
    'zero_std_count', zeros(1, num_grid), ...
    'infinite_flag', false(1, num_grid));
end

function ci = local_pack_ci(name, lower, upper, zero_std_count, true_values)
swap_mask = lower > upper;
tmp = lower(swap_mask);
lower(swap_mask) = upper(swap_mask);
upper(swap_mask) = tmp;

ci = local_empty_ci(numel(lower));
ci.CI_method = name;
ci.lower = lower;
ci.upper = upper;
ci.zero_std_count = double(zero_std_count);
ci.infinite_flag = ~(isfinite(lower) & isfinite(upper));
ci.length = upper - lower;
ci.length(ci.infinite_flag) = Inf;
valid_mask = ~(isnan(lower) | isnan(upper));
ci.contains_true(valid_mask) = ...
    (lower(valid_mask) <= true_values(valid_mask)) & ...
    (true_values(valid_mask) <= upper(valid_mask));
end

function [q_lower, q_upper] = local_two_sided_quantiles_matrix(x, alpha)
x_sorted = sort(x, 1);
n = size(x_sorted, 1);
lower_index = max(1, floor(n * alpha / 2));
upper_index = min(n, ceil(n * (1 - alpha / 2)));
q_lower = x_sorted(lower_index, :);
q_upper = x_sorted(upper_index, :);
end

function q = local_student_t_quantile(p, nu)
if p == 0.5
    q = 0;
    return;
end

tail_prob = 2 * min(p, 1 - p);
x = betaincinv(tail_prob, nu / 2, 0.5);
q = sign(p - 0.5) * sqrt(nu * (1 / x - 1));
end
