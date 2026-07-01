clc;
clear;
close all;

%% =========================================================
% Dual-plate MC reference-truth full-flow configuration
% =========================================================
sample_count = 10000;
sigma = 1;
w_list = linspace(-10, 1000, 800);

project_root = fileparts(mfilename('fullpath'));
cd(project_root);

mc_pipeline_config = struct();
mc_pipeline_config.sample_count = sample_count;
mc_pipeline_config.sigma = sigma;
mc_pipeline_config.w_list = w_list;
mc_pipeline_config.voronoi_aux_sample_count = 200000;
mc_pipeline_config.voronoi_train_aux_count = 20000;
mc_pipeline_config.voronoi_exact_block_size = 1000;
mc_pipeline_config.save_figures = false;
mc_pipeline_config.sample_fig_on = false;
mc_pipeline_config.draw_figures = false;
mc_pipeline_config.save_png_figures = false;
mc_pipeline_config.plot_gf_figures = false;
save('mc_pipeline_config.mat', 'mc_pipeline_config');

ensure_eval_centers_mc_file();

fprintf('Shared stage 1 / 5: gf_normal_final_mc\n');
gf_normal_final_mc;
fprintf('Shared stage 2 / 5: voronoi_pca_toolbox_final_blockwise_mc\n');
voronoi_pca_toolbox_final_blockwise_mc;
fprintf('Shared stage 3 / 5: KL2D_final_complete_3D_mc\n');
KL2D_final_complete_3D_mc;

fprintf('Plate stage 4a / 5: thin plate wc calculation\n');
jisuanpro_mc_quad_paperstyle_cloud;
load('wc_all_mc.mat', 'wc_all_mc');
wc_all_mc_thin = wc_all_mc;
save('wc_all_mc_thin.mat', 'wc_all_mc_thin');

fprintf('Plate stage 4b / 5: Mindlin plate wc calculation\n');
jisuanpro_mc_mindlin_quad_paperstyle_cloud;
load('wc_all_mc.mat', 'wc_all_mc');
wc_all_mc_mindlin = wc_all_mc;
save('wc_all_mc_mindlin.mat', 'wc_all_mc_mindlin');
load('wc_all_mc_thin.mat', 'wc_all_mc_thin');

cfg_data = load('mc_pipeline_config.mat', 'mc_pipeline_config');
sample_count = cfg_data.mc_pipeline_config.sample_count;
sigma = cfg_data.mc_pipeline_config.sigma;
w_list = cfg_data.mc_pipeline_config.w_list;
PV_mc = load('pexact_mc.txt');

fprintf('Shared stage 5 / 5: probability estimate curve assembly\n');
[p_mc_thin, cumulative_mc_thin] = compute_probability_estimate_curve(wc_all_mc_thin, PV_mc, sigma, w_list);
[p_mc_mindlin, cumulative_mc_mindlin] = compute_probability_estimate_curve(wc_all_mc_mindlin, PV_mc, sigma, w_list);

thin = struct();
thin.p_mc = p_mc_thin;
thin.cumulative_mc = cumulative_mc_thin;
thin.wc_all_mc = wc_all_mc_thin;

mindlin = struct();
mindlin.p_mc = p_mc_mindlin;
mindlin.cumulative_mc = cumulative_mc_mindlin;
mindlin.wc_all_mc = wc_all_mc_mindlin;

result = struct();
result.sample_count = sample_count;
result.sigma = sigma;
result.w_list = w_list;
result.PV_mc = PV_mc;
result.thin = thin;
result.mindlin = mindlin;

mat_name = sprintf('mc_reference_truth_dual_plate_n%d.mat', sample_count);
thin_curve_name = sprintf('mc_reference_truth_dual_plate_thin_n%d_curve.txt', sample_count);
mindlin_curve_name = sprintf('mc_reference_truth_dual_plate_mindlin_n%d_curve.txt', sample_count);
save(mat_name, '-struct', 'result');
writematrix([w_list(:), p_mc_thin(:)], thin_curve_name, 'Delimiter', 'tab');
writematrix([w_list(:), p_mc_mindlin(:)], mindlin_curve_name, 'Delimiter', 'tab');

if exist('mc_pipeline_config.mat', 'file') == 2
    delete('mc_pipeline_config.mat');
end

fprintf('Dual-plate MC reference truth finished: n = %d\n', sample_count);
fprintf('Saved MAT file           : %s\n', mat_name);
fprintf('Saved thin curve file    : %s\n', thin_curve_name);
fprintf('Saved Mindlin curve file : %s\n', mindlin_curve_name);

function ensure_eval_centers_mc_file()
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

    writematrix(nodes, 'eval_centers_mc.txt', 'Delimiter', 'tab');
    fprintf('Prepared eval_centers_mc.txt with %d nodes.\n', size(nodes, 1));
end
