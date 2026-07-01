if ~exist('rqmc_driver_state', 'var')
    error('rqmc_driver_state is required before running run_single_dual_plate_rqmc_outer_sample.m');
end

rqmc_pipeline_config = rqmc_driver_state.config;
rqmc_pipeline_config.rng_seed = rqmc_driver_state.seed_list(rqmc_driver_state.current_rep);
save('rqmc_pipeline_config.mat', 'rqmc_pipeline_config');

fprintf('  shared stage 1 / 6: gf_normal_final_rqmc\n');
gf_normal_final_rqmc;
fprintf('  shared stage 2 / 6: voronoi_pca_toolbox_final_blockwise_rqmc\n');
voronoi_pca_toolbox_final_blockwise_rqmc;
fprintf('  shared stage 3 / 6: KL2D_final_complete_3D_rqmc\n');
KL2D_final_complete_3D_rqmc;

cfg_data = load('rqmc_pipeline_config.mat', 'rqmc_pipeline_config');
PV_rqmc = load('pexact_rqmc.txt');

fprintf('  plate stage 4a / 6: thin plate wc calculation\n');
jisuanpro_rqmc_quad_paperstyle_cloud;
load('wc_all_rqmc.mat', 'wc_all_rqmc');
wc_all_rqmc_thin = wc_all_rqmc;

fprintf('  plate stage 4b / 6: Mindlin plate wc calculation\n');
jisuanpro_rqmc_mindlin_quad_paperstyle_cloud;
load('wc_all_rqmc.mat', 'wc_all_rqmc');
wc_all_rqmc_mindlin = wc_all_rqmc;

fprintf('  shared stage 5 / 6: thin probability estimate curve assembly\n');
[p_rqmc_thin, cumulative_rqmc_thin] = compute_probability_estimate_curve( ...
    wc_all_rqmc_thin, PV_rqmc, ...
    cfg_data.rqmc_pipeline_config.sigma, ...
    cfg_data.rqmc_pipeline_config.w_list);

fprintf('  shared stage 6 / 6: Mindlin probability estimate curve assembly\n');
[p_rqmc_mindlin, cumulative_rqmc_mindlin] = compute_probability_estimate_curve( ...
    wc_all_rqmc_mindlin, PV_rqmc, ...
    cfg_data.rqmc_pipeline_config.sigma, ...
    cfg_data.rqmc_pipeline_config.w_list);

run_id = rqmc_driver_state.current_rep;
rqmc_driver_state.curve_pool_thin(run_id, :) = p_rqmc_thin;
rqmc_driver_state.curve_pool_mindlin(run_id, :) = p_rqmc_mindlin;
rqmc_driver_state.wc_runs_thin{run_id, 1} = wc_all_rqmc_thin;
rqmc_driver_state.wc_runs_mindlin{run_id, 1} = wc_all_rqmc_mindlin;
rqmc_driver_state.weight_runs{run_id, 1} = PV_rqmc;
rqmc_driver_state.cumulative_runs_thin{run_id, 1} = cumulative_rqmc_thin;
rqmc_driver_state.cumulative_runs_mindlin{run_id, 1} = cumulative_rqmc_mindlin;

single_run_result = struct();
single_run_result.run_id = run_id;
single_run_result.seed = cfg_data.rqmc_pipeline_config.rng_seed;
single_run_result.w_list = cfg_data.rqmc_pipeline_config.w_list;
single_run_result.PV_rqmc = PV_rqmc;
single_run_result.thin = struct( ...
    'p_rqmc', p_rqmc_thin, ...
    'cumulative_rqmc', cumulative_rqmc_thin, ...
    'wc_all_rqmc', wc_all_rqmc_thin);
single_run_result.mindlin = struct( ...
    'p_rqmc', p_rqmc_mindlin, ...
    'cumulative_rqmc', cumulative_rqmc_mindlin, ...
    'wc_all_rqmc', wc_all_rqmc_mindlin);

save(sprintf('rqmc_outer_sample_dual_plate_%03d.mat', run_id), '-struct', 'single_run_result');

rqmc_driver_state.current_rep = rqmc_driver_state.current_rep + 1;
