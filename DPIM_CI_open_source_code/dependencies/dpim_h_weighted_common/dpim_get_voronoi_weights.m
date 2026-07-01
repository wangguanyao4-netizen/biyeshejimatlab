function weightData = dpim_get_voronoi_weights(n, d, methodName, problem, U, cfg, curveIndex)
%DPIM_GET_VORONOI_WEIGHTS Obtain Voronoi cell probability weights.
% This wrapper first calls the problem-specific provider supplied by the
% user.  If unavailable, it falls back to an internal exact blockwise
% Voronoi assignment using auxiliary target samples.  The fallback is only
% for portability; the intended path is the supplied provider.
if nargin < 7; curveIndex = 1; end

context = struct();
context.project_dir = cfg.project_root;
context.point_seed = cfg.seed + 100000*curveIndex + 97*double(sum(char(string(methodName))));
context.pool_index = curveIndex;
context.n = n;
context.d = d;
context.weighting_cfg = cfg.weighting_cfg;
context.weighting_cfg.voronoi_center_transform = problem.center_transform;
context.weighting_cfg.voronoi_target_distribution = problem.target_distribution;
context.weighting_cfg.voronoi_output_tag = [char(problem.short_name) '_' char(string(methodName))];
if strcmpi(problem.target_distribution, 'lognormal')
    context.weighting_cfg.theta_lognormal_mu = problem.lognormal_mu;
    context.weighting_cfg.theta_lognormal_sigma = problem.lognormal_sigma;
    context.weighting_cfg.voronoi_lognormal_mu = problem.lognormal_mu;
    context.weighting_cfg.voronoi_lognormal_sigma = problem.lognormal_sigma;
    context.weighting_cfg.theta_lognormal_nonlinear_mu = problem.lognormal_mu;
    context.weighting_cfg.theta_lognormal_nonlinear_sigma = problem.lognormal_sigma;
    context.weighting_cfg.voronoi_lognormal_nonlinear_mu = problem.lognormal_mu;
    context.weighting_cfg.voronoi_lognormal_nonlinear_sigma = problem.lognormal_sigma;
end

providerName = char(problem.provider);
try
    if exist(providerName, 'file') == 2
        fh = str2func(providerName);
        weightData = fh(n, d, char(string(methodName)), problem.name, U, context);
        weightData.weight_source = ['external:' providerName];
    else
        error('Provider %s not found on MATLAB path.', providerName);
    end
catch ME
    warning('Voronoi provider failed (%s). Using internal fallback.', ME.message);
    weightData = dpim_internal_voronoi_weights(n, d, methodName, problem, U, cfg, context);
    weightData.weight_source = 'internal_fallback';
end

w = double(weightData.weights(:));
if numel(w) ~= n
    error('Weight provider returned %d weights; expected %d.', numel(w), n);
end
s = sum(w);
if ~isfinite(s) || s <= 0
    error('Invalid Voronoi weights: sum=%g.', s);
end
w = max(w, 0);
w = w / sum(w);
weightData.weights = w;
weightData.probabilities = w;
if ~isfield(weightData, 'centers_target') || isempty(weightData.centers_target)
    weightData.centers_target = dpim_transform_unit_to_target(U, problem);
end
weightData.sum_weights = sum(w);
weightData.l1_from_equal = sum(abs(w - 1/n));
weightData.weight_cv = std(w) / mean(w);
weightData.weight_entropy = -sum(w(w>0).*log(w(w>0)));
weightData.weight_ess = 1 / sum(w.^2);
end
