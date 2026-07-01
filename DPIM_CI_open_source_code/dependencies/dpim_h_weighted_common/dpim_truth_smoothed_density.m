function truth = dpim_truth_smoothed_density(problem, y0, h, cfg, seedOffset)
%DPIM_TRUTH_SMOOTHED_DENSITY High-sample reference for E[K_h(y0-Y)].
if nargin < 5; seedOffset = 0; end
theta = dpim_sample_target(problem, cfg.truth_N, cfg.seed + 31001 + seedOffset);
y = problem.response_fun(theta);
truth = mean(dpim_gaussian_kernel(y0 - y(:), h));
end
