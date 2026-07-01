function truth = dpim_truth_failure_prob(problem, yThreshold, h, cfg, seedOffset)
%DPIM_TRUTH_FAILURE_PROB High-sample reference for smoothed failure probability.
if nargin < 5; seedOffset = 0; end
theta = dpim_sample_target(problem, cfg.truth_N, cfg.seed + 41001 + seedOffset);
y = problem.response_fun(theta);
z = (yThreshold - y(:)) ./ h;
truth = mean(0.5 * erfc(z ./ sqrt(2)));
truth = min(max(truth, 0), 1);
end
