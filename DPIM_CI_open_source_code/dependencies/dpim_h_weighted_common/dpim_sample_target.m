function theta = dpim_sample_target(problem, N, seed)
%DPIM_SAMPLE_TARGET Generate reference samples directly in target space.
if nargin < 3; seed = 1; end
oldState = rng;
rng(seed, 'twister');
cleanupObj = onCleanup(@() rng(oldState)); %#ok<NASGU>

d = problem.d;
switch lower(char(problem.target_distribution))
    case {'standard_normal','normal'}
        theta = randn(N, d);
    case 'lognormal'
        theta = exp(problem.lognormal_mu + problem.lognormal_sigma * randn(N, d));
    otherwise
        error('Unsupported target distribution: %s', problem.target_distribution);
end
end
