function theta = dpim_transform_unit_to_target(U, problem)
%DPIM_TRANSFORM_UNIT_TO_TARGET Map U in [0,1]^d into target probability space.
epsProb = 1e-15;
U = min(max(U, epsProb), 1-epsProb);
switch lower(char(problem.center_transform))
    case 'normal_icdf'
        theta = sqrt(2) * erfinv(2*U - 1);
    case 'lognormal_icdf'
        z = sqrt(2) * erfinv(2*U - 1);
        theta = exp(problem.lognormal_mu + problem.lognormal_sigma * z);
    case 'identity'
        theta = U;
    otherwise
        error('Unsupported center_transform: %s', problem.center_transform);
end
end
