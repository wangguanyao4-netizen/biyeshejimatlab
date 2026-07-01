function x = dpim_curve_failure_estimates_weighted(curves, yThreshold, h)
%DPIM_CURVE_FAILURE_ESTIMATES_WEIGHTED Compute P_f estimates for all curves.
R = numel(curves);
x = zeros(R,1);
for r = 1:R
    x(r) = dpim_weighted_failure_prob(curves(r), yThreshold, h);
end
end
