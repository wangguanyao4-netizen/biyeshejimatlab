function pf = dpim_weighted_failure_prob(curve, yThreshold, h)
%DPIM_WEIGHTED_FAILURE_PROB Smoothed failure probability.
% Integral of weighted Gaussian-smoothed density over y > yThreshold.
z = (yThreshold - curve.y(:)) ./ h;
tail = 0.5 * erfc(z ./ sqrt(2));
pf = sum(curve.weights(:) .* tail);
pf = min(max(pf, 0), 1);
end
