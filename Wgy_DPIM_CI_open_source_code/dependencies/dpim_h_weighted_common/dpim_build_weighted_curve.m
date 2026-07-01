function curve = dpim_build_weighted_curve(problem, methodName, n, cfg, curveIndex)
%DPIM_BUILD_WEIGHTED_CURVE Build one probability-weighted DPIM curve.
%
% The point set U is generated in unit probability space; the Voronoi weight
% provider maps the same centers to the target probability space and returns
% cell probabilities P_i.  The DPIM curve is therefore
%     Q_h(y) = sum_i P_i K_h(y - g(theta_i)),
% not an equal-weight KDE unless the Voronoi probabilities happen to be equal.

seed = cfg.seed + 100000 * curveIndex + 1009 * double(sum(char(string(methodName)))) + 17 * problem.d;
[U, pointInfo] = dpim_generate_point_set(n, problem.d, methodName, seed);

weightData = dpim_get_voronoi_weights(n, problem.d, methodName, problem, U, cfg, curveIndex);
theta = double(weightData.centers_target);
y = problem.response_fun(theta);
y = y(:);
w = double(weightData.weights(:));

curve = struct();
curve.problem = problem.name;
curve.method = char(string(methodName));
curve.curve_index = curveIndex;
curve.pointInfo = pointInfo;
curve.point_requested_method = char(pointInfo.requested_method);
curve.point_actual_method = char(pointInfo.actual_method);
curve.point_fallback_used = logical(pointInfo.fallback_used);
curve.point_message = char(pointInfo.message);
curve.U = U;
curve.theta = theta;
curve.y = y;
curve.weights = w;
curve.weightData = weightData;
curve.weight_source = weightData.weight_source;
curve.sum_weights = sum(w);
curve.weight_ess = 1 / sum(w.^2);
curve.weight_cv = std(w) / mean(w);
curve.l1_from_equal = sum(abs(w - 1/n));
end
