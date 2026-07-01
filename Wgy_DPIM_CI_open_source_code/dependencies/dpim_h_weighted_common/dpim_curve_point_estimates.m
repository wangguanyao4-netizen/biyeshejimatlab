function x = dpim_curve_point_estimates(curves, y0, h)
%DPIM_CURVE_POINT_ESTIMATES Compute Q_{r,h}(y0) for all curves.
R = numel(curves);
x = zeros(R,1);
for r = 1:R
    x(r) = dpim_weighted_kde_point(curves(r), y0, h);
end
end
