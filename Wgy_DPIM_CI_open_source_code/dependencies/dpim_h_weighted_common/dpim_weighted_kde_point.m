function q = dpim_weighted_kde_point(curve, y0, h)
%DPIM_WEIGHTED_KDE_POINT Weighted DPIM/KDE density estimate at y0.
q = sum(curve.weights(:) .* dpim_gaussian_kernel(y0 - curve.y(:), h));
end
