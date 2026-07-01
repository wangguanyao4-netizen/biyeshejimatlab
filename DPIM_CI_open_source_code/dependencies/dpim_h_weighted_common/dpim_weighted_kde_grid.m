function pdf = dpim_weighted_kde_grid(curve, yGrid, h)
%DPIM_WEIGHTED_KDE_GRID Weighted DPIM/KDE density estimate on a grid.
y = curve.y(:);
w = curve.weights(:);
yGrid = yGrid(:);
pdf = zeros(numel(yGrid),1);
for i = 1:numel(yGrid)
    pdf(i) = sum(w .* dpim_gaussian_kernel(yGrid(i) - y, h));
end
end
