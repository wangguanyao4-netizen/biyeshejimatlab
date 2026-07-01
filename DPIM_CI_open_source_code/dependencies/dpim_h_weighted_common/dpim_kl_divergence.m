function KL = dpim_kl_divergence(p, q, yGrid)
p = p(:); q = q(:); yGrid = yGrid(:);
p = max(p, 1e-300); q = max(q, 1e-300);
p = p / trapz(yGrid, p);
q = q / trapz(yGrid, q);
KL = trapz(yGrid, p .* log(p ./ q));
end
