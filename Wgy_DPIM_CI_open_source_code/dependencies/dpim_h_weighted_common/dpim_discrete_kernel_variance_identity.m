function out = dpim_discrete_kernel_variance_identity(curve, y0, h)
%DPIM_DISCRETE_KERNEL_VARIANCE_IDENTITY Weighted finite-measure identity.
% This is NOT the curve-to-curve variance.  It is the variance of the kernel
% response under the discrete weighted quadrature measure of one DPIM curve.
y = curve.y(:);
w = curve.weights(:);
kh = dpim_gaussian_kernel(y0 - y, h);
q = sum(w .* kh);
directSecond = sum(w .* kh.^2);
discVarDirect = directSecond - q.^2;
qHalf = sum(w .* dpim_gaussian_kernel(y0 - y, h/sqrt(2)));
secondIdentity = (1/(2*sqrt(pi)*h)) * qHalf;
discVarIdentity = secondIdentity - q.^2;

% Local normalized kernel influence weights; useful as a physical local
% effective number, not as a replacement for Voronoi probabilities.
unnorm = w .* exp(-0.5*((y0-y)./h).^2);
if sum(unnorm) > 0
    loc = unnorm / sum(unnorm);
    localNeff = 1 / sum(loc.^2);
else
    localNeff = 0;
end

out = struct();
out.q_h = q;
out.q_h_over_sqrt2 = qHalf;
out.second_direct = directSecond;
out.second_identity = secondIdentity;
out.disc_var_direct = max(discVarDirect, 0);
out.disc_var_identity = max(discVarIdentity, 0);
out.identity_abs_error = abs(discVarDirect - discVarIdentity);
out.identity_rel_error = abs(discVarDirect - discVarIdentity) / max(abs(discVarDirect), realmin);
out.local_neff_kernel = localNeff;
out.local_mass = sum(unnorm);
end
