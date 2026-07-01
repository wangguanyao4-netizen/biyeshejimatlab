function fit = dpim_fit_power_law(h, v, smallCount)
%DPIM_FIT_POWER_LAW Fit v = A h^{-beta} globally and on small-h points.
if nargin < 3; smallCount = min(4, numel(h)); end
h = h(:); v = v(:);
mask = isfinite(h) & isfinite(v) & h > 0 & v > 0;
h = h(mask); v = v(mask);
fit = struct('A_global',NaN,'beta_global',NaN,'R2_global',NaN, ...
             'A_small',NaN,'beta_small',NaN,'R2_small',NaN);
if numel(h) < 2; return; end
[hs, idx] = sort(h); vs = v(idx);
[Ag, bg, Rg] = local_fit(hs, vs);
fit.A_global = Ag; fit.beta_global = bg; fit.R2_global = Rg;
smallCount = min(smallCount, numel(hs));
[As, bs, Rs] = local_fit(hs(1:smallCount), vs(1:smallCount));
fit.A_small = As; fit.beta_small = bs; fit.R2_small = Rs;
end

function [A, beta, R2] = local_fit(h, v)
x = log(h(:)); y = log(v(:));
p = polyfit(x, y, 1);
yhat = polyval(p, x);
beta = -p(1);
A = exp(p(2));
SST = sum((y-mean(y)).^2);
SSE = sum((y-yhat).^2);
if SST > 0
    R2 = 1 - SSE/SST;
else
    R2 = NaN;
end
end
