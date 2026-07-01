function ci = dpim_bootstrap_t_ci(x, truth, alpha, B, epsReg)
%DPIM_BOOTSTRAP_T_CI Bootstrap-t CI for a scalar curve-pool sample.
% x is the R-vector of independently generated DPIM curve estimates.
if nargin < 5; epsReg = 0; end
x = x(:);
R = numel(x);
mu = mean(x);
s = std(x, 0);
se = s / sqrt(R);

idx = randi(R, B, R);
boot = x(idx);
bm = mean(boot, 2);
bs = std(boot, 0, 2);

if epsReg > 0
    bsEff = max(bs, epsReg);
    seEff = max(se, epsReg/sqrt(R));
else
    bsEff = bs;
    seEff = se;
end

T = sqrt(R) * (bm - mu) ./ bsEff;
zeroMask = bsEff <= 0;
T(zeroMask & bm == mu) = 0;
T(zeroMask & bm >  mu) = Inf;
T(zeroMask & bm <  mu) = -Inf;

tL = quantile(T, alpha/2);
tU = quantile(T, 1-alpha/2);
lo = mu - tU * seEff;
hi = mu - tL * seEff;
infFlag = ~(isfinite(lo) && isfinite(hi));
ci = struct();
ci.lower = lo;
ci.upper = hi;
ci.center = mu;
ci.se = se;
ci.width = hi - lo;
ci.infinite = infFlag;
ci.contains = ~infFlag && lo <= truth && truth <= hi;
ci.t_lower = tL;
ci.t_upper = tU;
end
