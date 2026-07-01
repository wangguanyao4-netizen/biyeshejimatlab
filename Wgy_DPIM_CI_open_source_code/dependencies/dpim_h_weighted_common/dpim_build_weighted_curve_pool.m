function curves = dpim_build_weighted_curve_pool(problem, methodName, n, cfg, Ncurves)
%DPIM_BUILD_WEIGHTED_CURVE_POOL Construct a sample pool of weighted DPIM curves.
%
% The curve pool is a struct array whose elements have identical top-level
% fields.  We deliberately initialize the array from the first generated
% curve, rather than using repmat(struct(),...), because the latter creates
% an empty-field struct array and causes MATLAB to throw
% "subscripted assignment between dissimilar structures".
if nargin < 5 || isempty(Ncurves)
    Ncurves = 1;
end
Ncurves = round(double(Ncurves));
if Ncurves <= 0
    curves = struct([]);
    return;
end

curves(1,1) = dpim_build_weighted_curve(problem, methodName, n, cfg, 1);
for r = 2:Ncurves
    curves(r,1) = dpim_build_weighted_curve(problem, methodName, n, cfg, r);
end
end
