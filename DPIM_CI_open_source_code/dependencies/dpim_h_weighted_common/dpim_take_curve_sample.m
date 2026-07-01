function x = dpim_take_curve_sample(poolValues, R, withReplacement)
%DPIM_TAKE_CURVE_SAMPLE Draw R curve values from a sample pool.
N = numel(poolValues);
if nargin < 3; withReplacement = true; end
if withReplacement || R > N
    idx = randi(N, R, 1);
else
    idx = randperm(N, R).';
end
x = poolValues(idx);
end
