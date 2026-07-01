function weightData = dpim_internal_voronoi_weights(n, d, methodName, problem, U, cfg, context)
%DPIM_INTERNAL_VORONOI_WEIGHTS Portable fallback for Voronoi probabilities.
Z = dpim_transform_unit_to_target(U, problem);
Naux = cfg.aux_sample_count;
Xaux = dpim_sample_target(problem, Naux, context.point_seed + 7919);
blockSize = min(1000, Naux);
[labels, minDist2] = local_assign_voronoi_exact_blockwise(Xaux, Z, blockSize);
counts = accumarray(labels, 1, [n, 1], @sum, 0);
prob = counts / Naux;
weightData = struct();
weightData.weights = prob;
weightData.probabilities = prob;
weightData.counts = counts;
weightData.centers_target = Z;
weightData.auxiliary_sample_count = Naux;
weightData.block_size = blockSize;
weightData.center_transform = problem.center_transform;
weightData.target_distribution = problem.target_distribution;
weightData.method_short_name = char(string(methodName));
weightData.min_cell_probability = min(prob);
weightData.max_cell_probability = max(prob);
weightData.mean_cell_probability = mean(prob);
weightData.empty_cell_count = sum(counts == 0);
weightData.mean_min_distance2 = mean(minDist2);
end

function [label, minDist2] = local_assign_voronoi_exact_blockwise(X, Z, blockSize)
[M, d] = size(X);
[N, d2] = size(Z);
if d ~= d2; error('X and Z dimensions differ.'); end
label = zeros(M,1);
minDist2 = zeros(M,1);
Z2 = sum(Z.^2, 2).';
nBlocks = ceil(M / blockSize);
for b = 1:nBlocks
    i1 = (b-1)*blockSize + 1;
    i2 = min(b*blockSize, M);
    Xb = X(i1:i2,:);
    D2 = sum(Xb.^2,2) + Z2 - 2*(Xb*Z.');
    D2(D2<0) = 0;
    [minDist2(i1:i2), label(i1:i2)] = min(D2, [], 2);
end
end
