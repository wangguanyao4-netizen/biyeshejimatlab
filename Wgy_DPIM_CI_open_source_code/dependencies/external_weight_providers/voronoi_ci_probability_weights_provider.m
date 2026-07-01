function weightData = voronoi_ci_probability_weights_provider(n, d, methodName, integrandName, U, context)
persistent cacheMap
if nargin < 6
    context = struct();
end

validateattributes(U, {'double'}, {'2d', 'nonempty'});

if size(U, 1) ~= n || size(U, 2) ~= d
    error("Input U must have size [%d, %d].", n, d);
end

opts = local_build_options(n, d, methodName, integrandName, context);

shortMethod = local_method_short_name(methodName);
opts.method_short = shortMethod;

cacheKey = local_cache_key(shortMethod, n, d, opts, context);

if isempty(cacheMap)
    cacheMap = containers.Map('KeyType', 'char', 'ValueType', 'any');
end

if opts.enable_cache && isKey(cacheMap, cacheKey)
    weightData = cacheMap(cacheKey);
    return;
end

Z = local_transform_to_target_space(U, opts.center_transform);

Xaux = local_generate_auxiliary_samples(opts, d, context);

[labels, minDist2, assignmentBackend] = local_assign_voronoi(Xaux, Z, opts);

counts = accumarray(labels, 1, [n, 1], @sum, 0);

probabilities = counts / opts.aux_sample_count;

weightData = struct();
weightData.probabilities = probabilities;
weightData.weights = probabilities;
weightData.counts = counts;
weightData.centers_target = Z;
weightData.auxiliary_sample_count = opts.aux_sample_count;
weightData.block_size = opts.block_size;
weightData.center_transform = char(opts.center_transform);
weightData.target_distribution = char(opts.target_distribution);
weightData.assignment_backend = char(assignmentBackend);
weightData.method_short_name = shortMethod;
weightData.note = sprintf([ ...
    'CI experiment Voronoi probabilities for method=%s, n=%d, d=%d. ', ...
    'Centers were mapped from U to target space using %s; target auxiliary ', ...
    'samples were drawn from %s.'], ...
    shortMethod, n, d, char(opts.center_transform), char(opts.target_distribution));
weightData.output_txt_path = "";
weightData.output_mat_path = "";
weightData.min_cell_probability = min(probabilities);
weightData.max_cell_probability = max(probabilities);
weightData.mean_cell_probability = mean(probabilities);
weightData.empty_cell_count = sum(counts == 0);
weightData.mean_min_distance2 = mean(minDist2);

if opts.save_outputs
    [txtPath, matPath] = local_save_outputs(weightData, probabilities, Z, Xaux, counts, minDist2, opts, context, shortMethod, n, d);
    weightData.output_txt_path = string(txtPath);
    weightData.output_mat_path = string(matPath);
    weightData.note = sprintf('%s Exported to %s and %s.', weightData.note, txtPath, matPath);
end

if opts.enable_cache
    cacheMap(cacheKey) = weightData;
end
end

function opts = local_build_options(n, d, methodName, integrandName, context)

opts = struct();

opts.center_transform = "normal_icdf";

opts.target_distribution = "standard_normal";

opts.enable_cache = true;

opts.save_outputs = false;

opts.export_txt = true;
opts.export_mat = true;

opts.output_tag = "ci_experiment";

opts.output_dir = "";

opts.assignment_backend = "auto";

if isfield(context, "project_dir") && strlength(string(context.project_dir)) > 0
    opts.output_dir = fullfile(char(string(context.project_dir)), "ci_probability_weights");
else
opts.output_dir = fullfile(pwd, "ci_probability_weights");
end

if isfield(context, "mode") && strcmpi(string(context.mode), "standalone_weight_export")
    opts.aux_sample_count = 10000;
else
opts.aux_sample_count = min(max(1000, 8 * n), 4000);
end

opts.block_size = min(1000, opts.aux_sample_count);

if isfield(context, "weighting_cfg") && isstruct(context.weighting_cfg)
    wcfg = context.weighting_cfg;

    if isfield(wcfg, "voronoi_aux_sample_count") && ~isempty(wcfg.voronoi_aux_sample_count)
        opts.aux_sample_count = double(wcfg.voronoi_aux_sample_count);
    end
    if isfield(wcfg, "voronoi_block_size") && ~isempty(wcfg.voronoi_block_size)
        opts.block_size = double(wcfg.voronoi_block_size);
    end
    if isfield(wcfg, "voronoi_center_transform") && ~isempty(wcfg.voronoi_center_transform)
        opts.center_transform = string(wcfg.voronoi_center_transform);
    end
    if isfield(wcfg, "voronoi_target_distribution") && ~isempty(wcfg.voronoi_target_distribution)
        opts.target_distribution = string(wcfg.voronoi_target_distribution);
    end
    if isfield(wcfg, "voronoi_save_outputs") && ~isempty(wcfg.voronoi_save_outputs)
        opts.save_outputs = logical(wcfg.voronoi_save_outputs);
    end
    if isfield(wcfg, "voronoi_output_dir") && strlength(string(wcfg.voronoi_output_dir)) > 0
        opts.output_dir = char(string(wcfg.voronoi_output_dir));
    end
    if isfield(wcfg, "voronoi_output_tag") && strlength(string(wcfg.voronoi_output_tag)) > 0
        opts.output_tag = string(wcfg.voronoi_output_tag);
    end
    if isfield(wcfg, "voronoi_enable_cache") && ~isempty(wcfg.voronoi_enable_cache)
        opts.enable_cache = logical(wcfg.voronoi_enable_cache);
    end
    if isfield(wcfg, "voronoi_assignment_backend") && ~isempty(wcfg.voronoi_assignment_backend)
        opts.assignment_backend = lower(string(wcfg.voronoi_assignment_backend));
    elseif isfield(wcfg, "voronoi_backend") && ~isempty(wcfg.voronoi_backend)
        opts.assignment_backend = lower(string(wcfg.voronoi_backend));
    elseif isfield(wcfg, "assignment_backend") && ~isempty(wcfg.assignment_backend)
        opts.assignment_backend = lower(string(wcfg.assignment_backend));
    end
end

opts.aux_sample_count = max(1, round(opts.aux_sample_count));

opts.block_size = max(1, round(min(opts.block_size, opts.aux_sample_count)));
if ~any(opts.assignment_backend == ["auto", "cpu", "gpu"])
    error("Unsupported Voronoi assignment backend '%s'.", opts.assignment_backend);
end

if nargin >= 4
end
end

function key = local_cache_key(shortMethod, n, d, opts, context)

if strcmpi(shortMethod, "qmc")
    key = sprintf('%s|n=%d|d=%d|aux=%d|blk=%d|tr=%s|tg=%s|backend=%s', ...
        shortMethod, n, d, opts.aux_sample_count, opts.block_size, ...
        char(opts.center_transform), char(opts.target_distribution), ...
        char(opts.assignment_backend));
else
if isfield(context, "point_seed")
        seedVal = double(context.point_seed);
    else
        seedVal = -1;
    end

    key = sprintf('%s|n=%d|d=%d|seed=%d|aux=%d|blk=%d|tr=%s|tg=%s|backend=%s', ...
        shortMethod, n, d, seedVal, opts.aux_sample_count, opts.block_size, ...
        char(opts.center_transform), char(opts.target_distribution), ...
        char(opts.assignment_backend));
end
end

function Z = local_transform_to_target_space(U, transformName)

switch lower(char(string(transformName)))
    case "normal_icdf"
epsProb = 1e-15;
        Uc = min(max(U, epsProb), 1 - epsProb);

Z = sqrt(2) * erfinv(2 * Uc - 1);

    case "identity"
Z = U;

    otherwise
        error("Unsupported center transform '%s'.", transformName);
end
end

function Xaux = local_generate_auxiliary_samples(opts, d, context)

seed = local_aux_seed(opts, context);

previousState = rng;

rng(seed, "twister");

cleanupObj = onCleanup(@() rng(previousState));

switch lower(char(string(opts.target_distribution)))
    case {"standard_normal", "normal"}
        Xaux = randn(opts.aux_sample_count, d);

    case "uniform"
        Xaux = rand(opts.aux_sample_count, d);

    otherwise
        error("Unsupported target distribution '%s'.", opts.target_distribution);
end
end

function seed = local_aux_seed(opts, context)

base = 730001 + 37 * opts.aux_sample_count + 19 * opts.block_size;

if isfield(context, "point_seed") && ~isempty(context.point_seed)
    pointSeed = double(context.point_seed);
else
    pointSeed = 0;
end

if isfield(context, "n") && ~isempty(context.n)
    nVal = double(context.n);
else
    nVal = 0;
end

if isfield(context, "d") && ~isempty(context.d)
    dVal = double(context.d);
else
    dVal = 0;
end

if isfield(context, "mode") && strcmpi(string(context.mode), "standalone_weight_export")
    seed = base + 101 * nVal + 13 * dVal;

elseif strcmpi(opts.method_short, "qmc")
    seed = base + 101 * nVal + 13 * dVal;

else
    seed = base + pointSeed + 101 * nVal + 13 * dVal;
end
end

function [label, minDist2, backendUsed] = local_assign_voronoi(X, Z, opts)
backend = lower(string(opts.assignment_backend));
if backend == "gpu" || (backend == "auto" && local_can_use_gpu())
    try
        [label, minDist2] = local_assign_voronoi_gpu(X, Z, opts.block_size);
        backendUsed = "gpu";
        return;
    catch ME
        warning("GPU Voronoi assignment failed (%s). Falling back to CPU.", ME.message);
        try
            reset(gpuDevice);
        catch
        end
    end
end

[label, minDist2] = local_assign_voronoi_exact_blockwise(X, Z, opts.block_size);
backendUsed = "cpu";
end

function tf = local_can_use_gpu()
tf = false;
try
    gpuDevice;
    tf = true;
catch
    tf = false;
end
end

function [label, minDist2] = local_assign_voronoi_exact_blockwise(X, Z, blockSize)

[M, d] = size(X);

[N, d2] = size(Z);

if d ~= d2
    error("X and Z must have the same dimension.");
end

label = zeros(M, 1);

minDist2 = zeros(M, 1);

Z2 = sum(Z.^2, 2).';

nBlocks = ceil(M / blockSize);

for b = 1:nBlocks
i1 = (b - 1) * blockSize + 1;

i2 = min(b * blockSize, M);

Xb = X(i1:i2, :);

Xb2 = sum(Xb.^2, 2);

D2 = Xb2 + Z2 - 2 * (Xb * Z.');

D2(D2 < 0) = 0;

[minDist2(i1:i2), label(i1:i2)] = min(D2, [], 2);
end

if any(label < 1) || any(label > N)
    error("Voronoi assignment produced out-of-range labels.");
end
end

function [label, minDist2] = local_assign_voronoi_gpu(X, Z, blockSize)
[M, d] = size(X);
[N, d2] = size(Z);
if d ~= d2
    error("X and Z must have the same dimension.");
end

dev = gpuDevice;
Xg = gpuArray(single(X));
Zg = gpuArray(single(Z));
Z2 = sum(Zg.^2, 2).';

labelG = gpuArray.zeros(M, 1, 'uint32');
minDist2G = gpuArray.zeros(M, 1, 'single');
nBlocks = ceil(M / blockSize);

for b = 1:nBlocks
    i1 = (b - 1) * blockSize + 1;
    i2 = min(b * blockSize, M);
    Xb = Xg(i1:i2, :);
    Xb2 = sum(Xb.^2, 2);
    D2 = Xb2 + Z2 - 2 * (Xb * Zg.');
    D2(D2 < 0) = 0;
    [minDist2G(i1:i2), labelG(i1:i2)] = min(D2, [], 2);
end
wait(dev);

label = double(gather(labelG));
minDist2 = double(gather(minDist2G));
if any(label < 1) || any(label > N)
    error("GPU Voronoi assignment produced out-of-range labels.");
end
end

function [txtPath, matPath] = local_save_outputs(weightData, probabilities, Z, Xaux, counts, minDist2, opts, context, shortMethod, n, d)

if ~exist(opts.output_dir, "dir")
    mkdir(opts.output_dir);
end

fileStem = local_output_stem(opts.output_tag, shortMethod, n, d, context);

txtPath = fullfile(opts.output_dir, fileStem + ".txt");
matPath = fullfile(opts.output_dir, fileStem + ".mat");

if opts.export_txt
    writematrix(probabilities, txtPath, "Delimiter", "space");
end

if opts.export_mat
    metadata = struct();
    metadata.method = shortMethod;
    metadata.n = n;
    metadata.d = d;
    metadata.output_tag = char(opts.output_tag);
    metadata.context = context;
    metadata.options = opts;
    save(matPath, "probabilities", "Z", "Xaux", "counts", "minDist2", "metadata", "weightData", "-v7.3");
end
end

function stem = local_output_stem(outputTag, shortMethod, n, d, context)

parts = {char(string(outputTag)), 'probweights', char(string(shortMethod)), sprintf('n%d', n), sprintf('d%d', d)};

if isfield(context, "point_seed") && ~isempty(context.point_seed)
    parts{end + 1} = sprintf('seed%d', double(context.point_seed));
end

if isfield(context, "pool_index") && ~isempty(context.pool_index)
    parts{end + 1} = sprintf('pool%04d', double(context.pool_index));
end

if isfield(context, "batch_label") && strlength(string(context.batch_label)) > 0
    parts{end + 1} = char(string(context.batch_label));
end

stem = string(strjoin(parts, "_"));
end

function shortMethod = local_method_short_name(methodName)

switch lower(char(string(methodName)))
    case "mc"
        shortMethod = "mc";
    case "sobol_qmc"
        shortMethod = "qmc";
    case "sobol_scrambled"
        shortMethod = "rqmc";
    otherwise
        shortMethod = char(string(methodName));
end
end
