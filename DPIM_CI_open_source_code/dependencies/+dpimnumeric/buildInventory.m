function outputs = buildInventory(projectRoot, reportsDir, runMode, entryName)
%buildInventory Scan MATLAB files and write code/data inventories.

dpimnumeric.ensureDir(reportsDir);

if nargin < 3 || strlength(string(runMode)) == 0
    runMode = "small";
end
if nargin < 4 || strlength(string(entryName)) == 0
    entryName = "run_dpim_ci_numeric_rebuild.m";
end

keywordList = ["DPIM", "RQMC", "confidence interval", "bootstrap", "bootstrap-t", ...
    "percentile", "t interval", "coverage", "probability weights", ...
    "Voronoi", "GF", "discrepancy", "plate", "FEM", "fracture", "fatigue", ...
    "nonlinear", "hybrid"];
paramList = ["n", "R", "B", "M", "h", "d", "alpha", "beta", "seed", ...
    "distribution", "bandwidth", "kernel", "method"];

mFiles = dir(fullfile(projectRoot, "**", "*.m"));
fileCount = numel(mFiles);

relPath = strings(fileCount, 1);
mainName = strings(fileCount, 1);
fileKind = strings(fileCount, 1);
inputSpec = strings(fileCount, 1);
outputSpec = strings(fileCount, 1);
keywordHits = strings(fileCount, 1);
paramHits = strings(fileCount, 1);
detectedCalls = strings(fileCount, 1);
isPlotting = false(fileCount, 1);
hasSeed = false(fileCount, 1);
hasParallel = false(fileCount, 1);
hasLogging = false(fileCount, 1);
hasCache = false(fileCount, 1);
readStatus = strings(fileCount, 1);
texts = cell(fileCount, 1);

for iFile = 1:fileCount
    absPath = fullfile(mFiles(iFile).folder, mFiles(iFile).name);
    relPath(iFile) = localRelativePath(projectRoot, absPath);

    try
        txt = string(fileread(absPath));
        texts{iFile} = txt;
        readStatus(iFile) = "ok";
    catch err
        txt = "";
        texts{iFile} = txt;
        readStatus(iFile) = "unreadable: " + string(err.message);
    end

    [nameValue, kindValue, inputsValue, outputsValue] = localParseSignature(txt, mFiles(iFile).name);
    mainName(iFile) = nameValue;
    fileKind(iFile) = kindValue;
    inputSpec(iFile) = inputsValue;
    outputSpec(iFile) = outputsValue;

    keywordHits(iFile) = localJoinMatches(txt, keywordList);
    paramHits(iFile) = localJoinMatches(txt, paramList);
    isPlotting(iFile) = localContainsAny(txt, ["figure", "plot(", "exportgraphics", "saveas", "print("]);
    hasSeed(iFile) = localContainsAny(txt, ["rng(", "RandStream", "seed"]);
    hasParallel(iFile) = localContainsAny(txt, ["parfor", "parpool", "UseParallel", "batch("]);
    hasLogging(iFile) = localContainsAny(txt, ["diary(", "fprintf(", ".log", "run_log"]);
    hasCache(iFile) = localContainsAny(txt, ["isfile(", "exist(", "load(", "save("]);
end

functionNames = unique(mainName(mainName ~= ""));
for iFile = 1:fileCount
    detectedCalls(iFile) = localDetectCalls(texts{iFile}, functionNames, mainName(iFile));
end

matlabInventory = table(relPath, fileKind, mainName, inputSpec, outputSpec, detectedCalls, ...
    keywordHits, paramHits, isPlotting, hasSeed, hasParallel, hasLogging, hasCache, readStatus, ...
    'VariableNames', {'relative_path', 'file_kind', 'main_function', 'inputs', 'outputs', ...
    'detected_calls', 'keyword_hits', 'parameter_hits', 'is_plotting_file', ...
    'has_seed_control', 'has_parallel_construct', 'has_logging', 'has_cache_or_save_load', 'read_status'});

artifactInventory = localBuildArtifactInventory(projectRoot);
reuseMap = localBuildReuseMap(matlabInventory);

writetable(matlabInventory, fullfile(reportsDir, "matlab_inventory.csv"));
writetable(artifactInventory, fullfile(reportsDir, "data_artifact_inventory.csv"));
writetable(reuseMap, fullfile(reportsDir, "reuse_map.csv"));

markdown = localBuildMarkdown(projectRoot, matlabInventory, artifactInventory, reuseMap, runMode, entryName);
dpimnumeric.writeText(fullfile(reportsDir, "code_inventory.md"), markdown);

conflictText = localConflictText();
dpimnumeric.writeText(fullfile(reportsDir, "conflict_resolution.md"), conflictText);

outputs = struct();
outputs.code_inventory_md = fullfile(reportsDir, "code_inventory.md");
outputs.matlab_inventory_csv = fullfile(reportsDir, "matlab_inventory.csv");
outputs.data_artifact_inventory_csv = fullfile(reportsDir, "data_artifact_inventory.csv");
outputs.reuse_map_csv = fullfile(reportsDir, "reuse_map.csv");
outputs.conflict_resolution_md = fullfile(reportsDir, "conflict_resolution.md");
end

function [nameValue, kindValue, inputsValue, outputsValue] = localParseSignature(txt, fileName)
nameValue = string(erase(fileName, ".m"));
kindValue = "script";
inputsValue = "";
outputsValue = "";

if strlength(txt) == 0
    return;
end

pattern = "(?m)^\s*function\s+(?:(?<outputs>\[[^\]]+\]|[A-Za-z]\w*)\s*=\s*)?(?<name>[A-Za-z]\w*)\s*(?:\((?<inputs>[^)]*)\))?";
match = regexp(char(txt), pattern, "names", "once");
if isempty(match)
    return;
end

kindValue = "function";
nameValue = string(match.name);
if isfield(match, "inputs")
    inputsValue = strtrim(string(match.inputs));
end
if isfield(match, "outputs")
    outputsValue = strtrim(string(match.outputs));
end
end

function out = localJoinMatches(txt, candidates)
hits = strings(0, 1);
lowerTxt = lower(txt);
for iCandidate = 1:numel(candidates)
    candidate = candidates(iCandidate);
    if contains(lowerTxt, lower(candidate))
        hits(end + 1, 1) = candidate; %#ok<AGROW>
    end
end
out = strjoin(unique(hits, "stable"), "; ");
end

function tf = localContainsAny(txt, candidates)
tf = false;
lowerTxt = lower(txt);
for iCandidate = 1:numel(candidates)
    if contains(lowerTxt, lower(candidates(iCandidate)))
        tf = true;
        return;
    end
end
end

function calls = localDetectCalls(txt, functionNames, ownName)
if strlength(txt) == 0
    calls = "";
    return;
end

hits = strings(0, 1);
charTxt = char(txt);
for iName = 1:numel(functionNames)
    fn = functionNames(iName);
    if fn == ownName || strlength(fn) < 4
        continue;
    end
    escapedName = regexptranslate('escape', char(fn));
    pattern = "(?<![\w\.])" + string(escapedName) + "\s*\(";
    if ~isempty(regexp(charTxt, char(pattern), "once"))
        hits(end + 1, 1) = fn; %#ok<AGROW>
    end
    if numel(hits) >= 40
        break;
    end
end

calls = strjoin(unique(hits, "stable"), "; ");
end

function artifactInventory = localBuildArtifactInventory(projectRoot)
patterns = ["*.mat", "*.csv", "*.xlsx", "*.txt", "*.json", "*.fig", "*.png", "*.pdf"];
allRows = table();

for iPattern = 1:numel(patterns)
    files = dir(fullfile(projectRoot, "**", patterns(iPattern)));
    n = numel(files);
    if n == 0
        continue;
    end
    relPath = strings(n, 1);
    extension = strings(n, 1);
    bytes = zeros(n, 1);
    modified = strings(n, 1);
    for iFile = 1:n
        absPath = fullfile(files(iFile).folder, files(iFile).name);
        relPath(iFile) = localRelativePath(projectRoot, absPath);
        [~, ~, ext] = fileparts(files(iFile).name);
        extension(iFile) = string(ext);
        bytes(iFile) = files(iFile).bytes;
        modified(iFile) = string(datetime(files(iFile).datenum, ConvertFrom="datenum"));
    end
    block = table(relPath, extension, bytes, modified, ...
        'VariableNames', {'relative_path', 'extension', 'bytes', 'modified'});
    allRows = [allRows; block]; %#ok<AGROW>
end

artifactInventory = allRows;
end

function reuseMap = localBuildReuseMap(matlabInventory)
experiments = [ ...
    "E1_finite_B"; ...
    "E2_R_order"; ...
    "E3_linear_beam_n_h"; ...
    "E4_nonlinear_tail"; ...
    "E5_bootstrap_t_instability"; ...
    "E6_weights_GF_RQMC"; ...
    "E7_plate_SFEM"; ...
    "E8_simultaneous_band"];

patterns = { ...
    ["finite_B", "finite B", "C0B", "smoke_finite_B", "run_exp1"]; ...
    ["run_exp2", "gaussian_theory", "Edgeworth", "coverage", "ci_methods"]; ...
    ["linear_gaussian", "beam", "green", "run_exp4_beam"]; ...
    ["nonlinear", "lognormal_nonlinear", "pilotfixed", "z^3"]; ...
    ["bootstrap-t", "hybrid", "small_denominator", "compute_hybrid"]; ...
    ["Voronoi", "GF", "RQMC", "weights", "discrepancy"]; ...
    ["plate", "FEM", "Mindlin", "Kirchhoff", "dual_plate"]; ...
    ["simultaneous", "band", "confidence band"]};

rows = table();
for iExp = 1:numel(experiments)
    mask = false(height(matlabInventory), 1);
    for iPattern = 1:numel(patterns{iExp})
        p = lower(patterns{iExp}(iPattern));
        mask = mask | contains(lower(matlabInventory.relative_path), p) ...
            | contains(lower(matlabInventory.keyword_hits), p) ...
            | contains(lower(matlabInventory.main_function), p);
    end
    selected = matlabInventory(mask, :);
    if height(selected) == 0
        pathText = "";
        noteText = "No direct existing implementation found; small wrapper provides an interface.";
    else
        pathText = strjoin(selected.relative_path(1:min(20, height(selected))), "; ");
        noteText = sprintf("Matched %d candidate MATLAB files. First 20 are listed.", height(selected));
    end
    block = table(experiments(iExp), string(pathText), string(noteText), ...
        'VariableNames', {'experiment', 'candidate_existing_files', 'reuse_note'});
    rows = [rows; block]; %#ok<AGROW>
end

reuseMap = rows;
end

function markdown = localBuildMarkdown(projectRoot, matlabInventory, artifactInventory, reuseMap, runMode, entryName)
numM = height(matlabInventory);
numArtifacts = height(artifactInventory);
numPlot = sum(matlabInventory.is_plotting_file);
numSeed = sum(matlabInventory.has_seed_control);
numParallel = sum(matlabInventory.has_parallel_construct);
numLogging = sum(matlabInventory.has_logging);
numCache = sum(matlabInventory.has_cache_or_save_load);

importantMask = matlabInventory.keyword_hits ~= "";
importantRows = matlabInventory(importantMask, :);
topN = min(80, height(importantRows));

lines = strings(0, 1);
lines(end + 1) = "# DPIM CI Code Inventory";
lines(end + 1) = "";
lines(end + 1) = "Generated by `" + string(entryName) + "` in `" + lower(string(runMode)) + "` mode.";
lines(end + 1) = "";
lines(end + 1) = "## Reading Plan";
lines(end + 1) = "";
lines(end + 1) = "- Treat the bundle TeX/PDF files as manuscript and checklist sources.";
lines(end + 1) = "- Treat this MATLAB tree as the implementation source to scan before adding wrappers.";
lines(end + 1) = "- Use the locked nonlinear model in `reports/model_lock.json` for every new nonlinear calculation.";
lines(end + 1) = "- Do not treat existing manuscript figures as numerical evidence unless they are traced to data/config metadata.";
lines(end + 1) = "";
lines(end + 1) = "## First-Stage Code Map Scheme";
lines(end + 1) = "";
lines(end + 1) = "- `reports/matlab_inventory.csv`: per-MATLAB-file signature, keywords, detected calls, and execution features.";
lines(end + 1) = "- `reports/data_artifact_inventory.csv`: existing MAT/CSV/TXT/JSON/FIG/PNG/PDF artifacts.";
lines(end + 1) = "- `reports/reuse_map.csv`: E1--E8 mapping to existing scripts/functions.";
lines(end + 1) = "- `reports/conflict_resolution.md`: explicit conflict and trust rules.";
lines(end + 1) = "";
lines(end + 1) = "## Summary";
lines(end + 1) = "";
lines(end + 1) = sprintf("- Project root: `%s`", projectRoot);
lines(end + 1) = sprintf("- MATLAB files scanned: `%d`", numM);
lines(end + 1) = sprintf("- Existing data/figure/report artifacts scanned: `%d`", numArtifacts);
lines(end + 1) = sprintf("- Plotting/export files: `%d`", numPlot);
lines(end + 1) = sprintf("- Files with seed controls: `%d`", numSeed);
lines(end + 1) = sprintf("- Files with parallel constructs: `%d`", numParallel);
lines(end + 1) = sprintf("- Files with logging calls: `%d`", numLogging);
lines(end + 1) = sprintf("- Files with cache/save/load behavior: `%d`", numCache);
lines(end + 1) = "";
lines(end + 1) = "## E1--E8 Reuse Map";
lines(end + 1) = "";
for iRow = 1:height(reuseMap)
    lines(end + 1) = sprintf("### %s", reuseMap.experiment(iRow));
    lines(end + 1) = "";
    lines(end + 1) = sprintf("- Reuse note: %s", reuseMap.reuse_note(iRow));
    lines(end + 1) = sprintf("- Candidate files: `%s`", reuseMap.candidate_existing_files(iRow));
    lines(end + 1) = "";
end
lines(end + 1) = "## Keyword-Relevant MATLAB Files";
lines(end + 1) = "";
lines(end + 1) = "| relative_path | main_function | keywords | detected_calls |";
lines(end + 1) = "| --- | --- | --- | --- |";
for iRow = 1:topN
    lines(end + 1) = sprintf("| `%s` | `%s` | %s | %s |", ...
        importantRows.relative_path(iRow), importantRows.main_function(iRow), ...
        importantRows.keyword_hits(iRow), importantRows.detected_calls(iRow));
end
if height(importantRows) > topN
    lines(end + 1) = sprintf("| ... | ... | %d additional keyword-relevant files omitted here; see CSV. | ... |", ...
        height(importantRows) - topN);
end
lines(end + 1) = "";
lines(end + 1) = "## Important Caveats";
lines(end + 1) = "";
lines(end + 1) = "- Call detection is conservative text matching, not a full MATLAB parser.";
lines(end + 1) = "- The wrapper creates small reproducibility bundles. It does not certify paper-scale numerical conclusions.";
lines(end + 1) = "- Nonlinear calculations in the new wrapper use `z=mean(theta,2)` exactly as locked in the prompt.";

markdown = strjoin(lines, newline) + newline;
end

function txt = localConflictText()
lines = [
    "# Conflict Resolution"
    ""
    "Rules used by the small rebuild wrapper:"
    ""
    "1. The prompt-level nonlinear model is authoritative for newly generated nonlinear results."
    "2. Existing MATLAB scripts are reuse evidence and implementation references, but older nonlinear definitions must not override `model_lock.json`."
    "3. Existing figures are treated as placeholders unless a data file, configuration, and metadata chain can be identified."
    "4. Small-run outputs are executable smoke outputs, not final paper-scale numerical evidence."
    "5. If existing code and manuscript text disagree, the final paper version must be decided from mathematical definition, code evidence, reproducible config, and regenerated results together."
    ""
];
txt = strjoin(lines, newline) + newline;
end

function rel = localRelativePath(rootPath, absPath)
rootPath = char(string(rootPath));
absPath = char(string(absPath));
if startsWith(absPath, rootPath)
    rel = string(extractAfter(string(absPath), strlength(string(rootPath)) + 1));
else
    rel = string(absPath);
end
rel = replace(rel, "\", "/");
end
