function rootDir = setup_matlab_paths()
%SETUP_MATLAB_PATHS Add the paper-code dependencies to the MATLAB path.
% Run this before calling the MATLAB entry points in this directory.

rootDir = fileparts(fileparts(mfilename("fullpath")));
addpath(fullfile(rootDir, "dependencies"));
addpath(fullfile(rootDir, "dependencies", "original_matlab"));
addpath(fullfile(rootDir, "dependencies", "dpim_h_weighted_common"));
addpath(fullfile(rootDir, "dependencies", "external_weight_providers"));
end

