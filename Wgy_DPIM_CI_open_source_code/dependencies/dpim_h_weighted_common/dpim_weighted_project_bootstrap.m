function dpim_weighted_project_bootstrap(projectRoot)
%DPIM_WEIGHTED_PROJECT_BOOTSTRAP Put the current weighted DPIM package first.
%
% This prevents stale helper functions from older generated packages from
% being used accidentally.  It intentionally does not call restoredefaultpath,
% so user/toolbox paths are preserved.

if nargin < 1 || strlength(string(projectRoot)) == 0
    projectRoot = fileparts(mfilename('fullpath'));
end
projectRoot = char(string(projectRoot));

% Put this package before any older package paths.
addpath(genpath(projectRoot), '-begin');
rehash;

% Clear only project helper functions that have changed across generated
% package versions.  Do not clear the currently running entry-point function.
clear('dpim_generate_point_set');
clear('dpim_build_weighted_curve');
clear('dpim_build_weighted_curve_pool');
clear('dpim_get_voronoi_weights');
clear('dpim_internal_voronoi_weights');
clear('dpim_weighted_h_config');
clear('dpim_weighted_define_problems');
clear('run_W0_weight_audit');
clear('run_W1_weighted_feature_scale');
clear('run_W2_weighted_point_ci_coverage');
clear('run_W3_weighted_failure_probability');
clear('run_W4_weighted_variance_h_formula');
clear('run_W5_weighted_selection_summary');
end
