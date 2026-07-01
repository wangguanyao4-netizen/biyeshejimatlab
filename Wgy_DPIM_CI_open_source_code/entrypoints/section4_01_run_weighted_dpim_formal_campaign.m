function outputs = section4_01_run_weighted_dpim_formal_campaign(runMode, resultsRoot)
%SECTION4_01_RUN_WEIGHTED_DPIM_FORMAL_CAMPAIGN Reproduce the main weighted DPIM campaign.
% Paper mapping: Section 4, formal probability-weighted DPIM experiments.

if nargin < 1 || strlength(string(runMode)) == 0
    runMode = "formal";
end
setup_matlab_paths();

if nargin < 2
    outputs = run_weighted_paper_formal_campaign(runMode);
else
    outputs = run_weighted_paper_formal_campaign(runMode, resultsRoot);
end
end
