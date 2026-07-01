function resultsRoot = section4_06_kirchhoff_plate_centered_coverage(runMode)
%SECTION4_06_KIRCHHOFF_PLATE_CENTERED_COVERAGE Reproduce the Kirchhoff plate study.
% Paper mapping: Section 4.3, stochastic finite-element Kirchhoff plate.

if nargin < 1 || strlength(string(runMode)) == 0
    runMode = "formal";
end
setup_matlab_paths();
resultsRoot = run_kirchhoff_plate_small_deflection_reanalysis(runMode);
end
