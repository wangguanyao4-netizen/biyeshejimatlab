# DPIM CI open-source code

This folder collects the scripts used for the paper
`Probability Error Estimation for Direct Probability Integral Method Based on Two-Sided Confidence Intervals`.

The folder is organized for paper readers:

- `entrypoints/` contains paper-readable scripts named by section, figure, or table.
- `dependencies/original_matlab/` contains the original MATLAB drivers copied from the working repository.
- `dependencies/original_python/` contains the original Python analysis and figure scripts copied from the working repository.
- `dependencies/+dpimnumeric/` contains MATLAB package helpers required by the drivers.
- `dependencies/dpim_h_weighted_common/` and `dependencies/external_weight_providers/` contain the probability-weighted DPIM helper functions and Voronoi probability-weight providers.

The entry-point names are intentionally different from the original internal names so that a reader can immediately map each script to the paper. The original function filenames are retained in `dependencies/` because MATLAB function files must keep their function names.

## Main entry points

- `entrypoints/section4_01_run_weighted_dpim_formal_campaign.m`:
  main probability-weighted DPIM formal campaign for Section 4.
- `entrypoints/02_section4_build_fig01_to_fig11.py`:
  rebuilds the numerical figures used in the paper.
- `entrypoints/03_section3_table1_direct_coverage_coefficients.py`:
  computes the direct coverage-expansion coefficients used in Section 3 and Table 1.
- `entrypoints/05_section4_internal_sample_size_coverage.py`:
  regenerates the internal-sample-size coverage analysis for \(n=192,384,768\).
- `entrypoints/section4_06_kirchhoff_plate_centered_coverage.m`:
  runs the stochastic finite-element Kirchhoff plate coverage study.
- `entrypoints/07_section4_plate_fixed_weight_formula_validation.py`:
  evaluates the fixed-probability-weight factorized formula for the plate example.

## MATLAB usage

From MATLAB, run:

```matlab
cd("<repository-root>/DPIM_CI_open_source_code/entrypoints")
setup_matlab_paths
section4_01_run_weighted_dpim_formal_campaign("diagnostic")
```

Use `"diagnostic"` first. Full formal runs can be time-consuming.

## Python usage

From PowerShell:

```powershell
cd <repository-root>/DPIM_CI_open_source_code
python .\entrypoints\03_section3_table1_direct_coverage_coefficients.py
```

Some scripts expect the same result folders used during paper preparation. If those large numerical result folders are not included, run the corresponding MATLAB/Python campaign first or adjust the input paths in the original script.

Note: the open-source folder includes lightweight scripts and code dependencies only. Large numerical result folders are intentionally excluded from version control.
