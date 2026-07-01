"""Compute direct coverage-expansion coefficients.

Paper mapping: Section 3 and Table 1, finite-B coefficient calculation for
Student t, Percentile Bootstrap, and Bootstrap-t intervals.
"""

from pathlib import Path
import runpy

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "dependencies" / "original_python" / "derive_direct_coverage_formulas.py"
runpy.run_path(str(SCRIPT), run_name="__main__")

