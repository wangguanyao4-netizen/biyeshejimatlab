"""Run fixed-probability-weight formula validation for the plate example.

Paper mapping: Section 4.3, comparison between centered stochastic-plate
coverage and the factorized fixed-probability-weight coverage formula.
"""

from pathlib import Path
import runpy

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "dependencies" / "original_python" / "random_weight_coverage_expansion__analyze_plate_original_formula.py"
runpy.run_path(str(SCRIPT), run_name="__main__")

