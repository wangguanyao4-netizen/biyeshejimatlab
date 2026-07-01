"""Regenerate internal-sample-size coverage results.

Paper mapping: Section 4.2 and Section 4.4, centered pointwise coverage for
n = 192, 384, 768.
"""

from pathlib import Path
import runpy

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "dependencies" / "original_python" / "docs__update_n_coverage_revision_20260630.py"
runpy.run_path(str(SCRIPT), run_name="__main__")

