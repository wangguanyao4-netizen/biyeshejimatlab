"""Build the paper numerical figures.

Paper mapping: Figures 1--11 in the numerical section.
This wrapper keeps a paper-readable name while executing the original script.
"""

from pathlib import Path
import runpy

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "dependencies" / "original_python" / "build_current_paper_numerical_figures.py"
runpy.run_path(str(SCRIPT), run_name="__main__")

