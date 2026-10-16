"""Regenerate every data file and figure in dependency order (about 15 minutes).

    python scripts/reproduce_all.py            # data + figures
    python scripts/reproduce_all.py --figures  # figures only, from the data already in figures/data/
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA_STEPS = ["scripts/linear_analysis.py", "scripts/global_sample.py", "scripts/nyquist_examples.py",
              "scripts/kinetic_control.py", "scripts/timescale_family.py", "scripts/growth_validation.py",
              "scripts/mode_dominance_search.py"]
FIGURE_STEPS = ["figures/plot_fig1.py", "figures/plot_fig2.py", "figures/plot_fig3.py", "figures/plot_fig4.py",
                "figures/plot_fig5.py", "figures/plot_supplementary.py"]


def main() -> None:
    steps = FIGURE_STEPS if "--figures" in sys.argv else DATA_STEPS + FIGURE_STEPS
    for step in steps:
        print(f"\n=== {step}", flush=True)
        subprocess.run([sys.executable, str(ROOT / step)], check=True, cwd=ROOT)


if __name__ == "__main__":
    main()
