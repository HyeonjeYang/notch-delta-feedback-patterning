"""Build the revised manuscript and the response letter as .docx and .pdf with pandoc.

    python manuscript/build.py

Requires pandoc on PATH. The .pdf outputs also need a LaTeX engine (xelatex). Figures are
read from ../figures/ and should be regenerated first (scripts/reproduce_all.py --figures).
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent
SOURCES = ["manuscript_revised.md", "response_to_reviewers.md", "robustness_supplement.md",
           "initial_condition_supplement.md", "timestep_supplement.md", "qss_validity_supplement.md",
           "recalibration_supplement.md", "geometry_supplement.md", "bifurcation_supplement.md",
           "pattern_stability_supplement.md", "growth_prior_supplement.md", "mechanism_supplement.md"]


def main() -> None:
    if shutil.which("pandoc") is None:
        raise SystemExit("pandoc not found on PATH")
    for name in SOURCES:
        src = HERE / name
        stem = src.with_suffix("")
        subprocess.run(["pandoc", str(src), "-o", f"{stem}.docx", "--resource-path", str(HERE)],
                       check=True, cwd=HERE)
        print(f"wrote {stem}.docx")
        if shutil.which("xelatex"):
            subprocess.run(["pandoc", str(src), "-o", f"{stem}.pdf", "--pdf-engine=xelatex",
                            "--resource-path", str(HERE), "-V", "mainfont=Cambria", "-V", "mathfont=Cambria Math"],
                           check=True, cwd=HERE)
            print(f"wrote {stem}.pdf")


if __name__ == "__main__":
    main()
