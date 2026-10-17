"""Regenerate reference analysis, sampling robustness, audits and figures.

    python scripts/reproduce_all.py            # data + figures
    python scripts/reproduce_all.py --figures  # figures only, from the data already in figures/data/
    python scripts/reproduce_all.py --robustness  # replicated sampling, tables and Figure S4
    python scripts/reproduce_all.py --initial-conditions  # paired random starts and Figure S5
    python scripts/reproduce_all.py --dt-audit  # step refinement and Figure S6
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.reproducibility import finish_manifest, new_manifest, write_manifest
DATA_STEPS = ["scripts/linear_analysis.py", "scripts/global_sample.py", "scripts/nyquist_examples.py",
              "scripts/kinetic_control.py", "scripts/timescale_family.py", "scripts/growth_validation.py",
              "scripts/mode_dominance_search.py"]
FIGURE_STEPS = ["figures/plot_fig1.py", "figures/plot_fig2.py", "figures/plot_fig3.py", "figures/plot_fig4.py",
                "figures/plot_fig5.py", "figures/plot_supplementary.py", "figures/plot_robustness.py",
                "figures/plot_initial_conditions.py", "figures/plot_timestep_audit.py",
                "figures/plot_reviewer_extensions.py"]
AUDIT_STEPS = ["scripts/audit_revision.py", "scripts/numerical_robustness.py"]
ROBUSTNESS_STEPS = ["scripts/robustness_study.py", "scripts/validate_robustness_cases.py", "scripts/write_robustness_report.py"]
INITIAL_STEPS = ["scripts/initial_condition_robustness.py", "scripts/mechanistic_predictions.py",
                 "scripts/write_initial_condition_report.py"]
TIMESTEP_STEPS = ["scripts/timestep_audit.py", "scripts/write_timestep_report.py"]
REVIEWER_STEPS = ["scripts/reviewer_linear.py", "scripts/reviewer_geometry.py",
                  "scripts/reviewer_pattern_stability.py", "scripts/write_reviewer_reports.py"]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--figures", action="store_true", help="regenerate figures from stored data")
    group.add_argument("--audit", action="store_true", help="audit stored data and numerical settings")
    group.add_argument("--robustness", action="store_true", help="run independent-seed, sample-size and wider-box robustness study")
    group.add_argument("--initial-conditions", action="store_true", help="run paired random-start lattice experiments and mechanistic predictions")
    group.add_argument("--dt-audit", action="store_true", help="audit time steps, tolerances, solvers and output sampling")
    group.add_argument("--reviewer-extensions", action="store_true", help="run seven additional matched controls and their paper supplements")
    parser.add_argument("--reuse", action="store_true", help="with --reviewer-extensions, reuse checksum-verified geometry and patterned-state cases")
    args = parser.parse_args()
    if args.reuse and not args.reviewer_extensions:
        parser.error("--reuse is available here only with --reviewer-extensions; other stages expose their own reuse flags")
    steps = (FIGURE_STEPS if args.figures else AUDIT_STEPS if args.audit else
             ROBUSTNESS_STEPS + ["figures/plot_robustness.py"] if args.robustness else
             INITIAL_STEPS + ["figures/plot_initial_conditions.py"] if args.initial_conditions else
             TIMESTEP_STEPS + ["figures/plot_timestep_audit.py"] if args.dt_audit else
             REVIEWER_STEPS + ["figures/plot_reviewer_extensions.py"] if args.reviewer_extensions else
             DATA_STEPS + AUDIT_STEPS + ROBUSTNESS_STEPS + INITIAL_STEPS + TIMESTEP_STEPS + REVIEWER_STEPS + FIGURE_STEPS)
    name = ("figure_manifest.json" if args.figures else "audit_manifest.json" if args.audit else
            "robustness_manifest.json" if args.robustness else
            "initial_conditions_manifest.json" if args.initial_conditions else
            "timestep_manifest.json" if args.dt_audit else
            "reviewer_manifest.json" if args.reviewer_extensions else "reproduction_manifest.json")
    manifest = new_manifest(steps)
    if args.reviewer_extensions:
        manifest["reuse_verified_geometry_and_pattern_cases"] = args.reuse
    write_manifest(manifest, name)
    try:
        for step in steps:
            print(f"\n=== {step}", flush=True)
            command = [sys.executable, str(ROOT / step)]
            if args.reuse and step in ("scripts/reviewer_geometry.py", "scripts/reviewer_pattern_stability.py"):
                command.append("--reuse")
            subprocess.run(command, check=True, cwd=ROOT)
    except BaseException as exc:
        manifest["error"] = str(exc)
        finish_manifest(manifest, "failed")
        write_manifest(manifest, name)
        raise
    finish_manifest(manifest, "complete")
    write_manifest(manifest, name)
    print(f"recorded figures/data/{name}")


if __name__ == "__main__":
    main()
