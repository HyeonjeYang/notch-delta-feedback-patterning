"""Generate Supplementary Note 3 from verified paired lattice experiments."""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import numpy as np
import pandas as pd
from scripts.initial_condition_robustness import (OUT, case_name, kernel_hash, make_pairs, specifications)


def verified_results():
    digest = kernel_hash()
    specs = list(specifications())
    summaries=[]
    for spec in specs:
        path = OUT/"cases"/(case_name(spec)+".npz")
        meta = json.loads(path.with_suffix(".json").read_text())
        if meta["kernel_sha256"] != digest or meta["npz_sha256"] != hashlib.sha256(path.read_bytes()).hexdigest():
            raise RuntimeError(f"stale or modified simulation: {path}")
        summaries.append(meta["summary"])
    frame = pd.read_csv(OUT/"runs.csv")
    if len(frame) != len(specs) or frame.duplicated(["point","amplitude","seed","model","downstream"]).any():
        raise RuntimeError("incomplete or duplicated run table")
    pd.testing.assert_frame_equal(frame,pd.DataFrame(summaries),check_dtype=False,rtol=1e-12,atol=1e-15)
    pairs=pd.read_csv(OUT/"pairs.csv")
    pd.testing.assert_frame_equal(pairs,make_pairs(frame),check_dtype=False,rtol=1e-12,atol=1e-15)
    return frame,pairs,json.loads((OUT/"summary.json").read_text())


def main():
    d,p,s = verified_results()
    text = ["# Supplementary Note 3. Robustness to initial conditions\n",
        "This experiment tests a timing comparison and the scope of the pure-mode endpoint example, "
        "rather than re-establishing that lateral inhibition can generate a pattern. The production "
        "design uses two operating points, P2 and P3, eight seeds (2026101901-2026101908), and "
        "relative amplitudes 1e-6, 1e-4 and 1e-2. Four independent uniform random fields are "
        "centred and each divided by its maximum absolute value. The same fields are reused across "
        "amplitudes and operating points, so these comparisons are paired. No noise is applied "
        "during integration. These are bounded local perturbations, not a global basin survey.\n",
        "For each main pair, D, N and A start from exactly the same fields x*(1+epsilon*xi). "
        "The full model starts with S=S_q(A); this is also the initial S reconstructed by the "
        "DNA model. The two dynamics thus agree on every physical concentration initially. "
        "A separate full-model control randomises S independently at amplitude 1e-4, keeping "
        "D, N and A matched to the main experiment. There are 96 main trajectories (48 pairs) "
        "and 16 downstream-initialisation controls. No clipping is applied.\n",
        "The design and endpoint rules were fixed before production. A DOP853 pilot completed "
        "one full-model trajectory; a stiff DNA pilot was interrupted before writing a result. "
        "For computational efficiency the production solver was changed to sparse BDF, retaining "
        "the design, thresholds and horizons. Both protocol versions are deposited. BDF uses "
        "rtol=1e-9, atol=1e-11, graph-based Jacobian sparsity, and recording interval 0.5. "
        "All trajectories end at t=600 (P2) or t=300 (P3), regardless of their apparent convergence.\n",
        "Onset is defined operationally as the first upward continuous event CV(D)=std(D)/mean(D)=0.1; "
        "CV=0.5 is a second threshold. An unreached threshold would be censored at the fixed "
        "horizon. These times differ from the K-mode-amplitude read-out in Figure 5. "
        "Early broad-spectrum evolution is compared with the exact linear matrix propagator "
        "for all 576 Fourier modes and all biochemical poles, on the initial prefix for which "
        "max|x/x*-1| remains at most 1e-3. No single exponential is fitted to a mixed-mode transient. "
        "Initial amplitude 1e-2 is outside this linear validation window.\n",
        f"All {len(d)} production trajectories crossed both CV thresholds and remained positive "
        f"at every recorded point (minimum {s['min_concentration']:.3g}). The full model crossed "
        f"CV=0.1 later in all {s['n_full_later']} paired cases. Table S2 gives paired onset-time "
        "ratios and translation-aligned endpoint correlations; ranges describe the eight seeds, "
        "not confidence intervals or biological variability.\n",
        "**Table S2.** Onset-time ratio is t_full/t_DNA at CV(D)=0.1. Correlation is maximised "
        "over all periodic translations, with no rotations or reflections.\n",
        "| point | amplitude | full onset, median [range] | DNA onset, median [range] | paired ratio, median [range] | endpoint correlation, median [range] |",
        "|---|---|---|---|---|---|"]
    for g in s["groups"]:
        def fmt(key):
            v=g[key]
            return f"{v['median']:.3f} [{v['min']:.3f}, {v['max']:.3f}]"
        text.append(f"| {g['point']} | {g['amplitude']:.0e} | {fmt('full_onset_01')} | {fmt('DNA_onset_01')} | {fmt('onset_ratio')} | {fmt('aligned_correlation')} |")
    frac=d.final_k_fraction
    text.extend(["\nThe endpoint fields are different in every main pair, even after periodic "
        f"translation: aligned relative L2 distances span {p.aligned_relative_L2.min():.3f}-"
        f"{p.aligned_relative_L2.max():.3f}. Exact K-mode power fractions span "
        f"{frac.min():.4f}-{frac.max():.4f} across all production endpoints. A fastest-growing "
        "K mode therefore does not imply a pure three-sublattice final field after a broad-spectrum "
        "start. Multiple competing spatial components yield heterogeneous endpoints. These "
        "results restrict the identical-pattern statement in Figure 5 to its symmetry-restricted "
        "initial conditions. They are not a proof of nonlinear pattern-selection mechanisms.\n",
        f"At the fixed horizons all {int(d.stationary_candidate.sum())} trajectories satisfy "
        f"max|dx/dt|/x* < 1e-7 (largest value {d.final_scaled_rhs.max():.3g}); the largest "
        f"last-ten-time-unit Delta change divided by D* is {d.last10_D_change.max():.3g}. "
        "They are numerically stationary candidates. A small endpoint residual does not prove "
        "attractor stability, basin equality or convergence for every initial state.\n",
        "**Table S3.** Full-model onset change when S is independently randomised, relative "
        "to starting on S_q(A), at amplitude 1e-4. D, N and A are identical within each pair.\n",
        "| point | free-S / QSS-S onset ratio, median [seed range] |",
        "|---|---|"])
    a=d[(d.model=="full") & (d.amplitude==1e-4)].pivot(index=["point","seed"],columns="downstream",values="onset_01")
    for point,g in (a["free"]/a["qss"]).groupby("point"):
        text.append(f"| {point} | {g.median():.3f} [{g.min():.3f}, {g.max():.3f}] |")
    tight=max(c["onset_abs_difference"] for c in s["tight_checks"])
    cross=max(max(c["onset_01_difference"],c["onset_05_difference"]) for c in s["solver_checks"])
    kend=max(c["k_fraction_abs_difference"] for c in s["tight_checks"])
    text.extend([f"\nOn {int((d.linear_n>=2).sum())} trajectories with a valid early window, "
        f"the largest relative L2 Delta error against the all-mode propagator is "
        f"{100*s['max_linear_relative_D_error']:.3f}%. Four representative cases (first seed, "
        "amplitude 1e-6, both models at both points) are rerun with BDF tolerances "
        "1e-11/1e-14 and interval 0.25. "
        f"CV=0.1 onset differences are at most {tight:.5f} time units and endpoint K-power "
        f"differences at most {kend:.2g}. Independent DOP853 runs with the tighter tolerances "
        "extend to five time units beyond the original CV=0.5 crossing; both onset thresholds "
        f"agree within {cross:.5f} time units. The cross-solver runs are not long-horizon "
        "endpoint tests. All 112 initial/final fields, trajectories, per-case hashes, seeds, "
        "solver checks and the protocol history are deposited under figures/data/initial_conditions/.\n",
        "![**Figure S5. Paired random initial conditions preserve the timing difference and change "
        "endpoint structure.** (A,B) First CV(D)=0.1 crossing at P2 and P3 versus perturbation "
        "amplitude. Dots: eight paired seeds per model; lines: medians. (C) Endpoint Delta "
        "correlation, maximised over periodic translations, for all 48 main pairs. "
        "(D) Exact K-mode power fraction at the fixed horizon, for the main paired runs. "
        "(E,F) Full and DNA Delta fields at P2, amplitude 1e-4, first seed; common colour scale. "
        "Fields are shown at t=600, after the stationary-candidate criterion is met. "
        "Polygons show the contact topology, not measured cell shape. No noise is applied "
        "during any trajectory.](../figures/supplementary/FigS5_initial_condition_robustness.png){width=100%}\n"])
    target=ROOT/"manuscript/initial_condition_supplement.md"
    target.write_text("\n".join(text)+"\n",encoding="utf-8")
    print(f"wrote {target}")


if __name__ == "__main__":
    main()
