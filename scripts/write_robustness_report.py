"""Generate numerical tables and the manuscript supplement from frozen results."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.validate_robustness_cases import load_samples

OUT = ROOT / "figures/data/robustness"


def main():
    samples = load_samples()
    report = json.loads((OUT / "summary.json").read_text())
    reps = pd.read_csv(OUT / "replicates.csv")
    hp = pd.read_csv(OUT / "high_precision_checks.csv")
    groups = pd.DataFrame(report["groups"])
    groups.to_csv(OUT / "group_summary.csv", index=False)
    primary = samples[samples.scenario == "primary"]
    wider = samples[samples.scenario == "wide"]
    checks = {"n_new_points": len(samples), "n_primary": len(primary), "n_wider": len(wider),
              "n_unstable_all_designs": int((samples.lambda_pattern > 0).sum()),
              "n_gain_mismatch_all_designs": int(((samples.lambda_pattern > 0) != (samples.G_minus3 > 1)).sum()),
              "n_DNA_missed": int(reps.DNA_missed.sum()), "n_DN_missed": int(reps.DN_missed.sum()),
              "n_DNA_not_faster": int(reps.DNA_not_faster.sum()), "n_DN_not_faster": int(reps.DN_not_faster.sum()),
              "n_other_mode_dominant": int(reps.n_other_mode_dominant.sum()),
              "n_other_unstable_K_stable": int(reps.n_other_unstable_K_stable.sum()),
              "n_gain_mismatch_not_two_rhp": int(reps.n_mismatch_not_two_rhp.sum()),
              "n_near_marginal_full": int(reps.n_near_marginal_full.sum()),
              "n_high_precision_distinct_cases": len(hp[["scenario","n","seed","sample_id"]].drop_duplicates()),
              "high_precision_signs_all_agree": bool(hp.signs_agree.all()),
              "primary_pooled_DNA_ratio_median": float((primary.loc[primary.lambda_pattern > 0,"lambda_DNA"] /
                                                          primary.loc[primary.lambda_pattern > 0,"lambda_pattern"]).median()),
              "wider_pooled_DNA_ratio_median": float((wider.loc[wider.lambda_pattern > 0,"lambda_DNA"] /
                                                        wider.loc[wider.lambda_pattern > 0,"lambda_pattern"]).median()),
              "n_uniform_only_full": int(((samples.uniform_re > 0) & (samples.lambda_pattern <= 0)).sum())}
    (OUT / "claim_checks.json").write_text(json.dumps(checks, indent=2) + "\n",encoding="utf-8")
    # This gate does not constrain the scientific computation; it prevents
    # publishing fixed narrative claims if future results contradict them.
    zero_claims = ["n_DNA_missed", "n_DN_missed", "n_DNA_not_faster", "n_DN_not_faster",
                   "n_other_mode_dominant", "n_other_unstable_K_stable",
                   "n_gain_mismatch_not_two_rhp", "n_uniform_only_full"]
    if any(checks[k] for k in zero_claims) or not checks["high_precision_signs_all_agree"]:
        raise RuntimeError("Qualitative conclusions changed; rewrite the supplement narrative before publication")
    main = groups[groups.scenario == "primary"].iloc[0]
    wide = groups[groups.scenario == "wide"].iloc[0]
    lines = [
        '---', 'title: "Sampling robustness of the Notch Delta dynamics comparison"',
        'author: Hyeonje Yang', 'date: October 2026', '---', '',
        '# Supplementary Note 2', '',
        'The main conclusions persist across independent parameter seeds and sample sizes. Their',
        'numerical magnitudes depend on the parameter measure. The original single-sample values',
        '(36.325% unstable; DNA/full median growth ratio 4.90) are retained as reference-sample',
        'results, rather than treated as universal circuit properties.', '',
        '## Design fixed before computation', '',
        'The protocol is deposited as `figures/data/robustness/protocol.json` and specifies 30',
        'separately generated scrambled Latin-hypercube designs, totalling 115000 new points.',
        'The primary analysis uses ten seeds, 2026101701 through 2026101710, each with 4000 points.',
        'Rates and constants span $[0.1,10]\\times p_0$ on a logarithmic scale; Hill coefficients',
        'are uniform on $[1,6]$. Five further designs at each of 1000, 2000 and 8000 points use',
        'the first five of those seeds. These designs are regenerated at each size, not nested',
        'prefixes. Five wider-box designs, seeds 2026101801 through 2026101805, use 4000 points',
        'and rate/constant bounds $[0.01,100]\\times p_0$, with unchanged Hill ranges.', '',
        'Every point is solved by scalar bracketing and scanned over all 54 distinct graph',
        'eigenvalues for the full, DNA and DN models. Actual parameters, equilibria and results',
        'are retained per replicate. No solver failure is labelled equilibrium nonexistence.',
        'The batched calculation reproduces the original 4000 reference growth rates exactly',
        'in this environment; its maximum static-gain difference is below $4\\times10^{-15}$.', '',
        '## Uncertainty and sample size', '',
        'For independent replicate fractions $p_1,\\ldots,p_R$, we report the approximate',
        '95% Monte Carlo interval $\\bar p\\pm t_{R-1,0.975}s_p/\\sqrt R$. The replication unit is',
        'the independently randomised LHS, not an individual point within it. These intervals',
        'describe numerical integration uncertainty under the specified parameter measure,',
        'not experimental or biological variability. Sample-size comparisons reuse seeds;',
        'they are not pooled into one inferential sample. The differing replication counts',
        'also preclude reading the interval widths as an empirical convergence rate.', '',
        f'The primary estimate is **{100*main.mean_unstable_fraction:.2f}%** (approximate 95% Monte Carlo',
        f'interval **{100*main.mc95_low:.2f}–{100*main.mc95_high:.2f}%**). Individual fractions range from',
        f'{100*main.min_unstable_fraction:.3f}% to {100*main.max_unstable_fraction:.3f}%. This is a more stable',
        'estimate of the chosen box measure than the original 36.325% single-sample value.', '',
        '| bounds | points per replicate | replicates | mean unstable (%) | approximate 95% MC interval (%) | DNA median ratio range across replicates |',
        '|---|---:|---:|---:|---|---|']
    for _, r in groups.iterrows():
        label = 'wider' if r.scenario == 'wide' else 'original'
        lines.append(f'| {label} | {int(r.n)} | {int(r.replicates)} | {100*r.mean_unstable_fraction:.2f} | '
                     f'{100*r.mc95_low:.2f}–{100*r.mc95_high:.2f} | '
                     f'{r.DNA_replicate_median_min:.2f}–{r.DNA_replicate_median_max:.2f} |')
    lines += ['', '## Dynamics and classification', '',
        f'Across the 40000 primary points, {int(main.gain_mismatches_total)} instabilities are missed by $G(-3)>1$.',
        f'The DNA reduction has {int(main.DNA_only_total)} reduced-only instabilities; DN has {int(main.DN_only_total)}.',
        'Neither reduction misses an instability of the full model in these samples. Both grow',
        'faster wherever the full model is unstable. The primary pooled DNA/full growth-ratio',
        f'median is {checks["primary_pooled_DNA_ratio_median"]:.2f}; the ten replicate medians range from',
        f'{main.DNA_replicate_median_min:.2f} to {main.DNA_replicate_median_max:.2f}. The corresponding DN range is',
        f'{main.DN_replicate_median_min:.0f}–{main.DN_replicate_median_max:.0f}. Medians are descriptive; a confidence',
        'interval for a mean replicate fraction is not an interval for these conditional medians.', '',
        f'The wider box changes the unstable fraction to {100*wide.mean_unstable_fraction:.2f}% and the DNA',
        f'replicate-median range to {wide.DNA_replicate_median_min:.2f}–{wide.DNA_replicate_median_max:.2f}; the DN range becomes',
        f'{wide.DN_replicate_median_min:.0f}–{wide.DN_replicate_median_max:.0f}. Thus the qualitative bias persists but its magnitude',
        'is not invariant to bounds. No empirical absence of an exception is treated as a theorem.', '',
        f'Over all 115000 analysed points there are {checks["n_unstable_all_designs"]} unstable samples and',
        f'{checks["n_gain_mismatch_all_designs"]} static-gain misses. Every miss has two RHP eigenvalues within',
        'the K-mode block. No faster competing nonuniform mode, full-only instability or reduced',
        'growth rate at or below the full-model rate was found. These totals count evaluations',
        'across separate designs; they are not used as an iid-binomial denominator.', '',
        '## Numerical edge cases and verification', '',
        'All primary and size-check classifications are unchanged at growth thresholds through',
        '$10^{-6}$. One wider-box point (seed 2026101803, ID 2229) has full growth',
        '$4.772659\\times10^{-9}$ and is explicitly flagged as near marginal. Its sign remains',
        'positive in an independent 80-digit calculation. Removing it with a $10^{-8}$ threshold',
        'changes the wider-box pooled fraction by 0.005 percentage points. It does not affect',
        'the primary estimate or qualitative conclusions.', '',
        f'Nyquist winding matches direct RHP counts at all {report["nyquist_checked_points"]} selected points: every static-gain',
        'miss and the 16 closest-to-unit-gain points per replicate, with overlaps removed.',
        'This is a targeted validation, not a claim to have run Nyquist at every new point.',
        'The five fastest complex-pair misses, five fastest real-pair misses and the marginal',
        'point were checked independently at 80-digit precision. All eleven agree in stability',
        'sign and K-mode RHP count. Equilibrium production/loss backward errors are below',
        '$8.2\\times10^{-15}$ throughout the new study.', '',
        '## Table S1', '',
        '| design | n | seed | unstable (%) | DC misses | DNA-only | DN-only | median DNA/full | median DN/full |',
        '|---|---:|---:|---:|---:|---:|---:|---:|---:|']
    for _, r in reps.iterrows():
        lines.append(f'| {r.scenario} | {int(r.n)} | {int(r.seed)} | {100*r.unstable_fraction:.3f} | '
                     f'{int(r.n_gain_mismatch)} | {int(r.DNA_only)} | {int(r.DN_only)} | '
                     f'{r.DNA_ratio_median:.3f} | {r.DN_ratio_median:.1f} |')
    lines += ['', '![**Figure S4. Robustness to sampling seed, sample size and parameter bounds.** '
              '(A) Nonuniform-instability fractions in the original and wider boxes. Each coloured point is an independent '
              '4000-point replicate; dark diamonds and bars show means and approximate 95% Monte Carlo intervals. '
              'Dashed line: original reference sample. (B) Separately generated sample-size checks in the original box. '
              'R is the number of independent replicates at that size; seeds are reused across sizes. '
              '(C) Per-replicate median reduced/full growth ratios among points unstable in both models. '
              '(D) Classification differences, with each coloured mark representing one replicate. Full-only means '
              'the explicit model is unstable but the named reduction is stable.](../figures/supplementary/FigS4_sampling_robustness.png){width=100%}', '']
    path = ROOT / "manuscript/robustness_supplement.md"
    path.write_text("\n".join(lines), encoding="utf-8")
    print(f"wrote {path}")
    print(json.dumps(checks, indent=2))


if __name__ == "__main__":
    main()
