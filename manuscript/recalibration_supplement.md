# Supplementary Note 6. Can a fitted reduced-model clock recover the dynamics?

## Clock calibration and prospective prediction

An obvious control for a faster reduction is to slow its entire clock. We multiply the DNA
vector field by one positive scalar $c$, preserving its equilibrium set. At P2, matching the
leading K-mode growth gives
$c=0.0645183724$.

We hold this scalar fixed across 161 equally spaced $\lambda_N$ values from 1.31 to 6.5,
augmented by exact P2/P3 values (163 conditions). All other parameters stay at baseline.
This is a prediction test for a **single global clock**, not a fit of all reduced-model
parameters or a claim that no three-variable model can be fitted to selected observations.

At P3, the P2-calibrated prediction is 16.51% below the full-model growth rate. A separate
least-squares clock fitted only on $1.5\leq\lambda_N\leq3$ is
$c=0.0718130664$; in the held-out range $3.5\leq\lambda_N\leq6.5$,
its growth errors range from -7.483% to -6.768%.

One clock therefore improves agreement but does not exactly transfer across operating points.

## Spatial spectrum and forced response

We also use the P2 clock to compare leading poles for every spatial mode at P2. Matching one
pole does not match the other poles. At stable P0, we compare the frequency response to a
small additive Delta-production input in mode $\mu=-3$. The full response is
$e_D^\top(i\omega I-J)^{-1}e_D$; the calibrated DNA response is
$e_D^\top(i\omega I-cR)^{-1}c e_D$. Scaling the input along with the clock is essential to
preserve the DC response. The responses share their DC value but differ dynamically.
We evaluate DC and 401 positive frequencies from $10^{-4}$ to $10^3$.

The conclusion concerns the precisely matched reduction and this specified one-parameter
recalibration. Re-estimating nonlinearities or fitting a frequency-dependent effective
memory would be a different model comparison. No experimental parameter-fitting claim is made.

## Reproduction

`scripts/reviewer_linear.py` writes `clock_recalibration.csv`, `clock_mode_predictions.csv`
and `clock_frequency_response.csv`; all are retained, including training and held-out labels
recoverable from the prespecified intervals. Use the report and figure commands in Note 5.

The design was recorded in `figures/data/reviewer_extensions/protocol.json` before the
production calculations. All selected cases and counterexamples are retained. Source hashes
and case checksums are verified by `scripts/write_reviewer_reports.py`. The primary package
versions and numerical reproduction contract are given in `docs/REPRODUCIBILITY.md`.

![**Figure S8. Clock calibration. (A) P2-fitted clock across Notch production. (B) Relative errors for the P2 fit and training-range least-squares fit; shading marks training and held-out ranges. (C) All P2 modes. (D) Stable-P0 Delta-production response; both input and vector field are clock-scaled.**](../figures/supplementary/FigS8_clock_recalibration.png){width=100%}
