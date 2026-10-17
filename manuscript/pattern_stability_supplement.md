# Supplementary Note 9. Local stability of selected nonlinear patterns

## Selection, refinement and complete spectra

A small ODE residual alone does not establish stability. We therefore select twelve distinct
24-by-24 periodic endpoints before calculating their spectra: the first-seed, amplitude
$10^{-4}$ full/DNA endpoint at P2/P3 (four); the largest recorded drift among remaining
main endpoints in each point/model group (four); the P2 free-S endpoint with largest K power
and P3 free-S endpoint with largest drift (two); and the two pure-mode P2 endpoints in main
Fig. 5 (two). Selection is deterministic and is stored in `pattern_selection.json`.

Only Delta fields were stored for the old pure-mode runs. We reconstruct their other fields
from the stationary equations for $N,A,S$ before refinement; the subsequent full-residual
test checks the reconstruction. For every endpoint, sparse Newton refinement with a positive
line search targets $\max|f_i|/x_i^*\leq10^{-11}$. We record the starting fields, every
residual and the change, rather than silently replacing a trajectory endpoint with a
different branch. The largest weighted relative refinement change is
$2.49e-10$, and the maximum final scaled residual is $5.49e-12$.


At each refined state we form the **complete patterned-state Jacobian**, with all local
derivatives evaluated at the heterogeneous concentrations and the original contact matrix.
No homogeneous mode decomposition is assumed here. Dense eigensolvers compute all 2304 full
or 1728 DNA poles. A separate shifted sparse solve recovers a leading eigenpair and checks
its normalised residual. The worst dense/sparse leading-pole difference is
$2.94e-14$ and the largest normalised sparse residual is $1.08e-17$.


All twelve spectra are strictly in the left half-plane: leading real parts range from
-0.10000000 to -0.09765733. These selected equilibria are therefore
locally asymptotically stable for their respective deterministic lattice equations.

This establishes local stability of these states, not a global attractor theorem or stability
of all 112 previously simulated endpoints.

## Finite perturbation and paired continuation

Each refined state receives a centred, max-normalised, relative $10^{-4}$ independent
species perturbation, using PCG64 seed 2026102201 plus its zero-based selection index.
Both the perturbed and **unperturbed** state are independently continued for 200 time units
using analytic-Jacobian BDF, rtol $10^{-10}$, atol $10^{-13}$, max step 0.5. We measure
the RMS difference between these continuations, weighted by homogeneous concentrations,
so residual baseline drift cannot masquerade as amplification. All perturbations recover;
the largest final/initial distance is
$1.83e-09$. This floor is controlled by floating-point/solver error;
we do not infer a physical decay rate from the final few digits.


**Table S11. Complete-spectrum and recovery results for all twelve selected states.**

| state index | point | model | equations | leading real part | RHP poles | final / initial distance |
|---|---|---|---|---|---|---|
| 0 | P2 | DNA | 1728 | -0.1 | 0 | 8.43862e-10 |
| 1 | P2 | full | 2304 | -0.0977229 | 0 | 1.35936e-09 |
| 2 | P3 | DNA | 1728 | -0.1 | 0 | 3.91081e-10 |
| 3 | P3 | full | 2304 | -0.0999987 | 0 | 1.83038e-09 |
| 4 | P2 | DNA | 1728 | -0.1 | 0 | 6.58076e-10 |
| 5 | P2 | full | 2304 | -0.0977129 | 0 | 1.57322e-09 |
| 6 | P3 | DNA | 1728 | -0.1 | 0 | 2.78792e-10 |
| 7 | P3 | full | 2304 | -0.0999982 | 0 | 1.67257e-09 |
| 8 | P2 | full | 2304 | -0.0976573 | 0 | 1.43434e-09 |
| 9 | P3 | full | 2304 | -0.0999985 | 0 | 1.78682e-09 |
| 10 | P2 | full | 2304 | -0.0990174 | 0 | 1.70987e-09 |
| 11 | P2 | DNA | 1728 | -0.1 | 0 | 6.53088e-10 |

## Reproduction

Run `python scripts/reviewer_pattern_stability.py --reuse` (omit `--reuse` to recompute).
All complete complex spectra, initial/refined fields and recovery curves are retained in
`figures/data/reviewer_extensions/pattern_stability/`. The selection table and
`pattern_stability.csv` expose every case rather than only the most stable example.

The design was recorded in `figures/data/reviewer_extensions/protocol.json` before the
production calculations. All selected cases and counterexamples are retained. Source hashes
and case checksums are verified by `scripts/write_reviewer_reports.py`. The primary package
versions and numerical reproduction contract are given in `docs/REPRODUCIBILITY.md`.

![**Figure S11. Patterned-state stability. (A) Leading poles from complete spectra, with table-order endpoint indices. (B) Paired perturbation recovery; blue full, orange DNA. (C) All poles of the pure-mode full P2 state; tiny imaginary parts are at numerical round-off scale, not evidence of oscillation. (D) Newton and independent sparse-eigenpair residuals.**](../figures/supplementary/FigS11_patterned_state_stability.png){width=100%}
