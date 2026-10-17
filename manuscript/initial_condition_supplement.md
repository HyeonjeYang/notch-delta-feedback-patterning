# Supplementary Note 3. Robustness to initial conditions

This experiment tests a timing comparison and the scope of the pure-mode endpoint example, rather than re-establishing that lateral inhibition can generate a pattern. The production design uses two operating points, P2 and P3, eight seeds (2026101901-2026101908), and relative amplitudes 1e-6, 1e-4 and 1e-2. Four independent uniform random fields are centred and each divided by its maximum absolute value. The same fields are reused across amplitudes and operating points, so these comparisons are paired. No noise is applied during integration. These are bounded local perturbations, not a global basin survey.

For each main pair, D, N and A start from exactly the same fields x*(1+epsilon*xi). The full model starts with S=S_q(A); this is also the initial S reconstructed by the DNA model. The two dynamics thus agree on every physical concentration initially. A separate full-model control randomises S independently at amplitude 1e-4, keeping D, N and A matched to the main experiment. There are 96 main trajectories (48 pairs) and 16 downstream-initialisation controls. No clipping is applied.

The design and endpoint rules were fixed before production. A DOP853 pilot completed one full-model trajectory; a stiff DNA pilot was interrupted before writing a result. For computational efficiency the production solver was changed to sparse BDF, retaining the design, thresholds and horizons. Both protocol versions are deposited. BDF uses rtol=1e-9, atol=1e-11, graph-based Jacobian sparsity, and recording interval 0.5. All trajectories end at t=600 (P2) or t=300 (P3), regardless of their apparent convergence.

Onset is defined operationally as the first upward continuous event CV(D)=std(D)/mean(D)=0.1; CV=0.5 is a second threshold. An unreached threshold would be censored at the fixed horizon. These times differ from the K-mode-amplitude read-out in Figure 5. Early broad-spectrum evolution is compared with the exact linear matrix propagator for all 576 Fourier modes and all biochemical poles, on the initial prefix for which max|x/x*-1| remains at most 1e-3. No single exponential is fitted to a mixed-mode transient. Initial amplitude 1e-2 is outside this linear validation window.

All 112 production trajectories crossed both CV thresholds and remained positive at every recorded point (minimum 2.19e-10). The full model crossed CV=0.1 later in all 48 paired cases. Table S2 gives paired onset-time ratios and translation-aligned endpoint correlations; ranges describe the eight seeds, not confidence intervals or biological variability.

**Table S2.** Onset-time ratio is t_full/t_DNA at CV(D)=0.1. Correlation is maximised over all periodic translations, with no rotations or reflections.

| point | amplitude | full onset, median [range] | DNA onset, median [range] | paired ratio, median [range] | endpoint correlation, median [range] |
|---|---|---|---|---|---|
| P2 | 1e-06 | 114.696 [112.559, 118.898] | 7.746 [7.632, 7.955] | 14.810 [14.697, 15.101] | 0.753 [0.596, 0.808] |
| P2 | 1e-04 | 67.432 [65.759, 70.164] | 4.707 [4.617, 4.824] | 14.436 [14.229, 14.782] | 0.742 [0.586, 0.838] |
| P2 | 1e-02 | 15.376 [14.522, 16.485] | 1.348 [1.325, 1.392] | 11.333 [10.960, 11.844] | 0.733 [0.584, 0.774] |
| P3 | 1e-06 | 54.885 [54.050, 56.475] | 4.273 [4.207, 4.402] | 12.844 [12.817, 12.895] | 0.927 [0.884, 0.953] |
| P3 | 1e-04 | 32.431 [31.830, 33.383] | 2.539 [2.488, 2.615] | 12.791 [12.751, 12.861] | 0.940 [0.866, 0.951] |
| P3 | 1e-02 | 7.719 [7.348, 8.119] | 0.619 [0.595, 0.648] | 12.432 [12.340, 12.529] | 0.919 [0.870, 0.939] |

The endpoint fields are different in every main pair, even after periodic translation: aligned relative L2 distances span 0.259-0.778. Exact K-mode power fractions span 0.0023-0.8014 across all production endpoints. A fastest-growing K mode therefore does not imply a pure three-sublattice final field after a broad-spectrum start. Multiple competing spatial components yield heterogeneous endpoints. These results restrict the identical-pattern statement in Figure 5 to its symmetry-restricted initial conditions. They are not a proof of nonlinear pattern-selection mechanisms.

At the fixed horizons all 112 trajectories satisfy max|dx/dt|/x* < 1e-7 (largest value 1.64e-08); the largest last-ten-time-unit Delta change divided by D* is 4.46e-08. They are numerically stationary candidates. A small endpoint residual does not prove attractor stability, basin equality or convergence for every initial state.

**Table S3.** Full-model onset change when S is independently randomised, relative to starting on S_q(A), at amplitude 1e-4. D, N and A are identical within each pair.

| point | free-S / QSS-S onset ratio, median [seed range] |
|---|---|
| P2 | 1.172 [1.101, 1.222] |
| P3 | 1.156 [1.100, 1.193] |

On 69 trajectories with a valid early window, the largest relative L2 Delta error against the all-mode propagator is 0.180%. Four representative cases (first seed, amplitude 1e-6, both models at both points) are rerun with BDF tolerances 1e-11/1e-14 and interval 0.25. CV=0.1 onset differences are at most 0.00582 time units and endpoint K-power differences at most 4e-13. Independent DOP853 runs with the tighter tolerances extend to five time units beyond the original CV=0.5 crossing; both onset thresholds agree within 0.00588 time units. The cross-solver runs are not long-horizon endpoint tests. All 112 initial/final fields, trajectories, per-case hashes, seeds, solver checks and the protocol history are deposited under figures/data/initial_conditions/.

![**Figure S5. Paired random initial conditions preserve the timing difference and change endpoint structure.** (A,B) First CV(D)=0.1 crossing at P2 and P3 versus perturbation amplitude. Dots: eight paired seeds per model; lines: medians. (C) Endpoint Delta correlation, maximised over periodic translations, for all 48 main pairs. (D) Exact K-mode power fraction at the fixed horizon, for the main paired runs. (E,F) Full and DNA Delta fields at P2, amplitude 1e-4, first seed; common colour scale. Fields are shown at t=600, after the stationary-candidate criterion is met. Polygons show the contact topology, not measured cell shape. No noise is applied during any trajectory.](../figures/supplementary/FigS5_initial_condition_robustness.png){width=100%}
