# Supplementary Note 4. Time-discretisation and solver audit

The reported stationary threshold, mode spectra and instability fractions are computed from scalar equilibrium bracketing and the continuous-time Jacobian, not from stepping a trajectory forward. A temporal integration step therefore cannot create those particular classifications. Floating-point eigensolver errors remain a separate issue: selected cases were previously checked at 80-digit precision (Supplementary Note 2). This note tests whether time-discrete dynamics reproduce the spectra and whether the nonlinear growth, emergence times and endpoint fields depend on solver choices.

**What dt means.** In the reference DOP853 runs and random-start BDF runs, the recording interval is an output sampling interval. It does not set the integrator's internal step. SciPy chooses adaptive internal steps using rtol and atol, while max_step imposes an upper bound. Event times are found from the solver's continuous interpolant. These distinctions follow the [SciPy solve_ivp documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html); the installed version used for the audit is 1.16.3. Four checks are kept separate below.

## 1. Genuine fixed-step refinement

Classical RK4 is implemented independently of SciPy's adaptive integrators, without clipping. The full 24x24 lattice, rather than a three-cell symmetry reduction, is integrated with fixed dt=0.08, 0.04, 0.02 and 0.01 at P2/P3 for full/DNA models (16 trajectories). Each starts on its own exact K eigenvector at Delta amplitude 1e-6 D*. Solutions are compared on the common 0.08 output grid to Radau at rtol=1e-12, atol=1e-15 and max_step=0.01. The fitting interval uses the same times for every integrator, ending at the last common time no greater than log(1000)/the predicted eigenvalue. Thus a shift of a sampled threshold cannot masquerade as a time-step effect. Finite-amplitude terms still cause a small growth-slope discrepancy from the exactly linear eigenvalue.

**Table S4.** Fixed-step nonlinear lattice checks. Delta-field errors are relative to the growing Delta perturbation, not to the much larger background concentration.

| point | model | dt | fitted growth | relative error vs eigenvalue | max relative Delta perturbation error vs Radau |
|---|---|---|---|---|---|
| P2 | full | 0.08 | 0.103077072 | 0.0001 | 2e-09 |
| P2 | full | 0.04 | 0.103077072 | 0.0001 | 2.49e-09 |
| P2 | full | 0.02 | 0.103077072 | 0.0001 | 6.12e-09 |
| P2 | full | 0.01 | 0.103077072 | 0.0001 | 1.76e-09 |
| P2 | DNA | 0.08 | 1.597565775 | 5.46e-05 | 1.53e-05 |
| P2 | DNA | 0.04 | 1.597568758 | 5.64e-05 | 1.01e-06 |
| P2 | DNA | 0.02 | 1.597568955 | 5.66e-05 | 6.42e-08 |
| P2 | DNA | 0.01 | 1.597568967 | 5.66e-05 | 2.56e-09 |
| P3 | full | 0.08 | 0.216315375 | 7.02e-05 | 7.07e-08 |
| P3 | full | 0.04 | 0.216315375 | 7.02e-05 | 2.62e-09 |
| P3 | full | 0.02 | 0.216315375 | 7.02e-05 | 2.29e-09 |
| P3 | full | 0.01 | 0.216315375 | 7.02e-05 | 4.51e-10 |
| P3 | DNA | 0.08 | 2.799099189 | 2.97e-05 | 0.000125 |
| P3 | DNA | 0.04 | 2.799143731 | 4.57e-05 | 8.87e-06 |
| P3 | DNA | 0.02 | 2.799146828 | 4.68e-05 | 5.83e-07 |
| P3 | DNA | 0.01 | 2.799147031 | 4.68e-05 | 3.72e-08 |

At dt=0.01, the largest Delta perturbation error is 3.72e-08; the largest fitted-growth error against the continuous eigenvalue is 0.0001. These are early-growth checks and do not validate fixed-step RK4 for the stiff late nonlinear regime.

An additional post-audit diagnostic reuses the stored high-accuracy reference amplitudes, without new trajectories or changing the prespecified assessments. Reducing the upper amplitude of the common theoretical window from 1e-3 to 1e-4 and 1e-5 reduces the largest relative growth-slope discrepancy from 0.0001 to 1.84e-05 and 4.05e-06. This supports finite-amplitude bias as the explanation for the residual slope difference once temporal error is small. The reference manuscript's reported fits and original window rule are retained; their last digits should not be interpreted as experimentally calibrated or exact accuracy.

A separate linear test uses every one of the 54 distinct graph eigenvalues at P0/P2/P3 for full and DNA models. At each fixed dt the RK4 stability polynomial R(hJ)=I+hJ+(hJ)^2/2+(hJ)^3/6+(hJ)^4/24 is compared with the exact matrix exponential at t=2. The number of unstable discrete eigenvalues is compared with the continuous RHP count; graph modes are not inferred from a single spatial start. There are 1296 mode/model/step comparisons and 216 count disagreements. This comparison is a check on the chosen steps at these operating points, not a statement that any explicit method at arbitrarily large dt preserves stability.

The 216 disagreements are all at the continuously stable P0: h=0.04 and 0.08 each make all 54 mode blocks unstable in both models. These are genuine numerical instabilities of the coarse RK4 discretisation. All 648 comparisons at h=0.02 and 0.01 agree with the continuous RHP counts. The largest fine-step propagator errors are 6.28e-07 and 3.84e-08, respectively. The failed coarse controls were not used in the adaptive production calculations.

## 2. Internal step caps, tolerance refinement and an independent stiff solver

Fourteen fixed initial states cover the first seed at every amplitude and both models/operating points (12 main cases), plus the two independent-S controls at amplitude 1e-4. For each case the existing BDF trajectory (rtol=1e-9, atol=1e-11, unbounded max_step) is compared with three step caps, 0.5, 0.25 and 0.125; BDF with cap 0.125 and tighter tolerances 1e-11/1e-14; and Radau with cap 0.25 and the same tighter tolerances. All 70 new trajectories run to the original fixed horizon, t=600 at P2 and t=300 at P3. Accepted internal step meshes are deposited. Cellwise endpoint differences use the original cell indices, without translation alignment, making this a stricter comparison than the full/DNA pattern comparison.

**Table S5.** Largest deviations across the five new configurations for each fixed initial state. Event change takes the maximum over CV(D)=0.1 and 0.5.

| point | model | initial S | amplitude | max relative event change | max relative endpoint L2, all species | max scaled endpoint drift |
|---|---|---|---|---|---|---|
| P2 | DNA | qss | 1e-06 | 6.29e-05 | 3.77e-12 | 2.08e-13 |
| P2 | DNA | qss | 1e-04 | 7.22e-07 | 1.65e-12 | 2.08e-13 |
| P2 | DNA | qss | 1e-02 | 2.25e-08 | 1.14e-12 | 2.08e-13 |
| P2 | full | free | 1e-04 | 1.98e-06 | 1.51e-13 | 2.08e-13 |
| P2 | full | qss | 1e-06 | 5.15e-05 | 5.11e-12 | 2.08e-13 |
| P2 | full | qss | 1e-04 | 8.32e-07 | 2.69e-13 | 2.08e-13 |
| P2 | full | qss | 1e-02 | 2.52e-08 | 1.92e-13 | 2.08e-13 |
| P3 | DNA | qss | 1e-06 | 0.000163 | 9.23e-12 | 2.78e-12 |
| P3 | DNA | qss | 1e-04 | 1.35e-06 | 2.34e-11 | 2.28e-12 |
| P3 | DNA | qss | 1e-02 | 2.8e-08 | 7.82e-12 | 1.88e-12 |
| P3 | full | free | 1e-04 | 1.56e-06 | 1.15e-11 | 6.01e-10 |
| P3 | full | qss | 1e-06 | 5.41e-05 | 1.11e-10 | 5.48e-09 |
| P3 | full | qss | 1e-04 | 6.52e-07 | 5.2e-12 | 5.7e-10 |
| P3 | full | qss | 1e-02 | 2.46e-08 | 1.27e-11 | 5.29e-11 |

Across these cases the largest relative event-time change is 0.01631%, the largest all-species endpoint relative L2 difference is 1.11e-10, and the largest absolute K-power-fraction change is 4.02e-13. The largest endpoint scaled drift is 5.48e-09. Endpoint agreement under refinement is numerical evidence of reproducibility for these starts; it does not prove attractor stability or equality of basins.

## 3. Output sampling without changing the integration

The identical dense Radau growth solution is resampled at recording intervals 0.5, 0.25 and 0.125. The original first-recorded-amplitude rule is applied separately to each sampled curve. Its endpoint may change with cadence, although the underlying numerical trajectory and its internal steps are identical.

**Table S6.** Growth fits under the original output-based window rule.

| point | model | output interval | window end | fitted growth |
|---|---|---|---|---|
| P2 | full | 0.5 | 67.5 | 0.103077650 |
| P2 | full | 0.25 | 67.25 | 0.103077363 |
| P2 | full | 0.125 | 67.125 | 0.103077222 |
| P2 | DNA | 0.5 | 4.5 | 1.597610723 |
| P2 | DNA | 0.25 | 4.5 | 1.597599461 |
| P2 | DNA | 0.125 | 4.375 | 1.597577098 |
| P3 | full | 0.5 | 32 | 0.216315949 |
| P3 | full | 0.25 | 32 | 0.216315729 |
| P3 | full | 0.125 | 32 | 0.216315618 |
| P3 | DNA | 0.5 | 2.5 | 2.799222287 |
| P3 | DNA | 0.25 | 2.5 | 2.799197517 |
| P3 | DNA | 0.125 | 2.5 | 2.799183755 |

The largest relative slope change from the finest recording interval is 2.1e-05. The output interval must therefore be reported for fitted slopes; it should not be described as the adaptive solver dt.

## 4. Prespecified assessments and limits

The protocol and numerical assessment thresholds were written before the new runs. Every assessment is reported below, including a failure if a coarse step exceeds the chosen accuracy target. A failed accuracy target need not mean an artificial instability; convergence at finer steps must be assessed separately.

| assessment | within prespecified target |
|---|---|
| adaptive_events | yes |
| adaptive_fields | yes |
| adaptive_k_fraction | yes |
| adaptive_stationarity | yes |
| fixed_growth | yes |
| fixed_fields | no: inspect the retained step-level results |
| output_cadence | yes |
| positive_to_threshold | yes |
| spectral_classifications | no: inspect the retained step-level results |

The combined all_assessments_pass flag is false because the sweep includes failed coarse RK4 controls. All four adaptive production-comparison targets pass. The fixed-field target is exceeded only by P3 DNA at h=0.08 (relative perturbation error 0.000125); it passes at every smaller tested step. The spectral failures are the coarse P0 cases described above. No failed control was discarded and no accuracy target was relaxed. This distinction is also recorded in assessment_interpretation.json. Across the 30 reintegrated matched model pairs, the full model remains later at CV=0.1 in all 30 comparisons.


The minimum accepted adaptive concentration or RK4 stage concentration is 2.46e-10. No positive concentration floor or post-step clipping is used. The fixed-step implementation is independently tested for fourth-order convergence against the nonlinear exact solution y'=y^2, y(0)=1. Tests supplement, and do not replace, the lattice comparisons.

This audit covers fourteen representative lattice initial states, not every seed or parameter vector. The 24x24 contact graph is a discrete biological modelling choice; this is not a spatial mesh-convergence study or a thermodynamic-limit result. Agreement among time integrators does not establish tissue validity, remove model misspecification or mathematically exclude all floating-point error. The deposited data support only the stated resolution and solver comparisons.

![**Figure S6. Time-step and solver checks.** (A) Number of mode blocks with an incorrect discrete RHP count, pooled over P0/P2/P3 and both models, versus fixed RK4 step. All failures occur at P0. (B) Maximum relative Delta perturbation error versus the high-accuracy Radau reference in the nonlinear early-growth checks. (C) Relative event-time changes from the original BDF run, taking the larger of the CV=0.1 and 0.5 differences; each point is one of fourteen matched initial states. (D) Raw all-species endpoint relative L2 difference from the original run, without spatial alignment, at the fixed horizons. Panel C uses a symmetric logarithmic scale with a linear region at zero. Horizontal dotted lines show prespecified accuracy targets, not biological confidence intervals.](../figures/supplementary/FigS6_timestep_audit.png){width=100%}
