# Supplementary material

The current manuscript sources are Markdown. Existing manuscript DOCX/PDF binaries predate
this revision and require rebuilding and visual inspection before submission. Main Figure 6
promotes the additional results on QSS accuracy, clock transfer, oscillatory entry and
nonlinear discrimination of linearly equivalent memory locations. Figure 1A stays blank
for the author's circuit illustration.

| Note | Question | Current source | Figures | Tables |
|---|---|---|---|---|
| 1 | Exact stationary identities and proofs | [THEORY.md](THEORY.md) | — | — |
| 2 | Sampling seed, sample size, parameter bounds | [Sampling](../manuscript/robustness_supplement.md) | S4 | S1 |
| 3 | Initial fields, competing modes and timing | [Initial conditions](../manuscript/initial_condition_supplement.md) | S5 | S2–S3 |
| 4 | Internal time steps, output cadence, solver accuracy | [Numerical audit](../manuscript/timestep_supplement.md) | S6 | S4–S6 |
| 5 | When is QSS growth accurate? | [QSS validity](../manuscript/qss_validity_supplement.md) | S7 | S7 |
| 6 | Does a calibrated reduced-model clock transfer? | [Recalibration](../manuscript/recalibration_supplement.md) | S8 | — |
| 7 | Lattice size, open edges, static heterogeneity | [Geometry](../manuscript/geometry_supplement.md) | S9 | S8 |
| 8 | How do gain-below-one instabilities arise? | [Oscillatory routes](../manuscript/bifurcation_supplement.md) | S10 | S9–S10 |
| 9 | Are selected nonlinear patterns stable? | [Pattern stability](../manuscript/pattern_stability_supplement.md) | S11 | S11 |
| 10 | Small growth denominators and alternative priors | [Growth and priors](../manuscript/growth_prior_supplement.md) | S12 | S12–S13 |
| 11 | Does linear evidence identify proneural memory? | [Mechanism](../manuscript/mechanism_supplement.md) | S13 | — |

Figures S1–S3 are described in the main source: the original versus corrected eigenvector
initialisation, one-parameter continuations, and subdominant unstable spatial modes.
Every plotted scientific figure has PNG and PDF versions; Figure 1A is intentionally blank.

The frozen additional protocol is
[`protocol.json`](../figures/data/reviewer_extensions/protocol.json). Data include all 52
geometry trajectories, all 12 selected complete patterned spectra, the paired recovery
curves, 52 QSS families, 12 constant-gain continuations and prior-weighted summaries.
Run everything with:

```bash
python scripts/reproduce_all.py --reviewer-extensions
```

For checksum-verified reuse of the expensive trajectories, run the stages separately:

```bash
python scripts/reviewer_linear.py
python scripts/reviewer_geometry.py --reuse
python scripts/reviewer_pattern_stability.py --reuse
python scripts/write_reviewer_reports.py
python figures/plot_reviewer_extensions.py
```

Removing `--reuse` recomputes trajectories. The linear stage recomputes directly from
verified stored parameter inputs. All modes are continuous-time spectra: no ODE time step
enters their stationary thresholds or instability fractions. The release provenance lists
kernel and artifact hashes. Finite sampling, local stability and linear observational
equivalence are distinguished from global theorems and biological identification throughout.
