# Computational reproducibility

The reference numerical environment observed during the 17 October 2026 audit is CPython
3.13.5 (Anaconda build), Windows x86_64, NumPy 2.2.6 and SciPy 1.16.3. `requirements.txt`
now matches the installed direct dependencies; the earlier version listed different versions
that were not the ones used for this audit. `requirements-lock.txt` additionally pins the
installed dependency closure. These are version pins, not a binary or container lock, and
installation in a newly created environment has not been verified in this session.

## Run the analysis

Use an isolated Python 3.13.5 environment, install `requirements-lock.txt`, then run:

```text
python -m pytest -q
python scripts/reproduce_all.py
```

`--figures` regenerates plots from saved data. `--audit` recalculates the global sample and
checks numerical robustness without repeating all data-generation steps. An unsupported
command-line argument is rejected rather than silently starting a complete reproduction.

`--robustness` runs the independent-seed, sample-size and wider-box study, high-precision
case validation, paper supplement generator and Figure S4. The full default run includes it.
The protocol was written before examining new results. It specifies ten independent
4000-point LHS replicates, five designs at each of 1000/2000/8000 points, and five wider-box
4000-point replicates, totalling 115000 new evaluations. Full details and seed-level results
are in `manuscript/robustness_supplement.md`.

Every replicate retains its inputs/results and checksum under `figures/data/robustness/`.
Cached results may be used explicitly with `python scripts/robustness_study.py --reuse`;
both kernel and file checksums must match. Default execution recalculates every point.
The computational-kernel hash normalises text line endings so a Git LF/CRLF conversion
does not prevent validation on a new checkout. File-byte hashes are also recorded separately.
The robustness script sets OpenBLAS, OMP and MKL thread counts to one before importing
NumPy to avoid oversubscription for small matrices. The top-level manifest's thread fields
describe the runner environment; this child-process override is explicit in the study source.

## Random inputs

The two parameter studies use different sampling designs. The parameter-column order is
`stability.sensitivity.PARAM_NAMES`, also recorded by the CSV headers.

| study | protocol | frozen inputs |
|---|---|---|
| 4000-point global sample | SciPy `qmc.LatinHypercube(d=13, scramble=True, strength=1, optimization=None, seed=20260817)`; one call to `.random(n=4000)` | `figures/data/global_sample_parameters.csv` |
| 40000-point additional mode search | `Generator(PCG64(20261016))`; independent uniform draws transformed into log-uniform rates and uniform Hill exponents; 20000 points in the base box, then 20000 in the wider box using the same continuing stream | `figures/data/mode_search_samples.csv.gz` |

The base bounds are $[0.1,10]\times p_0$ and wider bounds are $[0.01,100]\times p_0$
for rates/constants; Hill exponents stay in $[1,6]$.

Global-sample IDs are 0-based. Mode-search IDs are 0-based within each box and must be read
together with `box_widening`. The second stream is **not reset** between boxes. Reordering
columns, changing the draw shape, inserting random calls or resetting the stream changes the
sample even with the same integer seed.

Keep the legacy QMC `seed=` API when regenerating this existing study. In the validated
environment, replacing it with `rng=` and the same integer produces different points. This
API distinction is documented in the [SciPy LatinHypercube reference](https://docs.scipy.org/doc/scipy-1.17.0/reference/generated/scipy.stats.qmc.LatinHypercube.html).
NumPy also limits its stream-compatibility guarantee to specified generator, call-sequence
and environment conditions; see the [NumPy compatibility policy](https://numpy.org/doc/2.0/reference/random/compatibility.html).

Both generation scripts check regenerated parameters against the snapshots before replacing
results. Tests independently reconstruct both sampling streams. The relative comparison
tolerance is 5e-14, allowing final-bit differences in the log/exp transforms while detecting
a changed sample. The existing snapshots remain the authoritative inputs for the reported
counts; changing a seed creates a new study rather than a reproduction of the original one.

The reference nonlinear runs use deterministic spatial eigenvectors, without random initial conditions,
so they do not need a seed. The paired initial-condition study below uses explicit random
seeds for initial fields, then deterministic dynamics. The earlier stochastic experiment is archived and is not evidence
for the revised manuscript.

## Numerical settings and checks

Homogeneous equilibria use scalar Brent bracketing on `[0, lam_D/d_D]`, `xtol=float.tiny`,
`rtol=4*float.eps`, and at most 500 iterations. Residuals are normalised by
`max(1, ||x*||_infinity)`. Critical-point roots use the bracket `[1.2, 1.4]` and absolute
root tolerance 1e-13. Mode scans use 54 distinct graph eigenvalues, rounded to 12 decimal
places to merge numerical representations of the same eigenvalue.

The reported global classification uses `Re(s)>0`. The closest global-sample growth rate to
zero has magnitude 4.375e-5. The unstable count is unchanged at thresholds 0, 1e-10, 1e-8
and 1e-6. This is a check for this sample, not a prescription for classifying arbitrary
near-marginal parameter sets. A zero/imaginary-axis eigenvalue is a stability boundary.

The Nyquist count uses both frequency signs and zero. The positive grid is
`scale*logspace(-8, 4, 4000)`, where `scale=max(1,max(abs(eig(J_mu))))`. All 4000 K-mode
counts were checked against direct eigenvalues. `scripts/numerical_robustness.py` additionally
checks grids of 2000, 4000 and 8000 positive frequencies for the 12 gain mismatches and the
16 points closest to unit gain. A finite sampled winding is not assumed reliable at an
unresolved or marginal point.

The nonlinear reference runs use DOP853 with `rtol=1e-9`, `atol=1e-11` and recording interval
0.25. The initial Delta-mode amplitude is 1e-6 of its equilibrium for the corrected initial
states. The fitting rule ends at the **first recorded** amplitude at or above 1e-3 of the
equilibrium; the last point may overshoot that threshold. If it is not reached the script
now fails rather than silently fitting the entire trajectory. P2 full, P2 DNA and P3 full
are also checked with `rtol=1e-11`, `atol=1e-14` and recording interval 0.125. The largest
relative slope error against the linear prediction in these tighter runs is about 1.02e-4.
This is numerical agreement, not a confidence interval for a biological growth rate.

Central-difference sensitivity steps are also halved and doubled. The growth-sensitivity
vector changes by at most 0.14% of its largest component; the static vector changes by less
than 1e-7 of its largest component. Details are in `figures/data/numerical_robustness.json`.

## Provenance and limits

The runner writes `reproduction_manifest.json`, `figure_manifest.json`, `audit_manifest.json`
or `robustness_manifest.json` or `initial_conditions_manifest.json`
or `timestep_manifest.json`
under `figures/data/`. It records start/completion status, Python/platform/package versions,
BLAS details, relevant thread environment variables, Git HEAD, source-file SHA-256 hashes,
the global-input hash, random protocols and end-of-run artifact hashes. Git HEAD alone is
insufficient while the working tree contains edits; the source hashes identify the actual
files. A manifest reports whether source files changed during the run. Separate invocation
types have separate manifests so a figures-only run does not masquerade as full reproduction.

Numeric outputs are compared within stated tolerances. Different BLAS/LAPACK builds, platform
arithmetic and adaptive solver paths can change last digits and eigenvector signs. PDF hashes
can also change because of export metadata. Neither equal seeds nor dependency pins imply
byte-identical figures or trajectories across machines.

Independent-seed robustness is now reported separately from numerical reproducibility.
The primary ten-replicate estimate is 37.37% (approximate 95% Monte Carlo interval
37.07–37.67%); its pooled DNA/full growth-ratio median is 5.41. Replicate fractions, not
individual LHS points, form the interval calculation. The 36.325% and 4.9-fold median remain
the results of the original frozen reference sample. Wider bounds alter magnitudes, and
the size checks are not a formal convergence-rate claim. One wider-box growth rate below
1e-8 is flagged and confirmed positive at 80-digit precision. No universal dynamic ordering
is inferred solely from the absence of exceptions in these finite samples.
The Markdown manuscript/response are current; rebuild and visually verify the older checked-in
DOCX files and insert the author's circuit schematic into the intentionally empty Figure 1A
before submission.

## Paired initial-condition robustness

Run `python scripts/reproduce_all.py --initial-conditions` to regenerate the production
trajectories, stricter solver checks, mechanistic predictions, Supplementary Note 3 and
Figure S5. The numerical study supports `python scripts/initial_condition_robustness.py --reuse`
for per-case reuse after verifying its canonical numerical-kernel hash and NPZ SHA-256.
Report/figure generation independently checks the run tables against the case records and
recalculates translation-aligned pair metrics from the fields.

The frozen design is `figures/data/initial_conditions/protocol.json`. Seeds
2026101901–2026101908 initialise `Generator(PCG64(seed))`; four independent uniform
[-1,1] fields are centred separately and divided by their maximum absolute value.
The same random fields are reused across P2/P3 and amplitudes 1e-6, 1e-4 and 1e-2.
Within each pair, D, N and A match exactly and full-model S starts at S_q(A), also the
concentration reconstructed by DNA. Sixteen additional full runs randomise S independently
at amplitude 1e-4. Thus there are 96 main trajectories (48 pairs) and 16 controls, not 112
independent parameter samples. There is no stochastic forcing during integration.

Production uses BDF with graph Jacobian sparsity, rtol=1e-9, atol=1e-11, interval 0.5,
and fixed horizons t=600 (P2) and t=300 (P3). Four worker processes each set BLAS/OpenMP
thread counts to one before NumPy import. Runtime worker count does not alter the per-case
RNG. A DOP853 pilot completed one full P2 run; a slow DNA pilot was interrupted before
writing a result. Sparse BDF was selected before production; the earlier draft protocol
is retained as `protocol_draft_dop853.json`. The design, thresholds and horizons were
unchanged. This is a documented computational protocol, not a formal preregistration.

Onset is the first upward continuous crossing of spatial CV(D)=0.1; CV=0.5 is also
recorded. Unreached events would be censored at the fixed horizon. The single K-mode
amplitude in Figure 5 is a different read-out. Early small-amplitude trajectories are
compared with all Fourier-mode matrix propagators, retaining every biochemical pole,
on the initial prefix before max|x/x*-1| exceeds 1e-3. At amplitude 1e-2 no linear
window is claimed. The largest valid-window relative L2 Delta error is 0.1801%.

Endpoint comparisons use every cell of the periodic lattice, independent of the hexagonal
plot window. Correlation is maximised over all periodic translations, with no rotations;
relative L2 distance uses that alignment. The exact K fraction sums both conjugate K
points and divides by all nonuniform Delta power. All 112 endpoints have max|dx/dt|/x*
below 1e-7, but this is a numerical stationary-candidate criterion, not a stability theorem.
The two operating points and local perturbations do not survey global basins.

The first seed at amplitude 1e-6 is repeated for both models and points using BDF
rtol=1e-11, atol=1e-14 and interval 0.25 through the same fixed horizons. Four independent
DOP853 runs use these tighter settings through five time units beyond the original
CV=0.5 event. Onset differences remain below 0.006 dimensionless time units; endpoint
K-power differences in the tighter BDF runs are below 4.1e-13. The DOP853 checks do not
test long-horizon endpoints. All initial/final species fields and compact metric trajectories
are retained in per-case NPZ files; solver-check summaries and all seed-level results accompany them.

`scripts/mechanistic_predictions.py` separately reproduces the Discussion's alpha=0.25,
1 and 4 compensated-turnover rates, invariant static gains/critical values, clamped-input
relaxation times and half-power frequencies. These are model-based quantities, without
physical time calibration or experimental error bars.

## Time discretisation, independent methods and output cadence

`python scripts/reproduce_all.py --dt-audit` runs the frozen design in
`figures/data/timestep_audit/protocol.json`, generates Supplementary Note 4 and Figure S6,
and records `timestep_manifest.json`. Optional per-case reuse is available through
`python scripts/timestep_audit.py --reuse`; it verifies numerical-kernel, baseline-input
and NPZ hashes. The historical reference and initial-condition results are retained.

The output recording interval in DOP853/BDF is **not** a fixed integration dt.
The implementation chooses internal steps adaptively from `rtol` and `atol`;
`max_step` caps those internal steps. The distinction is documented in
[SciPy 1.16.3's implementation and API documentation](https://github.com/scipy/scipy/blob/v1.16.3/scipy/integrate/_ivp/ivp.py).
Changing only `t_eval` or a plot's sampling interval is not a time-step-convergence test.

The audit therefore combines four separate comparisons:

1. All 54 distinct mode blocks at P0/P2/P3, full and DNA: the RK4 one-step matrix
   and its t=2 propagator are compared with continuous eigenvalues and the matrix exponential
   for h=0.08, 0.04, 0.02 and 0.01 (1296 mode/model/step comparisons).
2. Sixteen full-lattice early-growth trajectories use genuine fixed-step RK4 at those
   steps, compared on the common 0.08 output grid with high-accuracy Radau. The same
   time window is used for all growth fits, ending at the last common time no greater
   than log(1000)/the predicted growth rate. No concentration clipping is applied;
   intermediate RK4 stage concentrations are monitored. This test does not license RK4
   at those steps in the much stiffer late nonlinear regime.
3. Fourteen first-seed initial states cover every amplitude/model/point in the main
   random-start design, plus two independent-S controls. Each original BDF run is compared
   with BDF caps of 0.5, 0.25 and 0.125, tighter BDF with cap 0.125, and independent
   tighter Radau with cap 0.25. All 70 new trajectories reach the same fixed horizons,
   t=600/P2 and 300/P3. Actual accepted internal step meshes, initial/final fields,
   event times and endpoint drift residuals are retained.
4. The identical dense Radau growth trajectory is resampled at intervals 0.5, 0.25
   and 0.125, reapplying the reference manuscript's first-recorded-amplitude window rule.
   Thus the effect of output cadence is isolated from integration error.

Assessment targets are fixed in the protocol before production. Every target, including
any missed coarse-step accuracy target, appears in the generated supplement and summary.
The time-continuous equilibrium/Jacobian calculations underlying the stationary threshold
and sampled instability fractions do not use a temporal integration dt; their separate
numerical checks concern equilibrium residuals, eigenvalue margins and selected high-precision
recomputations. No finite audit proves the absence of all numerical error. The tested
discrete graph is a model of cell contacts, not a spatial approximation whose continuum
limit is established here.

The coarse explicit controls expose, rather than hide, a numerical stability failure:
at P0, RK4 h=0.04 and 0.08 give incorrect RHP counts in all 54 full and all 54 DNA
blocks at each step, despite continuous stability. At h=0.02 and 0.01 all finer-step
counts agree. These fixed RK4 steps were not used in the production adaptive runs.
An accurate output cadence cannot repair an unstable integration step. Conversely,
the fine fixed steps tested near the homogeneous state are not claimed adequate for
the much stiffer late nonlinear trajectory or for arbitrary sampled parameters.

All adaptive production-comparison targets pass: maximum relative event-time change
0.01631%, raw all-species endpoint relative L2 difference 1.11331e-10, and endpoint
K-fraction change 4.02234e-13. In all 30 reintegrated matched full/DNA pairs the full
model remains later at CV=0.1. Output resampling changes fitted growth slopes by at most
0.002105%. At fixed RK4 h=0.01, early Delta perturbation errors relative to the high-accuracy
Radau reference stay below 3.8e-8.

The original `all_assessments_pass` flag remains **false** because coarse controls fail:
the 216 P0 spectral mismatches and the P3 DNA h=0.08 field error (1.25303e-4, just above
the prespecified 1e-4 accuracy target). The latter passes at every smaller tested step.
`assessment_interpretation.json` separates these diagnostic failures from the passing
adaptive-production comparisons without changing targets or discarding results.

After the prespecified audit, `scripts/write_timestep_report.py` additionally reuses the
stored high-accuracy reference amplitudes to shorten the theoretical fitting window.
This post-audit diagnostic performs no new trajectories and changes no assessment.
Reducing the upper amplitude from 1e-3 to 1e-4 and 1e-5 reduces the maximum growth-slope
discrepancy from 1.00451e-4 to 1.84378e-5 and 4.04517e-6, respectively, consistent with
finite-amplitude bias after temporal errors are small. Its source data are in
`amplitude_window_diagnostic.csv`. Final report/figure generation is recorded separately
in `timestep_report_manifest.json`; the original numerical run remains in `timestep_manifest.json`.


## Seven additional reviewer controls

`python scripts/reproduce_all.py --reviewer-extensions` runs all seven controls and writes
`reviewer_manifest.json`. The frozen protocol is
`figures/data/reviewer_extensions/protocol.json`; [the supplement index](SUPPLEMENTARY_INDEX.md)
links Notes 5-11, Tables S7-S13, Figures S7-S13 and individual commands. Main Figure 6
and Section 3.7 contain the promoted results. Existing numerical kernels remain unchanged.

The linear stage selects 50 original unstable cases by sorted growth rank plus P2/P3 for
52 QSS families, uses prespecified training/held-out clock ranges, follows all twelve
reference DC misses, and reweights the checksum-verified 40000 primary inputs. Prior
weights add no samples. The alternative production lag has an exact linear similarity
proved in Note 11; numerical equivalence checks verify implementation.

Variable-size controls use seeds 2026102001-2026102003 and a separate PCG64(seed+1) stream
for quenched heterogeneity. The 52 trajectories retain actual initial/final fields and
production multipliers. Open boundaries retain raw sums and are forced relative to the
periodic reference. Analytic sparse Jacobians are checked against patterned-state finite
differences and original graph-mode blocks.

Twelve endpoints are selected deterministically, refined by positive Newton iteration,
and analysed with complete dense spectra plus independent shifted sparse eigenpairs.
Recovery seeds are 2026102201-2026102212. Both perturbed and unperturbed states are
continued for 200 time units to distinguish recovery from baseline drift. Per-case reuse
verifies the kernel and NPZ checksums. Kernel digests canonicalise text line endings;
provenance additionally records actual working-file byte hashes.

Historical manifests describe their individual runs. Later report/plot changes do not
retroactively alter an earlier numerical run. The final extension manifest records the
final pipeline sources, actual packages and artifact hashes. Delivery bundles include
SHA256 checksums of their actual contents. Fresh runs remain available without --reuse.


The first combined extension run completed every scientific calculation but stopped at
an en-dash in the report script's console message under Windows cp949 output encoding.
`reviewer_initial_attempt_manifest.json` retains that execution record. The console
message was made ASCII; the final combined run with verified case reuse completed with
unchanged sources. A subsequent report/figure-only pass, after visual layout improvements,
is recorded in `reviewer_report_manifest.json`; no numerical kernel or trajectory changed.
The final test suite passes all 373 tests. Existing manuscript DOCX/PDF binaries are
historical drafts and are excluded from the delivery package; current manuscript and
supplement sources are Markdown and require typesetting and visual review for submission.
