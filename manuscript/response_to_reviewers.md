---
title: "Response to the referees"
subtitle: "Revised manuscript: Explicit downstream dynamics separate stationary onset from growth in Notch–Delta patterning"
author: Hyeonje Yang
date: October 2026
geometry: margin=2.5cm
fontsize: 11pt
---

I thank both referees for careful and specific reports. Each comment below is summarised in
italics and followed by my response and the corresponding changes. Section, figure and
equation numbers refer to the revised manuscript unless marked "v1".

Three things changed throughout:

* **Scope and claims.** The paper is reorganised around one question: which predictions require
  the explicit activated-Notch and proneural dynamics, and which follow from static responses
  alone? The Schur/static-gain identity is now presented as a standard tool, not a contribution.
  Material both referees regarded as known or unnecessary has been removed: the nonlinear
  saturation statistics, the "global rather than centred" discussion (v1 §3.3), and the
  stochastic check (v1 §3.4).
* **Numerics.** Every equilibrium-based statistic was recomputed with a corrected scalar
  solver. Growth-rate validation uses exact eigenvectors; a separate paired random-start
  study tests the scope of the resulting timing and endpoint comparisons.
* **Reproducibility.** Every number has a generating script in the public repository. A test
  suite numerically checks the analytic identities and limiting cases alongside
  the proofs. The v1 code and
  data are preserved unchanged in an archive folder of the repository.

# Referee 1

**1.1** *Steady-state failures distort the global statistics: 733 of 4000 samples were labelled
as having no positive fixed point, although the homogeneous equilibrium exists uniquely and can
be found by scalar bracketing. Recalculation gives 36.325% instead of 30.3%. Solver failure must
be distinguished from nonexistence.*

I agree on all points. The homogeneous steady state reduces to a scalar equation $F(D)=0$
(Eq. 3). $F$ is strictly decreasing, with $F(0)>0$ and $F(\lambda_D/d_D)<0$, so for every
positive parameter set the steady state exists, is unique, and lies in $(0,\lambda_D/d_D)$
(§2.3, Supplementary Note 1). It is now computed by Brent bracketing.

* All 733 v1 "no fixed point" samples have a steady state. Residuals are reported relative to
  $\max(1,\|x^*\|_\infty)$, rather than as unnormalised errors.
* The nonuniform-instability fraction is **1453/4000 = 36.325%**, which reproduces the
  referee's value (§3.2).
* The two-parameter maps (v1 Fig. 2C) had 115 missing grid points from the same cause. They
  have been recomputed (Fig. 3C,D).
* Every derived statistic (classification counts, gain-criterion accuracy, dominant-mode
  counts) has been recomputed.
* A test now checks all 733 previously failed samples.
* The v1 statement "multistable branches found: 0" is replaced by the uniqueness proof.

**1.2** *The explanation of the DC-gain misclassifications is incorrect. All seven originally
reported mismatches have dominant mode $\mu=-3$; e.g. sample 1730 has gain ≈0.366 but unstable
eigenvalues ≈0.159±0.262i. The failures should be explained within the tested mode, without
conflating the stationary zero-eigenvalue condition with general stability.*

The referee is right, and the v1 explanation ("another mode can dominate") was wrong. The
revised §3.1 and §3.3 explain the mismatches within the tested mode:

* **The plant is Hurwitz.** For every positive parameter set and every graph eigenvalue, the
  plant matrix $A_\mu$ is Hurwitz (§3.1(i), proved). Hence $\det(A_\mu)(-d_S)>0$, and by
  Eq. 5 $\operatorname{sign}(1-G)=\operatorname{sign}\det J_\mu$.
* **What $G$ counts.** $G(\mu)=1$ exactly at a zero eigenvalue, and $G(\mu)>1$ exactly when
  there is an odd number of positive real eigenvalues. $G<1$ is explicitly stated not to be a
  stability certificate.
* **The mismatches.** With the corrected sample there are 12, including the 7 of v1. All have
  dominant mode $\mu=-3$. Each has **exactly two** right-half-plane eigenvalues in the
  $\mu=-3$ block: a complex pair at 6 points (including 1730: $G=0.366$,
  $0.159\pm0.262i$) and two real eigenvalues at the other 6. An even count leaves
  $\det J_{-3}>0$, so these instabilities are invisible at zero frequency.
* **Complete test.** The complete per-mode test is the Nyquist criterion, which is valid here
  because the open loop has no unstable poles. For sample 1730, $L_{-3}(i\omega)$ crosses the
  real axis at 1.56 when $\omega=0.36$, and the curve encircles $+1$ twice (Fig. 4B).
  Encirclement counts equal the eigenvalue counts at all 4000 samples.
* **Text changes.** The v1 limitation sentence attributing the mismatches to another mode has
  been deleted. The stationary condition and general stability are now kept separate.

**1.3** *Componentwise scaling by equilibrium concentrations changes the eigenvector and
introduces mixed-mode transients. The late-window slope does match, but the initialization
should be corrected or justified and the actual validation window displayed.*

Corrected. Each simulation now starts from $x^*+\varepsilon v\cos\phi$, where the whole
eigenvector $v$ is multiplied by a single scalar (Methods). The fitting window was fixed in
advance: from $t=0$ until the three-sublattice amplitude of $D$ reaches $10^{-3}D^*$. It is
drawn in Fig. 5C. With the corrected start:

| point | predicted | fitted | window |
|---|---|---|---|
| P2 | 0.103067 | 0.103077 | $[0, 67.25]$ |
| P3 | 0.216300 | 0.216305 | $[0, 32]$ |

The corrected initial direction removes the linear mixing of biochemical eigenvectors;
finite-amplitude nonlinear terms and integration errors remain. Fig. S1 shows the v1 start for
comparison. It gives early slopes of 0.443 (P2) and 0.562 (P3), the mixed-mode transient the
referee identified. At P2 it also produces a long intermediate plateau before the final pattern
is reached. The v1 phrases "frozen validation window" and "the earliest fit contains initial
transient contamination" are removed.

**1.4** *The unit-DC-gain/zero-eigenvalue equivalence follows from the Schur-complement identity
under the relevant nonsingularity conditions; numerical agreement verifies consistency but does
not establish novelty. A substantive comparison with graph-based approaches and a matched
reduced model is needed to identify predictions that specifically require the explicit
downstream dynamics.*

I agree. The identity is now presented as standard (§1, §3.1, §4). It is related explicitly to
graph-based static-gain criteria [11, 12] through
$\det J_\mu=\det J_\mathrm{loc}\det(I-\mu K)$, where $K$ is the static input–output gain of an
isolated cell, when $J_\mathrm{loc}$ is nonsingular (Methods). This determinant identity does
not establish the additional monotonicity and cell-stability hypotheses of the graph-based
theorems, or subsume their results on patterned equilibria. The normalised graph convention
in those studies is also distinguished from our raw neighbour sums. The revision adds the two comparisons
requested (§3.4, Fig. 5, Tables 2–3). The Introduction and Discussion now distinguish
previous biological models, graph theorems and our matched-control question. We do not
claim that downstream timing has not been studied before, or that all reduced models fail
to fit a dynamic observable. The contribution is to isolate the error introduced by QSS
elimination when static responses and remaining kinetic parameters are held fixed.

*Matched reduced models.* Proneural activity, or both activated Notch and proneural activity,
are slaved by quasi-steady-state elimination (DNA and DN reductions). The DNA model is a
three-component control derived from our equations, with a composite repression law; it
is not identified with the exact equations of a published receptor–ligand–reporter model.
Both reductions have exactly the
same steady states, homogeneous and patterned, and their Jacobians are Schur complements of
$J_\mu$. They therefore share the stationary boundary, the critical $\lambda_N$, and the
selected mode exactly. These predictions do **not** require the explicit dynamics, and the
manuscript now says so.

*Predictions that do require the explicit dynamics:*

* **Growth rate.** Along a proneural-timescale family that leaves every steady state and every
  static gain fixed, $\lambda_N^c=1.304815343$ is unchanged for $\alpha\in[10^{-2},10^4]$. The
  P2 growth rate varies from $8\times10^{-4}$ to 1.59. When proneural turnover is slow, the rate
  is $\approx d_S[G-1]$. Over the global sample, the matched DNA reduction overestimates the
  growth rate at every unstable point (median 4.9-fold; 15.5-fold at P2). Nonlinear
  simulations confirm this (pattern at $t\approx111$ versus $\approx7$), and both models reach
  the same final pattern for these P2 initial states, to within $10^{-9}$. We do not infer
  general equivalence of attractor stability or basins from identical equilibrium sets.
* **Stability differences beyond the static threshold.** The DNA reduction predicts 26 instabilities where the
  explicit model is stable (DN reduction: 47), all with $G<1$. In the other direction, the
  static criterion misses the 12 two-eigenvalue instabilities of point 1.2.

A further analytic result has been added (§3.1(iii)): the static gain is maximised at the most
negative graph eigenvalue for **all** positive parameters, despite the two-channel trans
coupling and cis inactivation. The stationarily selected mode is therefore the three-sublattice
mode everywhere. This also answers Referee 2's point 2.3.

The Discussion adds two concrete prospective tests. Compensated changes in production and
degradation, $(\lambda_S,d_S)\mapsto\alpha(\lambda_S,d_S)$, preserve $S_q(A)$ and the
stationary onset while altering growth. A change in degradation alone does not implement
this control. A clamped activated-Notch step or oscillatory input tests the finite proneural
relaxation removed by DNA. Required static-response controls, reporter delays, falsifying
observations and the absence of tissue calibration are explicit. These are model-based
predictions, not experimental evidence or a claim of a uniquely identified mechanism.

# Referee 2

**2.1** *Variables are used before being defined (e.g. in the abstract and on p. 1, l. 57, not
defined until p. 5); $(m,n)$ is undefined on p. 5, l. 45. The paper needs serious editing.*

The manuscript has been rewritten with a definition-before-use rule:

* The abstract contains no symbols.
* The model, every variable and parameter (also tabulated in Table 1), the graph eigenvalue
  $\mu$, and the plane-wave indices $(m,n)$ with their wavevectors are all defined in §2,
  before any result.
* Fig. 1D now shows $\mu(m,n)$ over all indices.
* Symbols in captions are defined within each caption.

**2.2** *The novelty relative to previous Notch–Delta work is unclear; checking full-model
stability from the exact Jacobian under a pattern assumption is not a contribution.*

I agree that the per-mode Jacobian check is not a contribution. It has been reduced to one
sentence in Methods. The revised paper states its contributions in §4 ("What is and is not new"):

* uniqueness of the homogeneous steady state for all positive parameters;
* a proof that cis inactivation and two-channel trans coupling never change the stationarily
  selected mode;
* an exact separation of predictions that depend only on static responses, which every matched
  reduced model shares, from dynamic predictions that can differ under reduction;
* a quantitative account of where a static stability test fails.

**2.3** *There could be cases in which a pattern other than the assumed one is unstable, but no
such case is presented.*

This is now addressed directly (§3.1(iii) and §3.3).

* **Proof for stationary instabilities.** For stationary instabilities no other pattern can
  come first. $G(\mu)<G(\mu_{\min})$ for every other graph eigenvalue, for every positive
  parameter set and every regular contact graph.
* **Search for non-stationary cases.** The proof does not cover non-stationary instabilities,
  so 40,000 additional random parameter sets were searched, half in a box ten times wider than
  the original. Of these, 13,048 were unstable. None had a nonuniform mode other than
  $\mu=-3$ dominant, or unstable while $\mu=-3$ was stable.
* **Other unstable modes.** In 1451 of the 4000 points, subdominant nonuniform modes are also
  unstable. P2 has 21 distinct unstable nonuniform graph eigenvalues, illustrated in the added
  Fig. S3. Thus the K mode is dominant, not the only unstable spatial mode. A uniform complex
  pair is unstable at 45 points and always coexists with K instability.
* **Mismatches within $\mu=-3$.** The twelve static-gain mismatches (point 1.2) are a different
  kind of case. The tested mode itself is unstable through two eigenvalues.

**2.4** *That the system can be stable to uniform perturbations while unstable to patterned ones
is the same as in all other Notch–Delta models.*

Agreed. This is no longer presented as a finding. It appears once, as a one-line consequence of
§3.1(iii) ($G(6)<0$, so the uniform mode cannot lose stability through a zero eigenvalue). In
§4 it is explicitly described as a generic feature of lateral-inhibition models.

**2.5** *Putting equations as subfigures is not good practice.*

Removed. The v1 Fig. 1D (equation box) and Fig. 3A–B (block diagram with equations and the
Schur identity) no longer exist. All equations are in the text and are numbered.

**2.6** *Belabouring the growth of the linear instability into a stable nonlinear pattern is
unnecessary given earlier work.*

Agreed. The v1 Fig. 5 is removed: amplitude order parameter, snapshots, Delta histogram,
solver comparison and Moran's I. Nonlinear simulations now appear only where they test a new
claim (Fig. 5C,D): that the explicit model and its matched reduction grow at different, predicted
rates. The same final field is reached only for those symmetry-restricted starts. The new
§3.6 and Fig. S5 show distinct endpoint fields after paired random starts, despite the shared
equilibrium set. This limits the earlier example rather than promoting nonlinear pattern
formation itself as a contribution. The paragraph on fate commitment is cut to one sentence
in the Limitations.

**2.7** *Sections 3.1–3.2 (distinction from simpler three-component models) are the original
part; §3.3 seems trivial; the results of §3.4 are not shown in figures and their relevance to
stochasticity is unclear.*

Thank you. The revision is built around the part the referee found original. The comparison
with three-component models is now quantitative: an explicitly constructed matched
three-component (DNA) reduction, analysed in §3.4 and Fig. 5. v1 §3.3 has been deleted. v1 §3.4
(the stochastic check) has been removed from the manuscript, since it was not needed for any
conclusion. Its code remains in the repository archive.

# Additional robustness analysis

The revision now includes 115000 new parameter evaluations under a protocol fixed before
computation (§3.5, Supplementary Note 2, Table S1 and Fig. S4). Ten independent scrambled
4000-point Latin-hypercube designs in the original box give a mean nonuniform-instability
fraction of 37.37% (approximate 95% Monte Carlo interval 37.07–37.67%). The interval is based
on independent replicate fractions, rather than treating dependent LHS points as iid trials.
The referee's 36.325% is preserved exactly as the corrected result for the original reference
sample. It is not presented as a universal fraction.

The primary study has a pooled DNA/full growth-ratio median of 5.41, with replicate medians
5.25–5.61. DNA and DN miss no full-model instabilities and grow faster at every unstable
point. They have 216 and 445 reduced-only instabilities, respectively. The static threshold
misses 150 points; every miss has two RHP eigenvalues within the tested K-mode block.

Separate size checks and five wider-box designs retain these qualitative comparisons while
changing the quantitative magnitudes. In the wider box the mean unstable fraction is 28.89%
and DNA replicate medians are 2.99–3.64. Thus neither fraction nor median ratio is claimed
to be independent of the sampling measure. One near-marginal point is flagged explicitly.
Targeted Nyquist checks at 842 points and independent 80-digit recomputation of eleven cases
agree with the eigenvalue analysis. No sampled absence of an exception is treated as a proof.
All input vectors and per-point/per-replicate results are deposited with code and checksums.

# Initial-condition robustness and scientific scope

Section 3.6, Supplementary Note 3, Tables S2–S3 and Figure S5 report 112 production
trajectories: 48 full/DNA pairs with identical initial D, N, A and reconstructed S fields,
plus 16 independent-S controls, at two operating points and three amplitudes. In all
48 pairs the full model crosses CV(D)=0.1 later: onset ratios span 10.96–15.10 at P2
and 12.34–12.90 at P3. All 112 trajectories cross CV=0.5 as well. Mixed-mode early
transients are compared with the all-mode linear propagator, rather than fitted as a
single exponential. Four tighter BDF and four independent DOP853 checks support the
onset comparison; the documented solver change after a pilot is disclosed.

Endpoint fields differ in every main pair even after periodic translation. All fixed-horizon
endpoints meet the numerical stationary-candidate residual criterion, which is explicitly
distinguished from a proof of stability or a global basin survey. The identical fields
in Figure 5 therefore remain an initial-state-specific example. Timing robustness and
endpoint equality are no longer conflated. Prior studies and their assumptions are compared
in Table 3; biological predictions now specify the compensated intervention, transient
read-out, static controls and possible falsification. No claim of novelty is based solely
on retaining a fourth variable, a standard determinant identity or nonlinear saturation.

# Additional time-discretisation audit

Supplementary Note 4, Tables S4–S6 and Figure S6 distinguish output cadence from
the actual internal integration step. Fourteen representative initial states span both
operating points, all three amplitudes, both models and independent-S controls. Each
is reintegrated through the original horizon under five step-cap/tolerance/solver
settings, including independent Radau: 70 additional trajectories. Emergence times
change by at most 0.01631% and raw all-species endpoint relative L2 differences remain
below $1.12\times10^{-10}$. The checks are representative, not an exhaustive basin survey.

The audit retains a genuine coarse-step failure: fixed RK4 at h=0.04 or 0.08 falsely
destabilises the stable P0 mode blocks. These steps were not used for the production
adaptive runs. All finer-step spectral counts at h=0.02 and 0.01 agree with the
continuous-time Jacobian. Separate fixed-step nonlinear growth checks and output-cadence
resampling are reported with every assessment target, including any missed coarse-step
target. The stationary threshold and sampled instability fractions are obtained directly
from equilibrium/Jacobian calculations and have no temporal integration dt.

# Additional matched controls and scientific interpretation

We added seven individually reproducible supplementary notes (5–11), Tables S7–S13 and
Figures S7–S13, indexed in `docs/SUPPLEMENTARY_INDEX.md`. Main Figure 6 and Section 3.7
promote the results that affect the scientific interpretation:

* **QSS validity:** 52 compensated-turnover families show family-specific accuracy
  requirements. A small leading-pole/downstream-rate ratio alone does not bound relative
  growth error near onset. The complete pole diagnostics and every family are retained.
* **Reduced-model recalibration:** a single clock fitted at P2 underpredicts P3 growth by
  16.5%. A separate training-range fit retains held-out errors. This does not exclude
  fitting a more flexible reduced model to selected observations.
* **DC misses:** all twelve reference misses enter instability through nonzero-frequency
  K-mode crossings on constant-gain turnover paths. Six later become real unstable pairs.
  We report transversality but do not claim nonlinear Hopf nondegeneracy or stable cycles.
* **Mechanistic inference:** an explicit alternative with a finite Delta-production lag
  is exactly similar to the full homogeneous Jacobian and shares linear D,N,A responses.
  Thus these observables cannot uniquely identify proneural memory. A specified nonlinear
  clamped-A production response distinguishes the alternatives; independent proneural
  measurements supply another test. This standard similarity is used as a model-specific
  control, not presented as a new general mathematical identity.

Further controls test periodic sizes 23/24/25, raw open edges and quenched production
heterogeneity (52 trajectories). Periodic onset shifts by at most 0.253%, with slower full
emergence in every one of 18 matched pairs. Open edges create deterministic forcing; their
contrast is explicitly not called a homogeneous instability. Twelve selected patterned
states have strictly stable complete spectra and recover from relative 1e-4 perturbations,
establishing local stability only for these states. Growth differences persist away from
small denominators, while alternative prior weights alter fractions and effect sizes.
Full protocols, case selection, kernels, initial/final fields, complete spectra and
checksums are deposited. None of these finite controls is a global basin theorem or tissue
validation. The manuscript, Discussion and abstract now distinguish finite memory from
the unique biochemical localisation of that memory.

# Summary of changes

| v1 | revised |
|---|---|
| 30.3% unstable (733 solver failures excluded) | 36.3% unstable; steady state proved unique |
| mismatches attributed to "another mode" | two RHP eigenvalues in the $\mu=-3$ block; Nyquist count |
| componentwise-scaled initial state, "frozen" window | exact eigenvector; prespecified amplitude rule and displayed windows |
| unity-gain identity presented as a key result | presented as standard; equivalence to graph static-gain criteria |
| no reduced-model comparison | matched DNA and DN reductions; timescale family; Table 2 |
| equations inside Figs. 1D, 3A–B | equations only in text |
| nonlinear saturation figure (v1 Fig. 5) | removed; one targeted comparison kept |
| v1 §3.3 and §3.4 (stochastic) | removed |
| — | proof that stationary onset always selects the three-sublattice mode |
