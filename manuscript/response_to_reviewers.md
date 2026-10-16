---
title: "Response to the referees"
subtitle: "Mode-resolved feedback gain controls spatial symmetry breaking in an explicit Notch–Delta–proneural circuit (revised as: Explicit Notch–proneural dynamics set the tempo, not the onset, of lateral-inhibition patterning)"
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
* **Numerics.** Every statistic was recomputed with a corrected steady-state solver and an
  exact-eigenvector initial state.
* **Reproducibility.** Every number has a generating script in the public repository. A test
  suite (341 checks) verifies each analytic statement on random parameter sets. The v1 code and
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

* All 733 v1 "no fixed point" samples have a steady state, with residual below $3\times10^{-15}$.
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

The growth is exponential from $t=0$, with no transient. Fig. S1 shows the v1 start for
comparison. It gives early slopes of 0.443 (P2) and 0.568 (P3), the mixed-mode transient the
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
isolated cell (Methods; verified to $10^{-15}$). The revision adds the two comparisons
requested (§3.4, Fig. 5, Table 2).

*Matched reduced models.* Proneural activity, or both activated Notch and proneural activity,
are slaved by quasi-steady-state elimination (DNA and DN reductions). The DNA reduction is a
Sprinzak-type three-component Notch–Delta–reporter model. Both reductions have exactly the
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
  the same final pattern, to within $7\times10^{-11}$.
* **Non-stationary instabilities.** The DNA reduction predicts 26 instabilities where the
  explicit model is stable (DN reduction: 47), all with $G<1$. In the other direction, the
  static criterion misses the 12 two-eigenvalue instabilities of point 1.2.

A further analytic result has been added (§3.1(iii)): the static gain is maximised at the most
negative graph eigenvalue for **all** positive parameters, despite the two-channel trans
coupling and cis inactivation. The stationarily selected mode is therefore the three-sublattice
mode everywhere. This also answers Referee 2's point 2.3.

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

* uniqueness of the steady state for all positive parameters;
* a proof that cis inactivation and two-channel trans coupling never change the stationarily
  selected mode;
* an exact separation of predictions that depend only on static responses, which every matched
  reduced model shares, from those that depend on downstream timescales, which reduced models
  get wrong;
* a quantitative account of where a static stability test fails.

**2.3** *There could be cases in which a pattern other than the assumed one is unstable, but no
such case is presented.*

This is now addressed directly (§3.1(iii) and §3.3).

* **Proof for stationary instabilities.** For stationary instabilities no other pattern can
  come first. $G(\mu)<G(\mu_{\min})$ for every other graph eigenvalue, for every positive
  parameter set and every regular contact graph.
* **Search for non-stationary cases.** The proof does not cover non-stationary instabilities,
  so 40,000 additional random parameter sets were searched, half in a box ten times wider than
  the original. Of these, 13,047 were unstable. None had a nonuniform mode other than
  $\mu=-3$ dominant, or unstable while $\mu=-3$ was stable.
* **The one other unstable mode.** The only other unstable mode found is the uniform mode,
  through a complex pair (synchronous oscillation), at 45 of the 4000 sampled points. It always
  coexists with an unstable three-sublattice mode, and it is reported as such.
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
rates yet reach the same final pattern. The paragraph on fate commitment is cut to one sentence
in the Limitations.

**2.7** *Sections 3.1–3.2 (distinction from simpler three-component models) are the original
part; §3.3 seems trivial; the results of §3.4 are not shown in figures and their relevance to
stochasticity is unclear.*

Thank you. The revision is built around the part the referee found original. The comparison
with three-component models is now quantitative: an explicitly constructed matched
three-component (DNA) reduction, analysed in §3.4 and Fig. 5. v1 §3.3 has been deleted. v1 §3.4
(the stochastic check) has been removed from the manuscript, since it was not needed for any
conclusion. Its code remains in the repository archive.

# Summary of changes

| v1 | revised |
|---|---|
| 30.3% unstable (733 solver failures excluded) | 36.3% unstable; steady state proved unique |
| mismatches attributed to "another mode" | two RHP eigenvalues in the $\mu=-3$ block; Nyquist count |
| componentwise-scaled initial state, "frozen" window | exact eigenvector; window fixed in advance and drawn |
| unity-gain identity presented as a key result | presented as standard; equivalence to graph static-gain criteria |
| no reduced-model comparison | matched DNA and DN reductions; timescale family; Table 2 |
| equations inside Figs. 1D, 3A–B | equations only in text |
| nonlinear saturation figure (v1 Fig. 5) | removed; one targeted comparison kept |
| v1 §3.3 and §3.4 (stochastic) | removed |
| — | proof that stationary onset always selects the three-sublattice mode |
