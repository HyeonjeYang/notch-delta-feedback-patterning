# Supplementary Note 11. What linear data identify about downstream memory

## A finite-memory control at another location in the circuit

Does a growth discrepancy establish that memory must reside in the proneural variable?
To distinguish finite feedback memory from its biochemical location, retain $D,N,A$ and
replace proneural $S$ by a Delta-production variable $U$:

$$
\dot D_i=U_i-D_i[d_D+f_Ds_{N,i}+k_cN_i],\qquad
\dot U_i=d_S\{g(A_i)-U_i\},\quad
g(A)=\lambda_Dh_D[S_q(A)],\quad S_q(A)=\lambda_Sh_S(A)/d_S.
$$

The $N,A$ equations are unchanged. Every equilibrium maps exactly via
$U_i^*=g(A_i^*)=\lambda_Dh_D(S_i^*)$. Around a homogeneous equilibrium define
$\kappa_D=\lambda_Dh_D'(S^*)>0$ and
$T=\operatorname{diag}(1,1,1,\kappa_D)$. Direct differentiation gives

$$J^{\mathrm{prod\ lag}}_\mu=T J^{\mathrm{full}}_\mu T^{-1}.$$

Indeed, the first equation receives $\delta U$ with unit coefficient, while
$\partial\dot U/\partial A=d_Sg'(A^*)=
\kappa_D\lambda_Sh_S'(A^*)$. The fourth diagonal element remains $-d_S$.
This is a standard state-scaling similarity applied to this specified alternative, not a
new general control-theory theorem. The hidden perturbation coordinates map as
$\delta U=\kappa_D\delta S$.

Because the first three coordinates are unchanged, all homogeneous-mode linear $D,N,A$
responses to inputs in those coordinates are identical, including the growth spectrum.
Verification on 4000 reference inputs and 54 modes (216000 blocks) gives maximum normalised
similarity residual $1.12\times10^{-16}$. Stable-P0 forced $D,N,A$ responses at
$\mu=-3,0,6$ and 201 frequencies agree within $4.51\times10^{-16}$ in relative norm.
The identity, rather than these numerical agreements alone, establishes equivalence.

**Inference limit:** these linear observables cannot uniquely localise the memory to
proneural turnover. The full-versus-QSS comparison establishes the relevance of finite
feedback relaxation under fixed static responses. It does not uniquely identify the
biochemical implementation of that memory.

## A nonlinear response that separates the specified alternatives

At P2, clamp $A$ to 1.001, 1.2 or 2 times its equilibrium value and start the hidden states
at their old equilibria. In the explicit model,
$S(t)=S_1+(S_0-S_1)e^{-d_St}$ and Delta production is $\lambda_Dh_D[S(t)]$.
For the alternative, $U(t)=U_1+(U_0-U_1)e^{-d_St}$.
Both share the same initial/final production and linearised response, but the Hill
nonlinearity acts after relaxation in the first and before relaxation in the second.
We normalise each response by its own total signed production change, so it runs from zero
to one even when production decreases.

Maximum differences as a fraction of that total response are 0.0631%, 11.9854% and 31.9231%
for the three step sizes. The small-step agreement and large-step separation provide a
specified discrimination experiment. We predict a **production** response to externally
clamped A, not a whole-tissue trajectory or validated experimental measurement.
Independent time-resolved proneural measurements or a controlled nonlinear response can
test the implementation; unrestricted refitting of either alternative would require a new
identifiability analysis. No tissue data or assay feasibility is claimed here.

## Separate Hill thresholds and concentration units

Using the same threshold in both Hill functions can look like an extra constraint. Suppose
the activation threshold of $h_D(S)$ is $K_S$ and repression threshold of $h_S(A)$ is $K_A$.
For unrestricted positive proneural production and an uncalibrated S concentration scale,
define $S'=(K_A/K_S)S$, $\lambda'_S=(K_A/K_S)\lambda_S$ and $K_d=K_A$.
Then $h_D(S;K_S)=h_D(S';K_d)$, while $h_S(A;K_A)=h_S(A;K_d)$, and the S equation retains
the stated form. Thus the two-threshold deterministic model is represented by a rescaling
of S and its production rate. This does **not** preserve a prior on original concentrations
or production, and independently calibrated S units would constrain the rescaling.
One hundred positive random-state checks give maximum relative vector-field discrepancy
$7.63\times10^{-16}$. The algebra supplies the result; the checks confirm implementation.

## Reproduction

`scripts/reviewer_linear.py` and `scripts/reviewer_models.py` implement the alternative.
Data are `lag_similarity_checks.csv`, `mechanism_frequency_checks.csv`,
`mechanism_step_responses.csv` and `separate_Hill_threshold_checks.csv`. All original
parameters, step amplitudes and frequencies are retained. Use the report and figure commands
in Note 5.

The design was recorded in `figures/data/reviewer_extensions/protocol.json` before the
production calculations. All selected cases and counterexamples are retained. Source hashes
and case checksums are verified by `scripts/write_reviewer_reports.py`. The primary package
versions and numerical reproduction contract are given in `docs/REPRODUCIBILITY.md`.

![**Figure S13. Mechanism discrimination. (A–C) Normalised production after clamped-A steps of 1.001, 1.2 and 2 times baseline. Full and alternative models share equilibrium and linear D,N,A response but differ at finite amplitude. (D) Stable-P0 frequency-response differences are at floating-point precision; the exact similarity identity explains their agreement.**](../figures/supplementary/FigS13_mechanism_identifiability.png){width=100%}
