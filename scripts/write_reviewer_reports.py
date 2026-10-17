"""Write seven separate supplementary notes from verified extension outputs."""
from __future__ import annotations
import json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
import numpy as np
import pandas as pd
from scripts.reviewer_protocol import OUT,digest,sha,freeze


def table(frame,columns=None,formats=None):
    frame=frame if columns is None else frame[columns]
    formats=formats or {}
    labels={"alpha_required":r"$\alpha$ required on grid","target":"growth-error tolerance",
        "size":"periodic side length","min_mu":r"$\mu_{\min}$","min_mu_multiplicity":"multiplicity",
        "stationary_lam_N":r"$\lambda_N^c$","sample_id":"reference ID","alpha":r"$\alpha$",
        "mu":r"$\mu$","frequency":r"$|\mathrm{Im}\,s|$","G_K":r"$G(-3)$",
        "transverse_logalpha_slope":r"$d\mathrm{Re}\,s/d\log\alpha$","pair_real_mean":"collision real part",
        "pair_split":"numerical pair splitting","index":"state index","dimension":"equations",
        "leading_real":"leading real part","n_RHP":"RHP poles","recovery_final_over_initial":"final / initial distance",
        "cut_kind":"quantity cut","cut":"minimum value","n":"points","ratio_median":"growth ratio median",
        "gap_median":"absolute gap median","symmetric_error_median":"symmetric error median",
        "ratio_q25":"ratio Q25","ratio_q75":"ratio Q75","fraction":"unstable fraction",
        "mc95_low":"MC 95% lower","mc95_high":"MC 95% upper","weight_mean":"mean weight",
        "unstable_growth_ratio_median":"weighted growth-ratio median"}
    rows=["| "+" | ".join(labels.get(c,c) for c in frame.columns)+" |","|"+"|".join(["---"]*len(frame.columns))+"|"]
    for _,row in frame.iterrows():
        cells=[]
        for key,val in row.items():
            if pd.isna(val):value="not crossed / not applicable"
            elif key in formats:value=formats[key](val)
            elif isinstance(val,(float,np.floating)):value=f"{val:.6g}"
            else:value=str(val)
            cells.append(value)
        rows.append("| "+" | ".join(cells)+" |")
    return "\n".join(rows)


def read(name):return pd.read_csv(OUT/(name+".csv"))
def write(name,n,title,text,stem,caption):
    body=f"# Supplementary Note {n}. {title}\n\n"+text.strip()+f"\n\n![**Figure S{n+2}. {caption}**](../figures/supplementary/{stem}.png){{width=100%}}\n"
    (ROOT/"manuscript"/name).write_text(body.rstrip()+"\n",encoding="utf-8")


def main():
    freeze()
    summaries={}
    for stage,filename in (("linear","linear_summary.json"),("geometry","geometry_summary.json"),("pattern_stability","pattern_summary.json")):
        s=json.loads((OUT/filename).read_text())
        if s["kernel_sha256"]!=digest(ROOT/f"scripts/reviewer_{stage}.py"):raise RuntimeError("stale "+stage+" results")
        summaries[stage]=s
    for folder in ("geometry","pattern_stability"):
        for meta_path in (OUT/folder).glob("*.json"):
            meta=json.loads(meta_path.read_text())
            if sha(meta_path.with_suffix(".npz"))!=meta["npz_sha256"]:raise RuntimeError("modified "+str(meta_path))
    common="""
The design was recorded in `figures/data/reviewer_extensions/protocol.json` before the
production calculations. All selected cases and counterexamples are retained. Source hashes
and case checksums are verified by `scripts/write_reviewer_reports.py`. The primary package
versions and numerical reproduction contract are given in `docs/REPRODUCIBILITY.md`.
"""
    q=read("qss_validity");t=read("qss_accuracy_thresholds");sample=q[(q.family=="sample_3471")&np.isclose(q.alpha,.01)].iloc[0]
    qtext=r"""
## Question and matched comparison

When does eliminating proneural relaxation preserve growth accurately? We use 50 distinct
full-model-unstable cases selected at evenly spaced ranks of the original sample's sorted
growth rates, plus P2 and P3. This is a coverage design, not a new random sample. In each
family, $(\lambda_S,d_S)\mapsto\alpha(\lambda_S,d_S)$ preserves every equilibrium and every
static loop gain. We evaluate 61 logarithmically spaced values $10^{-2}\leq\alpha\leq10^4$
and all 54 distinct graph modes, including the uniform mode: 3172 family evaluations.
The matched DNA model is unchanged by this compensation.

The relative growth error is $e=|\lambda_{DNA}/\lambda_F-1|$ where $\lambda_F>10^{-8}$;
we retain a scaled absolute error separately. The diagnostic
$q=|s_{\rm lead}|/(\alpha d_S)$ measures the leading pole relative to downstream relaxation.
At that pole, replacing the exact downstream factor
$\alpha d_S/(s+\alpha d_S)$ by one gives the DC-normalised difference
$|1-\alpha d_S/(s+\alpha d_S)|=|s/(s+\alpha d_S)|$. These are mode-specific diagnostics,
not error bounds for an eigenvalue.
The sensitivity of a root of the complete characteristic equation matters as well.

## Results and practical use

For each family and tolerance, we report the first sampled $\alpha$ after which **every later
sampled value** meets the tolerance. These are grid-based empirical thresholds, not certified
continuous-interval bounds. All 52 families reach both tolerances by $10^4$.
The 10% thresholds range from 1 to 3162.28 (median 125.89); the 5% thresholds range from
1.995 to 6309.57 (median 251.19). At P2, the respective thresholds are 199.53 and 398.11;
at P3 they are 398.11 and 794.33. There are no full-faster-than-DNA cases in this sweep.

A small $q$ alone is insufficient for relative growth accuracy: 27 of 677 evaluated cases
with $q\leq0.01$ still exceed 10% growth error. These are dependent parameter-family
evaluations, not independent Bernoulli trials. At original sample 3471 and $\alpha=0.01$,
"""+f"$q={sample.pole_over_S_rate:.6g}$ and feedback substitution error is {sample.feedback_substitution_error:.6g},\nwhile the relative growth error is {sample.relative_growth_error:.4g}.\n"+r"""
Its gain is close to the stationary boundary. A small change in feedback can produce a large
*relative* change in a near-zero growth rate. The eigenvalue condition recorded in the data
is a descriptive conditioning measure, not a rigorous nonlinear QSS bound.

For quantitative growth prediction, compare the matched full and reduced poles over the
intended parameter range and apply the desired error tolerance there. If only an asymptotic
scale separation is available, state that assumption explicitly; $q\ll1$ does not by itself
guarantee a prescribed relative growth error near onset. No universal cutoff is inferred.

**Table S7. P2/P3 turnover factors required on the evaluated grid.**

"""+table(t[t.family.isin(["P2","P3"])])+"\n\n"+common+r"""

## Reproduction

Run `python scripts/reviewer_linear.py`, followed by
`python scripts/write_reviewer_reports.py` and `python figures/plot_reviewer_extensions.py`.
Full inputs and pole diagnostics are in `qss_families.csv`, `qss_validity.csv` and
`qss_accuracy_thresholds.csv` under `figures/data/reviewer_extensions/`.
"""
    write("qss_validity_supplement.md",5,"Accuracy of the matched quasi-steady-state reduction",qtext,"FigS7_qss_validity",
          "QSS validity. (A) P2/P3 relative growth error. (B) Error against the pole-to-relaxation ratio, coloured by eigenvalue condition. (C) Grid thresholds for all 52 families, sorted separately by tolerance. (D) Every family; repeated evaluations within a family are dependent.")
    s=summaries["linear"]["recalibration"]
    text=r"""
## Clock calibration and prospective prediction

An obvious control for a faster reduction is to slow its entire clock. We multiply the DNA
vector field by one positive scalar $c$, preserving its equilibrium set. At P2, matching the
leading K-mode growth gives
"""+f"$c={s['clock_at_P2']:.10f}$.\n"+r"""
We hold this scalar fixed across 161 equally spaced $\lambda_N$ values from 1.31 to 6.5,
augmented by exact P2/P3 values (163 conditions). All other parameters stay at baseline.
This is a prediction test for a **single global clock**, not a fit of all reduced-model
parameters or a claim that no three-variable model can be fitted to selected observations.

At P3, the P2-calibrated prediction is 16.51% below the full-model growth rate. A separate
least-squares clock fitted only on $1.5\leq\lambda_N\leq3$ is
"""+f"$c={s['training_clock']:.10f}$; in the held-out range $3.5\\leq\\lambda_N\\leq6.5$,\nits growth errors range from {100*s['heldout_LS_error_min']:.3f}% to {100*s['heldout_LS_error_max']:.3f}%.\n"+r"""
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
"""+common
    write("recalibration_supplement.md",6,"Can a fitted reduced-model clock recover the dynamics?",text,"FigS8_clock_recalibration",
          "Clock calibration. (A) P2-fitted clock across Notch production. (B) Relative errors for the P2 fit and training-range least-squares fit; shading marks training and held-out ranges. (C) All P2 modes. (D) Stable-P0 Delta-production response; both input and vector field are clock-scaled.")
    geom=read("geometry_runs");pairs=read("geometry_pairs");crit=read("geometry_stationary_thresholds")
    mainpairs=pairs[(pairs.condition=="periodic")&(pairs.point!="P0")]
    text=r"""
## Periodic size and commensurability

The original 24-cell side is divisible by three. To test dependence on this choice, we use
23, 24 and 25 cells per side, P2/P3, three PCG64 seeds (2026102001–2026102003), and full/DNA
models: 36 runs, 18 matched pairs. Each pair has identical initial $D,N,A$ fields with
relative amplitude $10^{-4}$, and full-model $S=S_q(A)$. Random fields are spatially centred
and normalised as in Note 3. Periodic graph sums retain six neighbours.

**Table S8. Actual minimum graph eigenvalue and stationary onset.**

"""+table(crit)+r"""

The maximum relative threshold change from 24 is approximately 0.253%. Thus exact
$\mu=-3$ commensurability changes onset slightly but is not required for instability or
the timing separation at these operating points. Across all 18 pairs, the full/DNA time
ratio to Delta CV 0.1 is
"""+f"{mainpairs.full_DNA_onset_ratio.min():.4f}–{mainpairs.full_DNA_onset_ratio.max():.4f}.\n"+r"""
Endpoint structures vary. To compare spectra across sizes, we report both exact-minimum-shell
power and power in the fixed band $\mu\leq-2.8$. Neither statistic implies a perfectly
three-sublattice final pattern. P0 is a stable periodic control in both models.

## Open edges and quenched heterogeneity

Six further runs use a 24-by-24 **open** graph at P0/P2/P3 in both models. Neighbours outside
the array are dropped; degrees range from two to six. Raw sums are retained without degree
normalisation. The periodic homogeneous reference is therefore **not an equilibrium** of
this open graph. P0 approaches a nonuniform boundary-forced state (Delta CV about 0.0698)
even though its periodic control returns to homogeneity. This is a useful negative control
against calling every nonzero CV a spontaneous instability. We report three-cell-rim and
interior CV separately. A degree-normalised open graph would be a different coupling model.

Eight runs add time-independent, mean-one lognormal Notch-production multipliers with
target CV 0.05 or 0.1 at P2/P3, using a separate PCG64(seed+1) stream; the realised multipliers
are stored and have sample mean exactly one. These are **quenched parameter differences**,
not dynamical noise. Full-model CV crossings remain later in all four pairs, but deterministic
forcing advances pattern contrast and reduces the time ratio. Neither open nor heterogeneous
crossing times are interpreted as homogeneous-instability thresholds.

All 52 trajectories use BDF with an analytic sparse Jacobian, rtol $10^{-10}$,
atol $10^{-13}$, maximum internal step 0.5, and recording interval 0.5. Horizons are
600 at P0/P2 and 300 at P3. All recorded concentrations are positive; the maximum final
scaled drift is $7.55\times10^{-10}$. These finite controls do not constitute a survey of
arbitrary contact graphs, tissue growth or stochastic dynamics.

## Reproduction

Run `python scripts/reviewer_geometry.py --reuse`, then the report and figure commands.
Omit `--reuse` for fresh trajectories. `geometry_runs.csv`, `geometry_pairs.csv`,
`geometry_stationary_thresholds.csv`, and all 52 initial/final field files are retained.
"""+common
    write("geometry_supplement.md",7,"Finite lattice, open boundaries and static heterogeneity",text,"FigS9_geometry_robustness",
          "Geometry controls. (A) Stationary thresholds. (B) All 18 paired periodic onset-time ratios. (C) Fixed near-K-band endpoint power at both operating points. (D) Stable P0: periodic decay versus open-edge forcing; open-graph contrast is not a homogeneous-instability test.")
    b=read("oscillatory_boundaries");c=read("real_pair_collisions")
    text=r"""
## Constant-gain continuation of every reference DC miss

The corrected reference sample has twelve unstable cases with $G(-3)<1$: six have an
unstable complex pair and six have two positive real poles at their original parameters.
All have dominant mode $\mu=-3$. To identify a route into these failures **within the tested
mode**, we follow $(\lambda_S,d_S)\mapsto\alpha(\lambda_S,d_S)$ at every one of the twelve
points. This holds the complete equilibrium and $G(\mu)$ fixed. We evaluate all 54 mode
blocks, including the uniform mode, over 401 logarithmically spaced $\alpha$ from
$10^{-6}$ to $10^2$, bracket the first stable-to-unstable crossing, and refine in $\log\alpha$.

All twelve first crossings occur in the K mode at nonzero frequency. The real-part residual
at a refined crossing is at most $1.57\times10^{-15}$; the derivative of the leading real
part with respect to $\log\alpha$ is positive in every case. Before and after crossing,
right-half-plane counts are computed over the 54 distinct blocks (without spatial
multiplicity). **We call these oscillatory linear boundaries or Hopf candidates.** We have
not computed a first Lyapunov coefficient or proved a nondegenerate nonlinear Hopf
bifurcation, the direction of a cycle branch, or cycle stability.

**Table S9. All first oscillatory boundaries on constant-gain paths.**

"""+table(b,["sample_id","alpha","mu","frequency","transverse_logalpha_slope","G_K"])+r"""

Six unstable complex pairs subsequently collide on the positive real axis before $\alpha=1$.
We locate this collision by the squared splitting of the leading pair. For example, sample
679 first crosses at $\alpha=0.1789253$ and frequency 0.1062911, then becomes a real pair at
$\alpha=0.7808506$, with mean real part 0.165840. Its $G(-3)=0.563539$ throughout.
Thus the real-pair misses need not enter instability through a zero eigenvalue on this path.
The Schur stationary condition and Nyquist unstable-pole count remain consistent.

**Table S10. Later unstable-pair collisions.**

"""+table(c)+r"""

At a collision, repeated eigenvalues are ill-conditioned. Stored tiny pair splittings are
numerical diagnostics, not a claim of exact root equality to that many digits. These paths
are constructed controls; no assertion is made that a biological tissue follows them or
that every possible route to the same endpoint must be oscillatory.

## Reproduction

`scripts/reviewer_linear.py` writes all K poles along the paths to
`bifurcation_poles.csv.gz`, the refined boundaries to `oscillatory_boundaries.csv` and
collisions to `real_pair_collisions.csv`. Use the report and figure commands in Note 5.
"""+common
    write("bifurcation_supplement.md",8,"Oscillatory routes to gain-below-one instability",text,"FigS10_oscillatory_routes",
          "Constant-gain routes. (A,B) K poles for samples 1730 and 679; blue shows real parts, orange absolute imaginary parts. Vertical lines mark the first oscillatory boundary and, where present, later real-pair collision. (C) Crossing frequency for every miss. (D) Boundaries and collisions at unchanged G below one.")
    ps=read("pattern_stability");s=summaries["pattern_stability"]
    text=r"""
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
"""+f"${s['max_newton_change']:.3g}$, and the maximum final scaled residual is ${s['max_refined_residual']:.3g}$.\n"+r"""

At each refined state we form the **complete patterned-state Jacobian**, with all local
derivatives evaluated at the heterogeneous concentrations and the original contact matrix.
No homogeneous mode decomposition is assumed here. Dense eigensolvers compute all 2304 full
or 1728 DNA poles. A separate shifted sparse solve recovers a leading eigenpair and checks
its normalised residual. The worst dense/sparse leading-pole difference is
"""+f"${s['max_sparse_dense_error']:.3g}$ and the largest normalised sparse residual is ${s['max_eigenpair_residual']:.3g}$.\n"+r"""

All twelve spectra are strictly in the left half-plane: leading real parts range from
"""+f"{s['leading_real_min']:.8f} to {s['leading_real_max']:.8f}. These selected equilibria are therefore\nlocally asymptotically stable for their respective deterministic lattice equations.\n"+r"""
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
"""+f"${s['recovery_ratio_max']:.3g}$. This floor is controlled by floating-point/solver error;\nwe do not infer a physical decay rate from the final few digits.\n"+r"""

**Table S11. Complete-spectrum and recovery results for all twelve selected states.**

"""+table(ps.assign(index=np.arange(len(ps))),["index","point","model","dimension","leading_real","n_RHP","recovery_final_over_initial"])+r"""

## Reproduction

Run `python scripts/reviewer_pattern_stability.py --reuse` (omit `--reuse` to recompute).
All complete complex spectra, initial/refined fields and recovery curves are retained in
`figures/data/reviewer_extensions/pattern_stability/`. The selection table and
`pattern_stability.csv` expose every case rather than only the most stable example.
"""+common
    write("pattern_stability_supplement.md",9,"Local stability of selected nonlinear patterns",text,"FigS11_patterned_state_stability",
          "Patterned-state stability. (A) Leading poles from complete spectra, with table-order endpoint indices. (B) Paired perturbation recovery; blue full, orange DNA. (C) All poles of the pure-mode full P2 state; tiny imaginary parts are at numerical round-off scale, not evidence of oscillation. (D) Newton and independent sparse-eigenpair residuals.")
    pri=read("prior_sensitivity");strata=read("growth_boundary_strata")
    text=r"""
## Small denominators and absolute growth differences

The 5.4-fold median ratio in the primary robustness study could be inflated by tiny full-model
growth rates. We therefore reuse all 40000 frozen primary inputs and examine the 14949
full-model-unstable points with three complementary statistics: the growth ratio,
the absolute gap $\lambda_{DNA}-\lambda_F$, and the bounded symmetric difference
$2(\lambda_{DNA}-\lambda_F)/(\lambda_{DNA}+\lambda_F)$. We also define a relative growth
margin $\eta=\lambda_F/\max_j|s_j(J_{-3})|$ and apply prespecified cuts in $\eta$ and raw
$\lambda_F$. The largest K-pole magnitude supplies a local characteristic rate, not a
new instability criterion.

**Table S12. Growth differences after excluding points near the boundary.**

"""+table(strata)+r"""

The ratio median decreases from 5.4119 to 3.4532 when $\eta\geq0.05$ (4048 points).
The median absolute gap in that subset is 0.8231 in the model's inverse dimensionless time.
Thus proximity to onset affects the magnitude but does not by itself explain the sampled
growth discrepancy. A ratio is never used as evidence of experimental fold accuracy.

## Prior sensitivity by importance reweighting

We retain the original support and all stored inputs; no new draws are made. In addition to
the reference independent log-uniform rates and uniform Hill exponents, we apply two
illustrative density ratios: (i) independent triangular Hill densities on [1,6] with mode 3;
(ii) a Gaussian copula with correlation 0.75 between log proneural production and log
proneural loss, preserving their original log-uniform marginals. The latter makes a strong
production/turnover relationship explicit without changing marginal ranges. Neither prior
is estimated from biology.

Each seed's instability fraction is a self-normalised weighted estimate. We average the ten
seed estimates and use a Student-t interval with nine degrees of freedom for Monte Carlo
variation across these designs. This is not uncertainty about biological prevalence.
The pooled effective sample size is $(\sum w)^2/\sum w^2$; we report it because reweighting
reduces information. Growth medians use pooled weights only among full-unstable points.

**Table S13. Prior-sensitive instability and growth summaries.**

"""+table(pri,formats={k:lambda x:f"{100*x:.3f}%" for k in ("fraction","mc95_low","mc95_high")})+r"""

Instability estimates vary from 37.37% to 39.13% and the weighted median growth ratio from
5.245 to 5.961 for these measures. The bounds study in Note 2 produces a larger fraction
change (to 28.89%). These percentages are properties of specified measures; they must not
be called the fraction of all possible biological systems. Importance sampling does not
test behaviour outside the reference support, and these two alternatives do not exhaust
reasonable priors.

## Reproduction

`scripts/reviewer_linear.py` writes `growth_error_points.csv.gz`,
`growth_boundary_strata.csv`, `prior_replicates.csv` and `prior_sensitivity.csv` from the
checksum-verified primary robustness cases. Use the report and figure commands in Note 5.
"""+common
    write("growth_prior_supplement.md",10,"Boundary proximity and parameter-measure sensitivity",text,"FigS12_growth_and_prior_sensitivity",
          "Growth and prior sensitivity. (A) Reweighted fractions with ten-design Monte Carlo intervals. (B) Ratio medians and IQRs after relative-margin cuts. (C) Absolute growth gap. (D) Bounded symmetric growth difference; all unstable primary points are retained.")
    text=r"""
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
"""+common
    write("mechanism_supplement.md",11,"What linear data identify about downstream memory",text,"FigS13_mechanism_identifiability",
          "Mechanism discrimination. (A–C) Normalised production after clamped-A steps of 1.001, 1.2 and 2 times baseline. Full and alternative models share equilibrium and linear D,N,A response but differ at finite amplitude. (D) Stable-P0 frequency-response differences are at floating-point precision; the exact similarity identity explains their agreement.")
    print("wrote seven verified supplementary notes (5-11)")


if __name__=="__main__":main()
