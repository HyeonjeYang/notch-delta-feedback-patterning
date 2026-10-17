---
title: "Explicit downstream dynamics separate stationary onset from growth in Notch–Delta patterning"
author: "Hyeonje Yang"
date: |
  Department of Biological Sciences, Seoul National University, Seoul, Republic of Korea
  (yanghj0206@snu.ac.kr). Revised manuscript draft, October 2026.
geometry: margin=2.5cm
fontsize: 11pt
linestretch: 1.4
header-includes:
  - \usepackage{lineno}
  - \linenumbers
  - \usepackage{float}
  - \floatplacement{figure}{H}
  - \usepackage{caption}
  - \captionsetup[figure]{labelformat=empty}
---

# Abstract

Lateral-inhibition models often compress the steps downstream of Notch activation. We ask
which predictions change when activated Notch and proneural activity remain separate
dynamical variables. We analyse a four-variable Delta–Notch–activated-Notch–proneural model
on a periodic six-neighbour contact lattice. The homogeneous equilibrium is unique for all
positive parameters. We prove that the static feedback gain is largest for the three-sublattice
mode, which therefore selects the first stationary instability. Matched quasi-steady-state
reductions preserve every equilibrium and the zero-eigenvalue boundary exactly. Their
dynamic predictions differ. Rescaling proneural turnover leaves the stationary threshold fixed
while growth varies over three decades. In ten independent 4000-point Latin-hypercube samples,
the three-variable reduction grows faster at every full-model-unstable point, with a pooled
median growth-rate ratio of 5.4. It also predicts 216 instabilities absent from the full model.
The mean nonuniform-instability fraction is 37.37%, with an approximate 95% Monte Carlo
interval of 37.07–37.67% under the specified parameter measure. Checks cover 115000 new points;
effect sizes depend on the sampling measure. Paired random starts retain slower full-model
emergence but yield different endpoints. A reduced-model clock fitted at one operating point
underpredicts growth by 16.5% at another. Constant-gain continuations explain all twelve
reference gain-below-one instabilities through oscillatory linear boundaries. A matched
Delta-production lag is exactly equivalent for linear Delta, Notch and activated-Notch
responses, yet differs under large inputs. Thus finite feedback memory affects growth and
non-stationary stability, but these linear observables cannot uniquely identify its
biochemical location. Nonlinear responses or proneural measurements can distinguish the
specified alternatives.

# 1. Introduction

In lateral inhibition, Delta on one cell activates Notch on its neighbours, and Notch signalling
represses the proneural activity that drives Delta expression. Small differences between
adjacent cells are amplified into a spatial pattern of signal-sending and signal-receiving
states [1, 2]. In *Drosophila* sensory-organ development this antagonism acts through
Achaete–Scute proneural self-stimulation and Notch/Enhancer-of-split repression [3–6].

Mathematical models resolve different parts of this process. Classical feedback models use
effective Notch and Delta variables [1], whereas cis/trans models distinguish receptor, ligand
and a downstream reporter [2, 7]. Tissue-scale work has already connected Notch dynamics and
proneural patterning to live observations [6]. We therefore do not claim that either downstream
regulation or its developmental timing is a new modelling topic. The coupled system is
decomposed into eigenmodes of the cell-contact graph, and each mode is tested separately
[1, 10, 13]. Under the cell-stability and monotonicity assumptions of graph-based analyses,
a static input–output characteristic also gives instability and patterned-state results
[11, 12]. These assumptions are substantive; a determinant factorisation alone does not
establish them. On the triangular contact graph of a hexagonal
tissue, the most negative graph eigenvalue belongs to a three-sublattice pattern, not a
checkerboard [13]. Two-dimensional contact coupling has also been analysed directly [14].

These results leave a question open. When the activated-Notch and proneural steps are kept as
separate dynamical variables, which predictions actually differ from those of a reduced model?
Retaining more variables gives more parameters to vary. It does not by itself show that the
extra dynamics matter for patterning.

We address this question by holding the entire steady-state response fixed while changing
the downstream relaxation time, and by comparing with reductions derived from the same
equations. This controls a confound in comparing unrelated four- and three-variable models:
different static nonlinearities could otherwise explain their different predictions. The
mode decomposition and Schur-complement factorisation are standard tools. Our main result
is a quantitative account of the dynamic information lost in an exactly matched reduction;
the following model-specific results establish the comparison:

1. **Uniqueness.** The homogeneous steady state is unique for every positive parameter set.
2. **Mode selection.** The static loop gain is largest at the most negative graph eigenvalue for
   every parameter set. Stationary onset therefore always selects the three-sublattice mode,
   even though cis inhibition and trans depletion of ligand make the gain a non-trivial function
   of the mode.
3. **Shared onset.** Any matched reduction that slaves the downstream layer to quasi-steady state
   has the same stationary boundary as the explicit model.
4. **Different dynamics.** The explicit layer changes the pattern growth rate, often by large
   factors, and changes the set of instabilities that do not pass through a zero eigenvalue.
   We identify the sampled parameter sets at which a static test fails, and explain each one.

# 2. Model and spatial modes

## 2.1 Equations

Each cell $i$ carries four concentrations: Delta ligand $D_i$, Notch receptor $N_i$, activated
Notch (NICD-like signalling activity) $A_i$, and proneural activity $S_i$. All are
dimensionless, as is time $t$. The dynamics are

$$
\begin{aligned}
\dot D_i &= \lambda_D\,h_D(S_i) - d_D D_i - f_D D_i\,s_{N,i} - k_c N_i D_i,\\
\dot N_i &= \lambda_N - d_N N_i - f_N N_i\,s_{D,i} - k_c N_i D_i,\\
\dot A_i &= -d_A A_i + f_N N_i\,s_{D,i},\\
\dot S_i &= -d_S S_i + \lambda_S\,h_S(A_i),
\end{aligned}
\qquad (1)
$$

where $s_{D,i}=\sum_{j\sim i}D_j$ and $s_{N,i}=\sum_{j\sim i}N_j$ are sums over the six
contact neighbours of cell $i$. The production terms use Hill functions:
$h_D(S)=1/[1+(K_d/S)^{n_D}]$ is increasing (proneural activity drives Delta) and
$h_S(A)=1/[1+(A/K_d)^{n_S}]$ is decreasing (activated Notch represses proneural activity).

The 13 positive parameters are:

* $\lambda_D,\lambda_N,\lambda_S$: maximal production rates of Delta, Notch and proneural activity;
* $d_D,d_N,d_A,d_S$: first-order loss rates;
* $f_N$: trans activation of Notch by neighbouring Delta, which consumes the receptor and produces $A$;
* $f_D$: trans depletion of Delta by neighbouring Notch;
* $k_c$: cis mutual inactivation of Notch and Delta in the same cell;
* $K_d$: the half-saturation constant shared by both Hill functions;
* $n_D,n_S$: the Hill coefficients.

The baseline parameter set $p_0$ (Table 1) is a literature-inspired operating point that is not
fitted to a particular tissue. The main continuation parameter is the Notch production rate
$\lambda_N$, with all other parameters held at $p_0$. We refer to three reference points:
$p_0$ ($\lambda_N=1$), P2 ($\lambda_N=2$) and P3 ($\lambda_N=5.623$).

## 2.2 Contact graph and spatial modes

Cells occupy a $24\times24$ periodic array. Cell $(r,c)$ contacts the six cells at offsets
$(\pm1,0)$, $(0,\pm1)$, $(-1,+1)$ and $(+1,-1)$, which is the contact graph of a hexagonal
tiling (Fig. 1B). Let $W$ be the adjacency matrix of this graph ($W_{ij}=1$ if cells $i$ and $j$
touch). Because $W$ is symmetric and every cell has $z=6$ neighbours, its eigenvalues $\mu$ lie in
$[-6,6]$. Its eigenvectors are plane waves indexed by integers $(m,n)$, $0\le m,n<24$, with
wavevector components $k_1=2\pi m/24$ and $k_2=2\pi n/24$ and

$$
\mu(m,n)=2\cos k_1+2\cos k_2+2\cos(k_1-k_2) \qquad (2)
$$

(Fig. 1D). The uniform mode $(m,n)=(0,0)$ has $\mu=6$. The minimum $\mu=-3$ occurs only at
$(m,n)=(8,16)$ and $(16,8)$. The real combinations of these two plane waves divide the lattice
into three interpenetrating sublattices (Fig. 1C). On this non-bipartite graph, this
three-sublattice pattern plays the role that a checkerboard plays on a square lattice [13].

Let $x^*=(D^*,N^*,A^*,S^*)$ be a homogeneous steady state of (1). A perturbation of the form
$\delta x_i(t)=v_i\,\xi(t)$, where $v$ is an eigenvector of $W$ with eigenvalue $\mu$, obeys
$\dot\xi=J_\mu\xi$ with a $4\times4$ matrix $J_\mu$ (Methods, Eq. 6). The neighbour sums enter
$J_\mu$ only through $\mu$, so the stability of the 2304-dimensional lattice linearisation
reduces to the eigenvalues of one $4\times4$ block per distinct graph eigenvalue. There are 54
distinct values of $\mu$ on the $24\times24$ lattice. We write

$$
\lambda_\mathrm{pattern}=\max_{\mu\neq 6}\ \max\,\mathrm{Re}\,\mathrm{eig}(J_\mu)
$$

for the growth rate of the fastest nonuniform mode. The homogeneous state is unstable to
spatial patterning when $\lambda_\mathrm{pattern}>0$.

## 2.3 The homogeneous steady state is unique

At a homogeneous state, $s_D=6D$ and $s_N=6N$. Solving the $N$, $A$ and $S$ equations in turn
expresses $N$, $A$ and $S$ as monotone functions of $D$, and the steady state is a root of

$$
\begin{aligned}
F(D)&=\lambda_D h_D\!\big(S(D)\big)-D\,[\,d_D+(6f_D+k_c)N(D)\,],\\
N(D)&=\frac{\lambda_N}{d_N+(6f_N+k_c)D},\qquad
S(D)=\frac{\lambda_S}{d_S}\,h_S\!\Big(\frac{6f_N N(D)\,D}{d_A}\Big).
\end{aligned}
\qquad (3)
$$

$F$ is strictly decreasing, with $F(0)>0$ and $F(\lambda_D/d_D)<0$ (Supplementary Note 1).
**For every positive parameter set the homogeneous steady state therefore exists, is unique, and
lies in $0<D^*<\lambda_D/d_D$.** We compute it by Brent bracketing on that interval.

This corrects the earlier version of this analysis, which used a multistart nonlinear solver and
counted solver failures as "no steady state". In the 4000-point parameter sample (Section 3.2),
733 points had been discarded in this way. All of them have a steady state. Over all 4000 points,
the maximum relative residual is below $3.2\times10^{-15}$, using
$\|f(x^*)\|_\infty/\max(1,\|x^*\|_\infty)$.

# 3. Results

## 3.1 The static loop gain fixes the stationary onset and the selected mode

To identify which part of the circuit controls stability, we split each mode's perturbation
into a "plant" $\zeta=(\delta D,\delta N,\delta A)$ and the proneural perturbation $\delta S$:

$$
\dot\zeta = A_\mu\zeta + b\,\delta S,\qquad \delta\dot S = c^{\top}\zeta - d_S\delta S. \qquad (4)
$$

Here $A_\mu$ is the upper-left $3\times3$ block of $J_\mu$. The vector $b$ carries proneural
activity into Delta production, $b=(\lambda_D h_D'(S^*),0,0)^{\top}$. The vector $c$ carries
activated Notch into proneural repression, $c^{\top}=(0,0,-\lambda_S|h_S'(A^*)|)$. The
open-loop transfer function around the $S\to D\to A\dashv S$ loop is
$L_\mu(s)=c^{\top}(sI-A_\mu)^{-1}b/(s+d_S)$, and its zero-frequency value

$$
G(\mu)=L_\mu(0)
$$

is the **static loop gain** of mode $\mu$. The Schur complement gives the standard identity

$$
\det(sI-J_\mu)=\det(sI-A_\mu)\,(s+d_S)\,[\,1-L_\mu(s)\,], \qquad
\det J_\mu=\det(A_\mu)\,(-d_S)\,[\,1-G(\mu)\,]. \qquad (5)
$$

Three properties make this identity informative here (proofs in Supplementary Note 1):

**(i) The plant is always stable.** $A_\mu$ is Hurwitz (all its eigenvalues have negative real part) for every $\mu\in[-6,6]$ and every
positive parameter set. Cis inactivation and trans depletion couple $D$ and $N$, but this
subsystem cannot destabilise any mode on its own when Delta production is held fixed. Every
instability involves the proneural loop.

**(ii) The static gain counts eigenvalue parity.** Because $\det(A_\mu)(-d_S)>0$, the sign of
$1-G(\mu)$ is the sign of $\det J_\mu$. So $G(\mu)=1$ exactly when $J_\mu$ has a zero
eigenvalue, and $G(\mu)>1$ exactly when $J_\mu$ has an odd number of positive real
eigenvalues. $G(\mu)<1$ is *not* a stability certificate. It is compatible with a complex pair,
or two real eigenvalues, in the right half-plane (Section 3.3). The complete test for a mode is
the Nyquist criterion [15], away from imaginary-axis eigenvalues. Since the open loop has no unstable poles, the number of unstable
eigenvalues of $J_\mu$ equals the number of clockwise encirclements of $+1$ by $L_\mu(i\omega)$.

**(iii) The static gain is largest at the most negative graph eigenvalue.** Writing out the
cofactors gives $G(\mu)=K\,n(\mu)/\Delta(\mu)$. Here $K>0$, $n(\mu)=6k_cD^*-\mu(d_N+k_cD^*)$
is linear and decreasing in $\mu$, and $\Delta(\mu)>0$ is the determinant of the $D$–$N$ block,
which is quadratic in $\mu$. A concavity argument (Supplementary Note 1) shows that $G(\mu)<G(\mu_{\min})$ for every
$\mu>\mu_{\min}$, and that $G(6)<0$. This holds for every positive parameter set and every
regular contact graph. Two consequences follow. Along any parameter path that starts from a
stable homogeneous state, **the first stationary instability is always in the
three-sublattice mode** ($\mu=-3$). The uniform mode can never lose stability through a zero
eigenvalue, which agrees with uniqueness of the steady state. In single-channel
lateral-inhibition models $G(\mu)$ is simply proportional to $-\mu$ and this ordering is
immediate [1, 11]. Here the neighbour coupling has two channels, trans activation and trans
depletion, and with cis inactivation $G(\mu)$ becomes a ratio of linear and quadratic
functions of $\mu$. Part (iii) shows that these extra channels never change the selected mode.

Fig. 2A shows $G(\mu)$ across the full spectrum at four values of $\lambda_N$. The maximum is
always at $\mu=-3$, and $G$ is negative near $\mu=6$. As $\lambda_N$ increases, the leading
eigenvalue of $J_{-3}$ crosses zero (Fig. 2B) at the same point where $G(-3)$ crosses one
(Fig. 2C). The two crossings, located by separate root-finding, give
$\lambda_N^c=1.304815343$ and agree within $10^{-12}$, as identity (5) requires. The
uniform mode keeps a leading eigenvalue of $-0.100$ over $\lambda_N\in[0.5,6.5]$. At P2,
$G(-3)=1.809$ and the three-sublattice mode grows at rate $0.103067$.

## 3.2 Which kinetic steps set the static margin

Because the stationary boundary is $G(-3)=1$, the parameters that move it are those that change
the static gain. $G(\mu)$ factorises exactly into three stage gains (Methods, Eq. 7):

* $\kappa_D=\lambda_D h_D'(S^*)$, from proneural activity to Delta production;
* $T_\mu$, the static response of activated Notch in the neighbouring cells to Delta
  production in mode $\mu$, including the $N$-mediated cis and trans paths;
* $\kappa_S=\lambda_S|h_S'(A^*)|/d_S$, from activated Notch to proneural repression.

Each logarithmic sensitivity $\partial\ln G(-3)/\partial\ln p$ splits accordingly. The split
includes each parameter's indirect effect through the shift of the steady state (Fig. 3A).

At $p_0$ the largest static sensitivities belong to the downstream parameters:

| parameter | $\partial\ln G(-3)/\partial\ln p$ |
|---|---|
| $n_S$ | $+6.9$ |
| $d_A$ | $-5.3$ |
| $\lambda_N$ | $+5.3$ |
| $K_d$ | $-3.5$ |
| $\lambda_S$ | $-1.8$ |
| $d_S$ | $+1.8$ |

The cis/trans parameters $f_N$, $k_c$ and $f_D$ are below 0.5 in magnitude. The stage
decomposition shows why the downstream parameters dominate. Each acts on more than one stage at
once, through the steady state. For example, raising $n_S$ steepens proneural repression and
also moves $S^*$ onto a steeper part of $h_D$. One-parameter continuations over a 100-fold range
cross the boundary only for $K_d,d_A,d_S,\lambda_S,\lambda_N$ and $n_S$ (Fig. S2). On
three two-parameter maps, $G(-3)>1$ and $\lambda_\mathrm{pattern}>0$ give identical
classifications at all 1875 grid points; two maps are shown in Fig. 3C,D and all three are
in the deposited data. This ranking describes the present operating point. It does not claim
that cis/trans interactions are unimportant in Notch biology [2, 7, 8].

The comparison in Fig. 3B anticipates Section 3.4. Static quantities depend on $\lambda_S$ and
$d_S$ only through their ratio $\lambda_S/d_S$, so their static sensitivities are equal and
opposite ($\mp1.80$). Their effects on the growth rate $\lambda_\mathrm{pattern}$ are not
($-0.060$ and $+0.044$), because $d_S$ also sets how fast the proneural stage relaxes.

**Reference sample.** Over the original 4000-point Latin-hypercube sample of all 13 parameters (each rate
within $[0.1,10]\times p_0$, log-uniform; Hill coefficients in $[1,6]$), **1453 points (36.3%)
are unstable to a spatial mode.** At every one of them the dominant mode is $\mu=-3$. The
earlier figure of 30.3% excluded the 733 points discussed in Section 2.3. The fraction depends
on the sampling measure and is not a biological probability. This is the reference sample
used in Figs. 4 and 5B; independent replication is reported in Section 3.5.

## 3.3 Where the static gain is not enough

$G(-3)>1$ and $\lambda_\mathrm{pattern}>0$ agree at 3988 of the 4000 sampled points (Fig. 4A). At
each of the 12 exceptions, $G(-3)<1$ while the $\mu=-3$ block has **exactly two**
right-half-plane eigenvalues. Six are complex pairs and six are pairs of real eigenvalues (Fig. 4C).
An even number of unstable eigenvalues leaves $\det J_{-3}>0$, so by (5) the static gain stays
below one. The mismatches are therefore not cases in which "another mode dominates". In all
12, as at every unstable point, the dominant mode is $\mu=-3$. They are instabilities of the
tested mode that are invisible at zero frequency.

Point 1730 is an example. There $G(-3)=0.366$, but $J_{-3}$ has eigenvalues
$0.159\pm0.262i$. Its Nyquist curve crosses the real axis at $L_{-3}(i\omega)=1.56$ when
$\omega=0.36$, and together with its mirror image it encircles $+1$ twice (Fig. 4B). The
stationary instability at P2 encircles $+1$ once, and $p_0$ not at all. Counting
encirclements reproduces the eigenvalue count at all 4000 points. The parameter sets of these
12 points are listed in the deposited data.

**Is another pattern ever selected?** By property (iii), no other nonuniform mode can be the
first to undergo a *stationary* instability. Non-stationary instabilities are not covered by
that argument, so we searched 40 000 further random parameter sets: 20 000 from the same box
and 20 000 with rate/constant bounds extended to $[0.01,100]\times p_0$
(Hill coefficients still in $[1,6]$). Of these, 13 048 were unstable to a spatial
mode. In none was a mode other than $\mu=-3$ dominant, or unstable while $\mu=-3$ was stable.
Dominance does not imply that other modes are stable. In the original 4000-point sample,
1451 points also have unstable nonuniform modes with $\mu>-3$. At P2, 21 distinct nonuniform
graph eigenvalues have positive growth rates, with the maximum at $\mu=-3$ (Fig. S3).
The uniform mode additionally has an unstable complex pair at 45 points, corresponding to
synchronous oscillation. All 45 points are also unstable to the three-sublattice mode.

## 3.4 What the explicit downstream layer adds

To isolate the role of the downstream dynamics, we compared the full model with two matched
reductions (Methods):

* **DNA:** proneural activity is slaved to activated Notch, $S=(\lambda_S/d_S)h_S(A)$. This
  gives three variables with receptor, ligand and activated-Notch roles, comparable in level
  of description to [2, 7]. Its composite Delta-production law is derived here, rather than
  copied from those models; it is an exactly matched control, not a claim of equation-level
  equivalence to a published three-component model.
* **DN:** $A$ is also slaved, giving two variables.

Both reductions only impose $\dot S=0$ (and $\dot A=0$) cell by cell, so they have exactly the
same steady states as (1), homogeneous or patterned. Their mode-resolved Jacobians are Schur
complements of $J_\mu$, and $\det J_\mu$ equals $-d_S\det J^\mathrm{DNA}_\mu$ and
$d_Ad_S\det J^\mathrm{DN}_\mu$. **The full model and both reductions therefore share the
stationary boundary, the critical $\lambda_N$, and the selected mode.** Section 3.1 and the
boundary in Section 3.2 do not need the explicit downstream dynamics. A reduced model with the
same static input–output response predicts them exactly.

Here "stationary boundary" means the zero-eigenvalue set. It is an actual boundary of the
stable region only where an oscillatory instability has not already intervened.

The explicit dynamics matter for timing. Consider the one-parameter family
$(\lambda_S,d_S)\to\alpha(\lambda_S,d_S)$. It leaves $S=(\lambda_S/d_S)h_S(A)$, every steady
state and every $G(\mu)$ unchanged, and only rescales how fast proneural activity relaxes. As
$\alpha\to\infty$ the model approaches the DNA reduction. As $\alpha\to0$ the plant equilibrates
quickly relative to $S$, and the slow eigenvalue is asymptotic to
$\alpha d_S^{(0)}[G(\mu)-1]$, where $d_S^{(0)}$ is the unscaled turnover rate
(Supplementary Note 1). This approximation describes the leading eigenvalue in that limit;
at finite timescales, $G<1$ does not guarantee stability. Along this family at P2,
$\lambda_N^c=1.304815343$ for every $\alpha$
from $10^{-2}$ to $10^{4}$, while the growth rate of the three-sublattice mode ranges from
$8\times10^{-4}$ to $1.59$ (Fig. 5A). At the baseline timescales ($\alpha=1$) the explicit model
lies close to the slow-proneural limit: it grows at $0.103$, against $0.081$ from
$d_S[G-1]$. The DNA and DN reductions predict $1.597$ (15.5-fold faster) and $11.28$
(109-fold).

The same direction of bias occurs in the reference sample (Fig. 5B). The DNA reduction is unstable wherever the
full model is unstable, and grows faster at every one of those points (median 4.9-fold,
interquartile range 2.6–13.0). It is also unstable at 26 points where the full model is stable,
all with $G(-3)<1$. These are instabilities not certified by the static-gain threshold that finite proneural
relaxation suppresses in the full model. For the DN reduction the median overestimate is 216-fold, with 47 spurious
instabilities.

Nonlinear simulations at P2 test the timing difference and compare the final pattern (Fig. 5C,D).
We started each model at P2 on the
exact $\mu=-3$ eigenvector of its own Jacobian, at relative Delta amplitude $10^{-6}$ (Methods).
The amplitude of the three-sublattice component then grows at the predicted rate from $t=0$:

| model | predicted | fitted | fitting window |
|---|---|---|---|
| full | 0.103067 | 0.103077 | $[0, 67.25]$ |
| DNA reduction | 1.597479 | 1.597498 | $[0, 4.5]$ |

The amplitude-based rule for the fitting window was specified before fitting; its endpoint
is determined by each trajectory and drawn in the figure. The pattern reaches 10% relative
amplitude at $t=110.8$ in the full model and at $t=7.3$ in the reduction. The two models
then settle onto the same three-sublattice steady state, with a maximum difference of
$9\times10^{-11}$ in the final Delta field. This agreement is specific to these P2 simulations
and their initial states. Identical equilibrium sets do not imply identical attractor stability
or basins of attraction for the full and reduced models.

## 3.5 Robustness to sampling seed, sample size and parameter bounds

To distinguish reproducibility of one sample from robustness to sampling, we fixed a further
protocol before computation (Methods; Supplementary Note 2). Ten independent scrambled LHS
designs of 4000 points use the original parameter bounds. They give 14949 nonuniform-unstable
points among 40000 evaluations: a mean fraction of **37.37%**, with an approximate 95% Monte
Carlo interval of **37.07–37.67%** across independent replicates (Fig. S4A). Individual fractions
range from 36.650% to 37.925%. The earlier 36.325% is the result of its particular reference
sample, not an invariant fraction. Five designs at each of 1000, 2000 and 8000 points give
mean fractions of 37.10%, 37.69% and 37.39%, respectively (Fig. S4B). These size checks are
consistent with the primary estimate but do not establish a formal convergence rate.

The dynamic comparison persists. The pooled DNA/full growth-ratio median is **5.41** in the
40000-point primary study; the ten replicate medians range from 5.25 to 5.61 (Fig. S4C).
The DNA reduction grows faster at every full-model-unstable point and misses none. It is
unstable at 216 points where the full model is stable, against 445 reduced-only points for
DN. The static-gain test misses 150 instabilities in the primary study (Fig. S4D).

Expanding rate/constant bounds from $[0.1,10]$ to $[0.01,100]\times p_0$, while retaining Hill
coefficients in $[1,6]$, changes the sampled question. Five independent 4000-point designs
then give a mean unstable fraction of 28.89% (approximate 95% Monte Carlo interval
28.58–29.20%). DNA replicate-median growth ratios fall to 2.99–3.64, whereas DN ratios rise
to 1493–2019. The direction of the dynamic bias persists, but its magnitude and the unstable
fraction depend on the sampling measure.

Across all 115000 new evaluations, 41307 points are unstable and 392 are missed by the
static threshold. Every miss has two RHP eigenvalues in the tested K-mode block. No faster
competing nonuniform mode, full-model instability missed by a reduction, or reduced growth
rate at or below the full rate was found. These are empirical findings over the stated
ranges, not global theorems about dynamic spectra. One wider-box point has growth
$4.77\times10^{-9}$ and is explicitly flagged as near marginal; its positive sign persists
at 80-digit precision. Numerical threshold choices do not affect the primary results.

## 3.6 Random initial states preserve the timing difference, not endpoint equality

The pure-mode comparison in Figure 5 suppresses competing spatial modes. We therefore
ran 48 full/DNA pairs at P2 and P3: eight seeds at each of three relative perturbation
amplitudes, $10^{-6},10^{-4},10^{-2}$ (Methods; Supplementary Note 3). The models receive
identical broad-spectrum D, N and A fields; full-model S initially equals $S_q(A)$, the
same concentration reconstructed by the reduction. These pairs begin with identical
physical fields, rather than each model's own eigenvector. A further 16 full-model runs
randomise S independently at amplitude $10^{-4}$.

We measure spatial contrast by the Delta coefficient of variation,
CV$(D)=\operatorname{std}(D)/\operatorname{mean}(D)$. All 112 trajectories cross both
spatial-contrast thresholds, CV$(D)=0.1$ and 0.5.
The first CV=0.1 crossing occurs later in the full model in every main pair (Fig. S5A,B).
The paired full/DNA onset-time ratio spans 10.96–15.10 at P2 and 12.34–12.90 at P3.
The same ordering holds in all 48 pairs at CV=0.5 (ratios 11.87–14.84 at P2 and
12.08–12.77 at P3). These are operational nonlinear emergence times, not estimates of an eigenvalue ratio.
Independent initial S variation changes full-model onset times by 10–22% relative to the
QSS-initialised controls; timing therefore depends on initial biochemical state as well as
the asymptotic growth rate. On the 69 trajectories with at least two observations within
the linear-window rule, the largest relative Delta error against the all-mode propagator
is 0.180%; no single exponential is fitted to these mixed-mode transients.

The endpoint comparison changes the scope of the pure-mode example. At fixed horizons
$t=600$ (P2) and $t=300$ (P3), full and DNA Delta fields differ in all 48 pairs even after
optimising over periodic translations (Fig. S5C,E,F). Their aligned correlations span
0.584–0.838 at P2 and 0.866–0.953 at P3. Exact K-mode power fractions across the 96
main-pair endpoints span 0.0023–0.424 (Fig. S5D): a fastest linear K mode does
not force a pure three-sublattice nonlinear pattern. All endpoints satisfy the numerical
stationary-candidate criterion $\max|\dot x|/x^*<10^{-7}$; this residual test does not
establish attractor stability or a global basin classification. Thus the timing difference
survives these random starts, while identical final fields in Figure 5 remain specific to
its symmetry-restricted initial conditions.

Representative time-discretisation checks support this comparison (Supplementary Note 4,
Fig. S6). Fourteen initial states are reintegrated under five step-cap/tolerance/solver
settings, for 70 full-horizon trajectories. The largest relative change in either emergence
time is 0.01631%, and the largest raw endpoint relative L2 difference across all species is
$1.12\times10^{-10}$. These checks cover the stated initial states, rather than proving
solver independence for every possible basin or parameter vector.

## 3.7 Accuracy, recalibration and the location of feedback memory

The dynamic difference raises three further questions: when is QSS accurate, can a simple
recalibration remove the error, and do these results uniquely identify proneural memory?
We address them with matched controls (Fig. 6; Supplementary Notes 5–11).

In 52 compensated-turnover families, the DNA growth error falls as proneural relaxation
accelerates. On the evaluated grid, P2 requires $\alpha\geq199.53$ for less than 10% error
and P3 requires $\alpha\geq398.11$ (Fig. 6A). These are family-specific grid thresholds.
The ratio of the leading pole to the proneural loss rate is not by itself a relative-error
bound: 27 of 677 evaluations with $|s_{\rm lead}|/(\alpha d_S)\leq0.01$ still have more
than 10% growth error. Near onset, a small feedback substitution error can produce a large
relative growth error (Note 5). A single clock multiplier fitted to the DNA growth at P2
improves agreement but predicts P3 growth 16.51% below the full value (Fig. 6B).
A least-squares clock trained on $\lambda_N\in[1.5,3]$ retains 6.77–7.48% underprediction
in the held-out range $[3.5,6.5]$ (Note 6). This tests one global clock; it does not rule out
refitting a more flexible reduced model.

Constant-gain turnover paths also explain how the reference DC misses enter instability.
All twelve first lose stability through a nonzero-frequency K-mode crossing, with positive
transversality. Six subsequently turn into positive real pairs (Fig. 6C; Note 8). These are
oscillatory linear boundaries, or Hopf candidates; no nonlinear Hopf nondegeneracy or cycle
stability is established. Their gain remains below one along the entire path.

To test biochemical specificity, replace $S$ by a Delta-production variable satisfying
$\dot U=d_S\{\lambda_Dh_D[S_q(A)]-U\}$ and use $U$ in the Delta equation. This alternative
shares every equilibrium. Its homogeneous Jacobian satisfies the exact similarity
$J_\mu^{\rm prod\ lag}=T J_\mu^{\rm full}T^{-1}$ with
$T=\operatorname{diag}(1,1,1,\lambda_Dh_D'(S^*))$ (Note 11).
Consequently, **linear D, N and A responses cannot distinguish these two memory locations**.
At P2, however, doubling externally clamped A produces a maximum difference of 31.92% of
the total Delta-production response (Fig. 6D). The difference tends to zero with step size.
This suggests a nonlinear-response control or an independent proneural measurement, rather
than identifying proneural turnover from a growth rate alone.

Additional controls constrain the scope of these conclusions. Changing the periodic side
from 24 to 23 or 25 shifts stationary onset by at most 0.253%; full-model emergence stays
later in all 18 matched pairs (Note 7). Complete patterned-state spectra and perturbation
recovery establish local stability of twelve selected equilibria (Note 9). Excluding small
growth margins reduces, but does not remove, the growth discrepancy: the median ratio is
3.45 for 4048 points with relative growth margin at least 0.05 (Note 10). Prior reweighting
changes both instability fractions and effect sizes; they remain measure-specific.

# 4. Discussion

**Relation to existing approaches.** Table 2 summarises what each level of description predicts.
Our zero-frequency determinant condition has the same input–output structure used in
graph-based lateral-inhibition analyses [11, 12]. The identity
$\det J_\mu=\det J_\mathrm{loc}\det(I-\mu K)$ applies when $J_\mathrm{loc}$ is nonsingular,
with $K$ the local static input–output gain (Methods). This algebraic equivalence does not
transfer the full theorems of [11, 12] to this model. Their monotonicity and cell-stability
assumptions provide additional conclusions, including patterned-state existence and stability;
these assumptions have not been established here. Those studies use a degree-normalised
contact matrix; on our regular graph its eigenvalue is $\mu/6$, with the factor six absorbed
into the corresponding input gain. Here the unconditional result is the determinant identity,
and Section 3.1 proves the ordering of stationary thresholds for this model.

A matched reduced model reproduces the stationary boundary because it has the same static
responses. It also predicts growth rates and non-stationary instabilities, but these can differ
from the full model. Static gain alone predicts neither quantity. In the sampled region, collapsing the proneural layer
always overestimates the leading linear growth rate and sometimes predicts linear instability
where the full homogeneous state is linearly stable. Independent sampling strengthens this directional conclusion, while the
range checks show why its numerical magnitude must be tied to a stated parameter measure.

**Comparison with prior work.** Classical lateral-inhibition models establish symmetry
breaking [1]; cis/trans models establish how receptor–ligand interactions alter signalling
and patterning [2, 7]. Graph analyses provide general results under their stated cellular
assumptions [11, 12], and tissue-level models already address developmental timing [6].
The present comparison asks a narrower question: which errors remain after the full and
reduced models have been made exactly indistinguishable at equilibrium? Table 3 separates
this question from those earlier results. We have not shown that only a four-variable model
can reproduce a measured growth rate: a separately calibrated reduced model could fit one
dynamic observable. The single-clock control shows why matching that observable need not
transfer across conditions. The production-lag alternative further separates evidence for
finite feedback memory from evidence for its biochemical location. Our control holds the
static responses fixed and identifies what the stated QSS elimination discards.

**Testable predictions and required controls.** The compensated family
$(\lambda_S,d_S)\mapsto\alpha(\lambda_S,d_S)$ provides a specific intervention within the
model. It preserves the full steady-state curve $S_q(A)$, all patterned and homogeneous
equilibria and the stationary threshold. It changes both the relaxation of $S$ and the
growth of spatial perturbations. At P2, $\alpha=0.25,1,4$ give K-mode rates
$0.0220,0.1031,0.3831$, respectively, while $\lambda_N^c=1.304815343$ throughout.
These dimensionless operating-point predictions distinguish turnover from abundance.
Changing $d_S$ alone also changes $\lambda_S/d_S$ and does not implement this control.

A prospective experiment would first measure the proneural steady-state response to
activated Notch, then compare conditions with matched $S_q(A)$ and different measured
relaxation times. Surface receptor/ligand response and contact geometry would also need to
remain comparable. The prediction is a change in early spatial amplification and pattern
emergence time with unchanged stationary onset. If the measured static responses match but
the stationary onset shifts, the present model or the intervention assumptions are
insufficient. If growth is unchanged despite altered proneural relaxation, the claimed
rate-limiting role does not hold at that operating point. Long-time field equality is not
predicted by this intervention, because attractor stability and basins need not be preserved.

A second test measures a transient that the DNA reduction removes altogether. With $A$
externally held at a small step from its baseline, $S$ relaxes with time constant $1/d_S$;
under a small sinusoidal input of angular frequency $\omega$, its response relative to
the static repression has magnitude $[1+(\omega/d_S)^2]^{-1/2}$ and additional lag
$\arctan(\omega/d_S)$. The compensated intervention changes $d_S$ to $\alpha d_S$ in
both expressions. The DNA reduction has an instantaneous $S_q(A)$ response and zero
additional lag. Such a test would require accounting for reporter maturation and independently
measuring activated Notch; a lag alone would not uniquely identify proneural turnover.

Steady-state measurements cannot identify $\lambda_S$ and $d_S$ separately along this
family. Linear D, N and A measurements also cannot distinguish the explicitly constructed
production-lag alternative. The nonlinear clamped-A response in Figure 6D or independent
proneural measurements discriminate these specified alternatives; time resolution alone
does not guarantee mechanistic identification.
All proposed tests are model-based predictions. No tissue-calibrated rate, physical time,
experimental validation or uniquely preferred biological mechanism is inferred here.

**What is and is not new.** The mode-by-mode linear analysis and the Schur/static-gain identity
are standard [1, 10–12, 15]. Stability to uniform perturbations combined with instability to a
patterned mode is a generic feature of lateral-inhibition models. So is the growth of a linear
instability into a finite-amplitude pattern. We report neither as a finding. The central
contribution is the matched comparison and compensated-turnover prediction. Its analytic
support consists of:

* uniqueness of the homogeneous steady state for all positive parameters;
* the proof that cis inactivation and two-channel trans coupling cannot change the stationarily
  selected mode;
* the separation of static predictions, which reductions share exactly, from dynamic ones,
  which they do not;
* a quantitative account of where a static test fails and why.

**Limitations.** Parameters are literature-inspired operating points, not fits. The main
periodic graph removes edge effects and nucleation sites. Supplementary open-edge controls
retain raw neighbour sums: their degree variation creates deterministic forcing, so their
nonuniform contrast is not a homogeneous-instability criterion. Quenched production
heterogeneity is tested separately and does not represent dynamical noise. Property (iii) covers stationary onset
only; non-stationary selection was examined numerically, not proved. The model describes
continuous signalling states. "Sender-like" and "receiver-like" refer to the two Delta levels
of the final pattern, not to irreversible fate commitment, which would need an additional
switch. The identical final patterns in Fig. 5 are not a general equivalence of nonlinear
dynamics or patterned-state stability under QSS reduction. The absence of dynamic
counterexamples in sampled parameter boxes does not prove their global absence.

# 5. Methods

**Mode-resolved Jacobian.** At a homogeneous steady state $x^*$,

$$
J_\mu=\begin{pmatrix}
-d_D-(6f_D+k_c)N^* & -D^*(k_c+f_D\mu) & 0 & \lambda_D h_D'(S^*)\\
-N^*(k_c+f_N\mu) & -d_N-(6f_N+k_c)D^* & 0 & 0\\
f_N N^*\mu & 6f_N D^* & -d_A & 0\\
0 & 0 & -\lambda_S|h_S'(A^*)| & -d_S
\end{pmatrix}. \qquad (6)
$$

The Hill functions are evaluated in log space without a small-concentration floor or
large-concentration cap, so the drift and its analytic derivatives describe the same equations.
This block was checked against a finite-difference Jacobian of the full lattice. The
characteristic polynomial from (5) matches $\det(sI-J_\mu)$ within $2\times10^{-13}$ in its coefficients.

**Stage gains.** With $T_\mu=[A_\mu^{-1}]_{A,D}=f_NN^*n(\mu)/[d_A\Delta(\mu)]$,

$$
G(\mu)=\kappa_D\,T_\mu\,\kappa_S. \qquad (7)
$$

Logarithmic sensitivities of $G(-3)$ and of $\lambda_\mathrm{pattern}$ were computed by central
differences in $\ln p$, with steps $10^{-4}$ and $10^{-2}$ respectively.

**Graph input–output form.** Writing $J_\mu=J_\mathrm{loc}+\mu BC$, where $B$ carries the
neighbour sums $(s_D,s_N)$ into a cell and $C$ reads out $(D,N)$, Sylvester's identity gives
$\det J_\mu=\det J_\mathrm{loc}\det(I-\mu K)$ with $K=-CJ_\mathrm{loc}^{-1}B$, provided
$J_\mathrm{loc}$ is nonsingular. We verified this to
$2\times10^{-15}$ relative error at the reference points.

**Critical point.** The value of $\lambda_N$ at which the leading eigenvalue of $J_{-3}$ crosses
zero, and the value at which $G(-3)=1$, were each found by Brent's method on $[1.2,1.4]$ with
tolerance $10^{-13}$.

**Global sample.** 4000 points were drawn by Latin-hypercube sampling [16] (SciPy
`qmc.LatinHypercube`, seed 20260817) [17] over the ranges given in Table 1. For each point we
computed:

* the unique steady state (Section 2.3);
* the leading eigenvalue of $J_\mu$ for all 54 distinct $\mu$;
* $G(-3)$;
* the number of right-half-plane eigenvalues of $J_{-3}$, both directly and from the winding
  number of $1-L_{-3}(i\omega)$ on a logarithmic frequency grid.

The additional search of Section 3.3 used independent draws, rather than Latin-hypercube
sampling: a single NumPy PCG64 generator with seed 20261016 drew 20000 points from the original
box followed by 20000 with rate/constant bounds $[0.01,100]\times p_0$, without resetting the generator. Parameter order
was fixed as in Table 1. Both studies use log-uniform rates/constants and uniform Hill
coefficients. Actual sampled parameter vectors are deposited, including the additional search; sample IDs
are zero-based within each study/box.
The reference numerical environment was CPython 3.13.5, NumPy 2.2.6 and SciPy 1.16.3 on
Windows x86_64; remaining dependency versions and run provenance are recorded in the repository.

**Sampling robustness.** The protocol was recorded before new results were examined. The
primary study uses seeds 2026101701–2026101710 with 4000 points each. Sample-size checks
regenerate non-nested designs at 1000, 2000 and 8000 points using the first five seeds;
the wider-box study uses seeds 2026101801–2026101805 with 4000 points each. All use the
same scrambled LHS specification and parameter ordering as the reference sample. The
approximate 95% Monte Carlo interval for the mean fraction is
$\bar p\pm t_{R-1,0.975}s_p/\sqrt R$, where $R$ counts independent LHS replicates and $s_p$
is their sample standard deviation. Individual LHS points are not treated as independent
Bernoulli trials. Sizes reuse seeds and are not pooled for inference; changing bounds changes
the parameter measure. All inputs, per-point results and per-replicate summaries are deposited.

Every point is scanned over all 54 distinct graph eigenvalues for the full and both reduced
models. Growth thresholds of zero, $10^{-10}$, $10^{-8}$ and $10^{-6}$ are reported. Nyquist
counts were checked at 842 selected points: every static-gain miss and the 16 closest-to-unit-gain
points per replicate, with overlaps removed. Eleven selected cases, including the near-marginal
point, were independently recomputed at 80-digit precision, with agreement in stability signs
and K-mode RHP counts. Full protocol details and all seed-level results appear in Supplementary
Note 2 and Table S1.

**Reductions.** The DNA reduction replaces $S$ by $(\lambda_S/d_S)h_S(A)$ in the Delta
production term. The DN reduction also replaces $A$ by $f_N N s_D/d_A$. Their Jacobians are the
Schur complements of $J_\mu$ with respect to the eliminated variables.

**Nonlinear simulations.** The $24\times24$ lattice ($2304$ equations for the full model) was
integrated with the adaptive DOP853 method [17] at relative and absolute tolerances $10^{-9}$
and $10^{-11}$. The initial state was $x_i(0)=x^*+\varepsilon\,v\cos\phi_i$, with:

* $v$ the real eigenvector of the leading eigenvalue of $J_{-3}$ (or of its DNA reduction);
* $\phi_i=2\pi(8r_i+16c_i)/24$;
* $\varepsilon$ chosen so that the Delta component has relative amplitude $10^{-6}$.

The whole vector $v$ is scaled by one number. Scaling each component by its steady-state value,
as in the earlier version, does not give an eigenvector, and it adds a transient with an early
slope of 0.44 instead of 0.103 at P2 (Fig. S1). The read-out is the amplitude of the
$(m,n)=(8,16)$ Fourier component of $D-\bar D$. Its logarithm was fitted from $t=0$ to the first
recorded time at or above $10^{-3}D^*$, a rule fixed before the runs. Tightening integration
tolerances to $(10^{-11},10^{-14})$ and halving the recording interval to 0.125 retains
relative slope errors below $1.1\times10^{-4}$ at full-model P2/P3 and reduced-model P2.
The global instability count is unchanged for absolute growth thresholds from zero to
$10^{-6}$. Nyquist counts were additionally checked at 2000 and 8000 positive-frequency
points for the 12 gain mismatches and the 16 samples closest to unit gain.

**Random initial-condition experiments.** At P2 and P3, seeds 2026101901–2026101908
generate four independent $24\times24$ fields with NumPy PCG64 and uniform draws on
$[-1,1]$. Each field is centred and divided by its maximum absolute value, so
$x_i=x^*[1+\epsilon\xi_i]$ remains positive for the tested amplitudes. D, N and A are
shared between the full and DNA models, with full-model $S=S_q(A)$ in the main comparison;
independent S perturbations are a separate control. The same seed reuses the random fields
across amplitudes and points. There is no dynamical noise. The 112 production trajectories
use sparse BDF [17] with tolerances $10^{-9}/10^{-11}$ and interval 0.5, at fixed horizons
600 (P2) and 300 (P3). A solver change after an initial DOP853 pilot is documented in the
deposited protocol history; endpoints and design were retained.

Continuous upward events detect CV$(D)=0.1$ and 0.5; missed events would be censored at the
horizon. These thresholds differ from the single-mode amplitude used in Figure 5.
Early evolution is compared with $\widehat{\delta x}_{mn}(t)=\exp[J_{\mu(m,n)}t]
\widehat{\delta x}_{mn}(0)$ (or its DNA block), retaining all spatial modes and biochemical
poles, until any normalised concentration deviation exceeds $10^{-3}$. Endpoint metrics are
Delta CV, power at the two exact K points divided by all nonuniform Fourier power, scaled
drift residual and last-ten-time-unit Delta change. Pairwise field correlation is maximised
over periodic translations only. Ranges describe the seeds; no population confidence interval
is attached to this two-point deterministic study. Four first-seed cases are checked with
tighter BDF tolerances $10^{-11}/10^{-14}$ and interval 0.25, and independently with
DOP853 through five time units after the CV=0.5 event. Onset differences are below 0.006
dimensionless time units; only the BDF checks repeat the full endpoint horizon. Full protocol,
fields and all results are in Supplementary Note 3 and the deposited data.

**Time-discretisation audit.** The output recording interval is distinct from the internal
step chosen by the adaptive solver. Supplementary Note 4 tests these separately. Genuine
fixed-step RK4 trajectories use steps 0.08, 0.04, 0.02 and 0.01 in the early-growth regime
at P2/P3 for full and DNA models, with a common fitting grid and a high-accuracy Radau
reference. Separately, fourteen first-seed random initial states span both operating points,
both models, all three amplitudes and the independent-S controls. Each is reintegrated to
its original fixed horizon with BDF step caps 0.5, 0.25 and 0.125, tighter BDF tolerances,
and independent Radau integration. Continuous event times and raw cellwise endpoint fields
are compared with the original results. A separate resampling of the identical dense growth
solution isolates the effect of output cadence. Every distinct graph mode is also checked
against the fixed-step linear propagator at P0/P2/P3. These checks address temporal numerical
accuracy on the tested graph; they do not establish a spatial continuum limit or biological
validity. The design and numerical assessment thresholds are deposited with all results.

The stationary threshold and sampled instability fractions are obtained from the
continuous-time Jacobian and do not involve a temporal integration step. The fixed-step
controls nevertheless expose a numerical stability failure: at stable P0, RK4 with
$h=0.04$ or 0.08 incorrectly destabilises all 54 graph-mode blocks in both full and DNA
models. The fastest baseline decay rate has magnitude about 113, placing these steps
outside RK4's stability region. These failures account for 216 of 1296 spectral comparisons;
all 648 comparisons at $h=0.02$ and 0.01 agree. Those coarse fixed steps were not used in
the production adaptive calculations. Full step-level results are retained in Supplementary
Note 4; numerical verification is not presented as a separate scientific contribution.

**Additional matched controls.** Supplementary Notes 5–11 give the frozen designs, full
results and commands for QSS accuracy, clock calibration, geometry, constant-gain paths,
patterned-state stability, prior sensitivity and alternative memory location. All twelve
patterned-state Jacobians are analysed without homogeneous symmetry assumptions, using
complete dense spectra and independent sparse leading-eigenpair residuals. Paired perturbed
and unperturbed continuations isolate recovery from baseline drift. These checks establish
local stability only of the selected states. The alternative production-lag similarity and
the concentration rescaling for separate Hill thresholds are derived in Note 11.

**Code and data.** All code, data and figure scripts are at
https://github.com/HyeonjeYang/notch-delta-feedback-patterning, together with an automated test suite. These numerically check the analytic identities and specified limiting cases;
they supplement the proofs and do not replace them. `docs/MANUSCRIPT_NUMBERS.md` gives the script
and file behind every number. The version
that accompanied the first submission is preserved in the repository's `archive/` folder.

# Supplementary material

**Supplementary Note 1** (proofs of uniqueness, of properties (i)–(iii), of the reduction
identities and of the slow-proneural limit) corresponds to `docs/THEORY.md` in the repository
and is to be typeset as a supplement. **Supplementary Note 2**, sampling robustness and
Table S1, is provided in `manuscript/robustness_supplement.md`; Fig. S4 accompanies it.
**Supplementary Note 3**, paired initial-condition robustness, Tables S2–S3 and Fig. S5,
is provided in `manuscript/initial_condition_supplement.md`.
**Supplementary Note 4**, time-discretisation and solver checks, Tables S4–S6 and Fig. S6,
is provided in `manuscript/timestep_supplement.md`.

**Supplementary Notes 5–11**, Tables S7–S13 and Figures S7–S13, address QSS validity,
recalibration, graph geometry, oscillatory routes, patterned-state stability, growth/prior
sensitivity and mechanism identification, respectively. Their individual sources and
reproduction commands are indexed in `docs/SUPPLEMENTARY_INDEX.md`.

# Acknowledgements

I sincerely thank Professor Philip Nelson for helpful discussions and thoughtful feedback on this
work.

# References

1. J. R. Collier, N. A. M. Monk, P. K. Maini, J. H. Lewis. Pattern formation by lateral inhibition with feedback: a mathematical model of Delta–Notch intercellular signalling. *J. Theor. Biol.* 183, 429–446 (1996). doi:10.1006/jtbi.1996.0233
2. D. Sprinzak et al. Cis-interactions between Notch and Delta generate mutually exclusive signalling states. *Nature* 465, 86–90 (2010). doi:10.1038/nature08959
3. C. V. Cabrera. Lateral inhibition and cell fate during neurogenesis in *Drosophila*: the interactions between *scute*, *Notch* and *Delta*. *Development* 109, 733–742 (1990). doi:10.1242/dev.109.3.733
4. J. Culí, J. Modolell. Proneural gene self-stimulation in neural precursors: an essential mechanism for sense organ development that is regulated by Notch signaling. *Genes Dev.* 12, 2036–2047 (1998). doi:10.1101/gad.12.13.2036
5. B. Jennings, J. de Celis, C. Delidakis, A. Preiss, S. Bray. Role of Notch and achaete-scute complex in the expression of Enhancer of split bHLH proteins. *Development* 121, 3745–3752 (1995). doi:10.1242/dev.121.11.3745
6. F. Corson, L. Couturier, H. Rouault, K. Mazouni, F. Schweisguth. Self-organized Notch dynamics generate stereotyped sensory organ patterns in *Drosophila*. *Science* 356, eaai7407 (2017). doi:10.1126/science.aai7407
7. D. Sprinzak, A. Lakhanpal, L. LeBon, J. Garcia-Ojalvo, M. B. Elowitz. Mutual inactivation of Notch receptors and ligands facilitates developmental patterning. *PLoS Comput. Biol.* 7, e1002069 (2011). doi:10.1371/journal.pcbi.1002069
8. P. Formosa-Jordan, M. Ibañes. Competition in Notch signaling with cis enriches cell fate decisions. *PLoS One* 9, e95744 (2014). doi:10.1371/journal.pone.0095744
9. U. Binshtok, D. Sprinzak. Modeling the Notch response. *Adv. Exp. Med. Biol.* 1066, 79–98 (2018). doi:10.1007/978-3-319-89512-3_5
10. H. J. Wearing, M. R. Owen, J. A. Sherratt. Mathematical modelling of juxtacrine patterning. *Bull. Math. Biol.* 62, 293–320 (2000). doi:10.1006/bulm.1999.0152
11. M. Arcak. Pattern formation by lateral inhibition in large-scale networks of cells. *IEEE Trans. Autom. Control* 58, 1250–1262 (2013).
12. A. S. Rufino Ferreira, M. Arcak. A graph partitioning approach to predicting patterns in lateral inhibition systems. *SIAM J. Appl. Dyn. Syst.* 12, 2012–2031 (2013). doi:10.1137/130910142
13. E. Teomy, D. A. Kessler, H. Levine. Ordered hexagonal patterns via Notch–Delta signaling. *Phys. Biol.* 18 (6) (2021). doi:10.1088/1478-3975/ac28a4
14. F. Berkemeier, K. M. Page. Coupling dynamics of 2D Notch–Delta signalling. *Math. Biosci.* 360, 109012 (2023). doi:10.1016/j.mbs.2023.109012
15. K. J. Åström, R. M. Murray. *Feedback Systems: An Introduction for Scientists and Engineers*, 2nd ed. Princeton University Press (2021).
16. M. D. McKay, R. J. Beckman, W. J. Conover. A comparison of three methods for selecting values of input variables in the analysis of output from a computer code. *Technometrics* 21, 239–245 (1979). doi:10.1080/00401706.1979.10489755
17. P. Virtanen et al. SciPy 1.0: fundamental algorithms for scientific computing in Python. *Nat. Methods* 17, 261–272 (2020). doi:10.1038/s41592-019-0686-2

\clearpage

# Tables

**Table 1.** Parameters. Baseline $p_0$ and global-sample ranges.

| parameter | meaning | $p_0$ | sampled range |
|---|---|---|---|
| $\lambda_D$ | max. Delta production | 5.0 | 0.5–50 (log) |
| $\lambda_N$ | Notch production | 1.0 | 0.1–10 (log) |
| $\lambda_S$ | max. proneural production | 1.0 | 0.1–10 (log) |
| $K_d$ | half-saturation (both Hill functions) | 0.3 | 0.03–3 (log) |
| $n_D$ | Hill coefficient, $S\to D$ | 3 | 1–6 (linear) |
| $n_S$ | Hill coefficient, $A\dashv S$ | 3 | 1–6 (linear) |
| $d_D$, $d_N$ | Delta, Notch loss | 0.1, 0.1 | 0.01–1 (log) |
| $d_A$ | activated-Notch loss | 1.0 | 0.1–10 (log) |
| $d_S$ | proneural loss | 0.1 | 0.01–1 (log) |
| $f_D$ | trans depletion of Delta | 0.05 | 0.005–0.5 (log) |
| $f_N$ | trans activation of Notch | 0.6 | 0.06–6 (log) |
| $k_c$ | cis mutual inactivation | 0.35 | 0.035–3.5 (log) |

**Table 2.** Comparison with the full model. Exact agreement refers to the stated equations,
not to experimentally validated prediction accuracy. Growth-rate ratios refer to the sample.

| prediction | static graph criterion / $G(\mu)$ | matched reduction (DNA, DN) | explicit D–N–A–S |
|---|---|---|---|
| stationary onset ($\lambda_N^c$) | yes | yes (identical) | yes |
| mode selected at stationary onset | yes (property iii) | yes (identical) | yes |
| growth rate / time to pattern | not determined | differs (DNA: 5.4-fold pooled median in primary study) | computed from full dynamics |
| stability with $G(-3)<1$ | not determined (150 primary-study misses) | differs (DNA: 216 primary reduced-only points) | computed from full spectrum |
| set of patterned steady states | no | yes (identical set) | yes |

**Table 3.** Scope of the comparison with prior approaches. Existing dynamic models and
graph results retain their own predictions; this study tests the loss caused by a particular
matched QSS elimination.

| approach | established scope | question isolated here |
|---|---|---|
| feedback lateral inhibition [1] | spatial symmetry breaking from cell-to-cell inhibition | not claimed as new |
| receptor–ligand cis/trans models [2, 7] | signalling and patterning with explicit receptor and ligand interactions | retain these interactions while eliminating only the proneural relaxation |
| graph input–output analysis [11] | mode conditions under stable-cell and monotonicity assumptions | prove stationary-gain ordering for this model without assuming the complete graph theorem applies |
| equitable graph partitions [12] | candidate patterned equilibria and stability conditions | do not replace that theory or infer global basins from a linear mode |
| tissue-level Notch dynamics [6] | temporally evolving proneural patterns with experimental comparison | propose compensated-turnover and transient-response tests; no tissue fit is supplied |
| present matched control | same equilibrium set and zero-eigenvalue set, different dynamic spectra | quantify growth and stability errors, then test timing and endpoint structure with paired random starts |

\clearpage

# Figures

![**Figure 1. Circuit and contact geometry.** (A) The intracellular circuit of cell $i$ and its contacts with a neighbour $j$. Proneural activity $S$ drives Delta $D$. Neighbouring Delta activates Notch $N$, producing activated Notch $A$. $A$ represses $S$. Cis inactivation couples $D$ and $N$ within a cell. (B) Each cell contacts six neighbours (hexagonal tiling; polygons show contact topology, not measured cell shape). (C) The three-sublattice pattern spanned by the $\mu=-3$ eigenvectors. (D) Eigenvalue $\mu(m,n)$ of the contact adjacency matrix over all $24\times24$ plane-wave indices $(m,n)$ (Eq. 2). Circle: uniform mode, $\mu=6$. Squares: the two three-sublattice modes, $\mu=-3$.](../figures/main/Fig1_circuit_and_geometry.png){width=100%}

![**Figure 2. The static loop gain fixes stationary onset and the selected mode.** (A) Static loop gain $G(\mu)$ over every distinct graph eigenvalue at $\lambda_N=1$ ($p_0$), $\lambda_N^c$, 2 (P2) and 5.623 (P3). Markers show the maximum, which is at $\mu=-3$ in every case. Dashed line: $G=1$. (B) Leading eigenvalue of $J_\mu$ for $\mu=-3$ and $\mu=6$ as $\lambda_N$ varies, all other parameters at $p_0$. Dotted line: $\lambda_N^c$. (C) $G(-3)$ versus $\lambda_N$. It crosses one at $\lambda_N^c=1.304815343$, where the eigenvalue in (B) crosses zero.](../figures/main/Fig2_static_loop_gain.png){width=100%}

![**Figure 3. Kinetic control of the static margin.** (A) Logarithmic sensitivity of $G(-3)$ to each parameter at $p_0$ (bars), split into contributions through the three stage gains $\kappa_D$, $T_{-3}$ and $\kappa_S$ (markers). (B) Static sensitivity, $\partial\ln G(-3)/\partial\ln p$, against growth-rate sensitivity, $\partial\lambda_\mathrm{pattern}/\partial\ln p$. $\lambda_S$ and $d_S$ have equal and opposite static effects but unequal effects on the growth rate. (C, D) $\lambda_\mathrm{pattern}$ on the $(\lambda_N,d_A)$ and $(n_S,K_d)$ planes, other parameters at $p_0$. Black contour: $\lambda_\mathrm{pattern}=0$. Dashed white contour: $G(-3)=1$. The two coincide.](../figures/main/Fig3_kinetic_control.png){width=100%}

![**Figure 4. Instabilities that the static gain misses.** (A) All 4000 reference-sample points, plotted as $G(-3)$ against $\lambda_\mathrm{pattern}$ (values of $G<10^{-4}$, all stable, are drawn at $10^{-4}$). Points are coloured by the number of right-half-plane (RHP) eigenvalues of $J_{-3}$. Twelve points are unstable with $G<1$; each has exactly two RHP eigenvalues. (B) Nyquist curves $L_{-3}(i\omega)$ (bold: $\omega\ge0$; light: mirror image) for $p_0$, P2, and point 1730 ($G=0.366$). Dots: $G=L_{-3}(0)$. Each curve encircles $+1$ once per RHP eigenvalue. (C) RHP eigenvalues of the twelve points: six complex pairs and six real pairs.](../figures/main/Fig4_beyond_static_gain.png){width=100%}

![**Figure 5. Downstream dynamics preserve the stationary threshold and change growth.** (A) Growth rate of the $\mu=-3$ mode at P2 along the family $(\lambda_S,d_S)\to\alpha(\lambda_S,d_S)$, which leaves every steady state, every $G(\mu)$ and $\lambda_N^c$ unchanged. Dashed: slow-proneural limit $d_S[G(-3)-1]$. Horizontal lines: matched DNA and DN reductions. (B) Growth rate of the DNA reduction against the full model at all 1453 reference-sample points where both are unstable. Dashed line: equal rates. (C) Nonlinear simulations at P2. Amplitude of the three-sublattice component of $D$ (relative to $D^*$), for the full model and the DNA reduction, each started on its own exact eigenvector. Separate time axes make both fitting windows visible; time is dimensionless. Shaded spans: windows selected by the prespecified amplitude rule. Dashed lines: predicted growth. (D) Final Delta fields of the two simulations, shown in a hexagonal window of 397 of the 576 cells of the periodic lattice. The two fields coincide.](../figures/main/Fig5_explicit_dynamics.png){width=100%}

![**Figure 6. Dynamic information beyond a matched static response.** (A) Relative growth error of the DNA reduction in the compensated-turnover families at P2/P3; dashed line marks 10%. The required turnover factors are family-specific grid results, not universal QSS bounds. (B) A single DNA clock fitted at P2 underpredicts P3 growth by 16.5%. (C) Sample 679: real parts (blue) and absolute imaginary parts (orange) of K-mode poles along a constant-gain path. Green marks an oscillatory linear boundary; ochre marks the later unstable-pair collision on the real axis. Far-left stable poles are outside the display. G remains 0.564 throughout. No nonlinear Hopf or cycle-stability claim is made. (D) Normalised Delta-production response after doubling externally clamped A at P2. The explicit proneural model and matched production-lag alternative share every equilibrium and linear D,N,A response, but differ by up to 31.92% of the total production change at this amplitude. This is a specified mechanism-discrimination prediction, not tissue validation.](../figures/main/Fig6_dynamic_information.png){width=100%}

![**Figure S1. Initial-state comparison.** Growth of the three-sublattice amplitude at P2 (A) and P3 (B), started on the exact eigenvector or on the eigenvector scaled componentwise by the steady state (earlier protocol). Insets: $t\le15$, normalised to the initial amplitude. The componentwise start excites the other poles of $J_{-3}$. Its early slope (0.44 at P2, 0.56 at P3) differs from the predicted 0.103 and 0.216.](../figures/supplementary/FigS1_initial_state_comparison.png){width=100%}

![**Figure S2. One-parameter continuations.** $\lambda_\mathrm{pattern}$ as each parameter is varied about $p_0$ (rates $0.1$–$10\times$; Hill coefficients 1–6). Dotted: $p_0$ value. Dashed: $G(-3)=1$, which coincides with $\lambda_\mathrm{pattern}=0$ wherever a crossing occurs.](../figures/supplementary/FigS2_one_parameter_continuations.png){width=100%}

![**Figure S3. Dominance does not mean exclusive instability.** Leading eigenvalue of each mode block at P2, plotted against the graph eigenvalue $\mu$. Twenty-one distinct nonuniform graph eigenvalues have positive growth rates. The largest rate occurs at $\mu=-3$; other spatial modes are unstable but subdominant. Dashed line: zero growth.](../figures/supplementary/FigS3_subdominant_modes.png){width=85%}
