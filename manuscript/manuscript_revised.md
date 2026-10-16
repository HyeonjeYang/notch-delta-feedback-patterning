---
title: "Explicit Notch–proneural dynamics set the tempo, not the onset, of lateral-inhibition patterning"
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

Lateral inhibition through Notch–Delta signalling turns small differences between neighbouring
cells into alternating cell states. Most models compress the steps downstream of Notch
activation into one or two effective variables. We ask which predictions change when the
activated-Notch and proneural steps are kept as separate dynamical variables. Our model has four
variables per cell: Delta, Notch receptor, activated Notch, and proneural activity. Cells sit on a
periodic six-neighbour contact lattice. We prove that the homogeneous steady state is unique for
all positive parameters. We also prove that the static (zero-frequency) gain of the feedback loop
is always largest for the three-sublattice spatial mode, so this mode is the only one that can
lose stability first through a zero eigenvalue. These onset properties depend only on static
input–output responses. Collapsing the proneural layer, or both downstream layers, by
quasi-steady-state elimination therefore leaves the stationary instability boundary exactly
unchanged. What the explicit layer changes is the dynamics. Rescaling the proneural turnover
rate leaves the critical Notch production rate fixed while the pattern growth rate varies over
three decades. Across a 4000-point Latin-hypercube parameter sample, the reduced model
overestimates the growth rate at every unstable point (median 4.9-fold) and predicts
instabilities absent from the full model. Twelve sampled points are unstable even though the
static gain is below one. Each has two eigenvalues in the right half-plane, which a Nyquist test
detects and a static test cannot. Patterning onset is thus a static property of the circuit,
while its tempo is set by the slowest downstream step.

# 1. Introduction

In lateral inhibition, Delta on one cell activates Notch on its neighbours, and Notch signalling
represses the proneural activity that drives Delta expression. Small differences between
adjacent cells are amplified into a spatial pattern of signal-sending and signal-receiving
states [1, 2]. In *Drosophila* sensory-organ development this antagonism acts through
Achaete–Scute proneural self-stimulation and Notch/Enhancer-of-split repression [3–6].

Mathematical models of this process usually compress everything downstream of
receptor–ligand binding into one or two effective variables, such as a "Notch activity" that
directly represses Delta [1, 7–9]. For such models, the conditions under which the homogeneous
state loses stability are well understood. The coupled system is decomposed into eigenmodes of
the cell-contact graph, and each mode is tested separately [1, 10, 13]. Equivalently, instability
of a mode follows from a static (zero-frequency) input–output gain of a single cell multiplied by
the corresponding graph eigenvalue [11, 12]. On the triangular contact graph of a hexagonal
tissue, the most negative graph eigenvalue belongs to a three-sublattice pattern, not a
checkerboard [13]. Two-dimensional contact coupling has also been analysed directly [14].

These results leave a question open. When the activated-Notch and proneural steps are kept as
separate dynamical variables, which predictions actually differ from those of a reduced model?
Retaining more variables gives more parameters to vary. It does not by itself show that the
extra dynamics matter for patterning.

We answer this question for a four-variable Delta–Notch–activated-Notch–proneural model. We do
not claim the stability machinery as new. The factorisation of each mode's characteristic
polynomial into a static loop gain (Section 3.1) is a standard Schur-complement identity, and it
is equivalent to the graph-based static-gain conditions of [11, 12]. We use it to separate
predictions that depend only on static responses from predictions that depend on downstream
timescales. Four results follow:

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
733 points had been discarded in this way. All of them have a steady state, with residual below
$3\times10^{-15}$.

# 3. Results

## 3.1 The static loop gain fixes the stationary onset and the selected mode

To identify which part of the circuit controls stability, we split the four variables of each
mode into a "plant" $z=(D,N,A)$ and the proneural feedback state $S$:

$$
\dot z = A_\mu z + b\,S,\qquad \dot S = c^{\top}z - d_S S. \qquad (4)
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
the Nyquist criterion [15]. Since the open loop has no unstable poles, the number of unstable
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
$\lambda_N^c=1.304815343$ and agree to $1\times10^{-15}$, as identity (5) requires. The
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
two-parameter maps the contour $\lambda_\mathrm{pattern}=0$ coincides with $G(-3)=1$ at all 1875
grid points (Fig. 3C,D). This ranking describes the present operating point. It does not claim
that cis/trans interactions are unimportant in Notch biology [2, 7, 8].

The comparison in Fig. 3B anticipates Section 3.4. Static quantities depend on $\lambda_S$ and
$d_S$ only through their ratio $\lambda_S/d_S$, so their static sensitivities are equal and
opposite ($\mp1.80$). Their effects on the growth rate $\lambda_\mathrm{pattern}$ are not
($-0.060$ and $+0.044$), because $d_S$ also sets how fast the proneural stage relaxes.

**Global sample.** Over a 4000-point Latin-hypercube sample of all 13 parameters (each rate
within $[0.1,10]\times p_0$, log-uniform; Hill coefficients in $[1,6]$), **1453 points (36.3%)
are unstable to a spatial mode.** At every one of them the dominant mode is $\mu=-3$. The
earlier figure of 30.3% excluded the 733 points discussed in Section 2.3. The fraction depends
on the sampling measure and is not a biological probability.

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
and 20 000 from a box ten times wider in every rate. Of these, 13 047 were unstable to a spatial
mode. In none was a mode other than $\mu=-3$ dominant, or unstable while $\mu=-3$ was stable.
The only other mode we found unstable is the uniform mode, which at 45 of the 4000 sampled
points has a complex pair in the right half-plane, corresponding to synchronous oscillation.
All 45 points are also unstable to the three-sublattice mode.

## 3.4 What the explicit downstream layer adds

To isolate the role of the downstream dynamics, we compared the full model with two matched
reductions (Methods):

* **DNA:** proneural activity is slaved to activated Notch, $S=(\lambda_S/d_S)h_S(A)$. This
  gives a three-variable Notch–Delta–reporter model of the type studied in [2, 7].
* **DN:** $A$ is also slaved, giving two variables.

Both reductions only impose $\dot S=0$ (and $\dot A=0$) cell by cell, so they have exactly the
same steady states as (1), homogeneous or patterned. Their mode-resolved Jacobians are Schur
complements of $J_\mu$, and $\det J_\mu$ equals $-d_S\det J^\mathrm{DNA}_\mu$ and
$d_Ad_S\det J^\mathrm{DN}_\mu$. **The full model and both reductions therefore share the
stationary boundary, the critical $\lambda_N$, and the selected mode.** Section 3.1 and the
boundary in Section 3.2 do not need the explicit downstream dynamics. A reduced model with the
same static input–output response predicts them exactly.

The explicit dynamics matter for timing. Consider the one-parameter family
$(\lambda_S,d_S)\to\alpha(\lambda_S,d_S)$. It leaves $S=(\lambda_S/d_S)h_S(A)$, every steady
state and every $G(\mu)$ unchanged, and only rescales how fast proneural activity relaxes. As
$\alpha\to\infty$ the model approaches the DNA reduction. As $\alpha\to0$ the plant equilibrates
quickly relative to $S$, and the growth rate of mode $\mu$ approaches $d_S[G(\mu)-1]$
(Supplementary Note 1). The static gain sets whether a mode grows, and the slowest downstream
step sets how fast. Along this family at P2, $\lambda_N^c=1.304815343$ for every $\alpha$
from $10^{-2}$ to $10^{4}$, while the growth rate of the three-sublattice mode ranges from
$8\times10^{-4}$ to $1.59$ (Fig. 5A). At the baseline timescales ($\alpha=1$) the explicit model
lies close to the slow-proneural limit: it grows at $0.103$, against $0.081$ from
$d_S[G-1]$. The DNA and DN reductions predict $1.597$ (15.5-fold faster) and $11.28$
(109-fold).

The same holds across the global sample (Fig. 5B). The DNA reduction is unstable wherever the
full model is unstable, and grows faster at every one of those points (median 4.9-fold,
interquartile range 2.6–13.0). It is also unstable at 26 points where the full model is stable,
all with $G(-3)<1$. These are non-stationary instabilities that the finite proneural relaxation
removes. For the DN reduction the median overestimate is 216-fold, with 47 spurious
instabilities.

Nonlinear simulations confirm both statements (Fig. 5C,D). We started each model at P2 on the
exact $\mu=-3$ eigenvector of its own Jacobian, at relative Delta amplitude $10^{-6}$ (Methods).
The amplitude of the three-sublattice component then grows at the predicted rate from $t=0$:

| model | predicted | fitted | fitting window |
|---|---|---|---|
| full | 0.103067 | 0.103077 | $[0, 67.25]$ |
| DNA reduction | 1.597479 | 1.597504 | $[0, 4.5]$ |

The fitting window is fixed in advance and drawn in the figure. The pattern reaches 10% relative
amplitude at $t=110.8$ in the full model and at $t=7.3$ in the reduction. The two models
then settle onto the same three-sublattice steady state, with a maximum difference of
$7\times10^{-11}$ in the final Delta field. The explicit proneural layer changes when the pattern
appears, but not whether or which pattern appears.

# 4. Discussion

**Relation to existing approaches.** Table 2 summarises what each level of description predicts.
Graph-based static-gain criteria [11, 12] and our $G(\mu)$ test the same zero-frequency loop;
the identity $\det J_\mu=\det J_\mathrm{loc}\det(I-\mu K)$, with $K$ the static input–output gain of
an isolated cell, links them directly (Methods). They locate stationary onset and its mode. A
matched reduced model reproduces both exactly, because it has the same static responses. What
neither captures is (a) the growth rate, which depends on downstream timescales that the
static responses do not constrain, and (b) instabilities that do not pass through a zero
eigenvalue. These are the predictions that specifically require the explicit
activated-Notch and proneural dynamics. In the sampled region, collapsing the proneural layer
always overestimates the speed of patterning and sometimes predicts patterning where the full
model has none.

**Biological reading.** In this model the proneural stage is often the slowest step, and its
turnover sets the tempo of symmetry breaking. When it is slow, the growth rate is
approximately $d_S[G(-3)-1]$: the excess static gain times the proneural relaxation rate. A
perturbation that changes proneural protein stability while preserving the steady-state
dose–response of the circuit should therefore change the timing of fine-grained patterning
without shifting the threshold at which it occurs. A perturbation that changes a static gain,
such as the steepness or position of proneural repression ($n_S$, $K_d$), should shift the
threshold itself. These are qualitative, model-based expectations. The parameter values are
not fitted to a tissue.

**What is and is not new.** The mode-by-mode linear analysis and the Schur/static-gain identity
are standard [1, 10–12, 15]. Stability to uniform perturbations combined with instability to a
patterned mode is a generic feature of lateral-inhibition models. So is the growth of a linear
instability into a finite-amplitude pattern. We report neither as a finding. The contributions
here are:

* uniqueness of the steady state for all positive parameters;
* the proof that cis inactivation and two-channel trans coupling cannot change the stationarily
  selected mode;
* the separation of static predictions, which reductions share exactly, from dynamic ones,
  which they do not;
* a quantitative account of where a static test fails and why.

**Limitations.** Parameters are literature-inspired operating points, not fits. Periodic
boundaries remove edge effects and nucleation sites. Property (iii) covers stationary onset
only; non-stationary selection was examined numerically, not proved. The model describes
continuous signalling states. "Sender-like" and "receiver-like" refer to the two Delta levels
of the final pattern, not to irreversible fate commitment, which would need an additional
switch.

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

This block was checked against a finite-difference Jacobian of the full lattice. The
characteristic polynomial from (5) matches $\det(sI-J_\mu)$ to $10^{-13}$ in its coefficients.

**Stage gains.** With $T_\mu=[A_\mu^{-1}]_{A,D}=f_NN^*n(\mu)/[d_A\Delta(\mu)]$,

$$
G(\mu)=\kappa_D\,T_\mu\,\kappa_S. \qquad (7)
$$

Logarithmic sensitivities of $G(-3)$ and of $\lambda_\mathrm{pattern}$ were computed by central
differences in $\ln p$, with steps $10^{-4}$ and $10^{-2}$ respectively.

**Graph input–output form.** Writing $J_\mu=J_\mathrm{loc}+\mu BC$, where $B$ carries the
neighbour sums $(s_D,s_N)$ into a cell and $C$ reads out $(D,N)$, Sylvester's identity gives
$\det J_\mu=\det J_\mathrm{loc}\det(I-\mu K)$ with $K=-CJ_\mathrm{loc}^{-1}B$. We verified this to
$10^{-15}$ relative error.

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

The random search of Section 3.3 used the same procedure with a fixed seed.

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
time the amplitude reached $10^{-3}D^*$, a rule fixed before the runs.

**Code and data.** All code, data and figure scripts are at
https://github.com/HyeonjeYang/notch-delta-feedback-patterning, together with a test suite of
341 checks. The checks include every analytic statement above, evaluated on random parameter
sets. `docs/MANUSCRIPT_NUMBERS.md` gives the script and file behind every number. The version
that accompanied the first submission is preserved in the repository's `archive/` folder.

# Supplementary material

**Supplementary Note 1** (proofs of uniqueness, of properties (i)–(iii), of the reduction
identities and of the slow-proneural limit) corresponds to `docs/THEORY.md` in the repository
and is to be typeset as a supplement.

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

**Table 2.** What each level of description predicts (yes: predicted exactly; no: not predicted).

| prediction | static graph criterion / $G(\mu)$ | matched reduction (DNA, DN) | explicit D–N–A–S |
|---|---|---|---|
| stationary onset ($\lambda_N^c$) | yes | yes (identical) | yes |
| mode selected at stationary onset | yes (property iii) | yes (identical) | yes |
| growth rate / time to pattern | no | no (4.9-fold median overestimate) | yes |
| non-stationary instabilities | no (12 missed) | no (26 spurious, DNA) | yes |
| set of patterned steady states | no | yes (identical set) | yes |

\clearpage

# Figures

![**Figure 1. Circuit and contact geometry.** (A) The intracellular circuit of cell $i$ and its contacts with a neighbour $j$. Proneural activity $S$ drives Delta $D$. Neighbouring Delta activates Notch $N$, producing activated Notch $A$. $A$ represses $S$. Cis inactivation couples $D$ and $N$ within a cell. (B) Each cell contacts six neighbours (hexagonal tiling; polygons show contact topology, not measured cell shape). (C) The three-sublattice pattern spanned by the $\mu=-3$ eigenvectors. (D) Eigenvalue $\mu(m,n)$ of the contact adjacency matrix over all $24\times24$ plane-wave indices $(m,n)$ (Eq. 2). Circle: uniform mode, $\mu=6$. Squares: the two three-sublattice modes, $\mu=-3$.](../figures/main/Fig1_circuit_and_geometry.png){width=100%}

![**Figure 2. The static loop gain fixes stationary onset and the selected mode.** (A) Static loop gain $G(\mu)$ over every distinct graph eigenvalue at $\lambda_N=1$ ($p_0$), $\lambda_N^c$, 2 (P2) and 5.623 (P3). Markers show the maximum, which is at $\mu=-3$ in every case. Dashed line: $G=1$. (B) Leading eigenvalue of $J_\mu$ for $\mu=-3$ and $\mu=6$ as $\lambda_N$ varies, all other parameters at $p_0$. Dotted line: $\lambda_N^c$. (C) $G(-3)$ versus $\lambda_N$. It crosses one at $\lambda_N^c=1.304815343$, where the eigenvalue in (B) crosses zero.](../figures/main/Fig2_static_loop_gain.png){width=100%}

![**Figure 3. Kinetic control of the static margin.** (A) Logarithmic sensitivity of $G(-3)$ to each parameter at $p_0$ (bars), split into contributions through the three stage gains $\kappa_D$, $T_{-3}$ and $\kappa_S$ (markers). (B) Static sensitivity, $\partial\ln G(-3)/\partial\ln p$, against growth-rate sensitivity, $\partial\lambda_\mathrm{pattern}/\partial\ln p$. $\lambda_S$ and $d_S$ have equal and opposite static effects but unequal effects on the growth rate. (C, D) $\lambda_\mathrm{pattern}$ on the $(\lambda_N,d_A)$ and $(n_S,K_d)$ planes, other parameters at $p_0$. Black contour: $\lambda_\mathrm{pattern}=0$. Dashed white contour: $G(-3)=1$. The two coincide.](../figures/main/Fig3_kinetic_control.png){width=100%}

![**Figure 4. Instabilities that the static gain misses.** (A) All 4000 global-sample points, plotted as $G(-3)$ against $\lambda_\mathrm{pattern}$ (values of $G<10^{-4}$, all stable, are drawn at $10^{-4}$). Points are coloured by the number of right-half-plane (RHP) eigenvalues of $J_{-3}$. Twelve points are unstable with $G<1$; each has exactly two RHP eigenvalues. (B) Nyquist curves $L_{-3}(i\omega)$ (bold: $\omega\ge0$; light: mirror image) for $p_0$, P2, and point 1730 ($G=0.366$). Dots: $G=L_{-3}(0)$. Each curve encircles $+1$ once per RHP eigenvalue. (C) RHP eigenvalues of the twelve points: six complex pairs and six real pairs.](../figures/main/Fig4_beyond_static_gain.png){width=100%}

![**Figure 5. The explicit downstream layer sets the tempo, not the onset.** (A) Growth rate of the $\mu=-3$ mode at P2 along the family $(\lambda_S,d_S)\to\alpha(\lambda_S,d_S)$, which leaves every steady state, every $G(\mu)$ and $\lambda_N^c$ unchanged. Dashed: slow-proneural limit $d_S[G(-3)-1]$. Horizontal lines: matched DNA and DN reductions. (B) Growth rate of the DNA reduction against the full model at all 1453 global-sample points where both are unstable. Dashed line: equal rates. (C) Nonlinear simulations at P2. Amplitude of the three-sublattice component of $D$ (relative to $D^*$), for the full model and the DNA reduction, each started on its own exact eigenvector. Shaded spans: prespecified fitting windows. Dashed lines: predicted growth. (D) Final Delta fields of the two simulations, shown in a hexagonal window of 397 of the 576 cells of the periodic lattice. The two fields coincide.](../figures/main/Fig5_explicit_dynamics.png){width=100%}

![**Figure S1. Initial-state comparison.** Growth of the three-sublattice amplitude at P2 (A) and P3 (B), started on the exact eigenvector or on the eigenvector scaled componentwise by the steady state (earlier protocol). Insets: $t\le15$, normalised to the initial amplitude. The componentwise start excites the other poles of $J_{-3}$. Its early slope (0.44 at P2, 0.57 at P3) differs from the predicted 0.103 and 0.216.](../figures/supplementary/FigS1_initial_state_comparison.png){width=100%}

![**Figure S2. One-parameter continuations.** $\lambda_\mathrm{pattern}$ as each parameter is varied about $p_0$ (rates $0.1$–$10\times$; Hill coefficients 1–6). Dotted: $p_0$ value. Dashed: $G(-3)=1$, which coincides with $\lambda_\mathrm{pattern}=0$ wherever a crossing occurs.](../figures/supplementary/FigS2_one_parameter_continuations.png){width=100%}
