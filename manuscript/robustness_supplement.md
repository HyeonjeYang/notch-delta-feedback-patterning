---
title: "Sampling robustness of the Notch Delta dynamics comparison"
author: Hyeonje Yang
date: October 2026
---

# Supplementary Note 2

The main conclusions persist across independent parameter seeds and sample sizes. Their
numerical magnitudes depend on the parameter measure. The original single-sample values
(36.325% unstable; DNA/full median growth ratio 4.90) are retained as reference-sample
results, rather than treated as universal circuit properties.

## Design fixed before computation

The protocol is deposited as `figures/data/robustness/protocol.json` and specifies 30
separately generated scrambled Latin-hypercube designs, totalling 115000 new points.
The primary analysis uses ten seeds, 2026101701 through 2026101710, each with 4000 points.
Rates and constants span $[0.1,10]\times p_0$ on a logarithmic scale; Hill coefficients
are uniform on $[1,6]$. Five further designs at each of 1000, 2000 and 8000 points use
the first five of those seeds. These designs are regenerated at each size, not nested
prefixes. Five wider-box designs, seeds 2026101801 through 2026101805, use 4000 points
and rate/constant bounds $[0.01,100]\times p_0$, with unchanged Hill ranges.

Every point is solved by scalar bracketing and scanned over all 54 distinct graph
eigenvalues for the full, DNA and DN models. Actual parameters, equilibria and results
are retained per replicate. No solver failure is labelled equilibrium nonexistence.
The batched calculation reproduces the original 4000 reference growth rates exactly
in this environment; its maximum static-gain difference is below $4\times10^{-15}$.

## Uncertainty and sample size

For independent replicate fractions $p_1,\ldots,p_R$, we report the approximate
95% Monte Carlo interval $\bar p\pm t_{R-1,0.975}s_p/\sqrt R$. The replication unit is
the independently randomised LHS, not an individual point within it. These intervals
describe numerical integration uncertainty under the specified parameter measure,
not experimental or biological variability. Sample-size comparisons reuse seeds;
they are not pooled into one inferential sample. The differing replication counts
also preclude reading the interval widths as an empirical convergence rate.

The primary estimate is **37.37%** (approximate 95% Monte Carlo
interval **37.07–37.67%**). Individual fractions range from
36.650% to 37.925%. This is a more stable
estimate of the chosen box measure than the original 36.325% single-sample value.

| bounds | points per replicate | replicates | mean unstable (%) | approximate 95% MC interval (%) | DNA median ratio range across replicates |
|---|---:|---:|---:|---|---|
| original | 4000 | 10 | 37.37 | 37.07–37.67 | 5.25–5.61 |
| original | 1000 | 5 | 37.10 | 36.00–38.20 | 4.98–5.44 |
| original | 2000 | 5 | 37.69 | 36.90–38.48 | 5.28–5.41 |
| original | 8000 | 5 | 37.39 | 36.91–37.87 | 5.03–5.48 |
| wider | 4000 | 5 | 28.89 | 28.58–29.20 | 2.99–3.64 |

## Dynamics and classification

Across the 40000 primary points, 150 instabilities are missed by $G(-3)>1$.
The DNA reduction has 216 reduced-only instabilities; DN has 445.
Neither reduction misses an instability of the full model in these samples. Both grow
faster wherever the full model is unstable. The primary pooled DNA/full growth-ratio
median is 5.41; the ten replicate medians range from
5.25 to 5.61. The corresponding DN range is
202–240. Medians are descriptive; a confidence
interval for a mean replicate fraction is not an interval for these conditional medians.

The wider box changes the unstable fraction to 28.89% and the DNA
replicate-median range to 2.99–3.64; the DN range becomes
1493–2019. Thus the qualitative bias persists but its magnitude
is not invariant to bounds. No empirical absence of an exception is treated as a theorem.

Over all 115000 analysed points there are 41307 unstable samples and
392 static-gain misses. Every miss has two RHP eigenvalues within
the K-mode block. No faster competing nonuniform mode, full-only instability or reduced
growth rate at or below the full-model rate was found. These totals count evaluations
across separate designs; they are not used as an iid-binomial denominator.

## Numerical edge cases and verification

All primary and size-check classifications are unchanged at growth thresholds through
$10^{-6}$. One wider-box point (seed 2026101803, ID 2229) has full growth
$4.772659\times10^{-9}$ and is explicitly flagged as near marginal. Its sign remains
positive in an independent 80-digit calculation. Removing it with a $10^{-8}$ threshold
changes the wider-box pooled fraction by 0.005 percentage points. It does not affect
the primary estimate or qualitative conclusions.

Nyquist winding matches direct RHP counts at all 842 selected points: every static-gain
miss and the 16 closest-to-unit-gain points per replicate, with overlaps removed.
This is a targeted validation, not a claim to have run Nyquist at every new point.
The five fastest complex-pair misses, five fastest real-pair misses and the marginal
point were checked independently at 80-digit precision. All eleven agree in stability
sign and K-mode RHP count. Equilibrium production/loss backward errors are below
$8.2\times10^{-15}$ throughout the new study.

## Table S1

| design | n | seed | unstable (%) | DC misses | DNA-only | DN-only | median DNA/full | median DN/full |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| primary | 4000 | 2026101701 | 37.250 | 10 | 23 | 44 | 5.362 | 204.7 |
| primary | 4000 | 2026101702 | 37.300 | 16 | 18 | 32 | 5.542 | 225.1 |
| primary | 4000 | 2026101703 | 37.925 | 15 | 19 | 41 | 5.529 | 201.9 |
| primary | 4000 | 2026101704 | 37.525 | 11 | 17 | 37 | 5.305 | 235.4 |
| primary | 4000 | 2026101705 | 37.500 | 18 | 22 | 47 | 5.372 | 226.2 |
| primary | 4000 | 2026101706 | 36.650 | 14 | 29 | 50 | 5.316 | 217.0 |
| primary | 4000 | 2026101707 | 37.750 | 23 | 19 | 47 | 5.609 | 240.3 |
| primary | 4000 | 2026101708 | 37.375 | 16 | 24 | 46 | 5.573 | 233.4 |
| primary | 4000 | 2026101709 | 36.725 | 13 | 24 | 44 | 5.250 | 239.9 |
| primary | 4000 | 2026101710 | 37.725 | 14 | 21 | 57 | 5.347 | 240.0 |
| size | 1000 | 2026101701 | 38.000 | 6 | 4 | 8 | 5.441 | 233.8 |
| size | 1000 | 2026101702 | 37.600 | 3 | 10 | 13 | 4.976 | 229.6 |
| size | 1000 | 2026101703 | 36.300 | 2 | 5 | 11 | 5.200 | 232.9 |
| size | 1000 | 2026101704 | 37.600 | 4 | 2 | 8 | 5.081 | 225.3 |
| size | 1000 | 2026101705 | 36.000 | 3 | 7 | 13 | 5.330 | 183.1 |
| size | 2000 | 2026101701 | 36.850 | 7 | 10 | 22 | 5.366 | 240.6 |
| size | 2000 | 2026101702 | 37.550 | 5 | 14 | 30 | 5.395 | 224.7 |
| size | 2000 | 2026101703 | 38.250 | 5 | 6 | 14 | 5.410 | 215.8 |
| size | 2000 | 2026101704 | 37.400 | 9 | 14 | 25 | 5.279 | 194.7 |
| size | 2000 | 2026101705 | 38.400 | 7 | 10 | 24 | 5.393 | 190.9 |
| size | 8000 | 2026101701 | 37.388 | 30 | 46 | 97 | 5.397 | 218.4 |
| size | 8000 | 2026101702 | 37.375 | 29 | 48 | 86 | 5.032 | 188.4 |
| size | 8000 | 2026101703 | 38.025 | 30 | 44 | 92 | 5.105 | 229.9 |
| size | 8000 | 2026101704 | 37.062 | 24 | 47 | 101 | 5.475 | 268.2 |
| size | 8000 | 2026101705 | 37.100 | 32 | 61 | 101 | 5.283 | 225.6 |
| wide | 4000 | 2026101801 | 29.100 | 8 | 6 | 15 | 3.322 | 1597.8 |
| wide | 4000 | 2026101802 | 28.900 | 10 | 7 | 23 | 2.994 | 1493.1 |
| wide | 4000 | 2026101803 | 28.825 | 13 | 13 | 29 | 3.149 | 2018.7 |
| wide | 4000 | 2026101804 | 29.125 | 7 | 10 | 25 | 3.641 | 1596.7 |
| wide | 4000 | 2026101805 | 28.500 | 8 | 13 | 32 | 3.205 | 1845.2 |

![**Figure S4. Robustness to sampling seed, sample size and parameter bounds.** (A) Nonuniform-instability fractions in the original and wider boxes. Each coloured point is an independent 4000-point replicate; dark diamonds and bars show means and approximate 95% Monte Carlo intervals. Dashed line: original reference sample. (B) Separately generated sample-size checks in the original box. R is the number of independent replicates at that size; seeds are reused across sizes. (C) Per-replicate median reduced/full growth ratios among points unstable in both models. (D) Classification differences, with each coloured mark representing one replicate. Full-only means the explicit model is unstable but the named reduction is stable.](../figures/supplementary/FigS4_sampling_robustness.png){width=100%}
