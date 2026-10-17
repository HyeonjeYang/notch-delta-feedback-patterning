"""Figure 5: what the explicit downstream layer adds.

A  proneural-timescale family at P2. The growth rate of the mu = -3 mode depends on alpha,
   while the equilibrium, G(-3) and lam_N^c do not. The DNA and DN
   DNA is the alpha -> infinity limit; DN additionally eliminates A.
B  global sample: growth rate of the matched DNA reduction against the
   explicit model, for every point unstable in both.
C  nonlinear simulations at P2 on the 24 x 24 lattice: K-mode amplitude of D
   for the explicit model and the DNA reduction, each started on its own
   exact eigenvector. Shaded spans are the prespecified fitting windows and
   dashed lines are the predicted slopes.
D  final Delta fields of the two models, which are identical.
Data: scripts/timescale_family.py, scripts/global_sample.py, scripts/growth_validation.py.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from figures.hex_plot import draw_hex_field
from figures.style import (AQUA, BLUE, DATA_DIR, INK, INK_2, ORANGE, SEQUENTIAL, TIME_LABEL, apply_style, panel_label,
                           save)


def panel_A(ax, fam, lam_c):
    ax.plot(fam.alpha, fam.P2_lambda, "-", color=BLUE, label="explicit D-N-A-S")
    slow = fam[fam.alpha <= 1]
    ax.plot(slow.alpha, slow.P2_slow_S_estimate, "--", color=INK_2, lw=1.0, label=r"slow-S limit $d_S[G(-3)-1]$")
    ax.axhline(fam.P2_lambda_DNA.iloc[0], color=ORANGE, lw=1.4, label="DNA reduction (S slaved)")
    ax.axhline(fam.P2_lambda_DN.iloc[0], color=AQUA, lw=1.4, ls=(0, (1, 1)), label="DN reduction (A, S slaved)")
    k = int(np.argmin(np.abs(np.log(fam.alpha))))
    ax.plot(fam.alpha.iloc[k], fam.P2_lambda.iloc[k], "o", color=BLUE, ms=5, mec="white", mew=0.8, zorder=5)
    ax.annotate(r"$\alpha=1$", (1, fam.P2_lambda.iloc[k]), xytext=(5, -11), textcoords="offset points", fontsize=7.5)
    ax.text(0.03, 0.68, rf"$\lambda_N^c={lam_c:.9f}$" "\nfor every " r"$\alpha$", transform=ax.transAxes,
            ha="left", va="top", fontsize=7, color=INK)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_ylim(4e-4, 3e1)
    ax.set_xlabel(r"proneural turnover factor $\alpha$")
    ax.set_ylabel(r"growth rate of $\mu=-3$ at P2")
    ax.legend(loc="lower right", fontsize=6.6)
    ax.set_title("Fixed stationary threshold", loc="left")


def panel_B(ax, r):
    both = r[(r.lambda_pattern > 0) & (r.lambda_DNA > 0)]
    ratio = both.lambda_DNA / both.lambda_pattern
    ax.scatter(both.lambda_pattern, both.lambda_DNA, s=7, c=BLUE, alpha=0.5, edgecolors="none", rasterized=True)
    lo = 10 ** np.floor(np.log10(min(both.lambda_pattern.min(), both.lambda_DNA.min())))
    ax.plot([lo, 1e2], [lo, 1e2], "--", color=INK_2, lw=0.9)
    ax.text(2e-4, 4e-4, "equal rates", rotation=38, fontsize=7, color=INK_2)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(lo, 10); ax.set_ylim(lo, 2e2)
    ax.set_xlabel(r"explicit model $\lambda_{\rm pattern}$")
    ax.set_ylabel(r"DNA reduction $\lambda_{\rm pattern}$")
    ax.text(0.97, 0.05, f"n = {len(both)}\nmedian ratio {ratio.median():.1f}\n"
            f"reduction faster in {(ratio > 1).mean():.0%}", transform=ax.transAxes, ha="right", va="bottom",
            fontsize=7, color=INK)
    ax.set_title("Reference parameter sample", loc="left")


def panel_C(ax, ts, summ, D_star):
    ax.axis("off")
    ax.set_title("Growth validation at P2", loc="left")
    # Separate time axes let the 4.5-unit DNA fitting window remain visible.
    for i, (run, colour, name, end) in enumerate([
        ("full_P2", BLUE, "explicit model", 175), ("DNA_P2", ORANGE, "DNA reduction", 11)
    ]):
        sub = ax.inset_axes([0.0 + 0.56 * i, 0.08, 0.44, 0.78])
        g = ts[ts.run == run]
        s = summ[summ.run == run].iloc[0]
        a0 = g.a_D.iloc[0] / D_star
        sub.axvspan(s.window_start, s.window_end, color=colour, alpha=0.10, lw=0)
        sub.plot(g.time, g.a_D / D_star, "-", color=colour, zorder=3)
        t = np.linspace(0, s.window_end, 50)
        sub.plot(t, a0 * np.exp(s.predicted_lambda * t), "--", color=INK, lw=1, zorder=4)
        sub.text(0.96, 0.06, f"fit {s.fitted_slope:.4f}\npole {s.predicted_lambda:.4f}\n"
                 f"[0, {s.window_end:g}]", transform=sub.transAxes, ha="right", fontsize=7, color=INK)
        sub.set_yscale("log"); sub.set_ylim(5e-7, 1e3); sub.set_xlim(0, end)
        sub.set_yticks([1e-6, 1e-3, 1, 1e3])
        sub.set_title(name, fontsize=8, pad=6, color=colour)
        sub.set_xlabel(r"time $t$", fontsize=8)
        sub.tick_params(labelsize=7)
        if i == 0:
            sub.set_ylabel(r"K-mode amplitude / $D^*$", fontsize=8)
        else:
            sub.set_yticklabels([])


def panel_D(ax, fields, diff):
    """ax is a frame occupying the same grid cell as the other panels; the two
    fields are drawn as hexagonal windows inside it."""
    ax.axis("off")
    ax.set_title("Final pattern at P2 ($t=400$)", loc="left")
    vmin = min(f.min() for f in fields.values()); vmax = max(f.max() for f in fields.values())
    for i, (name, field) in enumerate(fields.items()):
        sub = ax.inset_axes([0.0 + 0.52 * i, 0.26, 0.48, 0.68])
        sm = draw_hex_field(sub, field, vmin, vmax, cmap=SEQUENTIAL)
        sub.set_title(name, fontsize=8)
    cax = ax.inset_axes([0.1, 0.12, 0.8, 0.045])
    cb = plt.colorbar(sm, cax=cax, orientation="horizontal")
    cb.set_label(r"Delta $D$", fontsize=7.5)
    cb.ax.tick_params(labelsize=7); cb.ax.grid(False)
    mant, expo = f"{diff:.0e}".split("e")
    ax.text(0.5, -0.12, rf"max $|D_{{\mathrm{{explicit}}}}-D_{{\mathrm{{DNA}}}}|={mant}\times10^{{{int(expo)}}}$",
            transform=ax.transAxes, fontsize=7, color=INK, ha="center", va="top")


def main():
    apply_style()
    fam = pd.read_csv(DATA_DIR / "timescale_family.csv")
    lam_c = float(fam.critical_lam_N.dropna().iloc[0])
    r = pd.read_csv(DATA_DIR / "global_sample_results.csv")
    ts = pd.read_csv(DATA_DIR / "growth_timeseries.csv")
    summ = pd.read_csv(DATA_DIR / "growth_validation.csv")
    finals = np.load(DATA_DIR / "pattern_final_fields.npz")
    D_star = float(fam.P2_D_star.iloc[0])

    fig = plt.figure(figsize=(8.8, 7.5))
    gs = fig.add_gridspec(2, 2, wspace=0.42, hspace=0.55)
    ax_a, ax_b = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
    ax_c, ax_d = fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])
    panel_A(ax_a, fam, lam_c); panel_B(ax_b, r); panel_C(ax_c, ts, summ, D_star)
    panel_D(ax_d, {"explicit model": finals["D_full_P2"], "DNA reduction": finals["D_DNA_P2"]},
            float(summ.final_D_maxabs_diff_full_vs_DNA_P2.iloc[0]))
    for ax, letter in zip([ax_a, ax_b, ax_c, ax_d], "ABCD"):
        panel_label(ax, letter, x=-0.24)
    save(fig, "Fig5_explicit_dynamics")


if __name__ == "__main__":
    main()
