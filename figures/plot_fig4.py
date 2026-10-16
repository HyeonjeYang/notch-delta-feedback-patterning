"""Figure 4: where the static loop gain is not enough.

A  all 4000 global-sample points, G(-3) against lambda_pattern, coloured by
   how many right-half-plane eigenvalues the mu = -3 block has. The 12 points
   with G < 1 that are still unstable have exactly two such eigenvalues.
B  Nyquist curves L_-3(i omega) for a stable point (p0), a stationary
   instability (P2), and point 1730 (complex pair, G = 0.366). Each curve
   encircles +1 once per right-half-plane eigenvalue.
C  the right-half-plane eigenvalues of the 12 two-eigenvalue points.
Data: scripts/global_sample.py, scripts/nyquist_examples.py.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from figures.style import AQUA, BLUE, CONTEXT, DATA_DIR, INK, INK_2, MUTED, ORANGE, apply_style, panel_label, save

G_FLOOR = 1e-4
CLASSES = [("stable", CONTEXT, "o", 6, "stable (0 RHP)"),
           ("stationary", BLUE, "o", 7, "stationary (1 RHP, $G>1$)"),
           ("two_eig_complex", ORANGE, "s", 26, "complex pair (2 RHP, $G<1$)"),
           ("two_eig_real", AQUA, "^", 30, "real pair (2 RHP, $G<1$)")]


def panel_A(ax, r):
    assert (r.loc[r.G_minus3 < G_FLOOR, "lambda_pattern"] <= 0).all()
    for cls, colour, marker, size, name in CLASSES:
        sub = r[r.classification == cls]
        big = cls.startswith("two")
        ax.scatter(np.maximum(sub.G_minus3, G_FLOOR), sub.lambda_pattern, s=size, c=colour, marker=marker,
                   edgecolors="white" if big else "none", linewidths=0.6 if big else 0,
                   alpha=1.0 if big else 0.55, zorder=4 if big else 2, label=f"{name}: {len(sub)}",
                   rasterized=not big)
    ax.axvline(1.0, color=INK_2, lw=0.9, ls="--")
    ax.axhline(0.0, color=INK_2, lw=0.9, ls="--")
    ax.set_xscale("log"); ax.set_yscale("symlog", linthresh=0.02)
    ax.set_xlim(G_FLOOR * 0.7, 60)
    ax.set_xlabel(r"static loop gain $G(-3)$")
    ax.set_ylabel(r"$\lambda_{\rm pattern}$ (symlog)")
    ax.legend(loc="upper left", bbox_to_anchor=(-0.05, -0.25), fontsize=6.8, handletextpad=0.2, markerscale=1.1)
    ax.set_title("Global sample (n = 4000)", loc="left")


def panel_B(ax, ny):
    styles = {"p0 (stable)": (CONTEXT, "p0: 0 RHP"), "P2 (stationary)": (BLUE, "P2: 1 RHP"),
              "sample 1730 (complex pair)": (ORANGE, "#1730: 2 RHP")}
    for ex, (colour, name) in styles.items():
        g = ny[ny.example == ex]
        ax.plot(g.re_L, g.im_L, "-", color=colour, lw=1.5, label=name)
        ax.plot(g.re_L, -g.im_L, "-", color=colour, lw=0.8, alpha=0.55)
        ax.plot(g.re_L.iloc[0], 0, "o", color=colour, ms=4.5, mec="white", mew=0.6, zorder=5)
    ax.plot(1, 0, "x", color=INK, ms=7, mew=1.6, zorder=6)
    ax.annotate("+1", (1, 0), xytext=(1.04, 0.12), fontsize=8, color=INK)
    ax.axhline(0, color=MUTED, lw=0.6); ax.axvline(0, color=MUTED, lw=0.6)
    ax.set_xlim(-0.5, 2.0); ax.set_ylim(-1.35, 1.35)
    ax.set_xlabel(r"Re $L_{-3}(i\omega)$"); ax.set_ylabel(r"Im $L_{-3}(i\omega)$")
    ax.legend(loc="upper left", bbox_to_anchor=(-0.05, -0.25), fontsize=6.8)
    ax.set_title("Nyquist curves, $\\mu=-3$", loc="left")


def panel_C(ax, cases):
    for kind, colour, marker, name in [("complex", ORANGE, "s", "complex pairs"), ("real", AQUA, "^", "real pairs")]:
        sub = cases[cases.pair_type == kind]
        ax.plot(sub.re, sub.im, marker, color=colour, ms=5, mec="white", mew=0.6, ls="none",
                label=f"{name} ({sub.sample_id.nunique()} points)")
    c = cases[cases.sample_id == 1730]
    ax.annotate("#1730", (c.re.iloc[0], abs(c.im).max()), xytext=(-34, -3), textcoords="offset points",
                fontsize=7.5, color=INK)
    ax.axhline(0, color=MUTED, lw=0.6)
    ax.set_xscale("log")
    ax.set_xlabel(r"Re $s$ (right half-plane)"); ax.set_ylabel(r"Im $s$")
    ax.set_ylim(-0.32, 0.32)
    ax.legend(loc="upper left", bbox_to_anchor=(-0.05, -0.25), fontsize=6.8)
    ax.set_title(r"Unstable points with $G<1$", loc="left")


def main():
    apply_style()
    r = pd.read_csv(DATA_DIR / "global_sample_results.csv")
    ny = pd.read_csv(DATA_DIR / "nyquist_examples.csv")
    cases = pd.read_csv(DATA_DIR / "two_eigenvalue_cases.csv")
    fig, axes = plt.subplots(1, 3, figsize=(7.4, 2.6), gridspec_kw={"wspace": 0.45, "width_ratios": [1.15, 1, 1]})
    panel_A(axes[0], r); panel_B(axes[1], ny); panel_C(axes[2], cases)
    for ax, letter in zip(axes, "ABC"):
        panel_label(ax, letter, x=-0.26)
    save(fig, "Fig4_beyond_static_gain")


if __name__ == "__main__":
    main()
