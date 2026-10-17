"""Figure 2: the static loop gain G(mu) sets the stationary onset and the
selected mode.

A  G(mu) over every distinct graph eigenvalue at four lam_N values. The
   maximum is always at mu = -3, and G(6) < 0.
B  leading pole of J_mu versus lam_N for mu = -3 and mu = 6.
C  G(-3) and G(6) versus lam_N. G(-3) = 1 exactly where the pole in B
   crosses zero.
Data: scripts/linear_analysis.py.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.ticker import FixedLocator, NullFormatter, ScalarFormatter

from figures.style import BLUE, DATA_DIR, INK_2, MUTED, ORANGE, SERIES, apply_style, panel_label, save

MARKERS = ["o", "s", "^", "D"]
LINESTYLES = ["-", "--", "-.", ":"]
LABELS = {"baseline": r"$\lambda_N=1$ (p0)", "critical": r"$\lambda_N=\lambda_N^c$",
          "P2": r"$\lambda_N=2$ (P2)", "P3": r"$\lambda_N=5.623$ (P3)"}


def panel_A(ax, spec):
    for (label, colour, marker, ls) in zip(LABELS, SERIES, MARKERS, LINESTYLES):
        sub = spec[spec.label == label].sort_values("mu")
        ax.plot(sub.mu, sub.G, ls=ls, color=colour, lw=1.7, label=LABELS[label])
        k = sub.G.idxmax()
        ax.plot(sub.mu[k], sub.G[k], marker, color=colour, ms=5, mec="white", mew=0.8, zorder=4)
    ax.axhline(1.0, color=INK_2, lw=0.9, ls="--")
    ax.text(5.9, 1.05, "$G=1$", ha="right", va="bottom", fontsize=7.5, color=INK_2)
    ax.axhline(0.0, color=MUTED, lw=0.6)
    ax.set_xlabel(r"graph eigenvalue $\mu$")
    ax.set_ylabel(r"static loop gain $G(\mu)=L_\mu(0)$")
    ax.set_xlim(-3.3, 6.3)
    ax.legend(loc="lower left")
    ax.set_title("Spatial gain spectrum", loc="left")


def panel_B(ax, poles, lam_c):
    for mu, colour, marker, name in [(-3.0, BLUE, "o", r"$\mu=-3$"), (6.0, ORANGE, "s", r"$\mu=6$ (uniform)")]:
        lead = poles[(poles.mu == mu) & (poles.pole_index == 0)].sort_values("lam_N")
        ax.plot(lead.lam_N, lead.re, ls="-" if mu == -3 else "--", color=colour, label=name)
    ax.axhline(0.0, color=INK_2, lw=0.9, ls="--")
    ax.axvline(lam_c, color=MUTED, lw=0.8, ls=":")
    ax.set_xscale("log")
    ax.set_xlabel(r"Notch production $\lambda_N$")
    ax.set_ylabel(r"leading Re$\,s$ of $J_\mu$")
    ax.legend(loc="upper left")
    ax.set_title("Pole crossing", loc="left")


def panel_C(ax, poles, lam_c):
    sub = poles[(poles.mu == -3.0) & (poles.pole_index == 0)].sort_values("lam_N")
    ax.plot(sub.lam_N, sub.G_mu, "-", color=BLUE)
    ax.axhline(1.0, color=INK_2, lw=0.9, ls="--")
    ax.axvline(lam_c, color=MUTED, lw=0.8, ls=":")
    ax.annotate(rf"$\lambda_N^c={lam_c:.4f}$", (lam_c, 1.0), xytext=(2.0, 0.08), fontsize=7.5,
                arrowprops=dict(arrowstyle="-", color=INK_2, lw=0.7))
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel(r"Notch production $\lambda_N$")
    ax.set_ylabel(r"$G(-3)$")
    ax.set_title("Stationary threshold", loc="left")


def log_ticks(ax):
    ax.xaxis.set_major_locator(FixedLocator([0.5, 1, 2, 5]))
    ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.xaxis.set_minor_formatter(NullFormatter())


def main():
    apply_style()
    spec = pd.read_csv(DATA_DIR / "loop_gain_spectrum.csv")
    poles = pd.read_csv(DATA_DIR / "pole_continuation.csv")
    lam_c = float(pd.read_csv(DATA_DIR / "critical_point.csv").lam_N_pole.iloc[0])
    fig, axes = plt.subplots(1, 3, figsize=(8.6, 2.9), gridspec_kw={"wspace": 0.55})
    panel_A(axes[0], spec); panel_B(axes[1], poles, lam_c); panel_C(axes[2], poles, lam_c)
    log_ticks(axes[1]); log_ticks(axes[2])
    for ax, letter in zip(axes, "ABC"):
        panel_label(ax, letter, x=-0.24)
    save(fig, "Fig2_static_loop_gain")


if __name__ == "__main__":
    main()
