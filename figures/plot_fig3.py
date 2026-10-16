"""Figure 3: which kinetic steps move the stationary margin.

A  d ln G(-3) / d ln p at p0, total (bar), split into the three stage gains
   kappa_D (S -> D), T_-3 (D -> neighbour A) and kappa_S (A -| S).
B  static margin (d ln G) against the growth-rate sensitivity
   (d lambda_pattern / d ln p). lam_S and d_S act on G only through their
   ratio, so their static effects are equal and opposite. Their dynamic
   effects are not.
C, D  lambda_pattern on two parameter planes. The black contour is
   lambda_pattern = 0 and the dashed white contour is G(-3) = 1; the two coincide.
Data: scripts/kinetic_control.py.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm

from figures.style import (AQUA, BLUE, DATA_DIR, DIVERGING, INK, INK_2, MUTED, ORANGE, PARAM_TEX,
                           apply_style, panel_label, save)
from model.params import LAM_N_P2, LAM_N_P3, P0


def panel_A(ax, g):
    g = g.sort_values("dlnG_total", key=np.abs).reset_index(drop=True)
    y = np.arange(len(g))
    ax.barh(y, g.dlnG_total, height=0.62, color="#d6d5ce", edgecolor="none", label="total")
    for col, colour, marker, name in [("via_kappa_D", BLUE, "o", r"via $\kappa_D$ (S$\to$D)"),
                                      ("via_T", ORANGE, "s", r"via $T_{-3}$ (D$\to$A)"),
                                      ("via_kappa_S", AQUA, "^", r"via $\kappa_S$ (A$\dashv$S)")]:
        ax.plot(g[col], y, marker, color=colour, ms=4.2, mec="white", mew=0.5, ls="none", label=name)
    ax.axvline(0, color=INK_2, lw=0.8)
    ax.set_yticks(y, [PARAM_TEX[p] for p in g.parameter])
    ax.grid(axis="y", visible=False)
    ax.set_xlabel(r"$\partial \ln G(-3)\,/\,\partial \ln p$ at p0")
    ax.legend(loc="lower right", handletextpad=0.3)
    ax.set_title("Static loop gain: stage contributions", loc="left")


def panel_B(ax, g, s):
    m = g.merge(s, on="parameter")
    ax.axhline(0, color=MUTED, lw=0.6); ax.axvline(0, color=MUTED, lw=0.6)
    ax.plot(m.dlnG_total, m.S_lambda_pattern, "o", color=BLUE, ms=4.5, mec="white", mew=0.6)
    offsets = {"lam_ASC": (6, -9), "d_ASC": (6, 3), "n_ASC": (-14, -10), "lam_N": (5, -9), "d_A": (5, 2), "Kd": (5, 2)}
    for _, r in m.iterrows():
        if abs(r.dlnG_total) > 0.5 or r.parameter in offsets:
            ax.annotate(PARAM_TEX[r.parameter], (r.dlnG_total, r.S_lambda_pattern),
                        xytext=offsets.get(r.parameter, (4, 2)), textcoords="offset points", fontsize=7.5, color=INK)
    ax.set_xlabel(r"static: $\partial \ln G(-3)/\partial \ln p$")
    ax.set_ylabel(r"dynamic: $\partial \lambda_{\rm pattern}/\partial \ln p$")
    ax.set_title(r"Static vs. growth-rate sensitivity", loc="left")


def phase_map(ax, df, pair, xlabel, ylabel, marks):
    sub = df[df.pair == pair]
    lam = sub.pivot(index="y", columns="x", values="nonuniform_max_re_lambda")
    G = sub.pivot(index="y", columns="x", values="G_minus3")
    X, Y = np.meshgrid(lam.columns.values, lam.index.values)
    vmax = float(np.nanmax(np.abs(lam.values)))
    pc = ax.pcolormesh(X, Y, lam.values, cmap=DIVERGING, norm=TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax),
                       shading="nearest", rasterized=True)
    ax.contour(X, Y, lam.values, levels=[0], colors=INK, linewidths=1.6)
    ax.contour(X, Y, G.values, levels=[1], colors="white", linewidths=1.0, linestyles="--")
    for (x, y, marker, name) in marks:
        ax.plot(x, y, marker, ms=6, mfc="white", mec=INK, mew=1.1, ls="none")
        ax.annotate(name, (x, y), xytext=(0, 6), textcoords="offset points", ha="center", va="bottom",
                    fontsize=7.5, color=INK, bbox=dict(fc="white", ec="none", alpha=0.85, pad=0.8))
    ax.grid(False)
    ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)
    return pc


def main():
    apply_style()
    g = pd.read_csv(DATA_DIR / "stage_gain_sensitivity.csv")
    s = pd.read_csv(DATA_DIR / "local_sensitivity.csv")
    maps = pd.read_csv(DATA_DIR / "two_dimensional_phase_maps.csv")
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 6.2), gridspec_kw={"wspace": 0.6, "hspace": 0.42})
    panel_A(axes[0, 0], g)
    panel_B(axes[0, 1], g, s)
    marks = [(1.0, P0.d_A, "o", "p0"), (LAM_N_P2, P0.d_A, "^", "P2"), (LAM_N_P3, P0.d_A, "s", "P3")]
    pc = phase_map(axes[1, 0], maps, "lam_N_x_d_A", r"$\lambda_N$", r"$d_A$", marks)
    axes[1, 0].set_xscale("log"); axes[1, 0].set_yscale("log")
    axes[1, 0].set_title(r"$\lambda_N \times d_A$", loc="left")
    fig.colorbar(pc, ax=axes[1, 0], fraction=0.046, pad=0.03, label=r"$\lambda_{\rm pattern}$")
    pc = phase_map(axes[1, 1], maps, "n_ASC_x_Kd", r"$n_S$", r"$K_d$", [(P0.n_ASC, P0.Kd, "o", "p0")])
    axes[1, 1].set_yscale("log")
    axes[1, 1].set_title(r"$n_S \times K_d$", loc="left")
    fig.colorbar(pc, ax=axes[1, 1], fraction=0.046, pad=0.03, label=r"$\lambda_{\rm pattern}$")
    for ax, letter in zip(axes.flat, "ABCD"):
        panel_label(ax, letter, x=-0.22)
    save(fig, "Fig3_kinetic_control")


if __name__ == "__main__":
    main()
