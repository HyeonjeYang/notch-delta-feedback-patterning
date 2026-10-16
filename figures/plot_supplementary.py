"""Supplementary figures.

S1  the growth-rate validation under the v1 initial state (eigenvector
    scaled componentwise by the equilibrium) and under the corrected
    exact-eigenvector state, at P2 and P3. The fitting window is shaded.
S2  one-parameter continuations of lambda_pattern for all 13 parameters
    about p0. The dashed line marks where G(-3) crosses 1.
Data: scripts/growth_validation.py, scripts/kinetic_control.py.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from figures.style import (BLUE, CONTEXT, DATA_DIR, INK, INK_2, MUTED, PARAM_TEX, TIME_LABEL, apply_style,
                           panel_label, save)
from model.params import P0
from stability.sensitivity import PARAM_NAMES


def figure_S1():
    ts = pd.read_csv(DATA_DIR / "growth_timeseries.csv")
    summ = pd.read_csv(DATA_DIR / "growth_validation.csv").set_index("run")
    fam = pd.read_csv(DATA_DIR / "timescale_family.csv")
    k = int(np.argmin(np.abs(np.log(fam.alpha))))
    D_star = {"P2": float(fam.P2_D_star.iloc[k]), "P3": float(fam.P3_D_star.iloc[k])}
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.8), gridspec_kw={"wspace": 0.3})
    for ax, point, t_max in [(axes[0], "P2", 600), (axes[1], "P3", 300)]:
        for run, colour, lw, name in [(f"full_{point}", BLUE, 1.6, "exact eigenvector (this work)"),
                                      (f"full_{point}_v1ic", CONTEXT, 1.4, "componentwise-scaled (v1)")]:
            g = ts[ts.run == run]
            s = summ.loc[run]
            ax.plot(g.time, g.a_D / D_star[point], "-", color=colour, lw=lw,
                    label=f"{name}: fit {s.fitted_slope:.4f}, early [0,8] {s.early_slope_t0_8:.3f}")
        s = summ.loc[f"full_{point}"]
        a0 = ts[ts.run == f"full_{point}"].a_D.iloc[0] / D_star[point]
        t = np.linspace(0, s.window_end, 40)
        ax.axvspan(0, s.window_end, color=BLUE, alpha=0.08, lw=0)
        ax.plot(t, a0 * np.exp(s.predicted_lambda * t), "--", color=INK, lw=0.9,
                label=f"predicted pole {s.predicted_lambda:.6f}")
        ax.set_yscale("log"); ax.set_xlim(0, t_max)
        ax.set_xlabel(TIME_LABEL); ax.set_ylabel(r"K-mode amplitude of $D$ / $D^*$")
        ins = ax.inset_axes([0.52, 0.1, 0.44, 0.38])
        for run, colour in [(f"full_{point}", BLUE), (f"full_{point}_v1ic", CONTEXT)]:
            g = ts[(ts.run == run) & (ts.time <= 15)]
            ins.plot(g.time, g.a_D / g.a_D.iloc[0], "-", color=colour, lw=1.3)
        tt = np.linspace(0, 15, 40)
        ins.plot(tt, np.exp(s.predicted_lambda * tt), "--", color=INK, lw=0.8)
        ins.set_yscale("log"); ins.set_xlim(0, 15)
        ins.tick_params(labelsize=6); ins.set_title(r"$a_D(t)/a_D(0)$, $t\leq15$", fontsize=6.8)
        ax.legend(loc="upper left", bbox_to_anchor=(-0.02, -0.22), fontsize=6.6)
        lam_N = 2.0 if point == "P2" else 5.623
        ax.set_title(rf"{point} ($\lambda_N={lam_N:g}$)", loc="left")
    for ax, letter in zip(axes, "AB"):
        panel_label(ax, letter, x=-0.2)
    save(fig, "FigS1_initial_state_comparison", supplementary=True)


def figure_S2():
    one = pd.read_csv(DATA_DIR / "one_dimensional_scans.csv")
    fig, axes = plt.subplots(3, 5, figsize=(7.4, 4.6), gridspec_kw={"wspace": 0.45, "hspace": 0.75})
    for ax, name in zip(axes.flat, PARAM_NAMES):
        g = one[one.parameter == name].sort_values("value")
        ax.plot(g.value, g.nonuniform_max_re_lambda, "-", color=BLUE, lw=1.3)
        ax.axhline(0, color=MUTED, lw=0.6)
        ax.axvline(getattr(P0, name), color=INK_2, lw=0.7, ls=":")
        cross = np.where(np.diff(np.sign(g.G_minus3.values - 1.0)) != 0)[0]
        for i in cross:
            ax.axvline(np.sqrt(g.value.values[i] * g.value.values[i + 1]), color=INK, lw=0.8, ls="--")
        if name not in ("n_d", "n_ASC"):
            ax.set_xscale("log")
        ax.set_title(PARAM_TEX[name], fontsize=8.5)
        ax.tick_params(labelsize=6.5)
    for ax in axes.flat[len(PARAM_NAMES):]:
        ax.axis("off")
    axes.flat[-1].text(0, 0.6, "dotted: p0 value\ndashed: $G(-3)=1$\ny: $\\lambda_{\\rm pattern}$",
                       fontsize=7.5, color=INK, transform=axes.flat[-1].transAxes)
    save(fig, "FigS2_one_parameter_continuations", supplementary=True)


if __name__ == "__main__":
    apply_style()
    figure_S1()
    figure_S2()
