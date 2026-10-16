"""Figure 3 -- Control-theoretic mechanism."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

DATA_DIR = Path(__file__).resolve().parent / "data"
OUT_DIR = Path(__file__).resolve().parent / "main"
OUT_DIR.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({"font.size": 10.5, "figure.facecolor": "white",
                      "savefig.facecolor": "white", "savefig.bbox": "tight", "savefig.dpi": 300})

MARK_COLOR = {"baseline": "black", "P2": "#D68910", "P3": "#C0392B"}


def box(ax, xy, text, color, w=1.7, h=0.9, fontsize=10):
    b = FancyBboxPatch((xy[0] - w / 2, xy[1] - h / 2), w, h,
                       boxstyle="round,pad=0.04,rounding_size=0.08", fc=color, ec="black", lw=1.2)
    ax.add_patch(b)
    ax.text(*xy, text, ha="center", va="center", fontsize=fontsize)


def arrow(ax, p0, p1, color="black", label=None, rad=0.0, label_dy=0.0, va="center"):
    a = FancyArrowPatch(p0, p1, arrowstyle="-|>", mutation_scale=16, color=color,
                        connectionstyle=f"arc3,rad={rad}", lw=1.6)
    ax.add_patch(a)
    if label:
        mx, my = (p0[0] + p1[0]) / 2, (p0[1] + p1[1]) / 2 + label_dy
        ax.text(mx, my, label, fontsize=9, ha="center", va=va, bbox=dict(fc="white", ec="none", pad=0.5))


def panel_A(ax):
    ax.set_xlim(-1, 8); ax.set_ylim(-0.4, 3.9); ax.axis("off")
    box(ax, (1, 1.5), r"$z=(D,N,A)$" "\n" r"plant $A_\mu$", "#AED6F1", w=2.6, h=1.5)
    box(ax, (6, 1.5), r"$S=\mathrm{ASC}$" "\n" r"feedback $a_S$", "#F9E79F", w=2.4, h=1.5)
    # Straight, vertically separated arrows with labels offset clear of the
    # arrow lines (above/below) so text never sits on top of the arrowhead.
    arrow(ax, (2.3, 2.55), (4.7, 2.55), label=r"$c$ ($A\to S$)", label_dy=0.25, va="bottom")
    arrow(ax, (4.7, 0.45), (2.3, 0.45), label=r"$b$ ($S\to D$)", label_dy=-0.25, va="top")
    ax.text(3.5, 3.5, r"$\dot z = A_\mu z + bS,\quad \dot S = c^{\!\top}z + a_S S$", fontsize=10,
           ha="center", va="center", bbox=dict(fc="white", ec="#888", pad=4))
    ax.set_title("A. Mode-resolved feedback block diagram", fontsize=11)


def panel_B(ax):
    ax.axis("off")
    text = (
        r"$P_\mu(s) = c^{\!\top}(sI-A_\mu)^{-1} b$" "\n\n"
        r"$L_\mu(s) = P_\mu(s)/(s-a_S)$" "\n\n"
        r"$\det(sI-J_\mu) = \det(sI-A_\mu)\,(s-a_S)\,[\,1-L_\mu(s)\,]$" "\n\n"
        r"Closed-loop poles $\equiv$ roots of $1-L_\mu(s)=0$" "\n"
        "(verified against eig($J_\\mu$) to machine precision," "\n"
        "max coeff. diff. 1.4e-13 -- control/schur_verification.py)"
    )
    ax.text(0.5, 0.5, text, ha="center", va="center", fontsize=10.5, transform=ax.transAxes,
           bbox=dict(fc="#FCF3CF", ec="#B7950B", pad=8))
    ax.set_title("B. Exact Schur-complement relation", fontsize=11)


def panel_C(ax):
    df = pd.read_csv(DATA_DIR / "pole_continuation.csv")
    for mu, marker in [(-3.0, "o"), (6.0, "s")]:
        sub = df[df.mu == mu]
        ax.scatter(sub["re"], sub["im"], s=6, alpha=0.35,
                  color="#C0392B" if mu == -3 else "#2E86AB", marker=marker, label=fr"$\mu={mu:g}$")
        for mark_lam, mark_name in [(1.0, "baseline"), (2.0, "P2"), (5.623, "P3")]:
            m = sub[(sub.lam_N - mark_lam).abs() < 0.01]
            if len(m):
                ax.scatter(m["re"], m["im"], s=90, color=MARK_COLOR[mark_name], marker=marker,
                          edgecolor="black", zorder=5)
    ax.axvline(0, color="#333", lw=1, ls="--")
    ax.set_xlabel(r"$\mathrm{Re}(s)$"); ax.set_ylabel(r"$\mathrm{Im}(s)$")
    ax.set_xlim(-3, 1.5)
    ax.legend(fontsize=8, frameon=False, loc="upper left")
    ax.set_title(r"C. Root locus vs. $\lambda_N$ ($\mu=-3$ circles, $\mu=6$ squares)", fontsize=10.5)


def panel_D(ax):
    df = pd.read_csv(DATA_DIR / "exact_loop_gain.csv")
    ax.plot(df["lam_N"], df["G_exact_minus3"], "-", color="#C0392B", lw=2, label=r"$G_{exact}(-3)=L_{-3}(0)$")
    ax.plot(df["lam_N"], df["G_eff"], "--", color="#7D3C98", lw=1.6, label=r"$G_{eff}$ (heuristic)")
    ax.axhline(1.0, color="#333", lw=1, ls=":", label="unity gain")
    for lam_N, color in [(1.0, "black"), (2.0, "#D68910"), (5.623, "#C0392B")]:
        ax.axvline(lam_N, color=color, lw=0.9, ls=":", alpha=0.6)
    ax.set_xscale("log")
    ax.set_xlabel(r"$\lambda_N$"); ax.set_ylabel("loop gain")
    ax.legend(fontsize=8, frameon=False)
    ax.set_title(r"D. Exact DC loop gain vs. $\lambda_N$, vs. $G_{eff}$", fontsize=11)


def main():
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    panel_A(axes[0, 0]); panel_B(axes[0, 1]); panel_C(axes[1, 0]); panel_D(axes[1, 1])
    fig.suptitle("Figure 3 -- Control-theoretic mechanism", y=1.02, fontsize=14)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "Fig3_control_mechanism.png", dpi=300)
    fig.savefig(OUT_DIR / "Fig3_control_mechanism.pdf")
    plt.close(fig)
    print(f"wrote {OUT_DIR / 'Fig3_control_mechanism.png'}")


if __name__ == "__main__":
    main()
