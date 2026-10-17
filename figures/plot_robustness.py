"""Figure S4: independent seeds, sample size and parameter-measure dependence."""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from figures.style import BLUE, ORANGE, INK_2, MUTED, DATA_DIR, apply_style, panel_label, save

NAMES = {"primary": "Original bounds", "wide": "Wider bounds"}


def main():
    apply_style()
    d = pd.read_csv(DATA_DIR / "robustness/replicates.csv")
    groups = json.loads((DATA_DIR / "robustness/summary.json").read_text())["groups"]
    gmap = {(g["scenario"], g["n"]): g for g in groups}
    fig, axes = plt.subplots(2, 2, figsize=(8.4, 7.0), gridspec_kw={"hspace": .58, "wspace": .45})
    ax = axes[0, 0]
    for x, scenario, color in ((0, "primary", BLUE), (1, "wide", ORANGE)):
        s = d[d.scenario == scenario]
        g = gmap[scenario, 4000]
        ax.plot(x + np.linspace(-.13, .13, len(s)), s.unstable_fraction * 100, "o", color=color, ms=4.5, alpha=.8)
        mean = 100 * g["mean_unstable_fraction"]
        ax.errorbar(x + .23, mean, yerr=[[mean - 100*g["mc95_low"]], [100*g["mc95_high"] - mean]],
                    fmt="D", color=INK_2, capsize=3, ms=4)
        ax.text(x, 100*s.unstable_fraction.max() + .9, f"{len(s)} seeds", ha="center", fontsize=8)
    ax.axhline(36.325, color=MUTED, lw=1, ls="--")
    ax.text(.98, .70, "original single sample: 36.325%", transform=ax.transAxes, ha="right", fontsize=7, color=INK_2)
    ax.set_xticks([0,1], [NAMES["primary"], NAMES["wide"]]); ax.set_xlim(-.4,1.6)
    ax.set_ylim(27, 40.5); ax.set_ylabel("Nonuniform-unstable points (%)")
    ax.set_title("Instability across seeds and bounds", loc="left")

    ax = axes[0, 1]
    for i, n in enumerate((1000, 2000, 4000, 8000)):
        scenario = "primary" if n == 4000 else "size"
        s = d[(d.scenario == scenario) & (d.n == n)]
        g = gmap[scenario, n]
        ax.plot(i+np.linspace(-.1,.1,len(s)),s.unstable_fraction*100,"o",color=BLUE,ms=3.5,alpha=.65)
        mean = 100*g["mean_unstable_fraction"]
        ax.errorbar(i+.17,mean,yerr=[[mean-100*g["mc95_low"]],[100*g["mc95_high"]-mean]],fmt="D",color=INK_2,ms=4,capsize=3)
    ax.set_xticks(range(4), ["1000\nR=5","2000\nR=5","4000\nR=10","8000\nR=5"])
    ax.set_xlabel("Points per LHS replicate")
    ax.set_ylabel("Nonuniform-unstable points (%)")
    ax.set_title("Sample size in original bounds", loc="left")

    ax = axes[1, 0]
    for offset, scenario, color, marker in ((-.12,"primary",BLUE,"o"),(.12,"wide",ORANGE,"s")):
        s = d[d.scenario == scenario]
        for x, model in enumerate(("DNA","DN")):
            ax.plot(x+offset+np.linspace(-.05,.05,len(s)),s[f"{model}_ratio_median"],marker,
                    color=color,ls="none",ms=4,label=NAMES[scenario] if model=="DNA" else None)
    ax.axhline(1, color=MUTED, lw=.9, ls="--")
    ax.set_yscale("log"); ax.set_xticks([0,1],["DNA reduction","DN reduction"]);ax.set_xlim(-.5,1.5)
    ax.set_ylabel("Median growth ratio, reduced / full")
    ax.set_title("Growth-rate bias across bounds", loc="left")
    ax.legend(loc="upper left", fontsize=7)

    ax = axes[1, 1]
    fields = ["n_gain_mismatch", "DNA_only", "DN_only", "DNA_missed", "DN_missed"]
    for offset, scenario, color, marker in ((-.12,"primary",BLUE,"o"),(.12,"wide",ORANGE,"s")):
        s = d[d.scenario == scenario]
        for x, field in enumerate(fields):
            ax.plot(x+offset+np.linspace(-.045,.045,len(s)),100*s[field]/s.n,marker,
                    color=color,ls="none",ms=4,label=NAMES[scenario] if x==0 else None)
    ax.set_xticks(range(5), ["DC misses", "DNA-only", "DN-only", "Full-only\nvs DNA", "Full-only\nvs DN"], rotation=22, ha="right")
    ax.set_xlim(-.5,4.5);ax.set_ylabel("Points per replicate (%)")
    ax.set_ylim(-.08, max(100*d.DN_only/d.n)*1.18)
    ax.set_title("Classification differences across seeds", loc="left")
    ax.legend(loc="upper right", fontsize=7)
    for ax, letter in zip(axes.flat,"ABCD"):
        panel_label(ax,letter,x=-.21)
    save(fig,"FigS4_sampling_robustness",supplementary=True)


if __name__ == "__main__":
    main()
