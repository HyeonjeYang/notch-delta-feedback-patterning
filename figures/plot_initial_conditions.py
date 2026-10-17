"""Figure S5: random starts, paired emergence times and endpoint fields."""
from __future__ import annotations

import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))

import matplotlib.pyplot as plt
import numpy as np

from figures.hex_plot import draw_hex_field
from figures.style import (AQUA,BLUE,ORANGE,MUTED,SEQUENTIAL,TIME_LABEL,
                           apply_style,panel_label,save)
from scripts.initial_condition_robustness import AMPLITUDES,OUT,SEEDS,case_name
from scripts.write_initial_condition_report import verified_results


def main():
    apply_style()
    runs,pairs,summary=verified_results()
    main=runs[runs.downstream=="qss"]
    fig=plt.figure(figsize=(10.7,7.1))
    gs=fig.add_gridspec(2,3,hspace=.6,wspace=.47,left=.075,right=.96,top=.91,bottom=.11)
    axes=[fig.add_subplot(gs[i,j]) for i in range(2) for j in range(3)]
    for ax,point in zip(axes[:2],("P2","P3")):
        for model,color,marker,offset in (("full",BLUE,"o",-.07),("DNA",ORANGE,"s",.07)):
            medians=[]
            for i,amp in enumerate(AMPLITUDES):
                group=main[(main.point==point)&(main.model==model)&(main.amplitude==amp)]
                ax.scatter(i+offset+np.linspace(-.025,.025,len(group)),group.onset_01,
                           color=color,marker=marker,s=14,alpha=.65)
                medians.append(group.onset_01.median())
            ax.plot(np.arange(3)+offset,medians,color=color,label="Full" if model=="full" else "DNA",lw=1.3)
        ax.set_xticks(range(3),[r"$10^{-6}$",r"$10^{-4}$",r"$10^{-2}$"])
        ax.set_xlabel("Relative initial amplitude")
        ax.set_yscale("log")
        ax.set_ylabel(r"Time to CV$(D)=0.1$")
        ax.set_title(f"{point}: paired onset times",loc="left")
        ax.legend(loc="upper right",fontsize=7)
    ax=axes[2]
    for point,color,marker,offset in (("P2",BLUE,"o",-.07),("P3",AQUA,"s",.07)):
        med=[]
        for i,amp in enumerate(AMPLITUDES):
            g=pairs[(pairs.point==point)&(pairs.amplitude==amp)]
            ax.scatter(i+offset+np.linspace(-.025,.025,len(g)),g.aligned_correlation,
                       color=color,marker=marker,s=15,alpha=.7)
            med.append(g.aligned_correlation.median())
        ax.plot(np.arange(3)+offset,med,color=color,label=point,lw=1.3)
    ax.axhline(1,color=MUTED,lw=.8,ls="--")
    ax.set_ylim(.5,1.04)
    ax.set_xticks(range(3),[r"$10^{-6}$",r"$10^{-4}$",r"$10^{-2}$"])
    ax.set_xlabel("Relative initial amplitude")
    ax.set_ylabel("Translation-aligned correlation")
    ax.set_title("Different endpoint fields",loc="left")
    ax.legend(loc="lower right",fontsize=7)
    ax=axes[3]
    for model,color,marker,offset in (("full",BLUE,"o",-.1),("DNA",ORANGE,"s",.1)):
        for i,(point,amp) in enumerate((p,a) for p in ("P2","P3") for a in AMPLITUDES):
            g=main[(main.point==point)&(main.model==model)&(main.amplitude==amp)]
            ax.scatter(i+offset+np.linspace(-.035,.035,len(g)),g.final_k_fraction,
                       color=color,marker=marker,s=16,alpha=.65,label=model if i==0 else None)
    ax.set_xticks(range(6),[r"$10^{-6}$",r"$10^{-4}$",r"$10^{-2}$"]*2)
    ax.text(1,-.26,"P2",transform=ax.get_xaxis_transform(),ha="center",fontsize=8)
    ax.text(4,-.26,"P3",transform=ax.get_xaxis_transform(),ha="center",fontsize=8)
    ax.axvline(2.5,color=MUTED,lw=.8,ls=":")
    ax.set_ylim(0,.5)
    ax.set_ylabel("K power / nonuniform power")
    ax.set_title("Broad-spectrum endpoints",loc="left")
    ax.legend(loc="upper right",fontsize=7)
    arrays=[np.load(OUT/"cases"/(case_name(("P2",1e-4,SEEDS[0],model,"qss"))+".npz"))["final"][0]
            for model in ("full","DNA")]
    vmin=min(a.min() for a in arrays); vmax=max(a.max() for a in arrays)
    for ax,field,title in zip(axes[4:],arrays,("Full: first-seed endpoint","DNA: matched initial state")):
        sm=draw_hex_field(ax,field,vmin,vmax,cmap=SEQUENTIAL,window=10)
        ax.set_title(title,loc="left")
    fig.colorbar(sm,ax=axes[4:],orientation="horizontal",fraction=.045,pad=.08,shrink=.74,label=r"Delta $D$")
    for ax,letter in zip(axes,"ABCDEF"):
        panel_label(ax,letter,x=-.22 if letter in "ABCD" else -.07)
    fig.suptitle("Random initial conditions: slower emergence, distinct endpoint patterns",fontsize=11,x=.075,ha="left",y=.995)
    save(fig,"FigS5_initial_condition_robustness",supplementary=True)


if __name__ == "__main__":
    main()
