#!/usr/bin/env python
"""
This script aims to cross reference the mentioned mutated genes by Lohr et al with my pseudobulk DE results and visualise as a dot plot
x = log2FC (Memory B vs Plasma)
y = MM mutation frequency % 
size = -log10(padj)
colour = direction
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# gene, log2FC, padj, MM_mut_freq_%
DATA = [
    ("PRDM1",  -5.59, 2.4e-185,  5),
    ("TP53",   +1.52, 2.9e-9,    8),
    ("RB1",    +1.02, 2.2e-22,   2),
    ("BRAF",   +1.02, 1.0e-9,    6),
    ("DIS3",   +0.97, 1.8e-7,   11),
    ("TRAF3",  +0.56, 1.6e-5,    5),
    ("CYLD",   +0.31, 0.03,      3),
    ("NRAS",   -0.32, 0.09,     20),
    ("KRAS",   -0.31, 0.13,     23),
]

BLUE  = "#1565C0"
RED   = "#C62828"
GREY  = "#B0BEC5"
SIG   = 0.05

fig, ax = plt.subplots(figsize=(9, 7), facecolor="white")

for gene, lfc, padj, freq in DATA:
    sig = padj < SIG
    if not sig:
        col = GREY
    elif lfc < 0:
        col = RED
    else:
        col = BLUE

    size = max(80, -np.log10(padj) * 12)
    size = min(size, 800)

    ax.scatter(lfc, freq, s=size, c=col,
               alpha=0.85, edgecolors="white",
               linewidths=1.5, zorder=3)

    # label offsets tuned per gene to avoid overlap
    offsets = {
        "PRDM1": (-0.25, 1.5),
        "KRAS":  ( 0.15, 1.2),
        "NRAS":  ( 0.15,-1.8),
        "TP53":  ( 0.15, 1.2),
        "BRAF":  ( 0.15,-1.8),
        "DIS3":  (-0.5,  1.2),
        "TRAF3": ( 0.15, 1.2),
        "RB1":   ( 0.15,-1.8),
        "CYLD":  ( 0.15, 1.2),
    }
    dx, dy = offsets.get(gene, (0.15, 1.0))
    ha = "right" if dx < 0 else "left"
    ax.text(lfc + dx, freq + dy, gene,
            fontsize=10, fontweight="bold",
            color=col if col != GREY else "#666666",
            ha=ha, va="center", zorder=4)

# reference lines
ax.axvline(0, color="black", lw=1.0, zorder=1)
ax.axvline(-1, color=GREY, lw=0.7, ls="--", alpha=0.5)
ax.axvline(+1, color=GREY, lw=0.7, ls="--", alpha=0.5)
ax.axhspan(0, 7,  alpha=0.04, color=GREY)
ax.axhspan(7, 25, alpha=0.04, color=GREY)

# axis labels
ax.set_xlabel("log₂ Fold Change  (Memory B → Plasma, positive = higher in Memory B)",
              fontsize=10, fontweight="bold")
ax.set_ylabel("MM mutation frequency %\n(Lohr et al. 2014, n=203 patients)",
              fontsize=10, fontweight="bold")

ax.set_xlim(-8, 4)
ax.set_ylim(-2, 27)

# direction labels
ax.text(-7.7, -1.6, "← Higher in Plasma",
        fontsize=9, color=RED, fontstyle="italic")
ax.text(3.7, -1.6, "Higher in Memory B →",
        fontsize=9, color=BLUE, fontstyle="italic", ha="right")

ax.set_title(
    "MM-associated genes (Lohr et al. 2014) vs Pseudobulk DE\n"
    "Dot size = statistical significance  ·  Colour = direction",
    fontsize=12, fontweight="bold", pad=12)

# legend
for label, col in [("Higher in Memory B (sig.)", BLUE),
                   ("Higher in Plasma (sig.)",   RED),
                   ("Not significant",            GREY)]:
    ax.scatter([], [], c=col, s=120, label=label,
               edgecolors="white", linewidths=1.5, alpha=0.85)
ax.legend(loc="upper left", fontsize=9,
          framealpha=0.9, edgecolor="lightgrey")

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.grid(axis="y", color="#ECEFF1", lw=0.8)

fig.tight_layout()
out = "/rds/projects/r/russdr-bb-data/wao597/lohr_dotplot.png"
fig.savefig(out, dpi=250, bbox_inches="tight", facecolor="white")
plt.close(fig)
print(f"Saved: {out}")
