import matplotlib.pyplot as plt
import numpy as np
from adjustText import adjust_text

from utils import create_ds

def create_volcano_plot(dds, variant, control):
    res = create_ds(dds, variant, control)

    # find up and down regulated genes by comparing against log-fold change and adjusted p-value
    lfc_cut, padj_cut = 1, 0.05
    up = (res["log2FoldChange"] >= lfc_cut) & (res["padj"] < padj_cut)
    down = (res["log2FoldChange"] <= -lfc_cut) & (res["padj"] < padj_cut)
    ns = ~(up | down)

    # create subplots for each regulation (up, down, and not significant)
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(res.loc[ns, "log2FoldChange"], res.loc[ns, "neglog10_padj"], s=8, c="lightgrey", label="NS")
    ax.scatter(res.loc[up, "log2FoldChange"], res.loc[up, "neglog10_padj"], s=10, c="tab:red", label=f"Up ({up.sum()})")
    ax.scatter(res.loc[down, "log2FoldChange"], res.loc[down, "neglog10_padj"], s=10, c="tab:blue",
               label=f"Down ({down.sum()})")

    # create separators between each plot
    ax.axvline(lfc_cut, ls="--", c="k", lw=0.8)
    ax.axvline(-lfc_cut, ls="--", c="k", lw=0.8)
    ax.axhline(-np.log10(padj_cut), ls="--", c="k", lw=0.8)

    # create labels and legend
    ax.set_xlabel(f"log2 fold change ({variant} vs {control})")
    ax.set_ylabel("-log10 adjusted p-value")
    ax.legend()

    # label most significant genes (up or down regulated)
    top = res.loc[up | down].nsmallest(10, "padj")
    texts = [ax.text(r.log2FoldChange, r.neglog10_padj, g, fontsize=8) for g, r in top.iterrows()]
    adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="-", lw=0.5))

    plt.title(f'{variant} Gene Expression Compared to {control}', fontsize=14, fontweight='bold')

    plt.tight_layout()
    plt.savefig(f"results/volcano_{variant}_{control}.png", dpi=300, bbox_inches="tight")