import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

from utils import create_ds, find_top_genes

def create_heatmap(dds, disease_types, control, metadata):
    # normalize dataframe
    norm = pd.DataFrame(dds.layers["normed_counts"], index=dds.obs_names, columns=dds.var_names)
    log_norm = np.log2(norm + 1)  # samples x genes
    top_genes = set()

    for type in disease_types:
        res = create_ds(dds, type, control)
        top_genes = top_genes.union(find_top_genes(res))

    mat = log_norm[list(top_genes)]

    types = metadata.loc[mat.index, "Type"]
    palette = dict(zip(types.unique(), sns.color_palette("Greens", types.nunique())))
    col_colors = types.map(palette)

    g = sns.clustermap(
        mat.T,  # genes x samples
        z_score=0,  # z-score each gene (row) across samples
        cmap="RdBu_r", center=0, vmin=-2, vmax=2,
        col_colors=col_colors,
        xticklabels=True, yticklabels=False,
        figsize=(8, 10),
    )

    # legend for the group colors
    for t, c in palette.items():
        g.ax_col_dendrogram.bar(0, 0, color=c, label=t, linewidth=0)
    g.ax_col_dendrogram.legend(title="Type", loc="center", ncol=len(palette))
    g.ax_cbar.set_title("Z-score", fontsize=10)

    g.fig.suptitle('Gene Expression by Sample and Type', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(f"results/heatmap.png", dpi=300)
