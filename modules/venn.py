import matplotlib.pyplot as plt
from matplotlib_venn import venn2

from utils import create_ds

def create_venn_diagram(dds, disease_types, control):
    # find stats for both types compared against control
    subsets = []
    for type in disease_types:
        res = create_ds(dds, type, control)
        sig_genes = set(res[(res["padj"] < 0.05) & (res["log2FoldChange"].abs() > 1)].index)
        subsets.append(sig_genes)

    # plot both subsets in a venn diagram
    plt.figure(figsize=(9, 7))
    venn2(
        subsets=subsets,
        set_labels=[f"{type} vs {control}" for type in disease_types],
        set_colors=("#1f77b4", "#ff7f0e"),
        alpha=0.6
    )

    plt.title(f'Common Differentially Expressed Genes (Compared to {control})', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(f"results/venn_diagram.png", dpi=300)
