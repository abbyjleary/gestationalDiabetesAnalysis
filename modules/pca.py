import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
import seaborn as sns

def create_pca(dds, metadata, labels):
    # create a variant stabilization transformation based on the deseq dataset
    dds.vst(use_design=False)
    vst = pd.DataFrame(
        dds.layers["vst_counts"],
        index=dds.obs_names,
        columns=dds.var_names,
    )

    # fit PCA to data
    pca = PCA(n_components=2)
    coords = pca.fit_transform(vst)

    pca_df = pd.DataFrame(coords, columns=["PC1", "PC2"], index=vst.index)

    pca_df["Type"] = metadata.loc[vst.index, "Type"].values

    plt.figure(figsize=(9, 7))

    # create plot for samples plotted on PCA
    pca_df["Type"] = labels
    ax = sns.scatterplot(
        data=pca_df,
        x='PC1',
        y='PC2',
        hue='Type',
        s=100,
        palette='Set1',
        alpha=0.6
    )

    # label each point with sample name
    for x, y, label in zip(pca_df["PC1"], pca_df["PC2"], pca_df.index):
        ax.annotate(
            label,
            (x, y),
            xytext=(5, -5),
            textcoords="offset points",
            va="center",
        )


    # Label PC axes with explained variance percentage
    var_exp = pca.explained_variance_ratio_ * 100
    plt.xlabel(f'PC1 ({var_exp[0]:.1f}% Variance Explained)', fontsize=12)
    plt.ylabel(f'PC2 ({var_exp[1]:.1f}% Variance Explained)', fontsize=12)

    plt.title(f'PCA of Gene Expression', fontsize=14, fontweight='bold')
    plt.grid(True, linestyle='--', alpha=0.5)

    plt.tight_layout()
    plt.savefig(f"results/pca.png", dpi=300, bbox_inches="tight")
