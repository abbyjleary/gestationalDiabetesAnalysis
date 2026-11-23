import sys
import pandas as pd
import numpy as np
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
from sklearn.decomposition import PCA
from adjustText import adjust_text
import gseapy as gp

def read_inputs(data_infile, metadata_infile):
    # validate data input file and read to dataframe
    try:
        data = pd.read_csv(data_infile, sep='\t', index_col=0, header=0)
    except Exception as e:
        print(f"Problem reading data file({data_infile}): {e}")
        exit(-1)

    # validate metadata input file and read to dataframe
    try:
        metadata = pd.read_csv(metadata_infile, sep='\t', index_col=0, header=0)
    except Exception as e:
        print(f"Problem reading metadata file ({metadata_infile}): {e}")
        exit(-1)

    return data, metadata

def create_pca_deseq(dds):
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
    plt.savefig(f"results/pca_deseq.png", dpi=300, bbox_inches="tight")

def create_pca_sklearn(data):
    # normalize data with log2 transformation
    log_df = np.log2(data + 1)

    # drop genes with no variance
    log_df = log_df.loc[:, log_df.std(axis=0) > 0]

    # z-score each gene across samples (number of stdevs away from mean)
    z = (log_df - log_df.mean(axis=0)) / log_df.std(axis=0)

    # project data into PCA components
    pca = PCA(n_components=2)
    pca_coords = pca.fit_transform(z)

    # create dataframe from PCA components, with labels from original data
    pca_df = pd.DataFrame(pca_coords, columns=["PC1", "PC2"], index=z.index)
    pca_df["Type"] = labels

    # find mean of each group of genes
    group_means = log_df.groupby(labels).mean()

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
    plt.title('PCA of Gene Expression', fontsize=14, fontweight='bold')
    plt.grid(True, linestyle='--', alpha=0.5)

    plt.tight_layout()
    plt.savefig(f"results/pca_sklearn.png", dpi=300, bbox_inches="tight")

def create_ds(dds, type_1, type_2):
    # save various statistics from comparison of the two types
    stats = DeseqStats(dds, contrast=["Type", type_1, type_2])
    stats.summary()
    res = stats.results_df

    # removes results with no adjusted p-value
    res = res.dropna(subset=["padj"]).copy()

    # normalize results with -log10 p-value
    res["neglog10_padj"] = -np.log10(res["padj"])

    return res

def find_top_genes(res):
    # find most significant genes from comparison between two types
    sig = res[(res["padj"] < 0.05) & (res["log2FoldChange"].abs() >= 1)]

    return sig.nsmallest(50, "padj").index

def create_volcano_plot(dds, type_1, type_2):
    res = create_ds(dds, type_1, type_2)

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
    ax.set_xlabel(f"log2 fold change ({type_1} vs {type_2})")
    ax.set_ylabel("-log10 adjusted p-value")
    ax.legend()

    # label most significant genes (up or down regulated)
    top = res.loc[up | down].nsmallest(10, "padj")
    texts = [ax.text(r.log2FoldChange, r.neglog10_padj, g, fontsize=8) for g, r in top.iterrows()]
    adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="-", lw=0.5))

    plt.title(f'{type_2} Gene Expression Compared to {type_1}', fontsize=14, fontweight='bold')

    plt.tight_layout()
    plt.savefig(f"results/volcano_{type_1}_{type_2}.png", dpi=300, bbox_inches="tight")

def create_heatmap(dds, type_1, type_2, type_3, metadata):
    # normalize dataframe
    norm = pd.DataFrame(dds.layers["normed_counts"], index=dds.obs_names, columns=dds.var_names)
    log_norm = np.log2(norm + 1)  # samples x genes

    # find statistics for each comparison and union together to compare all
    res1 = create_ds(dds, type_1, type_2)
    top1 = find_top_genes(res1)

    res2 = create_ds(dds, type_1, type_3)
    top2 = find_top_genes(res2)

    res3 = create_ds(dds, type_2, type_3)
    top3 = find_top_genes(res3)

    top_genes = top1.union(top2).union(top3)

    mat = log_norm[top_genes]

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

def create_gsea(dds, type_1, type_2, gene_sets):
    s = DeseqStats(dds, contrast=["Type", type_1, type_2])
    s.summary()
    res = s.results_df

    # remove empty stats, duplicates, then sort
    rnk = res.dropna(subset=["stat"])["stat"]
    rnk = rnk[~rnk.index.duplicated()]
    rnk = rnk.sort_values(ascending=False)

    if gene_sets == "MSigDB":
        gene_sets_name = "MSigDB_Hallmark_2020"
    elif gene_sets == "GO":
        gene_sets_name = "GO_Biological_Process_2023"

    # compare genes to hallmark set of annotations
    pre = gp.prerank(
        rnk=rnk,
        gene_sets=gene_sets_name,
        min_size=5, max_size=500,
        permutation_num=1000,
        seed=42, threads=4,
        outdir=None,
    )

    pre.res2d[["Term", "NES", "NOM p-val", "FDR q-val"]].head(10)

    # classic enrichment curve for one pathway
    term = pre.res2d.Term.iloc[0]
    gp.gseaplot(rank_metric=pre.ranking, term=term, **pre.results[term])

    # several pathways overlaid
    terms = pre.res2d.Term[:5].tolist()
    gp.gseaplot2(
        terms=terms,
        RESs=[pre.results[t]["RES"] for t in terms],
        hits=[pre.results[t]["hits"] for t in terms],
        rank_metric=pre.ranking,
        figsize=(6, 8),
    )

    # summary dotplot across pathways
    gp.dotplot(pre.res2d, column="FDR q-val", cutoff=0.25, title=f"{type_1} vs {type_2} ({gene_sets})", figsize=(4, 6), ofname=f"results/gsea_{type_1}_{type_2}_MSigDB.png")

def create_venn_diagram(dds, type_1, type_2):
    # find stats for both types compared against control
    ds = create_ds(dds, 'Control', type_1)
    sset1 = set(ds[(ds["padj"] < 0.05) & (ds["log2FoldChange"].abs() > 1)].index)

    ds = create_ds(dds, 'Control', type_2)
    sset2 = set(ds[(ds["padj"] < 0.05) & (ds["log2FoldChange"].abs() > 1)].index)

    # plot both subsets in a venn diagram
    plt.figure(figsize=(9, 7))
    venn2(
        subsets=(sset1, sset2),
        set_labels=(f"{type_1} vs Control", f"{type_2} vs Control"),
        set_colors=("#1f77b4", "#ff7f0e"),
        alpha=0.6
    )

    plt.title('Common Differentially Expressed Genes (Compared to Control)', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(f"results/venn_diagram.png", dpi=300)

# read data and metadata from input files
data, metadata = read_inputs(sys.argv[1], sys.argv[2])

data = data.astype(int)

# create series of labels from metadata
common = data.index.intersection(metadata.index)
data = data.loc[common]
labels = metadata.loc[common, "Type"]

# filter low-count genes
genes_to_keep = data.columns[data.sum(axis=0) >= 10]
data = data[genes_to_keep]

# create deseq dataset to store counts and metadata to be analyzed later
dds = DeseqDataSet(counts=data, metadata=metadata, design_factors="Type")
dds.deseq2()

# create_pca_deseq(dds)
# create_pca_sklearn(data)
#
# create_volcano_plot(dds, "Control", "GDMA1")
# create_volcano_plot(dds, "Control", "GDMA2")
# create_volcano_plot(dds, "GDMA1", "GDMA2")
#
create_heatmap(dds, "Control", "GDMA1", "GDMA2", metadata)

# create_gsea(dds, "Control", "GDMA1", "MSigDB")
# create_gsea(dds, "Control", "GDMA2", "MSigDB")
# create_gsea(dds, "GDMA1", "GDMA2", "MSigDB")
#
# create_gsea(dds, "Control", "GDMA1", "GO")
# create_gsea(dds, "Control", "GDMA2", "GO")
# create_gsea(dds, "GDMA1", "GDMA2", "GO")

# create_venn_diagram(dds, "GDMA1", "GDMA2")