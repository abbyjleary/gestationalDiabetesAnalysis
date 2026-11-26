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
import argparse

def get_args():
    parser = argparse.ArgumentParser(description="Differential GeneExpression Analysis")
    parser.add_argument("data_infile", help="Path to the input data file (tab-separated)")
    parser.add_argument("metadata_infile", help="Path to the metadata file (tab-separated)")
    parser.add_argument("control", help="Label in metadata tsv of control group")
    parser.add_argument("--run_all", action="store_true", help="Run all analyses (PCA, volcano plots, heatmap, GSEA, Venn diagram)")
    parser.add_argument("--run_pca", action="store_true", help="Run PCA analysis")
    parser.add_argument("--run_volcano", action="store_true", help="Run volcano plot analysis")
    parser.add_argument("--run_heatmap", action="store_true", help="Run heatmap analysis")
    parser.add_argument("--run_gsea", action="store_true", help="Run GSEA analysis")
    parser.add_argument("--run_venn", action="store_true", help="Run Venn diagram analysis")
    parser.add_argument("--gsea_gene_set", choices=["MSigDB", "GO-BP"], help="Gene set for GSEA analysis (MSigDB or GO-BP)")

    return parser.parse_args()

def read_inputs(data_infile, metadata_infile):
    # validate data input file and read to dataframe
    try:
        data = pd.read_csv(data_infile, sep='\t', index_col=0, header=0)
    except Exception as e:
        print(f"Problem reading data file({data_infile}): {e}")
        exit(-1)

    # validate metadata input file and read to dataframe
    try:
        metadata = pd.read_csv(metadata_infile, sep='\t', index_col=0, names=["Sample", "Type"])
    except Exception as e:
        print(f"Problem reading metadata file ({metadata_infile}): {e}")
        exit(-1)

    return data, metadata

def create_pca(dds):
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

def create_ds(dds, variant, control):
    # save various statistics from comparison of the two types
    stats = DeseqStats(dds, contrast=["Type", variant, control])
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

def create_gsea(dds, variant, control, gene_sets):
    s = DeseqStats(dds, contrast=["Type", variant, control])
    s.summary()
    res = s.results_df

    # remove empty stats, duplicates, then sort
    rnk = res.dropna(subset=["stat"])["stat"]
    rnk = rnk[~rnk.index.duplicated()]
    rnk = rnk.sort_values(ascending=False)

    if gene_sets == "MSigDB":
        gene_sets_name = "MSigDB_Hallmark_2020"
    elif gene_sets == "GO-BP":
        gene_sets_name = "GO_Biological_Process_2023"

    # compare genes to set of pathways
    pre = gp.prerank(
        rnk=rnk,
        gene_sets=gene_sets_name,
        min_size=5, max_size=500,
        permutation_num=1000,
        seed=42, threads=4,
        outdir=None,
    )

    # find the most significant up and down regulated pathways
    by_nes = pre.res2d.copy()
    by_nes["NES"] = by_nes["NES"].astype(float)
    by_nes = by_nes.sort_values("NES")
    plot_df = pd.concat([by_nes.head(10), by_nes.tail(10)]).drop_duplicates("Term")

    # plot the pathways
    for term in [by_nes.Term.iloc[0], by_nes.Term.iloc[-1]]:
        gp.gseaplot(rank_metric=pre.ranking, term=term, **pre.results[term])

    terms = by_nes.Term[:5].tolist() + by_nes.Term[-5:].tolist()
    gp.gseaplot2(
        terms=terms,
        RESs=[pre.results[t]["RES"] for t in terms],
        hits=[pre.results[t]["hits"] for t in terms],
        rank_metric=pre.ranking,
        figsize=(6, 8),
    )

    gp.dotplot(plot_df, column="FDR q-val", cutoff=1, top_term=len(plot_df),
               title=f"{variant} vs {control}", figsize=(4, 8),
               ofname=f"results/gsea_{variant}_{control}_{gene_sets}.png")
    
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

args = get_args()

if args.run_all:
    args.run_pca = args.run_volcano = args.run_heatmap = args.run_gsea = args.run_venn = True
    
# read data and metadata from input files
data, metadata = read_inputs(args.data_infile, args.metadata_infile)

control = args.control

disease_types = metadata["Type"].unique()
disease_types = [t for t in disease_types if t != control]

if args.run_gsea:
    if args.gsea_gene_set is None:
        print("Please specify a gene set for GSEA analysis using --gsea_gene_set (MSigDB or GO-BP).")
        sys.exit(1)

if args.run_venn:
    if len(disease_types) < 2:
        print("At least two disease types are required to create a Venn diagram.")
        sys.exit(1)

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


if args.run_pca:
    create_pca(dds)

if args.run_volcano:
    for type in disease_types:
        create_volcano_plot(dds, type, control)

if args.run_heatmap:
    create_heatmap(dds, disease_types, control, metadata)

if args.run_gsea:
    for type in disease_types:
        create_gsea(dds, type, control, args.gsea_gene_set)

if args.run_venn:
    create_venn_diagram(dds, disease_types, control)