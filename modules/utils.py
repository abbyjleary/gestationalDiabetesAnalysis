import argparse
import pandas as pd
import numpy as np
from pydeseq2.ds import DeseqStats

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
