import sys
from pydeseq2.dds import DeseqDataSet

from modules.pca import create_pca
from modules.heatmap import create_heatmap
from modules.volcano import create_volcano_plot
from modules.gsea import create_gsea
from modules.venn import create_venn_diagram
from modules.utils import read_inputs, get_args

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
    create_pca(dds, metadata, labels)

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