import sys
import pandas as pd

infile = sys.argv[1]
outfile = sys.argv[2]

metadata = pd.read_csv(infile, sep='\t', header=None)
df = pd.DataFrame()

series = {}
for sample in metadata.iloc[:,0]:
    df = pd.read_csv(('data/' + sample + '.tsv'), sep="\t", header=None, names=["ensembl", "gene", "count"])
    series[sample] = df.set_index("gene")["count"]  # sample name from filename

counts = pd.DataFrame(series).T  # rows = samples, columns = genes
counts.index.name = None
counts.columns.name = None

counts.to_csv(outfile, sep='\t', index=True)

print(counts.columns)
print(counts.index)