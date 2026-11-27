# Gestational Diabetes Data Analysis

## main.py Arguments
| Argument / Flag | Required | Description |
| --- | --- | --- |
| `<data_infile> `| Yes | Path to the input data file (tab-separated) |
| `<metadata_infile>` | Yes | Path to the input metadata file (tab-separated) |
| `<control>` | Yes | Label in metadata tsv of control group |
| `--run_all` | No | Run all analyses (PCA, volcano plots, heatmap, GSEA, Venn diagram) |
| `--run_pca` | No | Run PCA analysis |
| `--run_volcano` | No | Run volcano plot analysis |
| `--run_heatmap` | No | Run heatmap analysis |
| `--run_gsea` | No | Run GSEA analysis |
| `--run_venn` | No | Run Venn diagram analysis |
| `--gsea_gene_set` | Yes, if `--run_gsea` or `--run_all` | Gene set for GSEA analysis (MSigDB or GO-BP) |


## Repository Layout

```
gestationalDiabetesAnalysis/
├── main.py
├── data/
│   ├── metadata.tsv (disease types of each sample)
│   ├── data.tsv (master list of gene counts)
│   └── SRR*.tsv (gene counts per sample)
├── modules/
│   ├── gsea.py
│   ├── heatmap.py
│   ├── pca.py
│   ├── venn.py
│   ├── volcano.py
│   └── utils.py (functions for data input and manipulation)
└── results/
    └── images of results graphs
```

## Reference Research Paper
#### Discrete placental gene expression signatures of diabetes during pregnancy
https://pubmed.ncbi.nlm.nih.gov/38763341/

## Data
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE249311

| Type    | Sample | NCBI Link                                                      |
|---------|--------|----------------------------------------------------------------|
| Control | SRR27056481 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM7933405  |
| Control | SRR27056480 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM7933406  |
| Control | SRR27056455 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM7933407  |
| GDMA1   | SRR27056470 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM7933392  |
| GDMA1   | SRR27056469 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM7933393  |
| GDMA1   | SRR27056468 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM7933394  |
| GDMA2   | SRR27056487 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM7933399  |
| GDMA2   | SRR27056486 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM7933400  |
| GDMA2   | SRR27056485 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM7933401  |

## Pipeline

1. FastQC (v0.11.7)
2. Trimmomatic (v0.38)
3. FastQC (v0.11.7)
4. STAR (v2.7.9a) aligning against GRCh38-2020-A
5. Picard (v2.20.8)
6. HTSeq (v0.11.2)

### FastQC

Use for quality control of FastQ files.

```
fastqc *_{1,2}.fastq.gz
```

### Trimmomatic

Use to trim barcodes and adapters.

```    
trimmomatic PE -threads 4 -phred33 \\
        ${sample_id}_1.fastq ${sample_id}_2.fastq \\
        ${sample_id}_1.trimmed.fastq.gz ${sample_id}_1.fastq.gz \\
        ${sample_id}_2.trimmed.fastq.gz ${sample_id}_2.fastq.gz \\
        ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \\
        SLIDINGWINDOW:4:15 MINLEN:50 \\
        2> ${sample_id}.trimmomatic.log
```

### FastQC

Use to compare quality of FastQ files after trimming.

```
fastqc *_{1,2}.trimmed.fastq.gz
```

### STAR

Align sequences to the human genome (GRCh38).

```
mkdir star_index
STAR --runMode genomeGenerate \\
  --runThreadN 8 \\
  --genomeDir star_index \\
  --genomeFastaFiles refdata-gex-GRCh38-2020-A/fasta/genome.fa \\
  --sjdbGTFfile refdata-gex-GRCh38-2020-A/genes/genes.gtf \\
  --sjdbOverhang 99

STAR --runMode alignReads \\
  --runThreadN 8 \\
  --genomeDir star_index \\
  --readFilesIn ${sample_id}_1.trimmed.fastq.gz ${sample_id}_2.trimmed.fastq.gz \\
  --readFilesCommand zcat \\
  --outFilterType BySJout \\
  --outFilterMultimapNmax 20 \\
  --alignSJoverhangMin 8 \\
  --alignSJDBoverhangMin 1 \\
  --outFilterMismatchNmax 999 \\
  --outFilterMismatchNoverReadLmax 0.04 \\
  --alignIntronMin 20 \\
  --alignIntronMax 1000000 \\
  --alignMatesGapMax 1000000 \\
  --outSAMtype BAM SortedByCoordinate \\
  --outFileNamePrefix ${sample_id}
```

### Picard

Remove duplicate reads from alignments.

```
picard MarkDuplicates \\
  INPUT=${sample_id}Aligned.sortedByCoord.out.bam \\
  OUTPUT=${sample_id}.dedup.bam \\
  METRICS_FILE=${sample_id}.dup_metrics.txt \\
  REMOVE_DUPLICATES=true \\
  CREATE_INDEX=true \\
  VALIDATION_STRINGENCY=SILENT
```

### HTSeq

Convert BAM files to tsv with gene counts for further analysis.

```
htseq-count -f bam -s no -r pos \\
  --additional-attr=gene_name \\
  ${sample_id}.dedup.bam null > ${sample_id}.htseq.counts.txt
```