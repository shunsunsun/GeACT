# scRNA-seq
A pipeline for single-cell RNA-seq data analysis, maily for the data from MALBAC-DT.

## Prerequisite
1. Perl (5.010 or later)  
2. R (3.0.1 or later)  

## Pipeline
### Dataset level
For the sequencing data in each dataset (e.g. datasets/13_D1_48):

1. Demultiplexing  
Reads (R1) from each multiplexed sample were split into 96 samples (cells) according to cell barcodes (first 8bp) in reads (R2).  
Unique molecular identifiers (UMIs) (20bp) in reads (R2) were extracted and added into the corresponding name in reads (R1).  
PolyA sequences were trimmed and reads were filtered based on reads length and quality.  
```
bash do_cleanFq.sh
```
2. Read alignment  
Reads (R1) were mapped to the human genome and ERCC sequences.  
Genome assembly: GRCh38.primary_assembly.genome.fa ([Download page](https://www.gencodegenes.org/human/release_26.html))  
Gene annotation: gencode.v26.primary_assembly.annotation.gtf ([Download page](https://www.gencodegenes.org/human/release_26.html))  
```
bash do_hisat2.sh
```
3. Expression calling  
UMI count for each gene was calculated.
```
bash do_htseq.sh
```

Note: `27_210303_07` is a special dataset, where the 5’ primer sequences (19bp) need to be removed in Read 2 before the standard pipeline.  

### Single-organ level (at one specific developmental stage)
For each organ (e.g. pooled_data/01_stomach):

1. Cell pooling (from different datasets)
```
Rscript do_poolData.r
```
2. Filtering (cells / genes)
```
Rscript do_filter.r
```
3. Cell clustering and cell type annotation
```
Rscript do_cluster.r
```

### Multiple-organ level (at one specific developmental stage)
The data from different organs were merged (e.g. pooled_data/All):

1. Data pooling (from different organs)
```
Rscript do_poolData.r
```

2. Dimension reduction
```
Rscript do_cluster.r
```
