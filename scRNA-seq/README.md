# scRNA-seq data analysis pipeline
A pipeline for single-cell RNA-seq data analysis, maily for the data from MALBAC-DT.

## Prerequisite
1. Perl (5.16.3)  
2. R (3.5.1)  
3. Slurm (18.08.7)

## Pipeline
### Dataset level
Dataset means the data from one sequencing experiment (e.g. datasets/13_D1_48).

Pipeline:
1. Demultiplexing  
Reads (R1) were assigned into one of 96 cells according to the [Cell barcode table](/scRNA-seq/experiments/cell_barcodes_0x96.txt) localized in the first 8bp in reads (R2).  
Unique molecular identifiers (UMIs) (20bp) in reads (R2) were extracted and added into the corresponding name in reads (R1).  
PolyA sequences were trimmed and reads (R1) were filtered based on reads length and quality.  
```
bash do_cleanFq.sh
```
2. Read alignment  
Reads (R1) were mapped to the human genome and ERCC sequences.

| Information | Version | Path |
| :------ | :------ | :------ |
| Genome assembly | GRCh38 (primary assembly) | [Download page](https://www.gencodegenes.org/human/release_26.html) |
| Gene annotation | Gencode v26 (primary assembly) | [Download page](https://www.gencodegenes.org/human/release_26.html) |
| ERCC | - | [Download](https://tools.thermofisher.com/content/sfs/manuals/ERCC92.zip) |

```
bash do_hisat2.sh
```
3. Expression calling  
First, in each cell, the UMI count for each gene was calculated, where similar UMIs (hamming distance <= 2) were counted only once.  
Then, the UMI counts in different cells were merged to produce a gene expression matrix (row: gene, column: cell).  
Gene IDs were transformed into gene names (e.g. ENSG00000223972.5 -> DDX11L1) based on a [Gene ID mapping table](Data/gene_ID2Name_fixed.txt?raw=true).  
Since some different genes contain the same gene name (e.g. ENSG00000278757.1 and ENSG00000283136.1 -> U6), in this Gene ID mapping table, the suffix was added into these names (ENSG00000278757.1 -> U6--1, ENSG00000283136.1 -> U6--2).  
```
bash do_htseq.sh
```

Note:  
`27_210303_07` is a special dataset, where the modified version of `do_cleanFq.sh` in `27_210303_07` folder should be used to remove the 5’ primer sequences (19bp) in Read 2.

### Pooled data level
Pooled data means the data pooled from multiple datasets. Pooled data was used for the analysis at organ or stage level.

| Stage | Organ | Path | Scripts |
| :------ | :------ | :------ | :------ |
| 19-22-week-old | Single (e.g. stomach) | pooled_data_19-22w/01_stomach | do_poolData.r, do_filter.r, do_cluster.r |
| 19-22-week-old | Multiple (e.g. all 17 organs) | pooled_data_19-22w/All | do_poolData.r, do_cluster.r |
| 11-14-week-old | Single (e.g. stomach) | pooled_data_11-14w/01_stomach | do_poolData.r, do_filter.r, do_cluster.r |
| 11-14-week-old | Multiple (e.g. all 17 organs) | pooled_data_11-14w/All | do_poolData.r, do_cluster.r |
| All | Single (e.g. stomach) | pooled_data_all/01_stomach | do_poolData.r, do_cluster.r |
| All | Multiple (e.g. all 17 organs) | pooled_data_all/All | do_poolData.r, do_cluster.r |

Pipeline:
1. Cell pooling
```
Rscript do_poolData.r
```
2. Filtering (cells / genes) (only used at single-stage / single-organ level)
```
Rscript do_filter.r
```
3. Cell clustering and cell type annotation
```
Rscript do_cluster.r
```
