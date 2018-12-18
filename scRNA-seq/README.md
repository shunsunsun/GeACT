# scRNA-seq
A pipeline for single-cell RNA-seq data analysis, tuning for MALBAC-DT protocol.

## Prerequisite
1. perl (5.010 or later)  
2. R (3.0.1 or later)  

## Protocol
1. demultiplex  
Reads (R1) from each multiplexed sample were splited into 96 samples (cells) according to cell barcodes in reads (R2).  
Unique molecular identifier (UMI) were extracted and added into the corresponding name in reads (R1).  
PolyA sequences were trimmed and reads were filtered based on read length and quality.  
```
bash do_cleanFq.sh
```
2. alignment  
Reads (R1) were mapped to genome.
```
bash do_hisat2.sh
```
3. expression calling  
UMI count for each gene was calculated.
```
bash do_htseq.sh
```
