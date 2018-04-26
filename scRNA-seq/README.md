# scRNA-seq
A pipeline for single-cell RNA-seq data analysis, tuning for MALBAC-DT protocol.

## Protocol
1. preparation
The sample list for raw data should be created for following analysis.
```
ls -1 00-rawdata/*/*.fastq.gz | tr '/' '\t' > sample_list_rawdata.txt
```
2. demultiplex
Reads (R1) from each multiplexed sample were splited into 96 samples (cells) according to cell barcodes in reads (R2). 
Unique molecular identifier (UMI) were extracted and added into the corresponding name in reads (R1). 
PolyA sequences were trimmed and reads were filtered based on read length and quality. 
```
bash do_cleanFq.sh
```
The sample list for clean data should be created for following analysis.
```
ls 02-cleandata/*/*/*.gz | tr '/' '\t' | sed 's/_read1\.fastq\.gz//' > sample_list_cleandata.txt
```
3. alignment
Reads (R1) were mapped to genome.
```
bash do_hisat2.sh
```
4. expression calling
UMI count for each gene was calculated.
```
bash do_htseq.sh
```
