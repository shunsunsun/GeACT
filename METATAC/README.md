# scATAC-seq data analysis pipeline
A pipeline for single-cell ATAC-seq data analysis, maily for the data from METATAC.

## Prerequisite
1. Perl (5.16.3)  
2. R (3.5.1)  
3. Slurm (18.08.7)

## Pipeline
### Dataset level
Dataset means the data from one sequencing experiment (e.g. datasets/20190611).

Pipeline:
1. Demultiplexing and reads quality control  
```
bash do_cleanFq.sh
```
2. Read alignment  
Reads (R1) were mapped to the human genome and ERCC sequences.  
The genome assembly and gene annotation used in this step is the same as in scRNA-seq data analysis.
```
bash do_mapping.sh
```
3. Remove contanmination  
```
bash run_decontamination.sh
```

### Pooled data level
Pooled data means the data pooled from multiple datasets. Pooled data was used for the analysis at organ or stage level.

Pipeline:  
The fragment files were used as input for downstream analysis using [ArchR](https://www.archrproject.com/).  
Working directory: `new_pipeline`.
