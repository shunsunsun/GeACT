# File format

## cell_metatable_RNA_global.txt
| Index | Name | Description |
| :------| :------ | :------ |
| 1 | cell | Unique cell ID |
| 2 | seqID | Sequencing batches |
| 3 | species | Species |
| 4 | tissue | Organs or tissues |
| 5 | samplingPos | Sampling positions |
| 6 | plate | 96-well plates |
| 7 | QC | Whether the cell passed quality control |
| 8 | ident | Cell types |
| 9| group | Cell groups |
| 10 | ABratio | The ratio of primer A in primer A and B |
| 11 | cleanReads | The number of reads after quality control |
| 12 | mpRatio | The ratio of the reads mapped to the corresponding genome |
| 13 | nGene | The number of detected genes |
| 14 | nUMI | The number of detected UMIs (after removing duplicates) |
| 15 | mitoRatio | The ratio of mitochondrial genes (UMI count) |
| 16 | ERCCratio | The ratio of ERCC spike-in (UMI count) |
| 17 | doublet | Whether cells were doublets |
| 18-21 | `[tSNE/UMAP]_[1/2]` | The embeddings of the cells from one organ in one stage |
| 22 | stage | The development stage of corresponding cells |
| 23-26 | `[tSNE/UMAP]_[1/2]_ali` | The embeddings of the cells from one organ in all stages |
| 27-30 | `[tSNE/UMAP]_[1/2]_glo` | The embeddings for all cells from all organs in all stages |
