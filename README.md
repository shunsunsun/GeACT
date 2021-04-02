# GeACT

The first draft of Genome Architecture of Cells in Tissues (GeACT) of humans. GeACT is an extension of the Human Genome Project, offering a comprehensive database for documenting and understanding the functional human genome based on high-precision single-cell omics technologies.

## Data
version: v1.0.2
| Class | Technology | Cell number | Organ number |
| :------| :------ | ------: | ------: |
| sc-RNA-seq | MALBAC-DT | 38,671 | 17 |
| sc-ATAC-seq | METATAC | 45,643 | 16 |

To download raw data, please request it in [the Genome Sequence Archive](https://bigd.big.ac.cn/gsa-human/) with the ID of HRA000330.

### Transcriptome landscape
Gene expression matrix is accessible at http://geact.gao-lab.org/  
Cell meta table: [Download](https://github.com/gao-lab/GeACT/blob/master/scRNA-seq/pooled_data_all/All/cell_metatable_filtered_aligned.txt)

### Chromatin accessibility landscape
Gene activity score matrix is accessible at http://geact.gao-lab.org/  
Cell meta table: [Download](https://github.com/gao-lab/GeACT/blob/master/METATAC/pooled_data_all/All/cell_metatable_ATAC_global.txt)

## Pipeline
Transcriptome: https://github.com/gao-lab/GeACT/blob/master/scRNA-seq/README.md  
Chromatin accessibility landscape: NA

## Citation
(in preparation)

## Contact
geact@mail.cbi.pku.edu.cn
