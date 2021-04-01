# Dataset
## Definition
A dataset consists of the data from one sequencing experiment.  
A dataset corresponds to several samples, each of which contain the multiplex data from one 96-well plate.  

## Naming rules
For the first five datasets, they were named as "order + (system abbreviation + organ number) + sample number" (e.g. 13_D1_48).  
As the system increased, other datasets were named as "order + (date) + sample number" (e.g. 18_190930_48).  

## Note
`27_210303_07` is a special dataset, where the 5’ primer sequences (19bp) need to be removed in Read 2 before the standard pipeline.  

## File format

### estimateCleanFq.txt

| Column | Description |
| ---------- | ----------- |
| 1 | Species |
| 2 | 96-well plate |
| 3 | Total reads number |
| 4 | Running time for spliting cell barcodes and reads quality control |
| 5 | Cell ID |
| 6 | Recognized read number (no mismatch in cell barcode) |
| 7 | Recognized read number (allow X mismtach in cell barcode) |
| 8 | Primer A number (less than 3 bases are inconsistent with "HBDV" pattern in UMIs) |
| 9 | Primer B number (less than 3 bases are inconsistent with "VDBH" pattern in UMIs) |
| 10 | Recognized UMI number (no mismatch) |
| 11 | Recognized UMI number (allow X mismatch) |
| 12 | Read number after quality control |
| 13 | Sequencing batch (optional column) |

### estimateMap.txt

| Column | Description |
| ---------- | ----------- |
| 1 | Species |
| 2 | Cell ID |
| 3 | Number of total reads |
| 4 | Number of mapped reads |
| 5 | Number of uniquely-mapped reads |
| 6 | Mapped ratio (extracted from Hisat2 log) |
| 7 | Mapped ratio (calculated) |
| 8 | Uniquely-mapped ratio |
| 9 | Strand (single or pair-end) |
| 10 | Running time for reads mapping |

### estimateExpr.txt

| Column | Description |
| ---------- | ----------- |
| 1 | Species |
| 2 | Cell ID |
| 3 | Number of uniquely-mapped reads (total) |
| 4 | Number of uniquely-mapped reads annotated as "no feature" (assigned to non-gene regions) |
| 5 | Number of uniquely-mapped reads annotated as "ambiguous" (assigned to multiple genes) |
| 6 | Number of uniquely-mapped reads (assigned to ERCC) |
| 7 | Number of uniquely-mapped (assigned to genes) |
| 8 | Number of non-redundant UMIs (total) |
| 9 | Number of non-redundant UMIs (ERCC) |
| 10 | Number of non-redundant UMIs (genes) |
| 11 | Number of detected ERCC |
| 12 | Number of detected gene |
| 13 | Running time for gene expression calculation |
