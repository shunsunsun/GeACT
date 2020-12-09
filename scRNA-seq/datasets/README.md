# Dataset
## Definition
A dataset consists of the data from one sequencing experiment.  
A dataset corresponds to several samples, each of which contain the multiplex data from one 96-well plate.  

## Naming rules
For the first five datasets, they were named as "order + (system abbreviation + organ number) + sample number" (e.g. 13_D1_48).  
As the system increased, other datasets were named as "order + (date) + sample number" (e.g. 18_190930_48).

# File format

## estimateCleanFq.txt

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

## estimateMap.txt

| Column | Description |
| ---------- | ----------- |
| 1 | Cell ID |
| 2 | Number of total reads |
| 3 | Number of mapped reads |
| 4 | Number of uniquely-mapped reads |
| 5 | Mapped ratio (extracted from Hisat2 log) |
| 6 | Mapped ratio (calculated) |
| 7 | Uniquely-mapped ratio |
| 8 | Strand (single or pair-end) |

## estimateExpr.txt

| Column | Description |
| ---------- | ----------- |
| 1 | Number of uniquely-mapped reads (total) |
| 2 | Number of uniquely-mapped reads annotated as "no feature" (assigned to non-gene regions) |
| 3 | Number of uniquely-mapped reads annotated as "ambiguous" (assigned to multiple genes) |
| 4 | Number of uniquely-mapped reads (assigned to ERCC) |
| 5 | Number of uniquely-mapped (assigned to genes) |
| 6 | Number of non-redundant UMIs (total) |
| 7 | Number of non-redundant UMIs (ERCC) |
| 8 | Number of non-redundant UMIs (genes) |
| 9 | Number of detected ERCC |
| 10 | Number of detected gene |
