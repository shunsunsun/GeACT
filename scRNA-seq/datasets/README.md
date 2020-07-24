# File format

## estimateCleanFq.txt

| Column | Description |
| ---------- | ----------- |
| 1 | Cell ID |
| 2 | Recognized read number (no mismatch in cell barcode) |
| 3 | Recognized read number (allow X mismtach in cell barcode) |
| 4 | Primer A number (less than 3 bases are inconsistent with "HBDV" pattern in UMIs) |
| 5 | Primer B number (less than 3 bases are inconsistent with "VDBH" pattern in UMIs) |
| 6 | Recognized UMI number (no mismatch) |
| 7 | Recognized UMI number (allow X mismatch) |
| 8 | Read number after quality control |

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
