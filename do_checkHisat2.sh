#!/bin/bash

cat do_hisat2_*.e.txt | grep -v 'of these' | grep -v 'aligned' | grep -v 'overall' | grep -v 'bam_sort_core' | wc -l
grep -l 'ERR' do_hisat2_*.e.txt | sed 's/\.e\.txt/.o.txt/' | while read fl; do grep -L 'Done' $fl; done | wc -l
grep -l 'CANCELLED' do_hisat2_*.e.txt | wc -l

grep -l -e 'ERR' -e 'CANCELLED' do_hisat2_*.e.txt | sed -e 's/^do_hisat2_//' -e 's/\.e\.txt//' -e 's/^\([a-zA-Z]*\)_/\1\t/' > hisat2_error.txt

wc -l hisat2_error.txt

echo "--------------"
grep -f <(cat hisat2_error.txt | tr '\t' '_' | sed 's/$/.tmp.sh/') <(squeue -u gaog_pkuhpc -o "%.7i %.9P %.120j %.8u %.2t %.10M") | awk '{print $1}' | paste -s -d ' '
echo "--------------"

