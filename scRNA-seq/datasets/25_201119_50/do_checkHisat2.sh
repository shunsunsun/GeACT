#!/bin/bash

cat do_hisat2_*.e.txt | grep -v -e '^HISAT2 summary' -e 'Total reads' -e 'Aligned 0 time' -e 'Aligned 1 time' -e 'Aligned >1 times' -e 'Overall alignment rate' -e 'bam_sort_core' -e '^processing' -e '^$' | wc -l
grep -l 'ERR' do_hisat2_*.e.txt | sed 's/\.e\.txt/.o.txt/' | while read fl; do grep -L 'Done' $fl; done | wc -l
grep -l 'CANCELLED' do_hisat2_*.e.txt | wc -l

grep -l -e 'ERR' -e 'CANCELLED' do_hisat2_*.e.txt | sed -e 's/^do_hisat2_//' -e 's/\.e\.txt//' -e 's/^\([a-zA-Z]*\)_/\1\t/' > hisat2_error.txt

wc -l hisat2_error.txt

echo "--------------"
grep -f <(cat hisat2_error.txt | tr '\t' '_' | sed 's/$/.tmp.sh/') <(squeue -u gaog_pkuhpc -o "%.7i %.9P %.120j %.8u %.2t %.10M") | awk '{print $1}' | paste -s -d ' '
echo "--------------"

if [[ $(cat hisat2_error.txt | wc -l) -eq 0 ]]
then
	rm -f hisat2_error.txt
fi

