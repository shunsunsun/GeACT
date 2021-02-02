#! /bin/bash

OUT=all_organ.bed.gz

# clear output
:>$OUT

for organ_bed in $(ls 19-21w/*/*.bed.gz)
do
	echo $organ_bed
	cat $organ_bed >> $OUT
done
