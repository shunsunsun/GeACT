#!/bin/bash

IN=/share/newdata4/lix/ATAC/datasets/20191018/01-cleandata_g
OUT=/share/newdata4/lix/ATAC/datasets/20191018/merge
mkdir -p $OUT

function do_go {
ls $IN |cut -d '_' -f 4 |uniq |while read ct
do
	mkdir -p ${OUT}/${ct}/ATAC
	output_file=${OUT}/${ct}/ATAC/cleanFqStat.txt
	echo "Batch	Total	Cell_id	Barcode	META" > $output_file
	ls $IN |grep ${ct} |while read Fsp
	do
		total=`grep '^Total' $IN/$Fsp/do_cleanFq_${Fsp}.o.txt |cut -d ' ' -f 2`
		grep ^${Fsp}_ $IN/$Fsp/do_cleanFq_${Fsp}.o.txt |awk 'BEGIN{OFS="\t"}{print "'$Fsp'\t'$total'",$0}' >> $output_file
	done
done
}
do_go; exit
