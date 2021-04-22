#!/bin/bash

IN=/share/newdata4/lix/ATAC/datasets
OUT=$IN/All
sp=human
output_file=$OUT/merge_bedpe.bed

if [[ -e $output_file ]]
then
    rm -f $output_file
fi

for organ in SI Stomach Kidney Lung Pan Eso Heart Liver BM Bladder Dia Spleen Thymus LI
do
    cells=`ls $IN/$organ/20*/merge/kept_cell.txt`
	sed 's/^M//g' $cells |while read id
	do
		cat $IN/*/02-mapping/*/${id}/${id}_${sp}_frag_decon.bed |cut -f 1-3 >> $output_file
	done
done
