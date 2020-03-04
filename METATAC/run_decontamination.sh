#!/bin/bash

IN=/share/newdata4/lix/ATAC/datasets/20191018/02-mapping
OUT=/share/newdata4/lix/ATAC/datasets/20191018/merge
sp=human
mkdir -p $OUT/con_freq
act=do_decon

function do_go {
ls $IN |python scripts/decon.py $IN $OUT $sp $OUT/con_freq/con_freq_${sp}.csv.gz

ls $IN |while read batch
do
    ls $IN/$batch |while read id
	do
	    frag_decon=$IN/$batch/$id/${id}_${sp}_frag_decon.bed
		insert_decon=$IN/$batch/$id/${id}_${sp}_insert_decon.bed
		cat ${frag_decon} |awk 'BEGIN{OFS="\t"}{print $1,$2-1,$2}{print $1,$3-1,$3}' |sort -V -u > ${insert_decon}
	done
done
}
do_go
