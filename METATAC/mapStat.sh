#!/bin/bash

IN=/share/newdata4/lix/ATAC/datasets/20191018/02-mapping
OUT=/share/newdata4/lix/ATAC/datasets/20191018/merge
mkdir -p $OUT
sp=human

function do_go {
ls $IN |cut -d '_' -f 4 |uniq |while read ct
do
	output_file=$OUT/${ct}/ATAC/mapStat_${sp}.txt
	echo "Cell	Reads	Paired	Aligned_1_time	Aligned	Aligned_1_time_ratio	Aligned_ratio	Num_frags	Num_frags_decon	Num_frags_con	Num_mito_frags" > ${output_file}
	ls $IN |grep $ct |while read Fsp
	do
		ls $IN/$Fsp |while read Fid
		do
			aligned_info=`cat $IN/$Fsp/${Fid}/do_mapping_${Fid}_${sp}.e.txt |tail -n 12 |grep -v $'\t' |grep -v 'Warning' |python scripts/mapStat.py`
			bedpe=$IN/$Fsp/${Fid}/${Fid}_${sp}_bedpe.bed
			frag=$IN/$Fsp/${Fid}/${Fid}_${sp}_frag.bed
			frag_decon=$IN/$Fsp/${Fid}/${Fid}_${sp}_frag_decon.bed
			frag_con=$IN/$Fsp/${Fid}/${Fid}_${sp}_frag_con.bed
			num_frags=`cat $frag |wc -l`
			num_frags_decon=`cat $frag_decon |wc -l`
			num_frags_con=`cat $frag_con |wc -l`
			num_mito_frags=`cat $bedpe |grep -w chrM |wc -l`
			echo "${Fid}	${aligned_info}	${num_frags}	${num_frags_decon}	${num_frags_con}	${num_mito_frags}" >> ${output_file}
		done
	done
done
}
do_go; exit
