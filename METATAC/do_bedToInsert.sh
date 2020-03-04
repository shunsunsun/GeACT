#!/bin/bash

act=do_bedToInsert
cpu=1
mem=1g
sp=human
organ=All

tss=~/references/human/gene_TSS.hg38.v26.bed
chrom_size=~/references/human/hg38.chrom.sizes

narrowPeaks=/share/newdata4/lix/ATAC/datasets/${organ}/nomodel_narrow/${organ}_peaks_sorted_filt.narrowPeak
narrowPeaks_split=/share/newdata4/lix/ATAC/datasets/${organ}/nomodel_narrow/${organ}_narrowPeak_split.bed

function do_main {
IN=$1
batch=$2
id=$3

insert=$IN/$batch/${id}/${id}_${sp}_insert_decon.bed
output1=$IN/$batch/${id}/${id}_${sp}_TSS1kb.txt
output2=$IN/$batch/${id}/${id}_${sp}_TSS1kb_pos.txt
output3=$IN/$batch/${id}/${id}_${organ}_NarrowPeaks.txt
output4=$IN/$batch/${id}/${id}_${organ}_NarrowPeaks_split.txt
output6=$IN/$batch/${id}/${id}_${organ}_NarrowPeaks_pos.txt
output7=$IN/$batch/${id}/${id}_${organ}_NarrowPeaks_split_pos.txt

cat << EOF > ${act}_${id}.tmp.sh
#!/bin/bash

echo All\$'\\t'`cat $insert | wc -l` > $output1
cat $tss |grep ^chr[0-9X] |bedtools slop -i - -g ${chrom_size} -b 1000 |bedtools intersect -c -a - -b $insert |cut -f 4,7 >> $output1

cat $tss |grep ^chr[0-9X] |awk 'BEGIN{OFS="\\t"}{\$7=\$2}{\$8=\$3}{print \$0}' |bedtools slop -i - -g ${chrom_size} -b 1000 |bedtools intersect -a - -b $insert -wa -wb |awk 'BEGIN{OFS="\\t"}{print \$1,\$7,\$8,\$4,\$5,\$6,\$10,\$11}' > $output2

echo All\$'\\t'`cat $insert | wc -l` > $output3
bedtools intersect -c -a $narrowPeaks -b $insert |cut -f 4,11 >> $output3

echo All\$'\\t'`cat $insert | wc -l` > $output4
bedtools intersect -c -a $narrowPeaks_split -b $insert |cut -f 4,5 >> $output4

bedtools intersect -a $narrowPeaks -b $insert -wa -wb |awk 'BEGIN{OFS="\t"}{print \$4,\$12-\$2+1}' > $output6

bedtools intersect -a $narrowPeaks_split -b $insert -wa -wb |awk 'BEGIN{OFS="\t"}{print \$4,\$6-\$2+1}' > $output7

EOF

chmod 750 ${act}_${id}.tmp.sh

cat << EOF > run_${act}_${id}.tmp.sh
#! /bin/bash

set -x
sbatch -J ${act}_${id} \\
-D . \\
--export=ALL \\
-c $cpu \\
--mem=$mem \\
--nodelist=node06 \\
-o $IN/$batch/${id}/${act}_${id}.o.txt \\
-e $IN/$batch/${id}/${act}_${id}.e.txt \\
/share/newdata4/lix/ATAC/METATAC_pipeline/${act}_${id}.tmp.sh
set +x
EOF

echo $id
bash run_${act}_${id}.tmp.sh
}

function do_go {
for Date in 20190611 20190621 20190705 20190823 20191018 20191025 20191122 20200117
do
    In=/share/newdata4/lix/ATAC/datasets/$Date/02-mapping
	for ct in SIU SIM SIL WD WT WDOU RC RM SY LC LP Pan SD YXF YXS ZXF ZXS S4 S5 S7623 S8 GS AC DC TC PG GJ PZ XoX
	do
		ls $In |grep _${ct}[0-9]*_ |while read Batch
		do
			echo $Batch
			ls $In/$Batch |while read Id
			do
				do_main $In $Batch $Id
			done
			sleep 8s
		done
	done
done
}

do_go; exit
