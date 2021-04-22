#!/bin/bash

# organ=Lung
# IN=/share/newdata4/lix/ATAC/datasets
# OUT=/share/newdata4/lix/ATAC/datasets/$organ/20191018_1025/narrow
# mkdir -p $OUT
# python output_insertcount_merge.py $IN $OUT $organ narrow

act=do_output_insertcount
cpu=1
mem=5g
sp=human

IN=/share/newdata4/lix/ATAC/datasets
OUT=/share/newdata4/lix/ATAC/datasets/All/20190611_20200117
mkdir -p $OUT

function do_main {
organ=$1
cells=`ls $IN/$organ/20*/cell_id_rename.txt`
mkdir -p $OUT/$organ/narrow_split

f1=$OUT/$organ/narrow_split/insertcount_All_NarrowPeaks_split_filter.tsv.gz
f2=$IN/All/nomodel_narrow/All_narrowPeak_split.bed
f3=$OUT/$organ/narrow_split/insertcount_All_NarrowPeaks_split_filter_coo.tsv.gz

cat << EOF > ${act}_${organ}.tmp.sh
#!/bin/bash
python output_insertcount_merge.py $IN $OUT/$organ/narrow_split All $cells narrow_split
python scripts/insertcount_tocoo.py $f1 $f2 $f3 False
EOF

chmod 750 ${act}_${organ}.tmp.sh

cat << EOF > run_${act}_${organ}.tmp.sh
#! /bin/bash

set -x
sbatch -J ${organ} \\
-D . \\
--export=ALL \\
-c $cpu \\
--mem=$mem \\
--nodelist=node08 \\
-o run_${act}_${organ}.o.txt \\
-e run_${act}_${organ}.e.txt \\
/share/newdata4/lix/ATAC/METATAC_pipeline/${act}_${organ}.tmp.sh
set +x
EOF

echo $id
bash run_${act}_${organ}.tmp.sh
}

function do_go {
for Organ in SI Stomach Kidney Lung Pan Eso Heart Liver BM Bladder Dia Spleen Thymus LI
do
    do_main $Organ
done
}
do_go
