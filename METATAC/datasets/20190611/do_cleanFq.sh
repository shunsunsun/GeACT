#!/bin/bash

IN=../datasets/20190407/00-rawdata
OUT=../datasets/20190407/01-cleandata

mkdir -p $OUT

cpu=1
mem=1G

act=do_cleanFq

function do_main {
batch=$1

read1=`ls $IN/$batch/${batch}*_R1_001.fastq.gz`
read2=`ls $IN/$batch/${batch}*_R2_001.fastq.gz`

cat << EOF > ${act}_${batch}.tmp.sh
#!/bin/bash
echo -ne ">:>\t"; date
echo $batch
echo "$read1 $read2"

python scripts/cleanFq.py $read1 $read2 ../../experiments/indexToSample.txt ../../experiments/META16_Sequence.txt $OUT $batch

echo -ne ">:>\t"; date
echo "Done..."
EOF
chmod 750 ${act}_${batch}.tmp.sh

cat << EOF > run_${act}_${batch}.tmp.sh
#! /bin/bash

set -x
mkdir -p $OUT/$batch
sbatch -J  ${act}_${batch} \\
-D . \\
--export=ALL \\
-c $cpu \\
--nodelist=node03 \\
--mem=$mem \\
-o $OUT/$batch/${act}_${batch}.o.txt \\
-e $OUT/$batch/${act}_${batch}.e.txt \\
${act}_${batch}.tmp.sh
set +x
EOF

chmod 750 run_${act}_${batch}.tmp.sh
bash run_${act}_${batch}.tmp.sh

}

ls $IN |while read Batch
do
    if [[ ! -e $OUT/${Batch}/${act}_${Batch}.o.txt ]]
	then
		mkdir -p $OUT/$Batch
		echo $Batch
		do_main $Fbatch
	else
	    echo "[skip] $Batch"
	fi
done
