#!/bin/bash

act=do_mapping
sp=human
IN=../datasets/20190407/01-cleandata
OUT=../datasets/20190407/02-mapping
mkdir -p $OUT
cpu=2
mem=8g

ref=/share/home/lix/genomes/human/GRCh38/index

function do_main {
batch=$1
id=$2

read1=$IN/$batch/${id}/${id}_R1.fastq.gz
read2=$IN/$batch/${id}/${id}_R2.fastq.gz
cut1=$OUT/$batch/${id}/${id}_cut_R1.fastq.gz
cut2=$OUT/$batch/${id}/${id}_cut_R2.fastq.gz

aln=$OUT/$batch/${id}/${id}_${sp}_aln.bam
dedup=$OUT/$batch/${id}/${id}_${sp}_dedup.bam
bedpe=$OUT/$batch/${id}/${id}_${sp}_bedpe.bed
frag=$OUT/$batch/${id}/${id}_${sp}_frag.bed
insert=$OUT/$batch/${id}/${id}_${sp}_insert.bed

cat << EOF > ${act}_${id}_${sp}.tmp.sh
#!/bin/bash

# Adapter trimming
cutadapt -e 0.22 -j 10 -a CTGTCTCTTATACACATCT -o - $read1 |cutadapt -e 0.22 -j 10 -g AGATGTGTATAAGAGACAG -o $cut1 -
cutadapt -e 0.22 -j 10 -a CTGTCTCTTATACACATCT -o - $read2 |cutadapt -e 0.22 -j 10 -g AGATGTGTATAAGAGACAG -o $cut2 -

# Mapping
bowtie2 -p $cpu -X 2000 --local --mm --no-discordant --no-mixed -x $ref -1 $cut1 -2 $cut2 |samtools view -bS -F 12 |samtools sort -n - > $aln

# Deduplicated
bamToBed -i $aln -bedpe |grep chr |awk 'BEGIN{OFS="\\t"}(\$8>=30){print \$0}' |sort -k 1,1 -k 2g,2 -k 6g,6 |python scripts/bedpe_process.py > $bedpe

# Generate fragments and insertions
cat $bedpe |grep -w chr[0-9X]* |awk 'BEGIN{OFS="\\t"}{print \$1,\$2+4,\$6-5,\$7,\$11}' > $frag
cat $frag |awk 'BEGIN{OFS="\\t"}{print \$1,\$2-1,\$2}{print \$1,\$3-1,\$3}' |sort -V -u > $insert

EOF

chmod 750 ${act}_${id}_${sp}.tmp.sh

cat << EOF > run_${act}_${id}_${sp}.tmp.sh
#! /bin/bash

set -x
sbatch -J ${act}_${id}_${sp} \\
-D . \\
--export=ALL \\
-c $cpu \\
--nodelist=node06 \\
--mem=$mem \\
-o $OUT/$batch/${id}/${act}_${id}_${sp}.o.txt \\
-e $OUT/$batch/${id}/${act}_${id}_${sp}.e.txt \\
${act}_${id}_${sp}.tmp.sh
set +x
EOF

echo $id
bash run_${act}_${id}_${sp}.tmp.sh
}

function do_go {
ls $IN |while read Batch
do
    echo $Batch
	ls $IN/$Batch |grep -v txt |while read Id
	do
		if [[ ! -e $OUT/$Batch/${Batch}/${Id}_${sp}_bedpe.bed ]]
		then
		    mkdir -p $OUT/$Batch/${Id}
			do_main $Batch $Id
		else
            echo "[skip] $Id"
        fi
	done
done
}

do_go; exit
