#!/bin/bash
IN="00-rawdata"
IN2="../../experiments"
OUT="01-cleandata"

cpu=20
mem=10g
account=gaog_g1
partition=cn-long       # fat4way cn-long fat4long
qos=$(grep -w "$partition" ~/users/tianf/.qos | cut -f 2)

act=`basename $0 .sh`

function do_main {
bpp=$1	# species
id=$2

### sbatch
cat << EOF > ${act}_${bpp}_${id}.tmp.sh
#!/bin/bash
echo -ne "[Start]\t"; date
export PATH=/home/gaog_pkuhpc/users/tianf/tools:\$PATH

echo $bpp $id
mkdir -p $OUT/$bpp/$id

read1=\`ls $IN/$id/${id}_*_R1_001.fastq.gz\`
read2=\`ls $IN/$id/${id}_*_R2_001.fastq.gz\`
echo "\$read1 \$read2"

# check output
if [[ -e $OUT/${bpp}/${id}.fastq.gz ]]
then
	echo "File existed, exit."
	exit
else
	echo "Going..."
fi

awk -F "\t" '{\$3="${id}_"\$3; print \$0}' OFS="\t" $IN2/cell_barcodes_0x96.txt | sed -e 's/P[A-Z]10_HCA_//' > $OUT/$bpp/$id/barcodeToSample.txt

perl scripts/cleanFq.pl \$read1 \$read2 $OUT/$bpp/$id/barcodeToSample.txt $OUT/$bpp/$id 0 151 > $OUT/$bpp/$id/info.txt

echo -ne "[End]\t"; date
echo "Done..."
EOF

#exit	# XXX

rm -f ${act}_${bpp}_${id}.[oe].txt
sbatch -x $(cat /tmp/node_BL.txt) -p $partition -A $account --qos=$qos --cpus-per-task=$cpu -o ${act}_${bpp}_${id}.o.txt -e ${act}_${bpp}_${id}.e.txt ${act}_${bpp}_${id}.tmp.sh
}

function do_go {
squeue -t R -h -p $partition -o "%.8i %.9P %.100j %.15u %.2t %.10M %.6D %C %R" | grep -v 'gaog_pkuhpc' | awk '{print $9}' | sort -u | paste -s -d ',' > /tmp/node_BL.txt

if [[ ! -e sample_list_rawdata.txt ]]
then
	ls -1 00-rawdata/*/*.fastq.gz | awk -F "/" '{idx=$1"\t"$2}{hs[idx]=(idx in hs)?hs[idx]"\t"$3:$3} END{for(i in hs) print i,hs[i]}' OFS="\t" | sort -V > sample_list_rawdata.txt
fi

Fsp="human"
cut -f 2 sample_list_rawdata.txt | while read Fid
do
	if [[ ! -e ${act}_${Fsp}_${Fid}.tmp.sh ]] && [[ $(ls $OUT/$Fsp/$Fid/*.fastq.gz 2>/dev/null | wc -l) -eq 0 ]]
	then
		echo "> $Fsp $Fid"
		do_main $Fsp $Fid
		#exit
	fi
done
}

function do_merge {
echo "[info] Merge results..."
# merge info
mkdir -p 01-cleandata/merged; find 01-cleandata -name "info.txt" | sort -V | while read fl; do id=$(echo $fl | cut -d '/' -f 3); tail -n +8 $fl | head -n -1 | awk -v ks=$id '{print ks,$0}' OFS="\t"; done > 01-cleandata/merged/cleanFqStat.txt0
find 01-cleandata -name "info.txt" | sort -V | while read fl; do sp=$(echo $fl | cut -d '/' -f 2); id=$(echo $fl | cut -d '/' -f 3); head -n 2 $fl | tail -n 1 | awk -v ks=$sp -v ki=$id '{print ks,ki,$0}' OFS="\t"; done > 01-cleandata/merged/rawFqStat.txt1
ls Logs/do_cleanFq/do_cleanFq_*.o.txt | sort -V | while read fl; do perl scripts/timeStat.pl $fl; done > 01-cleandata/merged/rawFqStat.txt2
myjoin -m -F 2 -f 2 01-cleandata/merged/rawFqStat.txt1 01-cleandata/merged/rawFqStat.txt2 | cut -f 1,2,4,7 > 01-cleandata/merged/rawFqStat.txt
rm -f 01-cleandata/merged/rawFqStat.txt[12]
myjoin -m -F 2 -f 1 01-cleandata/merged/rawFqStat.txt 01-cleandata/merged/cleanFqStat.txt0 | cut -f 1-4,6-13 > 01-cleandata/merged/cleanFqStat.txt
rm -f 01-cleandata/merged/rawFqStat.txt 01-cleandata/merged/cleanFqStat.txt0

targetIP="162.105.250.222"
targetPT="/rd/user/tianf/06-Human_cell_atlas/datasets"
targetDT=$(pwd -P | cut -d '/' -f 8)
for fl in 01-cleandata/merged/cleanFqStat.txt
do
	scp $fl tianf@$targetIP:$targetPT/$targetDT/01-cleandata/merged/
done
}

if [[ $1 != "merge" ]]
then
	do_go
else
	do_merge
fi

