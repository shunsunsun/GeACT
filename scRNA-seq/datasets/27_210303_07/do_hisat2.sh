#!/bin/bash
IN=01-cleandata
IN2=../../Genomes
OUT=02-alignment

cpu=5
mem=10g
account=gaog_g1	# gaog_g1
partition=cn_icg	# fat4way cn-long fat4long
qos=$(grep -w "$partition" ~/users/tianf/.qos | cut -f 2)

act=`basename $0 .sh`

function do_main {
bpp=$1	# species
id=$2	# GSM

### sbatch
cat << EOF > ${act}_${bpp}_${id}.tmp.sh
#!/bin/bash
echo -ne "[Start]\t"; date
export PATH=/home/gaog_pkuhpc/users/tianf/tools:\$PATH

echo $bpp $id

if [[ -e $OUT/$bpp/$id/${id}_sorted.bam ]]
then
	echo "File existed, but not exit."
	#exit
else
	echo "Going..."
fi

mkdir -p $OUT/$bpp/$id

read=\`ls $IN/$bpp/*/${id}_read1.fastq.gz\`
echo "Input file: \$read"

hisat2 -p $cpu -x $IN2/$bpp/genome -U \$read --new-summary --summary-file $OUT/$bpp/$id/${id}.stat | samtools view -@ $cpu -bh - | samtools sort -@ $cpu - $OUT/$bpp/$id/${id}_sorted

# read distri
samtools view -H -@ $cpu $OUT/$bpp/$id/${id}_sorted.bam > $OUT/$bpp/$id/${id}_unique.sam
samtools view -@ $cpu $OUT/$bpp/$id/${id}_sorted.bam | grep -w 'NH:i:1' >> $OUT/$bpp/$id/${id}_unique.sam

srcdir="/home/gaog_pkuhpc/users/tianf/Tools/Python-3.6.4/install/bin"
python \$srcdir/read_distribution.py -i $OUT/$bpp/$id/${id}_unique.sam -r $IN2/$bpp/genes.bed12 > $OUT/$bpp/$id/readsDistri.txt
rm -f $OUT/$bpp/$id/${id}_unique.sam

echo -ne "[End]\t"; date
echo "$OUT/$bpp/$id/${id}_sorted.bam"
echo "Done"
EOF

#exit	# XXX

rm -f ${act}_${bpp}_${id}.[oe].txt
sbatch -x "$(cat /tmp/node_BL.txt)" -p $partition -A $account --qos=$qos --cpus-per-task=$cpu -o ${act}_${bpp}_${id}.o.txt -e ${act}_${bpp}_${id}.e.txt ${act}_${bpp}_${id}.tmp.sh
}

function do_go {
squeue -t R -h -p cn-long -o "%.8i %.9P %.100j %.15u %.2t %.10M %.6D %C %R" | grep -v 'gaog_pkuhpc' | awk '{print $9}' | sort -u | paste -s -d ',' > /tmp/node_BL.txt

if [[ ! -e sample_list_cleandata.txt ]]
then
	ls 01-cleandata/*/*/*.gz | sort -V | tr '/' '\t' | sed 's/_read1\.fastq\.gz//' > sample_list_cleandata.txt
fi

function do_subgo {
Nto=600
Ned=$(squeue -p $partition -o "%.9i %.9P %.100j" -u gaog_pkuhpc | grep $act | wc -l)
Nsm=$(expr $Nto - $Ned)
echo -e "> [Total] $Nto; [Running] $Ned; [Submit] $Nsm"

Fic=0
Fcd=0
while read Fsp Fid
do
	if [[ $Fic -ge $Nsm ]]
	then
		Fcd=1
		break
	fi

	if [[ ! -e $OUT/$Fsp/$Fid/${Fid}_sorted.bam ]] && [[ ! -e ${act}_${Fsp}_${Fid}.tmp.sh ]]
	then
		Fic=$(expr $Fic + 1)
		do_main $Fsp $Fid
		#exit
	else
		continue
		#echo "[skip] $Fsp $Fid"
	fi
done < <(cat sample_list_cleandata.txt | cut -f 2,4)
return $Fcd
}

while [[ 1 ]]; do sleep 2; do_subgo; if [[ $? -eq 0 ]]; then echo "> [All submitted]"; break; fi; done
}

# rescue errors
function do_rescue {
cat hisat2_error.txt | while read Fsp Fid
do
	echo "$Fsp $Fid"
	do_main $Fsp $Fid
done
}

function do_merge {
echo "[info] Merge results..."
# merge info
mkdir -p 02-alignment/merged; find 02-alignment -name "*.stat" | sort -V | while read fl; do perl scripts/mapStat.pl $fl; done > 02-alignment/merged/mapStat.txt1
ls Logs/do_hisat2/do_hisat2_*.o.txt | sort -V | while read fl; do perl scripts/timeStat.pl $fl; done > 02-alignment/merged/mapStat.txt2
myjoin -m -F 2 -f 2 02-alignment/merged/mapStat.txt1 02-alignment/merged/mapStat.txt2 | cut -f 1-9,12 > 02-alignment/merged/mapStat.txt
rm -f 02-alignment/merged/mapStat.txt[12]

find 02-alignment -name "readsDistri.txt" | grep -v 'merged' | sort -V | while read fl; do N=$(grep '^Total Tags' $fl | awk '{print $3}'); n=$(grep '^Introns' $fl | awk '{print $3}'); ratio=$(echo "scale=4; $n / $N" | bc | sed 's/^\./0./'); echo -e "$fl\t$N\t$n\t$ratio"; done | tr '/' '\t' | cut -f 2-3,5-7 > 02-alignment/merged/readsDistri.txt

targetIP="162.105.250.222"
targetPT="/rd/user/tianf/06-Human_cell_atlas/datasets"
targetDT=$(pwd -P | cut -d '/' -f 8)
for fl in 02-alignment/merged/mapStat.txt 02-alignment/merged/readsDistri.txt
do
	scp $fl tianf@$targetIP:$targetPT/$targetDT/02-alignment/merged/
done
}

if [[ $1 != "merge" ]]
then
	do_go
else
	do_merge
fi

