#!/bin/bash
IN=01-cleandata
OUT=02-alignment
mkdir -p $OUT

cpu=4
mem=10g

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

hisat2 -p $cpu -x ../Genomes/$bpp/genome -U \$read --new-summary --summary-file $OUT/$bpp/$id/${id}.stat | samtools view -@ $cpu -bh - | samtools sort -@ $cpu - $OUT/$bpp/$id/${id}_sorted

echo -ne "[End]\t"; date
echo "$OUT/$bpp/$id/${id}_sorted.bam"
echo "Done"

EOF

#exit	# XXX

partition=cn-long	# fat4way cn-long fat4long
account=gaog_g1
qos=gaogcnl	# gaogf4w gaogcnl gaogf4l

rm -f ${act}_${bpp}_${id}.[oe].txt
sbatch -x $(cat node_blacklist.txt) -p $partition -A $account --qos=$qos --cpus-per-task=$cpu -o ${act}_${bpp}_${id}.o.txt -e ${act}_${bpp}_${id}.e.txt ${act}_${bpp}_${id}.tmp.sh

##############
}

function do_go {
### create clean data table
if [[ ! -e sample_list_cleandata.txt ]]
then
	ls 01-cleandata/*/*/*.gz | sort -V | tr '/' '\t' | sed 's/_read1\.fastq\.gz//' > sample_list_cleandata.txt
fi
###

function do_subgo {
Nto=500
Ned=$(squeue -o "%.9i %.9P %.100j" -u gaog_pkuhpc | grep $act | wc -l)
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

#do_go; exit
while [[ 1 ]]; do sleep 2; do_subgo; if [[ $? -eq 0 ]]; then echo "> [All submitted]"; break; fi; done; exit

# rescue errors
function do_rescue {
cat hisat2_error.txt | while read Fsp Fid
do
	echo "$Fsp $Fid"
	do_main $Fsp $Fid
done
}

#do_rescue; exit
}

function do_merge {
echo "[info] Merge results..."
# merge info
mkdir -p 02-alignment/merged; find 02-alignment -name "*.stat" | sort -V | while read fl; do perl scripts/mapStat.pl $fl; done > 02-alignment/merged/mapStat.txt1
ls Logs/do_hisat2/do_hisat2_*.o.txt | sort -V | while read fl; do perl scripts/timeStat.pl $fl; done > 02-alignment/merged/mapStat.txt2
myjoin -m -F 2 -f 2 02-alignment/merged/mapStat.txt1 02-alignment/merged/mapStat.txt2 | cut -f 1-9,12 > 02-alignment/merged/mapStat.txt
rm -f 02-alignment/merged/mapStat.txt[12]

workdir=$(pwd -P | cut -d '/' -f 7); echo "Directory: $workdir"
scp 02-alignment/merged/mapStat.txt tianf@162.105.250.222:/rd/user/tianf/06-Human_cell_atlas/$workdir/02-alignment/merged/
}

if [[ $1 != "merge" ]]
then
	do_go
else
	do_merge
fi

