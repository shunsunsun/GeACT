#!/bin/bash
IN=01-cleandata
OUT=02-alignment
mkdir -p $OUT

cpu=2
mem=10g

act=`basename $0 .sh`

function do_main {
bpp=$1	# species
id=$2	# GSM

### sbatch
cat << EOF > ${act}_${bpp}_${id}.tmp.sh
#!/bin/bash
echo -ne ">:>\t"; date
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

hisat2 -p $cpu -x Genomes/$bpp/genome -U \$read --summary-file $OUT/$bpp/$id/${id}.stat | samtools view -@ $cpu -bh - | samtools sort -@ $cpu - $OUT/$bpp/$id/${id}_sorted

echo -ne ">:>\t"; date
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
cat sample_list_cleandata.txt | cut -f 2,4 | while read Fsp Fid
do
	if [[ ! -e $OUT/$Fsp/$Fid/${Fid}_sorted.bam ]] && [[ ! -e do_hisat2_${Fsp}_${Fid}.tmp.sh ]]
	then
		do_main $Fsp $Fid
		#exit
	else
		echo "[skip] $Fsp $Fid"
	fi
done
}

#do_go; exit

# rescue errors
function do_rescue {
#grep -l 'ERR' do_hisat2_*.e.txt | sed -e 's/^do_hisat2_//' -e 's/\.e\.txt//' -e 's/_/\t/' > hisat2_error.txt
cat hisat2_error.txt | while read Fsp Fid
do
	echo "$Fsp $Fid"
	do_main $Fsp $Fid
done
}

do_rescue; exit

# stat
mkdir 02-alignment/merged; find 02-alignment -name "*.stat" | sort -V | while read fl; do perl scripts/do_mapStat.pl $fl; done > 02-alignment/merged/mapStat.txt

