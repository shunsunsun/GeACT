#!/bin/bash
IN=02-cleandata
OUT=STAR
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
export PATH=/home/gaog_pkuhpc/users/tianf/tools:\$PATH

echo $bpp $id

date
totid=\`echo $id | sed 's/_[0-9A-Z]*$//'\`
subid=\`echo $id | sed 's/.*_\([0-9A-Z]*\)$/\1/'\`

if [[ -e $OUT/$bpp/$id/sorted.bam ]]
then
	echo "File existed, but not exit."
	#exit
else
	echo "Going..."
fi

###
mkdir -p $OUT/$bpp/$id
rm -rf $OUT/$bpp/$id/STARtmp

STAR --runThreadN $cpu --genomeDir Genomes/$bpp/STARIndex --readFilesIn $IN/$bpp/\$totid/\${subid}_read1.fastq.gz --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $OUT/$bpp/$id/ --outTmpDir $OUT/$bpp/$id/STARtmp/ --outSAMunmapped Within

echo "Done"
echo "$OUT/$bpp/$id/sorted.bam"
date
###

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
cut -f 2 sample_list.txt | sort -u | while read Fsp
do
	cat barcode_ID.txt | while read Fbc
	do
		Fid=${Fsp}_${Fbc}
		if [[ ! -e $OUT/human/$Fsp/${Fbc}_sorted.bam ]] && [[ ! -e do_STAR_human_${Fid}.tmp.sh ]]
		then
			do_main human $Fid
			#exit
		else
			echo "[skip] $Fid"
		fi
	done
done
}

#do_go; exit
#do_main human PR10_XXL_501_1A; exit

# rescue errors
cat STAR_error.txt | while read Fsp Fid
do
	echo "$Fsp $Fid"
	do_main $Fsp $Fid
done

