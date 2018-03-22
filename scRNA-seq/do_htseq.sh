#!/bin/bash
IN=02-alignment
IN2=Genomes
OUT=03-expression
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

if [[ -e $OUT/$bpp/$id/${id}_UMIcount_allGenes.txt ]]
then
	echo "File existed, but not exit."
	#exit
else
	echo "Going..."
fi

mkdir -p $OUT/$bpp/$id

samtools view -h -@ $cpu $IN/$bpp/$id/${id}_sorted.bam | grep -w 'NH:i:1' > $OUT/$bpp/$id/${id}_sorted.sam
cat $OUT/$bpp/$id/${id}_sorted.sam | wc -l
htseq-count -f sam $OUT/$bpp/$id/${id}_sorted.sam $IN2/$bpp/genes.gtf -o $OUT/$bpp/$id/${id}_htseq.samout > $OUT/$bpp/$id/${id}_htseq.xls
rm -f $OUT/$bpp/$id/${id}_sorted.sam

perl scripts/htseq2count.pl $OUT/$bpp/$id/${id}_htseq.samout | sort -u | cut -f 1 | sort | uniq -c | awk '{print \$2,\$1}' OFS="\t" > $OUT/$bpp/$id/${id}_UMIcount.txt
myjoin -F 1 -f 1 $IN2/$bpp/genes.txt $OUT/$bpp/$id/${id}_UMIcount.txt | cut -f 4 | sed 's/^$/0/' > $OUT/$bpp/$id/${id}_UMIcount_allGenes.txt

echo -ne ">:>\t"; date
echo "$OUT/$bpp/$id/${id}_UMIcount_allGenes.txt"
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
	if [[ ! -e $OUT/$Fsp/$Fid/${Fid}_UMIcount_allGenes.txt ]] && [[ ! -e do_htseq_${Fsp}_${Fid}.tmp.sh ]]
	then
		do_main $Fsp $Fid
		#exit
	else
		echo "[skip] $Fsp $Fid"
	fi
done
}

do_go; exit

# rescue errors
cat htseq_error.txt | while read Fsp Fid
do
	echo "$Fsp $Fid"
	do_main $Fsp $Fid
done

# merge samples
mkdir -p 03-expression/merged; cat <(echo -ne "Gene\t";ls 03-expression/human/*/*_UMIcount_allGenes.txt | sort -V | cut -d '/' -f 4 | sed 's/_UMIcount_.*//' | paste -s) <(paste Genomes/human/genes.txt `ls 03-expression/human/*/*_UMIcount_allGenes.txt | sort -V | paste -s -d ' '`) > 03-expression/merged/UMIcount_ensemblGene.txt
cat <(head -n 1 03-expression/merged/UMIcount_ensemblGene.txt) <(myjoin <(myjoin Genomes/human/gene_ID2Name.txt <(tail -n +2 03-expression/merged/UMIcount_ensemblGene.txt) | cut -f 3,5-) Genomes/human/gene_dupName.txt | grep -v '^=' | cut -f 2-) > 03-expression/merged/UMIcount_allGenes.txt
#scp 03-expression/merged/UMIcount_allGenes.txt tianf@202.205.131.15:/home/tianf/lustre/06-Human_cell_atlas/test_samples/03-expression/merged/

