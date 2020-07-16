#!/bin/bash
IN=02-alignment
IN2=../../Genomes
OUT=03-expression

cpu=1
mem=10g
account=gaog_g1
partition=cn-long       # fat4way cn-long fat4long
qos=$(grep -w "$partition" ~/users/tianf/.qos | cut -f 2)

act=`basename $0 .sh`

function do_main {
bpp=$1	# species
id=$2	# GSM

### sbatch
cat << EOF > ${act}_${bpp}_${id}.tmp.sh
#!/bin/bash
echo -ne "[Start]\t"; date
echo "Node: \$SLURM_NODELIST"
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

htseq-count -q -f bam $IN/$bpp/$id/${id}_sorted.bam $IN2/$bpp/genes.gtf -o $OUT/$bpp/$id/${id}_htseq.samout > $OUT/$bpp/$id/${id}_htseq.xls
rm -f $OUT/$bpp/$id/${id}_htseq.xls

perl scripts/htseq2tab.pl $OUT/$bpp/$id/${id}_htseq.samout 1 > $OUT/$bpp/$id/${id}_htseq.tab
rm -f $OUT/$bpp/$id/${id}_htseq.samout
perl scripts/rmdup.pl $OUT/$bpp/$id/${id}_htseq.tab 2 | cut -f 1 | sort | uniq -c | awk '{print \$2,\$1}' OFS="\t" > $OUT/$bpp/$id/${id}_UMIcount.txt
myjoin -F 1 -f 1 $IN2/$bpp/genes.txt $OUT/$bpp/$id/${id}_UMIcount.txt | cut -f 4 | sed 's/^$/0/' > $OUT/$bpp/$id/${id}_UMIcount_allGenes.txt

# stat
echo -ne "$id\t" > $OUT/$bpp/$id/stat.txt
perl scripts/tabStat.pl $OUT/$bpp/$id/${id}_htseq.tab $OUT/$bpp/$id/${id}_UMIcount_allGenes.txt atxt >> $OUT/$bpp/$id/stat.txt

echo "$OUT/$bpp/$id/${id}_UMIcount_allGenes.txt"
echo -ne "[End]\t"; date
echo "Done."
EOF

#exit	# XXX

rm -f ${act}_${bpp}_${id}.[oe].txt
sbatch -x "$(cat /tmp/node_BL.txt)" -p $partition -A $account --qos=$qos --cpus-per-task=$cpu -o ${act}_${bpp}_${id}.o.txt -e ${act}_${bpp}_${id}.e.txt ${act}_${bpp}_${id}.tmp.sh
}

function do_go {
#squeue -t R -h -p cn-long -o "%.8i %.9P %.100j %.15u %.2t %.10M %.6D %C %R" | grep -v 'gaog_pkuhpc' | awk '{print $9}' | sort -u | paste -s -d ',' > /tmp/node_BL.txt

function do_subgo {
Nto=1500
Ned=$(squeue -o "%.9i %.9P %.100j" -u gaog_pkuhpc -p $partition | grep $act | wc -l)
Nsm=$(expr $Nto - $Ned)
echo -e "> [Total] $Nto; [Running] $Ned; [Submit] $Nsm"

Fic=0
cat sample_list_cleandata.txt | cut -f 2,4 | while read Fsp Fid
do
	if [[ $Fic -ge $Nsm ]]
	then
		return 1
	else
		if [[ ! -e $OUT/$Fsp/$Fid/${Fid}_UMIcount_allGenes.txt ]] && [[ ! -e ${act}_${Fsp}_${Fid}.tmp.sh ]]
		then
			Fic=$(expr $Fic + 1)
			do_main $Fsp $Fid
			#exit
		else
			continue
			#echo "[skip] $Fsp $Fid"
		fi
	fi
done
}

while [[ 1 ]]; do sleep 0.5; do_subgo; if [[ $? -eq 0 ]]; then echo "> [All submitted]"; break; fi; done
}

# rescue errors
function do_rescue {
cat htseq_error.txt | while read Fsp Fid
do
	echo "$Fsp $Fid"
	do_main $Fsp $Fid
done
}

# merge sample
function do_merge {
echo "[info] Merge results..."
mkdir -p 03-expression/merged; cat <(echo -ne "Gene\t"; ls 03-expression/human/*/*_UMIcount_allGenes.txt | sort -V | cut -d '/' -f 4 | sed 's/_UMIcount_.*//' | paste -s) <(paste ../../Genomes/human/genes.txt `ls 03-expression/human/*/*_UMIcount_allGenes.txt | sort -V | paste -s -d ' '`) > 03-expression/merged/UMIcount_ensemblGene.txt
cat <(head -n 1 03-expression/merged/UMIcount_ensemblGene.txt) <(myjoin -m <(cat ../../Genomes/human/gene_ID2Name_fixed.txt ../../Genomes/human/ERCC_ID2Name.txt) <(tail -n +2 03-expression/merged/UMIcount_ensemblGene.txt) | cut -f 2,4-) > 03-expression/merged/UMIcount_allGenes.txt
# ERCC
head -n 1 03-expression/merged/UMIcount_allGenes.txt > 03-expression/merged/UMIcount_ERCC.txt
grep '^ERCC-[0-9][0-9][0-9][0-9]' 03-expression/merged/UMIcount_allGenes.txt >> 03-expression/merged/UMIcount_ERCC.txt
grep -v '^ERCC-[0-9][0-9][0-9][0-9]' 03-expression/merged/UMIcount_allGenes.txt > 03-expression/merged/UMIcount_allGenes.txts
mv 03-expression/merged/UMIcount_allGenes.txts 03-expression/merged/UMIcount_allGenes.txt

gzip -c 03-expression/merged/UMIcount_allGenes.txt > 03-expression/merged/UMIcount_allGenes.txt.gz
gzip -c 03-expression/merged/UMIcount_ERCC.txt > 03-expression/merged/UMIcount_ERCC.txt.gz

workdir=$(pwd -P | cut -d '/' -f 8); echo "Directory: $workdir"
targdir="tianf@162.105.250.222:/rd/user/tianf/06-Human_cell_atlas/datasets/$workdir/03-expression/merged/"

scp 03-expression/merged/UMIcount_allGenes.txt $targdir
scp 03-expression/merged/UMIcount_allGenes.txt.gz $targdir
scp 03-expression/merged/UMIcount_ERCC.txt $targdir
scp 03-expression/merged/UMIcount_ERCC.txt.gz $targdir

# stat
ls 03-expression/human/*/stat.txt | sort -V | while read fl; do sp=$(echo $fl | cut -d '/' -f 2); cat $fl | awk -v ks=$sp '{print ks,$0}' OFS="\t"; done > 03-expression/merged/exprStat.txt1
ls Logs/do_htseq/do_htseq_*.o.txt | sort -V | while read fl; do perl scripts/timeStat.pl $fl; done > 03-expression/merged/exprStat.txt2
myjoin -m -F 2 -f 2 03-expression/merged/exprStat.txt1 03-expression/merged/exprStat.txt2 | cut -f 1-12,15 > 03-expression/merged/exprStat.txt
rm -f 03-expression/merged/exprStat.txt[12]

scp 03-expression/merged/exprStat.txt $targdir
}

if [[ $1 = "" ]]
then
	do_go
elif [[ $1 = "rescue" ]]
then
	do_rescue
elif [[ $1 = "merge" ]]
then
	do_merge
fi
