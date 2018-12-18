#!/bin/bash
IN=02-alignment
IN2=../Genomes
OUT=03-expression
mkdir -p $OUT

cpu=1
mem=10g
act=`basename $0 .sh`
partition=cn-long	# fat4way cn-long fat4long
account=gaog_g1
qos=gaogcnl	# gaogf4w gaogcnl gaogf4l

function do_main {
bpp=$1	# species
id=$2	# GSM

### sbatch
cat << EOF > ${act}_${bpp}_${id}.tmp.sh
#!/bin/bash
echo -ne "[Start]\t"; date
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

perl scripts/htseq2tab.pl $OUT/$bpp/$id/${id}_htseq.samout > $OUT/$bpp/$id/${id}_htseq.tab
rm -f $OUT/$bpp/$id/${id}_htseq.samout

perl scripts/UMIcollapse.pl $OUT/$bpp/$id/${id}_htseq.tab 4 > $OUT/$bpp/$id/${id}_UMIcount.txt
myjoin -F 1 -f 1 $IN2/$bpp/genes.txt $OUT/$bpp/$id/${id}_UMIcount.txt | cut -f 4 | sed 's/^$/0/' > $OUT/$bpp/$id/${id}_UMIcount_allGenes.txt

# stat
nd_tt=\$(cat $OUT/$bpp/$id/${id}_htseq.tab | grep -v -e '^__not_aligned' -e '^__alignment_not_unique' | wc -l)	# only count for uniquely-mapped reads
nd_nf=\$(cat $OUT/$bpp/$id/${id}_htseq.tab | grep '^__no_feature' | wc -l)
nd_ab=\$(cat $OUT/$bpp/$id/${id}_htseq.tab | grep '^__ambiguous' | wc -l)
nd_er=\$(cat $OUT/$bpp/$id/${id}_htseq.tab | grep -v '^__' | grep '^ERCC-' | wc -l)
nd_ge=\$(cat $OUT/$bpp/$id/${id}_htseq.tab | grep -v '^__' | grep '^ENSG' | wc -l)

nd_ch=\$(expr \$nd_nf + \$nd_ab + \$nd_er + \$nd_ge)
if [[ \$nd_tt -ne \$nd_ch ]]
then
	echo "Warning: inconsistent uniquely mapped reads: \$nd_tt (in tab) vs \$nd_ch (stat sum)." 1>&2
fi

nu_tt=\$(cat $OUT/$bpp/$id/${id}_UMIcount_allGenes.txt | awk '{sum+=\$1} END{print sum}')
nu_er=\$(paste $IN2/$bpp/genes.txt $OUT/$bpp/$id/${id}_UMIcount_allGenes.txt | grep '^ERCC-' | awk 'BEGIN{sum=0} {sum+=\$2} END{print sum}')
nu_ge=\$(paste $IN2/$bpp/genes.txt $OUT/$bpp/$id/${id}_UMIcount_allGenes.txt | grep -v '^ERCC-' | awk 'BEGIN{sum=0} {sum+=\$2} END{print sum}')

tp_er=\$(paste $IN2/$bpp/genes.txt $OUT/$bpp/$id/${id}_UMIcount_allGenes.txt | grep '^ERCC-' | awk '\$2!=0' | wc -l)
tp_ge=\$(paste $IN2/$bpp/genes.txt $OUT/$bpp/$id/${id}_UMIcount_allGenes.txt | grep -v '^ERCC-' | awk '\$2!=0' | wc -l)

echo -e "$id\t\$nd_tt\t\$nd_nf\t\$nd_ab\t\$nd_er\t\$nd_ge\t\$nu_tt\t\$nu_er\t\$nu_ge\t\$tp_er\t\$tp_ge" > $OUT/$bpp/$id/stat.txt

echo "$OUT/$bpp/$id/${id}_UMIcount_allGenes.txt"
echo -ne "[End]\t"; date
echo "Done."

EOF

#exit	# XXX

rm -f ${act}_${bpp}_${id}.[oe].txt
sbatch -x $(cat node_blacklist.txt) -p $partition -A $account --qos=$qos --cpus-per-task=$cpu -o ${act}_${bpp}_${id}.o.txt -e ${act}_${bpp}_${id}.e.txt ${act}_${bpp}_${id}.tmp.sh

##############
}

function do_go {

function do_subgo {
Nto=800
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
mkdir -p 03-expression/merged; cat <(echo -ne "Gene\t"; ls 03-expression/human/*/*_UMIcount_allGenes.txt | sort -V | cut -d '/' -f 4 | sed 's/_UMIcount_.*//' | paste -s) <(paste ../Genomes/human/genes.txt `ls 03-expression/human/*/*_UMIcount_allGenes.txt | sort -V | paste -s -d ' '`) > 03-expression/merged/UMIcount_ensemblGene.txt
cat <(head -n 1 03-expression/merged/UMIcount_ensemblGene.txt) <(myjoin <(myjoin <(cat ../Genomes/human/gene_ID2Name.txt ../Genomes/human/ERCC_ID2Name.txt) <(tail -n +2 03-expression/merged/UMIcount_ensemblGene.txt) | cut -f 3,5-) ../Genomes/human/gene_dupName.txt | grep -v '^=' | cut -f 2-) > 03-expression/merged/UMIcount_allGenes.txt
# ERCC
head -n 1 03-expression/merged/UMIcount_allGenes.txt > 03-expression/merged/UMIcount_ERCC.txt
grep '^ERCC-[0-9][0-9][0-9][0-9]' 03-expression/merged/UMIcount_allGenes.txt >> 03-expression/merged/UMIcount_ERCC.txt
grep -v '^ERCC-[0-9][0-9][0-9][0-9]' 03-expression/merged/UMIcount_allGenes.txt > 03-expression/merged/UMIcount_allGenes.txts
mv 03-expression/merged/UMIcount_allGenes.txts 03-expression/merged/UMIcount_allGenes.txt

gzip -c 03-expression/merged/UMIcount_allGenes.txt > 03-expression/merged/UMIcount_allGenes.txt.gz
gzip -c 03-expression/merged/UMIcount_ERCC.txt > 03-expression/merged/UMIcount_ERCC.txt.gz

workdir=$(pwd -P | cut -d '/' -f 7); echo "Directory: $workdir"
scp 03-expression/merged/UMIcount_allGenes.txt tianf@162.105.250.222:/rd/user/tianf/06-Human_cell_atlas/$workdir/03-expression/merged/
scp 03-expression/merged/UMIcount_allGenes.txt.gz tianf@162.105.250.222:/rd/user/tianf/06-Human_cell_atlas/$workdir/03-expression/merged/
scp 03-expression/merged/UMIcount_ERCC.txt tianf@162.105.250.222:/rd/user/tianf/06-Human_cell_atlas/$workdir/03-expression/merged/
scp 03-expression/merged/UMIcount_ERCC.txt.gz tianf@162.105.250.222:/rd/user/tianf/06-Human_cell_atlas/$workdir/03-expression/merged/

# stat
ls 03-expression/human/*/stat.txt | sort -V | while read fl; do sp=$(echo $fl | cut -d '/' -f 2); cat $fl | awk -v ks=$sp '{print ks,$0}' OFS="\t"; done > 03-expression/merged/exprStat.txt1
ls Logs/do_htseq/do_htseq_*.o.txt | sort -V | while read fl; do perl scripts/timeStat.pl $fl; done > 03-expression/merged/exprStat.txt2
myjoin -m -F 2 -f 2 03-expression/merged/exprStat.txt1 03-expression/merged/exprStat.txt2 | cut -f 1-12,15 > 03-expression/merged/exprStat.txt
rm -f 03-expression/merged/exprStat.txt[12]

scp 03-expression/merged/exprStat.txt tianf@162.105.250.222:/rd/user/tianf/06-Human_cell_atlas/$workdir/03-expression/merged/
}

if [[ $1 != "merge" ]]
then
	do_go
else
	do_merge
fi

