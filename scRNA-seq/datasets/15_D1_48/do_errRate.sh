#!/bin/bash
IN=01-cleandata
OUT=02-errRate

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
export PATH=/home/gaog_pkuhpc/users/tianf/tools:\$PATH

echo $bpp $id

if [[ -e $OUT/$bpp/$id/${id}_errRate.txt ]]
then
	echo "File existed, but not exit."
	#exit
else
	echo "Going..."
fi

mkdir -p $OUT/$bpp/$id

read=\`ls $IN/$bpp/*/${id}_read1.fastq.gz\`
echo "Input file: \$read"

zcat \$read | grep '^@' | cut -d '_' -f 2 | awk '{sum+=\$1} END{print sum/NR/20*3}' > $OUT/$bpp/$id/${id}_errRate.txt

echo -ne "[End]\t"; date
echo "$OUT/$bpp/$id/${id}_errRate.txt"
echo "Done"
EOF

#exit	# XXX

rm -f ${act}_${bpp}_${id}.[oe].txt
sbatch -x $(cat /tmp/node_BL.txt) -p $partition -A $account --qos=$qos --cpus-per-task=$cpu -o ${act}_${bpp}_${id}.o.txt -e ${act}_${bpp}_${id}.e.txt ${act}_${bpp}_${id}.tmp.sh
}

function do_go {
squeue -t R -h -p cn-long -o "%.8i %.9P %.100j %.15u %.2t %.10M %.6D %C %R" | grep -v 'gaog_pkuhpc' | awk '{print $9}' | sort -u | paste -s -d ',' > /tmp/node_BL.txt

function do_subgo {
Nto=1200
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

	if [[ ! -e $OUT/$Fsp/$Fid/${Fid}_errRate.txt ]] && [[ ! -e ${act}_${Fsp}_${Fid}.tmp.sh ]]
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
cat 02-errRate/human/*/*_errRate.txt > 02-errRate/merged/errRate.txt
paste <(ls 02-errRate/human/*/*_errRate.txt | cut -d '/' -f 3) 02-errRate/merged/errRate.txt > 02-errRate/merged/errRate.txts
mv 02-errRate/merged/errRate.txts 02-errRate/merged/errRate.txt

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

