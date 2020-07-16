#!/bin/bash
IN=00-rawdata
OUT=00-rawdata/QC

cpu=20
mem=10g
account=gaog_g1
partition=cn-long	# fat4way cn-long fat4long
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

reads=\`ls $IN/$id/${id}_*.fastq.gz\`
echo "\$reads"

# check output
if [[ -e $OUT/${bpp}/${id}.fastq.gz ]]
then
	echo "File existed, exit."
	exit
else
	echo "Going..."
fi

fastqc -t $cpu -O $OUT/$bpp/$id --extract \$reads

echo -ne "[End]\t"; date
echo "Done."
EOF

#exit	# XXX

rm -f ${act}_${bpp}_${id}.[oe].txt
sbatch -x $(cat /tmp/node_BL.txt) -p $partition -A $account --qos=$qos --cpus-per-task=$cpu -o ${act}_${bpp}_${id}.o.txt -e ${act}_${bpp}_${id}.e.txt ${act}_${bpp}_${id}.tmp.sh
}

function do_go {
squeue -t R -h -p $partition -o "%.8i %.9P %.100j %.15u %.2t %.10M %.6D %C %R" | grep -v 'gaog_pkuhpc' | awk '{print $9}' | sort -u | paste -s -d ',' > /tmp/node_BL.txt

Fsp="human"
cut -f 2 sample_list_rawdata.txt | while read Fid
do
	if [[ ! -e ${act}_${Fsp}_${Fid}.tmp.sh ]] && [[ ! -e $OUT/$Fsp/$Fid ]]
	then
		echo "> $Fsp $Fid"
		do_main $Fsp $Fid
		#exit
	fi
done
}

function do_merge {
echo "[info] Merge results..."
multiqc -f -o 00-rawdata/QC/merged 00-rawdata/QC

targetIP="162.105.250.222"
targetPT="/rd/user/tianf/06-Human_cell_atlas/datasets"
targetDT=$(pwd -P | cut -d '/' -f 8)

scp 00-rawdata/QC/merged/multiqc_report.html tianf@$targetIP:$targetPT/$targetDT/00-rawdata/QC/merged/
}

if [[ $1 != "merge" ]]
then
	do_go
else
	do_merge
fi

