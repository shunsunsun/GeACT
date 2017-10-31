#!/bin/bash
IN=00-rawdata
OUT=00-rawQC
mkdir -p $OUT

cpu=20
mem=10g

act=`basename $0 .sh`

function do_main {
bpp=$1	# species
id=$2

### sbatch
cat << EOF > ${act}_${bpp}_${id}.tmp.sh
#!/bin/bash
echo -ne ">:>\t"; date
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

echo -ne ">:>\t"; date
echo "Done..."
EOF

#exit	# XXX

partition=cn-long	# fat4way cn-long fat4long
account=gaog_g1
qos=gaogcnl	# gaogf4w gaogcnl gaogf4l

rm -f ${act}_${bpp}_${id}.[oe].txt
sbatch -x $(cat node_blacklist.txt) -p $partition -A $account --qos=$qos --cpus-per-task=$cpu -o ${act}_${bpp}_${id}.o.txt -e ${act}_${bpp}_${id}.e.txt ${act}_${bpp}_${id}.tmp.sh

##############
}

cut -f 2 sample_list_rawdata.txt | sort -u | while read Fid
do
	do_main human $Fid
done

# merge
multiqc -f -o 00-rawQC/merged 00-rawQC

