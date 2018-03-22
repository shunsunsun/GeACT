#!/bin/bash
IN=00-rawdata
OUT=01-cleandata
mkdir -p $OUT

cpu=2
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

read1=\`ls $IN/$id/${id}_R1_001.fastq.gz\`
read2=\`ls $IN/$id/${id}_R2_001.fastq.gz\`
echo "\$read1 \$read2"

# check output
if [[ -e $OUT/${bpp}/${id}.fastq.gz ]]
then
	echo "File existed, exit."
	exit
else
	echo "Going..."
fi

awk -F "\t" '\$1=="$id"' $IN/barcode/barcode_all.txt | cut -f 2- > $OUT/$bpp/$id/barcodeToSample.txt
wc -l $OUT/$bpp/$id/barcodeToSample.txt

perl scripts/cleanFq.pl \$read1 \$read2 $OUT/$bpp/$id/barcodeToSample.txt $OUT/$bpp/$id 2 131

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

