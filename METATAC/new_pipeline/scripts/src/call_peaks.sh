#!/usr/bin/bash
frag_file=""
blacklists=()
outdir=unset
force=0

usage()
{
	echo "Usage: call_peaks [ -f | --fragments FRAGMENTFILE ]
	[ -o | --outdir   OUTDIR ]
	[ -b | --blacklists BLACKLIST1 (can pass more than once) ]
	[ --force ]"
	exit 2
}

PARSED_ARGUMENTS=$(getopt -n call_peaks -o f:o:b: --long force,fragments:,outdir:,blacklists: -- "$@")
VALID_ARGUMENTS=$?
if [ "$VALID_ARGUMENTS" != "0" ]; then
	usage
fi

echo "PARSED_ARGUMENTS is $PARSED_ARGUMENTS"
eval set -- "$PARSED_ARGUMENTS"

while :
do
  case "$1" in
    -f | --fragments)   frag_file="$2"      ; shift 2 ;;
    -o | --outdir)    outdir="$2"       ; shift 2 ;;
    -b | --blacklists) blacklists+=("$2") ; shift 2 ;;
	--force) force=1 ; shift ;; 
    --) shift; break ;;
    *) echo "Unexpected option: $1 - this should not happen."
       usage ;;
  esac
done

blacklists+=$@

if [ ! -e $outdir ]; then
	mkdir $outdir
elif [ -f $outdir ]; then
	rm $outdir
	mkdir $outdir
else
	if [ $force -eq 1 ]; then
		rm -rf ${outdir}/*
	fi
fi

echo "Constructing insertion bed from fragments"
zcat $frag_file | awk 'BEGIN{FS="\t";OFS="\t"}{print $1, $2, $2+1; print $1, $3-1, $3}' > ${outdir}/insertion.bed
echo "Calling Peaks"
macs2 callpeak -t ${outdir}/insertion.bed -f BED -n organ --outdir ${outdir}/macs2  --nolambda --nomodel --keep-dup all -g hs -q 0.05 --extsize 150 --shift -75 --call-summits
echo "Filtering blacklists"
bedtools intersect -v -a ${outdir}/macs2/organ_peaks.narrowPeak -b ${blacklists[@]} > ${outdir}/macs2/organ_peaks_filtered.narrowPeak

