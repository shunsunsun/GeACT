#!/usr/bin/bash
frag_file=unset
outdir=unset
root=unset
force=0
bySummits=0
peakSetName=default

usage()
{
	echo "Usage: call_peaks [ -f | --fragments FRAGMENTFILE ]
	[ -o | --outdir   OUTDIR ]
	[ -r | --root ROOT ]
	[ -s | --bySummits ]
	[ -n | --peakSetName ]
	[ --force ]"
	exit 2
}

PARSED_ARGUMENTS=$(getopt -n call_peaks -o sf:o:r:n: --long bySummits,force,fragments:,outdir:,root:,peakSetName: -- "$@")
VALID_ARGUMENTS=$?
if [ "$VALID_ARGUMENTS" != "0" ]; then
	usage
fi

echo "PARSED_ARGUMENTS is $PARSED_ARGUMENTS"
eval set -- "$PARSED_ARGUMENTS"

while :
do
  case "$1" in
    -f | --fragments)   frag_file="$2";     shift 2;;
    -o | --outdir)    outdir="$2";          shift 2;;
    -r | --root) root="$2";                 shift 2;;
    -n | --peakSetName) peakSetName="$2";   shift 2;;
    -s | --bySummits) bySummits=1;          shift;;
	  --force) force=1;                       shift;; 
    --) shift;                              break;;
    *) echo "Unexpected option: $1 - this should not happen."
       usage ;;
  esac
done

mkdir -p ${outdir}
if [ ! -f ${outdir}/insertion.bed ] || [ $force -eq 1 ]; then
    echo "Constructing insertion bed from fragments"
    zcat $frag_file | awk 'BEGIN{FS="\t";OFS="\t"}{print $1, $2, $2+1; print $1, $3-1, $3}' > ${outdir}/insertion.bed
fi

if [ ! -f ${outdir}/macs2/organ_peaks.narrowPeak ] || [ $force -eq 1 ]; then
    echo "Calling Peaks"
    macs2 callpeak -t ${outdir}/insertion.bed -f BED -n organ --outdir ${outdir}/macs2  --nolambda --nomodel --keep-dup all -g hs -q 0.05 --extsize 150 --shift -75 --call-summits
fi

# echo "Filtering blacklists"
# bedtools intersect -v -a ${outdir}/macs2/organ_peaks.narrowPeak -b ${blacklists[@]} > ${outdir}/macs2/organ_peaks_filtered.narrowPeak

echo "Constructing peakSet(filtering and annotation)"
if [ $bySummits -eq 1 ]; then
  ${root}/scripts/src/construct_peakSet.R --inputdir ${outdir}/macs2 --output ${outdir}/summitPeaks --annodir ${root}/database/annotation --fromSummits --source $peakSetName
else
  ${root}/scripts/src/construct_peakSet.R --inputdir ${outdir}/macs2 --output ${outdir}/normalPeaks --annodir ${root}/database/annotation --source $peakSetName
fi
