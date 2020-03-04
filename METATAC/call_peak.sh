#!/bin/bash

prefix=Lung

IN=/share/newdata4/lix/ATAC/datasets/$prefix
mkdir -p $IN/nomodel_narrow
chrsz=/share/home/lix/references/human/hg38.chrom.sizes
blacklist=~/references/human/hg38.blacklist.bed.gz
gene_position=~/references/human/gene_position.hg38.v26.bed
genome=~/genomes/human/GRCh38/GRCh38.p10.genome.fa
peaks=$IN/nomodel_narrow/${prefix}_peaks_sorted_filt.narrowPeak

macs2 callpeak -t $IN/merge_bedpe.bed -f BEDPE -n ${prefix} --outdir $IN/nomodel_narrow --nomodel --nolambda --SPMR --keep-dup all
bedtools intersect -v -a $IN/nomodel_narrow/${prefix}_peaks.narrowPeak -b ${blacklist} \
|awk 'BEGIN{OFS="\t"}{if($5>1000) $5=1000;print $0}' > $IN/nomodel_narrow/${prefix}_peaks_filt.narrowPeak
sort -k 8gr,8gr $IN/nomodel_narrow/${prefix}_peaks_filt.narrowPeak |awk 'BEGIN{OFS="\t"}{$4="Peak-"NR;print $0}' > $peaks
rm $IN/nomodel_narrow/${prefix}_peaks_filt.narrowPeak
cat $gene_position |awk 'BEGIN{OFS="\t"}{$7=$2;$8=$3;print $0}' |bedtools slop -i - -g $chrsz -b 100000 \
|bedtools intersect -a - -b $peaks -wa -wb |awk 'BEGIN{OFS="\t"}{print $1,$7,$8,$4,$5,$6,$10,$11,$12}' \
|awk 'BEGIN{OFS="\t"}{m=($7+$8-1)/2}($6=="+"){$10=m-$2;$11=m-$3+1}($6=="-"){$10=$3-1-m;$11=$2-m}{print $0}' \
> $IN/nomodel_narrow/${prefix}_narrowPeak_gene_TSS100kb_TTS100kb_distance.bed

# fc_bedgraph=$IN/nomodel_narrow/${prefix}.fc.signal.bedgraph
# fc_bedgraph_srt=$IN/nomodel_narrow/${prefix}.fc.signal.sort.bedgraph
# fc_bigwig=$IN/nomodel_narrow/${prefix}.fc.signal.bigwig

# macs2 bdgcmp -t $IN/nomodel_narrow/${prefix}_treat_pileup.bdg -c $IN/nomodel_narrow/${prefix}_control_lambda.bdg --outdir $IN/nomodel_narrow --o-prefix ${prefix} -m FE
# slopBed -i $IN/nomodel_narrow/${prefix}_FE.bdg -g $chrsz -b 0 |bedClip stdin $chrsz $fc_bedgraph
# sort -k1,1 -k2,2n $fc_bedgraph > $fc_bedgraph_srt
# bedGraphToBigWig $fc_bedgraph_srt $chrsz $fc_bigwig
# rm -f $IN/nomodel_narrow/${prefix}_FE.bdg $fc_bedgraph $fc_bedgraph_srt

# pval_bedgraph=$IN/nomodel_narrow/${prefix}.pval.signal.bedgraph
# pval_bedgraph_srt=$IN/nomodel_narrow/${prefix}.pval.signal.srt.bedgraph
# pval_bigwig=$IN/nomodel_narrow/${prefix}_sig.pval.signal.bigwig

# sval=$(wc -l <(cat $IN/merge_bedpe.bed ) |awk '{printf "%f",$1/1000000}')
# macs2 bdgcmp -t $IN/nomodel_narrow/${prefix}_treat_pileup.bdg -c $IN/nomodel_narrow/${prefix}_control_lambda.bdg --outdir $IN/nomodel_narrow --o-prefix ${prefix} -m ppois -S "${sval}"
# slopBed -i $IN/nomodel_narrow/${prefix}_ppois.bdg -g $chrsz -b 0 |bedClip stdin $chrsz $pval_bedgraph
# rm -f $IN/nomodel_narrow/${prefix}_ppois.bdg
# sort -k1,1 -k2,2n $pval_bedgraph > $pval_bedgraph_srt
# bedGraphToBigWig $pval_bedgraph_srt $chrsz $pval_bigwig
# rm -f $IN/nomodel_narrow/${prefix}_ppois.bdg $pval_bedgraph $pval_bedgraph_srt
