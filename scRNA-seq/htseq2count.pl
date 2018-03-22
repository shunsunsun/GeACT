#!/usr/bin/perl
use warnings;
use strict;
use 5.010;

if (@ARGV!=1) {
	print "Usage:   perl $0 samout\n";
	print "Example: perl $0 htout/human/PR10_XXL_501/1A_htseq.samout\n";
	exit;
}

#TTAATTAATGTACTCGCACT_0T_FFJJJJJJJJJJJJFJJJJJ_K00261:92:HGVHFBBXX:4:1117:9008:12603_N:0:26_ATCGAGTCGCTTGGGTGTAGTGCGTTTGG_AA-FFJJ<F-AFA---<A<--A-7FJF<-   16      chr1    629890  60      50M  *        0       0       GAACCATAACCAATACTACCAATCAATACTCATCATTAATAATCATAATG      JJJJJFJJJJJJJJJJJFJJJJJJJFJJJJFJJJJJJJJJJJJJFFAF<A      AS:i:-6 XN:i:0  XM:i:1  XO:i:0  XG:i:0  NM:i:1  MD:Z:16C33    YT:Z:UU NH:i:1  XF:Z:ENSG00000225630.1

my $UMI_len_def = 20;
my ($UMI, $gene);
open my $FILE , $ARGV[0];
while(<$FILE>) {
	chomp;
	if(/^(\S+).*XF:Z:(\S+)/) {
		$UMI = (split(/_/, $1))[0];
#		print STDERR "Warning: [$.] Inconsistence UMI length\n" if length($UMI) != $UMI_len_def;
		$gene = $2;
		next if $gene =~ /^__/;
	} else {
		print STDERR "Warning: [$.] No pattern found.\n";
	}
	print "$gene\t$UMI\n";
}
close $FILE;

