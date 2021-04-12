#!/usr/bin/perl
use warnings;
use strict;
use 5.010;

if (@ARGV!=2) {
	print "Usage:   perl $0 samout remove_nm\n";
	print "Example: perl $0 htout/human/PR10_XXL_501/1A_htseq.samout 0\n";
	exit;
}

#CCGCATGACTTCTTGAATAA_0_5T_A00455:12:H7KTCDSXX:4:1103:1768:6183_GATATGA_GTCGCTTGGGTGTAGTGC_GTTGTT_0

my $rmnm = $ARGV[1];
my $UMI_len_def = 20;
my ($UMI, $gene);
my %stat;
open my $FILE , $ARGV[0];
while(<$FILE>) {
	chomp;
	my @lines = split "\t";
	my $nm = (split "_", $lines[0])[7];
	if($rmnm) {
		next if $nm > 0;
	}
	if(/^([^_]+)_.*XF:Z:(\S+)/) {
		$UMI = $1;
#		print STDERR "Warning: [$.] Inconsistence UMI length\n" if length($UMI) != $UMI_len_def;
		$gene = $2;
		#next if $gene =~ /^__/;
	} else {
		print STDERR "Warning: [$.] No pattern found.\n";
	}
	#print "$gene\t$UMI\n";
	$stat{$gene . "\t" . $UMI} ++;
}
close $FILE;

foreach my $i (sort keys %stat) {
	print "$i\t$stat{$i}\n";
}

