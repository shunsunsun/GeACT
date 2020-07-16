#!/usr/bin/perl
use warnings;
use strict;
use 5.010;

if (@ARGV!=2) {
	print "Remove the UMI duplication and adjust the UMI number.\n";
	print "Usage:   perl $0 foo.tab mismatch\n";
	print "Example: perl $0 03-expression/human/SD-1_cell1/SD-1_cell1_htseq.tab 4\n";
	exit;
}

# method from http://www.arrayserver.com/wiki/index.php?title=Ngs_ReportSingleCellCounts.pdf
my $cutoff = $ARGV[1];

# 1. read tab
my %expr;
open my $FILE, $ARGV[0];
while(<$FILE>) {
	chomp;
	next if /^__/;
	my ($gene, $umi, $num) = split(/\t/, $_);
	$expr{$gene}{$umi} = $num;
}
close $FILE;

# 2. remove duplication and adjust the UMI number
foreach my $gene (sort keys %expr) {	# for each gene
	#print "> $gene\n";
	next if keys %{$expr{$gene}} == 1;
	my @sorted_seqs = sort { $expr{$gene}{$b} <=> $expr{$gene}{$a} or $a cmp $b } keys %{$expr{$gene}};	# sort umi by num (reverse) and umi
	for(my $i=0; $i<@sorted_seqs; $i++) {
		next if $expr{$gene}{$sorted_seqs[$i]} == 0;
		for(my $j=($i+1); $j<@sorted_seqs; $j++) {
			next if $expr{$gene}{$sorted_seqs[$j]} == 0;
			my $dist = ($sorted_seqs[$i] ^ $sorted_seqs[$j]) =~ tr/\0//c;
			#print "[debug] $i ($sorted_seqs[$i]:$duplicate{$sorted_seqs[$i]})\t$j ($sorted_seqs[$j]:$duplicate{$sorted_seqs[$j]})\t$dist\n";
			if($dist <= $cutoff) {
				$expr{$gene}{$sorted_seqs[$i]} += $expr{$gene}{$sorted_seqs[$j]};
				$expr{$gene}{$sorted_seqs[$j]} = 0;
			}
		}
	}
}

# 3. write tab
foreach my $gene (sort keys %expr) {	# for each gene
	#print "> $gene\n";
	foreach my $umi (sort keys %{$expr{$gene}}) {	# for each umi
		next if $expr{$gene}{$umi} == 0;
		print "$gene\t$umi\t$expr{$gene}{$umi}\n";
	}
}

