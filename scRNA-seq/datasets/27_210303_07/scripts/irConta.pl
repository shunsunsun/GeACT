#!/usr/bin/perl
use strict;
use warnings;
use 5.010;

if(@ARGV!=1) {
	print "Remove the intra-cellular UMI contamination.\n";
	print "Usage:   perl $0 foo.tab\n";
	print "Example: perl $0 03-expression/human/SD-1_cell1/SD-1_cell1_htseq.tab\n";
	exit;
}

my $cutoff = 0.75;
my $rescue = 0;

# 1. read tab
my (%expr, %umi_stat);
open my $FILE, $ARGV[0];
while(<$FILE>) {
	chomp;
	next if /^__/;
	my ($gene, $umi, $num) = split(/\t/, $_);
	$expr{$umi}{$gene} = $num;
	$umi_stat{$gene} += $num;
}
close $FILE;

#my $umi_sum;
#my $umi_sum += $_ foreach values %umi_stat;
#print "Total umi number: $umi_sum\n";

# 2. remove contamination
my %conta_stat;
foreach my $umi (sort keys %expr) {
	#print "> $umi\n";
	next if keys %{$expr{$umi}} == 1;
	my $max_value = 0;
	my $max_gene;
	foreach my $gene (keys %{$expr{$umi}}) {
		#print ">> $umi\t$gene\t$expr{$umi}{$gene}\n";
		if($expr{$umi}{$gene} > $max_value) {
			$max_value = $expr{$umi}{$gene};
			$max_gene = $gene;
		}
	}
	my @max_gene_all;
	foreach my $gene (keys %{$expr{$umi}}) {
			push @max_gene_all, $gene if $expr{$umi}{$gene} == $max_value
	}
	#next if @max_gene_all > 1;
	foreach my $gene (keys %{$expr{$umi}}) {
		my $expr_ratio = $expr{$umi}{$gene} / $max_value;
		if( ($expr_ratio > 0) && ($expr_ratio < $cutoff) ) {
			$conta_stat{$gene} += $expr{$umi}{$gene};
			$expr{$umi}{$max_gene} += $expr{$umi}{$gene} if $rescue;
			$expr{$umi}{$gene} = 0;
		}
	}
}

#foreach my $gene (sort keys %conta_stat) {
#	print "$gene\t$umi_stat{$gene}\t$conta_stat{$gene}\t" . $conta_stat{$gene} / $umi_stat{$gene} . "\n";
#}

# 3. write tab
my @out_strs;
foreach my $umi (sort keys %expr) {
	foreach my $gene (sort keys %{$expr{$umi}}) {
		if($expr{$umi}{$gene} > 0) {
			push @out_strs, join("\t", $gene, $umi, $expr{$umi}{$gene})
		}
	}
}
print "$_\n" foreach (sort @out_strs);

