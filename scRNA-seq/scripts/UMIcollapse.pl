#!/usr/bin/perl
use warnings;
use strict;
use 5.010;

if (@ARGV!=2) {
	print "Usage:   perl $0 gene_UMI mismatch\n";
	print "Example: perl $0 test/test.geneUMI 1\n";
	exit;
}

sub cluster_UMI {
# method from http://www.arrayserver.com/wiki/index.php?title=Ngs_ReportSingleCellCounts.pdf
	my ($tmp, $cutoff) = @_;
	my @seqs = @{$tmp};
#print "@seqs\n";
	my %duplicate;
	$duplicate{$_}++ for @seqs;
	#map{print "$_=>$duplicate{$_}","\n"}keys %duplicate;

	my @sorted_seqs = sort { $duplicate{$b} <=> $duplicate{$a} or $a cmp $b } keys %duplicate;
	#foreach my $key (@sorted_seqs) { print "$key\t$duplicate{$key}\n"; }

	for(my $i=0; $i<@sorted_seqs; $i++) {
		next if $duplicate{$sorted_seqs[$i]} == 0;
		for(my $j=($i+1); $j<@sorted_seqs; $j++) {
			next if $duplicate{$sorted_seqs[$j]} == 0;
			my $dist = ($sorted_seqs[$i] ^ $sorted_seqs[$j]) =~ tr/\0//c;
#print "> $i ($sorted_seqs[$i]:$duplicate{$sorted_seqs[$i]})\t$j ($sorted_seqs[$j]:$duplicate{$sorted_seqs[$j]})\t$dist\n";
			if($dist <= $cutoff) {
				$duplicate{$sorted_seqs[$i]} += $duplicate{$sorted_seqs[$j]};
				$duplicate{$sorted_seqs[$j]} = 0;
			}
		}
	}
	#foreach my $key (@sorted_seqs) { print "$key\t$duplicate{$key}\n"; }
	my @outs = grep { $duplicate{$_} > 0 } (@sorted_seqs);
	return(scalar @outs);
}


#ENSG00000228463.9       ACAGTGGGTGTATGAATTTG

my $ctof = $ARGV[1];

my %gene;
open my $FILE , $ARGV[0];
while(<$FILE>) {
	chomp;
	my ($id, $seq) = split(/\t/, $_);
	if(! exists $gene{$id}) {
		$gene{$id} = $seq;
	} else {
		$gene{$id} .= ("-" . $seq);
	}
}
close $FILE;

foreach my $id (keys %gene) {	# for each gene
#print "> $id\n";
	my @seqs = split("-", $gene{$id});
	my @sequ = do { my %seen; grep { !$seen{$_}++ } @seqs };	# unique
	my $UMI_rmdup = &cluster_UMI(\@seqs, $ctof);

#print "> $id\t" . scalar @seqs . "\t->\t" . scalar @sequ . "\t->\t" . $UMI_rmdup . "\n";
	print "$id\t$UMI_rmdup\n";
}

