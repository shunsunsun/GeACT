#!/usr/bin/perl
use warnings;
use strict;
use 5.010;

if (@ARGV!=2) {
	print "Usage:   perl $0 gene_UMI mismatch\n";
	print "Example: perl $0 test/test.geneUMI 1\n";
	exit;
}

#ENSG00000228463.9       ACAGTGGGTGTATGAATTTG

my %gene;
my $cutoff = $ARGV[1];

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
#print "$id\t" . scalar @seqs . "\t->\t" . scalar @sequ . "\n";
	my @outs;
	foreach my $seq (@sequ) {
#print ">> $seq\n";
		my $res = 1;
		foreach my $ref (@outs) {
#print "$seq -> $ref\n";
			my $dist = ($seq ^ $ref) =~ tr/\0//c;
#print "$seq\t$ref\t$dist\n";
			if($dist <= $cutoff) {
				$res = 0;
				last;
			}
		}
		if($res == 1) {
			push @outs, $seq;
		}
	}
#print "$id\t" . scalar @seqs . "\t->\t" . scalar @sequ . "\t->\t" . scalar @outs . "\n";
	print "$id\t" . scalar @outs . "\n";
}

