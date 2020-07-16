#!/usr/bin/perl
use strict;
use warnings;
use 5.010;

if(@ARGV != 3) {
	print "Stat for the tab file.\n";
	print "Usage:   perl $0 foo.tab expr.txt [txt|atxt]\n";
	print "Example: perl $0 03-expression/human/293_B9-1_01/293_B9-1_01_htseq.tab0 03-expression/human/293_B9-1_01/293_B9-1_01_UMIcount.txt0 txt\n";
	exit;
}

my %stat = ("tt" => 0, "nf" => 0, "ab" => 0, "er" => 0, "ge" => 0);
open my $FILE1, $ARGV[0];
while(<$FILE1>) {
	chomp;
	my ($gene, $umi, $num) = split(/\t/, $_);

	if( ($gene !~ /^__not_aligned/) && ($gene !~ /^__alignment_not_unique/) ) {
		$stat{"tt"} += $num;
	}

	if($gene =~ /^__no_feature/) {
		$stat{"nf"} += $num;
	} elsif($gene =~ /^__ambiguous/) {
		$stat{"ab"} += $num;
	} elsif($gene =~ /^ERCC-/) {
		$stat{"er"} += $num;
	} elsif($gene =~ /^ENSG/) {
		$stat{"ge"} += $num;
	}
	
}
close $FILE1;

$stat{"ch"} = $stat{"nf"} + $stat{"ab"} + $stat{"er"} + $stat{"ge"};
if($stat{"tt"} != $stat{"ch"}) {
	print STDERR "Warning: inconsistent uniquely mapped reads: $stat{\"tt\"} (in tab) vs  $stat{\"ch\"} (stat sum).\n";
}

my %stat_count = ("tt" => 0, "er" => 0, "ge" => 0, "et" => 0, "gt" => 0);
if($ARGV[2] eq "txt") {
open my $FILE2, $ARGV[1];
	while(<$FILE2>) {
		chomp;
		my ($gene, $num) = split(/\t/, $_);
		next if $num == 0;
		$stat_count{"tt"} += $num;
		if($gene =~ /^ERCC-/) {
			$stat_count{"er"} += $num;
			$stat_count{"et"} += 1;
		} else {
			$stat_count{"ge"} += $num;
			$stat_count{"gt"} += 1;
		}
	}
close $FILE2;
} elsif($ARGV[2] eq "atxt") {
	open my $GENES, "../../Genomes/human/genes.txt" or die;
	my @genes;
	while(<$GENES>) {
		chomp;
		push @genes, $_;
	}
	close $GENES;

	my $i = 0;
	open my $FILE2, $ARGV[1];
	while(<$FILE2>) {
		chomp;
		my $num = $_;
		my $gene = $genes[$i];
		if($num == 0) {
			$i += 1;
			next;
		}
		$stat_count{"tt"} += $num;
		if($gene =~ /^ERCC-/) {
			$stat_count{"er"} += $num;
			$stat_count{"et"} += 1;
		} else {
			$stat_count{"ge"} += $num;
			$stat_count{"gt"} += 1;
		}
		$i += 1;
	}
	close $FILE2;
}
print join ("\t", $stat{"tt"}, $stat{"nf"}, $stat{"ab"}, $stat{"er"}, $stat{"ge"}, $stat_count{"tt"}, $stat_count{"er"}, $stat_count{"ge"}, $stat_count{"et"}, $stat_count{"gt"}) . "\n";

