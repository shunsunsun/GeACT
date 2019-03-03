#!/usr/bin/perl
use strict;
use warnings;
use 5.010;
use POSIX;

if(@ARGV!=2) {
	print "Remove the inter-cellular UMI contamination.\n";
	print "Usage:   perl $0 plate outpath\n";
	print "Example: perl $0 SD-1 03-expression/itConta\n";
	exit;
}

my $plate = $ARGV[0];
my $outpt = $ARGV[1];
my $cutoff = 0.75;
my $ncpu = 20;
my $rescue = 0;
my $verbose = 0;
my $writeTab = 1;
my $writeStat = 1;
my $writeExpr = 1;
my $writeMergedExpr = 1;

# 1. read tab files
my @files = glob("03-expression/human/" . $plate . '_cell*/*_htseq.tab');
@files = sort { $a=~s/.*_cell(\d+)_htseq.tab/$1/r <=> $b=~s/.*_cell(\d+)_htseq.tab/$1/r } @files;
#print "@files\n";
my @cells;
foreach my $file (@files) { push @cells, (split("/", $file))[2] };

my (%expr, %umi_stat, %umt_stat);
foreach my $file (@files) {
	my $cell = (split("/", $file))[2];
#print "> $cell\n";
	open my ($FILE), $file;
	while(<$FILE>) {
		chomp;
		next if /^__/;
		my ($gene, $umi, $num) = split(/\t/, $_);
	#print "$gene\t$umi\t$num\n";
		$expr{$gene}{$umi}{$cell} = $num;
		$umi_stat{$cell} += $num;
		$umt_stat{$cell} += 1;
	}
	close $FILE;
}

#foreach my $cell (@cells) {
#	print "$cell\t$umi_stat{$cell}\t$umt_stat{$cell}\n";
#}

# 2. de-contamination
my (%conta_stat, %contt_stat);
foreach my $gene (keys %expr) {
	#print "> $gene\n";
	foreach my $umi (keys $expr{$gene}) {
		next if keys $expr{$gene}{$umi} == 1;
		my $max_value = 0;
		my $max_cell;
		foreach my $cell (keys $expr{$gene}{$umi}) {
			#print ">>> $gene\t$umi\t$cell\t$expr{$gene}{$umi}{$cell}\n";
			if($expr{$gene}{$umi}{$cell} > $max_value) {
				$max_value = $expr{$gene}{$umi}{$cell};
				$max_cell = $cell;
			}
		}
		my @max_cell_all;
		foreach my $cell (keys $expr{$gene}{$umi}) {
			push @max_cell_all, $cell if $expr{$gene}{$umi}{$cell} == $max_value;
		}
		next if @max_cell_all > 1;
		my $max_cid = $max_cell=~s/.*cell//r;
		my $max_row = $max_cid % 8;
		$max_row = 8 if $max_row == 0;
		my $max_col = POSIX::ceil($max_cid / 8);
		#print "+++ Max cell: $max_cell\t$max_value\t$max_row\t$max_col\n";
		foreach my $cell (keys $expr{$gene}{$umi}) {
			my $expr_ratio = $expr{$gene}{$umi}{$cell} / $max_value;
			if( ($expr_ratio > 0) && ($expr_ratio < $cutoff) ) {
				my $cel_cid = $cell=~s/.*cell//r;
				my $cel_row = $cel_cid % 8;
				$cel_row = 8 if $cel_row == 0;
				my $cel_col = POSIX::ceil($cel_cid / 8);
				my $conta_type;
				if ( ($cel_row != $max_row) && ($cel_col != $max_col) ) {
					$conta_type = "contaOth";
				} elsif ($cel_row == $max_row) {
					$conta_type = "contaRow";
				} elsif ($cel_col == $max_col) {
					$conta_type = "contaCol";
				}
				if($verbose) {
					printf "[debug] %s\t%s\t%s\t%d\t%d\t%d\t->\t%s\t%d\t%d\t%d\t:\t%s\n", $gene, $umi, $cell, $cel_row, $cel_col, $expr{$gene}{$umi}{$cell}, $max_cell, $max_row, $max_col, $max_value, $conta_type;
				}
				if(! exists $conta_stat{$cell}) {
					$conta_stat{$cell} = {"contaRow" => 0, "contaCol" => 0, "contaOth" => 0};
				}
				$conta_stat{$cell}{$conta_type} += $expr{$gene}{$umi}{$cell};
				if(! exists $contt_stat{$cell}) {
					$contt_stat{$cell} = {"contaRow" => 0, "contaCol" => 0, "contaOth" => 0};
				}
				$contt_stat{$cell}{$conta_type} += 1;
				$expr{$gene}{$umi}{$max_cell} += $expr{$gene}{$umi}{$cell} if $rescue;
				$expr{$gene}{$umi}{$cell} = 0;
			}
		}
	}
}

#foreach my $cell (@cells) {
#	printf "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", $cell, $umi_stat{$cell}, $conta_stat{$cell}{"contaRow"}, $conta_stat{$cell}{"contaCol"}, $conta_stat{$cell}{"contaOth"}, $contt_stat{$cell}{"contaRow"}, $contt_stat{$cell}{"contaCol"}, $contt_stat{$cell}{"contaOth"};
#}

# 3. create expression table
if($writeTab) {
	foreach my $cell (@cells) {
		print "> $cell\n";
		my $outpath = $outpt . "/" . $cell;
		mkdir $outpath if ! -e $outpath;
		open my $FILE, '>', ($outpath . "/" . $cell . "_htseq.tab");
		foreach my $gene (sort keys %expr) {
			foreach my $umi (sort keys $expr{$gene}) {
				printf $FILE ("%s\t%s\t%d\n", $gene, $umi, $expr{$gene}{$umi}{$cell}) if exists $expr{$gene}{$umi}{$cell};
			}
		}
		close $FILE;
	}
}

# 4. write stat
if($writeStat) {
	my $outpath = "03-expression/merged";
	mkdir $outpath if ! -e $outpath;
	open my $FILE, '>', ($outpath . "/itConta.stat");
	for my $cell (@cells) {
		printf $FILE ("%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", $cell, $umi_stat{$cell}, $conta_stat{$cell}{"contaRow"}, $conta_stat{$cell}{"contaCol"}, $conta_stat{$cell}{"contaOth"}, $umt_stat{$cell}, $contt_stat{$cell}{"conttRow"}, $contt_stat{$cell}{"conttCol"}, $contt_stat{$cell}{"conttOth"});
	}
}

# 5. write expr
if( $writeExpr || $writeMergedExpr) {
	my @genes;
	open my $FILE, "../../Genomes/human/gene_ID2Name.txt";
	while(<$FILE>) {
		chomp;
		my ($gene, $symbol) = split("\t", $_);
		push @genes, $gene;
	}
	close $FILE;

	sub umi_count() {
		my ($cell) = @_;
		my %oxpr;
		for my $gene (@genes) {
			for my $umi (keys $expr{$gene}) {
				if( (exists $expr{$gene}{$umi}{$cell}) && ($expr{$gene}{$umi}{$cell} > 0) ) {
					$oxpr{$cell}{$gene} += 1;
				}
			}
		}
		return(%oxpr);
	}

	my %oxpr;
	for my $cell (@cells) {
		$oxpr{$cell} = &umi_count($cell);
	}
	print "processing finished\n";

	if($writeExpr) {
		for my $cell (@cells) {
			my $outpath = $outpt . "/" . $cell;
			open my $FILE, ($outpath . "/_UMIcount_allGenes.txt");
			for my $gene (@genes) {
				printf $FILE ("%d\n", $oxpr{$cell}{$gene});
			}
			close $FILE;
		}
	}
}

