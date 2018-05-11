#!/usr/bin/perl
use strict;
use 5.010;
use warnings;

if(@ARGV!=6) {
	print "Add barcode and UMI to read names, and clean reads.\n";
	print "Usage:   perl $0 fq1 fq2 barcode2sample.txt output_path trim_left trim_right\n";
	print "Example: perl $0 test/test_R1.fastq.gz test/test_R2.fastq.gz test/barcode2sample.txt cleanFq 2 131\n";
	exit;
}

my ($fq1, $fq2, $bcdToSp, $output_path, $trim_left, $trim_right) = @ARGV;
print "('$fq1', '$fq2', '$bcdToSp', '$output_path', '$trim_left', '$trim_right')\n";
#print "Start...\n";

# 0. functions
sub hamming_distances {
	my ($seq1, @seqks) = @_;
	my @seq1s = split(//, $seq1);
	my @dsts;
	for(my $i=0; $i<@seqks; $i++) {
		my $seq2 = $seqks[$i];
		my @seq2s = split(//, $seq2);
		my $dst = 0;
		for(my $j=0; $j<@seq1s; $j++) {	# for each base pair
			$dst ++ if $seq1s[$j] ne $seq2s[$j];
		}
		push @dsts, $dst;
	}
	return(@dsts);
}

sub hamming_match {
	my ($seq1, $seq2) = @_;
	my @seq1s = split(//, $seq1);
	my @seq2s = split(//, $seq2);
	my $match = 0;
	for(my $i=0; $i<@seq1s; $i++) {
		$match ++ if $seq1s[$i] eq $seq2s[$i];
	}
	return($match);
}

sub additional_dict {
	my (%input) = @_;
	my %dict_res;
	my @blacklists;
	foreach my $seq (keys %input) {	# for each barcode
		my $value = $input{$seq};
		my @seqs = split(//, $seq);
		for(my $i=0; $i<@seqs; $i++) {	# for each base pair
			my @seks = @seqs;
			foreach my $base (qw /A C G T N/) {
				$seks[$i] = $base;
				my $out = join "", @seks;
				next if $seq eq $out;
	#print "$out\t$value\n";
				if(exists $dict_res{$out}) {
					#print "[Warning] Sequence collision: $seq:$value ($i $base) $out -> $dict_res{$out}\n";
					if(! ($out ~~ @blacklists)) {
						push @blacklists, $out;
					}
				} else {
					$dict_res{$out} = $value;
				}
			}
		}
	}
	foreach my $blacklist (@blacklists) {
	#	print "> $blacklist\n";
		delete $dict_res{$blacklist};
	}
	return(%dict_res);
}

# 1. create_barcode_dict
#my @outer_barcodes = ("GATATG", "ATACG", "CCGTCTG", "TGCG", "GAACTCG", "ATGTAG", "CCCG", "TGTAG", "GAGTAAG", "ATCG", "CCTAG", "TGACCG");
#my @inner_barcodes = ("GTTGTT", "GTTAAA", "GTTTGG", "AGGGTT", "AGGAAA", "AGGTGG", "TAATGG", "GGAGAG");
my (@outer_barcodes, @inner_barcodes);
my $rt3            = ("AGTCGCTTGGGTGTAGTGC");

# read barcode info from file
my @barcode_file_ids;
my %bcdToSp;
open my $BCDTOSP , $ARGV[2];
while(<$BCDTOSP>) {
        chomp;
        my @lines = split(/\t/ , $_);
	push @outer_barcodes, $lines[0] unless $lines[0] ~~ @outer_barcodes;
	push @inner_barcodes, $lines[1] unless $lines[1] ~~ @inner_barcodes;
	push @barcode_file_ids , $lines[2];
        $bcdToSp{$lines[0] . "__" . $lines[1]} = $lines[2];
}
#print "@outer_barcodes" , "\n";
#print "@inner_barcodes" , "\n";

#foreach my $i (sort keys %bcdToSp) {
#	my @lines = split(/__/ , $i);
#	print "$lines[0]\t$lines[1]\t$bcdToSp{$i}\n";
#}
#exit;
#exit;

# extend the outer barcode, which will produce identical length
my $outer_longest_len;
foreach my $outer_barcode (@outer_barcodes) {
	my $outer_barcode_len = length($outer_barcode);
	if( (! defined $outer_longest_len) || ($outer_barcode_len > $outer_longest_len) ) { $outer_longest_len = $outer_barcode_len }
}
#print "Longest length: $outer_longest_len\n";
my(%outer_dict , %inner_dict);
# 1a. outer dict (basic)
for(my $i=0; $i<@outer_barcodes; $i++) {
	$outer_dict{substr( ($outer_barcodes[$i] . $rt3), 0, $outer_longest_len )} = ( $i . "_" . substr(($outer_barcodes[$i] . $rt3), $outer_longest_len) );	# outer_ext : id _ trimmed
}
#foreach my $key (keys %outer_dict) {
#	print "'$key' : ($outer_dict{$key})\n";
#}

# 1b. outer dict (additional)
my %outer_dict_mm = &additional_dict(%outer_dict);
#foreach my $key (sort keys %outer_dict_mm) {
#	print "$key\t$outer_dict_mm{$key}\n";
#}

# 2a. inner dict (basic)
for(my $i=0; $i<@inner_barcodes; $i++) {
	$inner_dict{$inner_barcodes[$i]} = $i;	# inner_barcode : id
}
#foreach my $key (keys %inner_dict) {
#	print "'$key' : $inner_dict{$key}\n";
#}

# 2b. inner dict (additional)
my %inner_dict_mm = &additional_dict(%inner_dict);
#foreach my $key (sort keys %inner_dict_mm) {
#	print "$key\t$inner_dict_mm{$key}\n";
#}

# the length of barcodes after extension
my $outer_len = length ( (keys %outer_dict)[0] );
my $inner_len = length ( (keys %inner_dict)[0] );
#print "$outer_len\t$inner_len\n";

# 2. split_fastq
sub mm_search {
	my ($query, @cands) = @_;
	my $out;
	#my $dst;
	my @distances = &hamming_distances($query, @cands);
	for(my $i=0; $i<@distances; $i++) {
		my $distance = $distances[$i];
		if($distance <= 1) {
			if(defined $out) {	# avoid multiple mapping
				$out = "NA";
				#$dst = -2;
				last;
			} else {
				$out = $cands[$i];
				#$dst = $distance;
			}
		}
	}
	if(! defined $out) {
		$out = "NA";
		#$dst = -1;
	}
#print "$query\t->\t$out\t$dst\n";
	return($out);
}

sub split_fastq {
	my ($seq) = @_;
	my ($outer_id, $spacer, $inner_id);
	my $mm_index = 0;	# whether rescued by mm
	my ($seq_left, $seq_right) = ( substr($seq,0,$outer_len), substr($seq,$outer_len) );	# split by outer barcode
	if(! exists $outer_dict{$seq_left}) {
#print "mm_search for ourter barcode...\n";
		if(! exists $outer_dict_mm{$seq_left}) {
			return ("NA", "NA", "-outer", $mm_index);
		} else {
			($outer_id, $spacer) = split (/_/ , $outer_dict_mm{$seq_left});
			$mm_index += 2;
		}
	} else {
		($outer_id, $spacer) = split (/_/ , $outer_dict{$seq_left});
	}
	my $spacer_len = length($spacer);
	#if( ($spacer_match eq "True") && ( (index $seq_right, $spacer) != 0) ) { return ($outer_id, -1, "-spacer") }
	# if False, the "$spacer" might not be that in the sequence
	my $seq_right_cut = substr($seq_right, $spacer_len, $inner_len);	# the inner barcode in sequences
	if(! exists $inner_dict{$seq_right_cut}) {
		if($mm_index > 0) {
			return ($outer_id, "NA", "-inner0", $mm_index);
		}
#print "mm_search for inner barcode...\n";
		if(! exists $inner_dict_mm{$seq_right_cut}) {
			return ($outer_id, "NA", "-inner", $mm_index);
		} else {
			$inner_id = $inner_dict_mm{$seq_right_cut};
			$mm_index += 1;
		}
	} else {
		$inner_id = $inner_dict{$seq_right_cut};
	}
	return ($outer_id, $inner_id, $seq_left . "_" . substr($seq_right, 0, $spacer_len) . "_" . $seq_right_cut, $mm_index);	# outer+spacer+inner
}

# 3. extract UMI
sub UMI_search {
	my ($seq) = @_;
	#my $umi1="HBDVHBDVHBDVHBDVHBDV";
	#my $umi2="VDBHVDBHVDBHVDBHVDBH";
	my $nec = "GACT" x 5;
	my $len = length $seq;
	my $comm = &hamming_match($seq, $nec);	# comm level
#print "$nec\t$seq\t$comm/$len\n";
	return($comm);
}

sub extract_umis {
	my ($seq, $umi_length, $min_t) = @_;
	#print STDERR "Warning: N in UMI region:\t[$seq]\n" if (substr($seq, 0, $umi_length)) =~ /[N]/;
	if ($seq =~ /^([ACGTN]{$umi_length})([T]{$min_t,})(.*)/) {
#		print ">>> $1 $2 $3\n";
		my ($UMI, $polyT, $cDNA) = ($1, $2, $3);
		my $mm_res = &UMI_search($UMI);
		if($mm_res < 3) {
			return ($UMI, $polyT, $cDNA, $mm_res);
		} else {
			return ("NA", "NA", $cDNA, $mm_res);
		}
	} else {
		return ("NA", "NA", "NA", -1);
	}
}

# 4. trim adapter

# 5. filter low quality and too many N reads
sub filter_by_quality {
	my ($seq, $quality) = @_;
	#print "$seq\n";
	#my $seq_len = length $seq;

	# trim polyA
	my $polyA = "A" x 7;
	my $polyA_r1 = index $seq, $polyA;
	if($polyA_r1 >= 50) {
		$seq = substr($seq, 0, $polyA_r1);
		$quality = substr($quality, 0, $polyA_r1);
	}

	# remove short reads
	my $seq_len = length $seq;
	if($seq_len < 40) { return ("NA", -1) };

	# remove low quality reads
	#my ($Q_20, $Q_30) = (53, 63);
	my $Q_min = 33 + 5;
	my $low_qual_num = 0;
	for(my $i=0; $i<$seq_len; $i++) {
		my $base_asc = substr($quality, $i, 1);
		my $base_quality = ord($base_asc);
		$low_qual_num ++ if($base_quality <= $Q_min);
	}
	if($low_qual_num >= ($seq_len * 0.5)) { return ("NA", -2) };

	# remove reads contain >=10% N
	my $N_num = ($seq =~ tr/N/N/) + 0;
	if($N_num >= ($seq_len * 0.1)) { return ("NA", -3) };

	return ($seq, $quality);
}


# X. main
system "mkdir -p $output_path";
my %read;
my %reap;
my $seqnum = 0;	# read number
my ($identified_num, $unidentified_num, $identified_noMm_num, $unidentified_noMm_num) = (0) x 4;
my ($UMIed_num, $unUMIed_num, $UMIed_noMm_num, $unUMIed_noMm_num, $qualified_num, $unqualified_num) = (0) x 6;
my (%identified_num, %identified_noMm_num, %UMIed_num, %UMIed_noMm_num, %qualified_num);
open my $READ1, "gzip -dc $ARGV[0] |" or die "Cannot open Read1 file.";
open my $READ2, "gzip -dc $ARGV[1] |" or die "Cannot open Read2 file.";

#open unidentified_read1_output , ">" , ($output_path . "/unidentified_read1.fastq") or die;
#open unidentified_read2_output , ">" , ($output_path . "/unidentified_read2.fastq") or die;
#open unUMIed_read1_output , ">" , ($output_path . "/unUMIed_read1.fastq") or die;
#open unUMIed_read2_output , ">" , ($output_path . "/unUMIed_read2.fastq") or die;

#my (%read1_outputs, %read2_outputs);	# barcode level
my (%UMIed_read1_outputs, %UMIed_read2_outputs);	# UMI level
#my @barcode_file_ids;
#for(my $i=0; $i<@outer_barcodes; $i++) {
#	foreach my $j (("A".."Z")[0..(@inner_barcodes-1)]) {
		#open $read1_outputs{"${i}${j}"} , ">" , ($output_path . "/" . ($i+1) . $j . "_read1.fastq") or die;
		#open $read2_outputs{"${i}${j}"} , ">" , ($output_path . "/" . ($i+1) . $j . "_read2.fastq") or die;
#		open $UMIed_read1_outputs{"${i}${j}"} , ">" , ($output_path . "/" . ($i+1) . $j . "_read1.fastq") or die;
#		open $UMIed_read2_outputs{"${i}${j}"} , ">" , ($output_path . "/" . ($i+1) . $j . "_read2.fastq") or die;
#		push @barcode_file_ids , (($i+1) . $j);
#	}
#}

foreach my $i (sort keys %bcdToSp) {
	my $output = $output_path . "/" . $bcdToSp{$i} . "_read1.fastq.gz";
	open $UMIed_read1_outputs{$i} , " | gzip -c > $output" or die;
}

while(1) {
	my $line1_1 = <$READ1>;
	my $line1_2 = <$READ1>;
	my $line1_3 = <$READ1>;
	my $line1_4 = <$READ1>;

	my $line2_1 = <$READ2>;
	my $line2_2 = <$READ2>;
	my $line2_3 = <$READ2>;
	my $line2_4 = <$READ2>;

	last if (! defined $line1_1);
	chomp ($line1_1, $line1_2, $line1_3, $line1_4, $line2_1, $line2_2, $line2_3, $line2_4);

	$read{name} = substr( (split(/ / , $line1_1))[0] , 1); $read{comment} = (split(/ / , $line1_1))[1];
	$reap{name} = substr( (split(/ / , $line2_1))[0] , 1); $reap{comment} = (split(/ / , $line2_1))[1];
	$read{sequence} = $line1_2;
	$reap{sequence} = $line2_2;
	$read{optional} = $line1_3;
	$reap{optional} = $line2_3;
	$read{quality} = $line1_4;
	$reap{quality} = $line2_4;

	$seqnum ++;
	
	#foreach my $key (qw /name comment sequence optional quality/ ) { print "[read1]\t$key\t$read{$key}\n"; }
	#foreach my $key (qw /name comment sequence optional quality/ ) { print "[read2]\t$key\t$reap{$key}\n"; }
	### name check for read pair
	if($read{name} ne $reap{name}) { print STDERR "Warning: inconsistence read name [$seqnum] $read{name} $reap{name}\n"; }
	my ($i1, $i2, $barcode_seq, $mm_idx) = &split_fastq($reap{sequence});	# read2
#print "$i1, $i2, $barcode_seq, $mm_idx\n";
	if($i2 eq "NA") {	# not allow mismatch in both outer and inner barcode (performed above)
		#print unidentified_read1_output "\@$read{name} $read{comment}\n$read{sequence}\n$read{optional}\n$read{quality}\n";
		#print unidentified_read2_output "\@$reap{name} $reap{comment}\n$reap{sequence}\n$reap{optional}\n$reap{quality}\n";
		$unidentified_num ++;
		$unidentified_noMm_num ++;
	} else {
		my $bcds = $outer_barcodes[$i1] . "__" . $inner_barcodes[$i2];
		my $sample = $bcdToSp{$bcds};
		#$read{name} .= ("_" . substr($read{comment},2) . "_" . $barcode_seq . "_" . substr( $reap{quality}, 0, length($barcode_seq)) );
		$read{name} .= ("_" . $barcode_seq . "_" . $mm_idx);
		# sequence no need trimming for R1
		# quality no need trimming for R1
		#$reap{name} .= ("_" . substr($reap{comment},2) . "_" . $barcode_seq . "_" . substr( $reap{quality}, 0, length($barcode_seq)) );
		$reap{name} .= ("_" . $barcode_seq . "_" . $mm_idx);
		$reap{sequence} = substr($reap{sequence}, length($barcode_seq)-2);
		$reap{quality} = substr($reap{quality}, length($barcode_seq)-2);
		#print { $read1_outputs{($i1 . ("A".."Z")[$i2])} } "\@$read{name} $read{comment}\n$read{sequence}\n$read{optional}\n$read{quality}\n";
		#print { $read2_outputs{($i1 . ("A".."Z")[$i2])} } "\@$reap{name} $reap{comment}\n$reap{sequence}\n$reap{optional}\n$reap{quality}\n";
		$identified_num{$sample} ++;
		$identified_num ++;
		if($mm_idx > 0) {	# should be missed without mm
			$unidentified_noMm_num ++;
		} else {
			$identified_noMm_num{$sample} ++;
			$identified_noMm_num ++;
		}

		my ($UMI, $polyT, $kept_seg, $mm_idx) = &extract_umis($reap{sequence}, 20, 0);	# read2
		if( ($UMI eq "NA") || length($kept_seg)<20 ) {
			$unUMIed_num ++;
			if( ($mm_idx >= 3) && (length($kept_seg)>=20) ) {	# should be kept without mm
				$UMIed_noMm_num{$sample} ++;
				$UMIed_noMm_num ++;
			}
			else {
				$unUMIed_noMm_num ++;
			}
			#print ">>> $reap{sequence} $UMI $polyT $kept_seg $mm_idx\n";
			#print unUMIed_read1_output "\@$read{name} $read{comment}\n$read{sequence}\n$read{optional}\n$read{quality}\n";
			#print unUMIed_read2_output "\@$reap{name} $reap{comment}\n$reap{sequence}\n$reap{optional}\n$reap{quality}\n";
		} else {
			$UMIed_num{$sample} ++;
			$UMIed_num ++;
			$UMIed_noMm_num{$sample} ++;
			$UMIed_noMm_num ++;

			my ($filter_result, $filter_quality) = &filter_by_quality($read{sequence}, $read{quality});	# read1
			if($filter_result eq "NA") {
				#print ">>>[$filter_result]\n\@$read{name} $read{comment}\n$read{sequence}\n$read{optional}\n$read{quality}\n";
				$unqualified_num ++;
			} else {
				#$read{name} = ($UMI . "_" . length($polyT) . "T_" . substr( $reap{quality}, 0, length($UMI) ) . "_" . $read{name});
				$read{name} = ($UMI . "_" . $mm_idx . "_" . length($polyT) . "T" . "_" . $read{name});
				$read{sequence} = $filter_result;
				$read{quality} = $filter_quality;
				if($trim_right > 0) {
					$read{sequence} = substr( $read{sequence}, $trim_left, ($trim_right - $trim_left) );
					$read{quality} = substr( $read{quality}, $trim_left, ($trim_right - $trim_left) );
				}
				#$reap{name} = $read{name};	# read2 will NOT be printed XXX
				#$reap{sequence} = substr( $kept_seg, 0, $trim_right);
				#$reap{quality} = substr( substr($reap{quality}, (length($UMI)+length($polyT)), ), 0, $trim_right);
				print { $UMIed_read1_outputs{$bcds} } "\@$read{name} $read{comment}\n$read{sequence}\n$read{optional}\n$read{quality}\n";
				#print { $UMIed_read2_outputs{$bcds} } "\@$reap{name} $reap{comment}\n$reap{sequence}\n$reap{optional}\n$reap{quality}\n";
				$qualified_num{$sample} ++;
				$qualified_num ++;
			}
		}
	}
}

sub do_stat {
# after all reads were processed
print "Total\t$seqnum\n";
print (("-" x 40) . "\n");
print "State\tBarcode_noMm\tBarcode\tUMI_noMm\tUMI\tQC\n";
print "Pass\t$identified_noMm_num\t$identified_num\t$UMIed_noMm_num\t$UMIed_num\t$qualified_num\n";
print "Fail\t$unidentified_noMm_num\t$unidentified_num\t$unUMIed_noMm_num\t$unUMIed_num\t$unqualified_num\n";
print (("-" x 40) . "\n");
for (@barcode_file_ids) { print "$_\t$identified_noMm_num{$_}\t$identified_num{$_}\t$UMIed_noMm_num{$_}\t$UMIed_num{$_}\t$qualified_num{$_}\n"; }
print (("-" x 40) . "\n");
#print "OK.\n";
}

do_stat;

