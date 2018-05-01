#!/usr/bin/perl
use strict;
use 5.010;
use warnings;
#use Text::Levenshtein qw(distance);

sub hamming_distances {
        my ($seq1, @seqks) = @_;
        my @dsts;
        for(my $i=0; $i<@seqks; $i++) {
                my $seq2 = $seqks[$i];
                my @seq1s = split(//, $seq1);
                my @seq2s = split(//, $seq2);
                my $dst = 0;
                for(my $j=0; $j<@seq1s; $j++) {
                        $dst ++ if $seq1s[$j] ne $seq2s[$j];
                }
                push @dsts, $dst;
        }
        return(@dsts);
}

sub hamming_distance {
        my ($seq1, $seq2) = @_;
        my @seq1s = split(//, $seq1);
        my @seq2s = split(//, $seq2);
        my $dst = 0;
        for(my $i=0; $i<@seq1s; $i++) {
                $dst ++ if $seq1s[$i] ne $seq2s[$i];
        }
        return($dst);
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

#my $distance = &hamming_distance($ARGV[0], $ARGV[1]);
#print "$distance\n";

my $input = "AAAAA";
my @seqks = qw (AAAAG AAAAC AAAAT GCAAT);
my @outputs = &hamming_distances($input, @seqks);
print "@outputs\n";

my $output = &hamming_match($input, "TTTTA");
print "$output\n";

