#!/usr/bin/perl
use warnings;
use strict;
use 5.010;
use Date::Parse;

if (@ARGV==0) {
print "perl $0 foo.o.txt\n";
exit;
}

my ($sp, $id);
if($ARGV[0]=~/.*do_[^_]+_([^_]+)_(.*)\.o\.txt/) {
	$sp = $1;
	$id = $2;
}

my ($time0, $time1);
open my $FILE , $ARGV[0] or die "Can't open o file.";
while(<$FILE>) {
	chomp;
	if(/^\[Start\]\t(.*)/) { $time0 = $1 }
	if(/^\[End\]\t(.*)/) { $time1 = $1 }
}
close $FILE;

#print "$time0\t$time1\n";
my $t0 = str2time($time0);
my $t1 = str2time($time1);
my $ts = ($t1 - $t0) / 60;
printf "%s\t%s\t%.4f\n", $sp, $id, $ts;

