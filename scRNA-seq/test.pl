use strict;
use Text::Levenshtein qw(distance);

my $out = distance("foo","four");
print "$out\n";

my @words     = qw/ four foo bar /;
my @distances = distance("foo", @words);

print "@distances\n";
