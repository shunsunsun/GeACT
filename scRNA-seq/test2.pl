use strict;
use Text::Levenshtein qw(distance);

sub mismatch_search {
        my ($query, @cands) = @_;
        my ($out, $dst);
        foreach my $cand (@cands) {
                my $distance = distance($query, $cand);
print "$query\t$cand\t$distance\n";
                if($distance <= 1) {
                        if(defined $out) {      # avoid multiple mapping
                                $out = "NA";
                                $dst = -2;
                                last;
                        } else {
                                $out = $cand;
                                $dst = $distance;
                        }
                }
        }
        if(! defined $out) {
                $out = "NA";
                $dst = -1;
        }
print "$query\t->\t$out\t$dst\n";
        return($out, $dst);
}

my $input = "AAAAAC";
my @inputs = qw /TAGAAC GGAAAG AATATC ATATAC AAACAC/;

&mismatch_search($input, @inputs);

