#!/usr/bin/perl

use strict;
use warnings;

my $samtools = "/home/eoz639/p30002/Applications/samtools-1.9/samtools";
my $primers = "/home/eoz639/Applications/artic-ncov2019/primer_schemes/nCoV-2019/V3/nCoV-2019.scheme.bed";

my $usage = "
artic_amplicon_coverage <primertrimmed.rg.sorted.bam>

Determines the total coverage across non-overlapping regions of each amplicon

Options:
  -n    Normalize results for each amplicon to median depth for the whole genome
        (default: output median depth per amplicon)

";

use Getopt::Std;
use vars qw( $opt_n );
getopts('n');

die $usage unless @ARGV;
my %amps;
## Read in the amplicon coords
open (my $pin, "<$primers") or die "ERROR: Can't open $primers: $!\n";
while (my $line = <$pin>){
    chomp $line;
    my @tmp = split("\t", $line);
    my ($start, $stop, $id) = @tmp[1..3];
    (my ($num, $pos)) = $id =~ m/nCoV-2019_(\d+)_([LR])/;
    @{$amps{$num}} = (0,30000) unless $amps{$num};
    if ($pos eq "L"){
        $amps{$num}[0] = ($stop + 1) if ($stop + 1) > $amps{$num}[0];
    } else {
        $amps{$num}[1] = ($start - 1) if ($start - 1) < $amps{$num}[1];
    }
}
close ($pin);

my @trim_amps;
my ($last_amp, $last_start, $last_stop);
foreach my $amp (sort{$a <=> $b} keys %amps){
    my ($start, $stop) = @{$amps{$amp}};
    #print "$amp $start-$stop\n";
    unless (@trim_amps){
        push @trim_amps, ([$amp, $start, $stop]);
        next;
    }
    my $last = $trim_amps[$#trim_amps][2];
    $trim_amps[$#trim_amps][2] = $start;
    push @trim_amps, ([$amp, $last, $stop]);
}

#foreach (@trim_amps){
#    print STDERR join(" ", @{$_}), "\n";
#}
#die;

my @results;
foreach my $file (@ARGV){
    (my $fid) = $file =~ m/([^\/]+)(\.primertrimmed.rg)*\.sorted.bam/;
    print STDERR "Reading $fid\n";
    my @tdepth;
    my @medians;
    foreach (@trim_amps){
        my ($amp, $start, $stop) = @{$_};
        open (my $din, "$samtools depth -aa -r MN908947.3:$start-$stop $file | ");
        my @deps;
        while (my $line = <$din>){
            chomp $line;
            my $dep = (split("\t", $line))[2];
            push @deps, $dep;
        }
        close ($din);
        push @tdepth, @deps;
        my ($avg, $stdev, $ste, $min, $negiqr, $q1, $med, $q3, $posiqr, $max) = stats(\@deps);
        push @medians, ([$amp, $avg, $med]);
        #print STDERR "\t$amp\t$avg\t$med\n";
    }
    my ($avg, $stdev, $ste, $min, $negiqr, $q1, $med, $q3, $posiqr, $max) = stats(\@tdepth);
    push @medians, (["tot", $avg, $med]);
    push @results, [$fid, [@medians]];
}

print "amplicon";
foreach my $slice (@results){
    my @tmp = @{$slice};
    print "\t$tmp[0]";
}
print "\tavg\tstdev\tste\tmin\tnegiqr\tq1\tmed\tq3\tposiqr\tmax\n";
for my $i (0 .. $#trim_amps){
    my ($amp, $start, $stop) = @{$trim_amps[$i]};
    print "$amp";
    my @vals;
    foreach my $slice (@results){
        my @tmp = @{$slice};
        my @tmp2 = @{$tmp[1]};
        my $avg_tdep = $tmp2[$#tmp2][1];
        my $med_tdep = $tmp2[$#tmp2][2];
        my $avg_dep = $tmp2[$i][1];
        my $med_dep = $tmp2[$i][2];

        my $norm = 0;
        $norm = $med_dep unless $opt_n;
        ## using average
        #if ($avg_tdep){
        #   $norm = sprintf("%.04f", $avg_dep / $avg_tdep);
        #}
        ## using median
        if ($med_tdep){
            $norm = sprintf("%.04f", $med_dep / $med_tdep) if $opt_n;
        }

        print "\t$norm";
        push @vals, $norm;
    }
    my @stat = stats(\@vals);
    print "\t", join("\t", @stat), "\n";
}
print "total";
my @vals;
foreach my $slice (@results){
    my @tmp = @{$slice};
    my @tmp2 = @{$tmp[1]};
    my $avg_tdep = $tmp2[$#tmp2][1];
    my $med_tdep = $tmp2[$#tmp2][2];

    my $val;
    ## using average
    #$val = $avg_tdep;
    ## using median
    $val = $med_tdep;

    print "\t$val";
    push @vals, $val;
}
my @stat = stats(\@vals);
print "\t", join("\t", @stat), "\n";



#------------------------------------------------------------------------------
sub stats {
    my @sizes = @{$_[0]};
    my @out = ("NA") x 10;
    #my ($avg, $stdev, $ste, $min, $negiqr, $q1, $med, $q3, $posiqr, $max) = ("NA") x 10;
    return (@out) unless scalar (@sizes) >= 2;
    @sizes = sort{$a <=> $b}@sizes;
    ($out[3], $out[9]) = ($sizes[0], $sizes[$#sizes]);
    my $num = scalar @sizes;
    my $sum = 0;
    foreach (@sizes){
        $sum += $_;
    }
    my $mean = $sum / $num;
    my $sqtotal = 0;
    foreach (@sizes){
        $sqtotal += ($mean - $_) ** 2;
    }
    my $stdev = sqrt($sqtotal / ($num - 1));
    my $ste = $stdev / sqrt($num);
    ($out[0], $out[1], $out[2]) = (sprintf("%.1f", $mean), sprintf("%.1f", $stdev), sprintf("%.1f", $ste));

    my $median = median(\@sizes);
    my ($q1, $q3) = ($median) x 2;
    if ($num > 1){
        my (@aq1, @aq3);
        if ($num % 2 == 0){ #even
            my $mid2 = $num / 2;
            my $mid1 = $mid2 - 1;
            @aq1 = @sizes[0 .. $mid1];
            @aq3 = @sizes[$mid2 .. $#sizes];
        } else { #odd
            my $mid = ($num / 2) - 0.5;
            @aq1 = @sizes[0 .. ($mid - 1)];
            @aq3 = @sizes[($mid + 1) .. $#sizes];
        }
        $q1 = median(\@aq1);
        $q3 = median(\@aq3);
    }
    my $iqr = $q3 - $q1;
    my $iqr2 = $iqr * 1.5;
    my $negiqr = $q1 - $iqr2;
    $negiqr = $out[3] if $negiqr < $out[3];
    my $posiqr = $q3 + $iqr2;
    $posiqr = $out[9] if $posiqr > $out[9];
    ($out[4], $out[5], $out[6], $out[7], $out[8]) = ($negiqr, $q1, $median, $q3, $posiqr);
    my @outliers;
    foreach my $val (@sizes){
        push @outliers, $val if ($val < $negiqr or $val > $posiqr);
    }
    push @out, @outliers;
    return(@out);
}

sub median {
    my @array = @{$_[0]};
    my $leng = scalar @array;
    my $med;
    if ($leng % 2 == 0){
        my $mid2 = $leng / 2;
        my $mid1 = $mid2 - 1;
        $med = ($array[$mid1] + $array[$mid2]) / 2;
    } else {
        my $mid = ($leng / 2) - 0.5;
        $med = $array[$mid];
    }
    return($med);
}
