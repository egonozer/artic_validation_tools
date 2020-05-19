#!/usr/bin/perl

use strict;
use warnings;

my $usage = "
Enter file prefix

";

die $usage unless @ARGV;

my $pref = $ARGV[0];

##read in the mask file
my @marray;
open (my $min, "<$pref.coverage_mask.txt") or die "ERROR: Can't open $pref.coverage_mask.txt: $!\n";
while (my $line = <$min>){
    chomp $line;
    my @tmp = split("\t", $line);
    push @marray, ([@tmp[1..2]]);
}
close ($min);

## read in the pass vcf file
print "Masked PASSING:\n";
open (my $pin, "zcat $pref.pass.vcf.gz | ");
while (my $line = <$pin>){
    chomp $line;
    next if $line =~ m/^#/;
    my @tmp = split("\t", $line);
    my $pos = $tmp[1];
    foreach my $slice (@marray){
        my ($start, $stop) = @{$slice};
        if ((sort{$a<=>$b}($start, $stop, $pos))[1] == $pos){
            print "$line\n";
            next;
        }
    }
}
close ($pin);

print "Masked FAILING:\n";
open (my $fin, "<$pref.fail.vcf");
while (my $line = <$fin>){
    chomp $line;
    next if $line =~ m/^#/;
    my @tmp = split("\t", $line);
    my $pos = $tmp[1];
    foreach my $slice (@marray){
        my ($start, $stop) = @{$slice};
        if ((sort{$a<=>$b}($start, $stop, $pos))[1] == $pos){
            print "$line\n";
            next;
        }
    }
}
close ($fin);
