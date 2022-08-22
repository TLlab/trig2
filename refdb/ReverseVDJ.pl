#!/usr/bin/perl -w
# usage : ReverseVDJ.pl vdj len

use strict;


# load VDJ info

my $ref;
my $rl = $ARGV[1] + 1;

my %gkv;

open IN, "<$ARGV[0]" || die "open $ARGV[0]: $!\n";

while (<IN>) {
    my @a = split "\t"; chomp $a[-1];

    $ref = $a[0];
    $gkv{$a[1]}{f} = $a[2];
    $gkv{$a[1]}{o} = $a[3] eq "-" ? "+" : "-";

    my @er = split ",", $a[4];
    foreach my $er (@er) {
	my ($s, $e) = $er =~ /(\d+)\.\.(\d+)/;
	my $rs = $rl - $e;
	my $re = $rl - $s;
	$er = "$rs..$re";
    }
    $gkv{$a[1]}{er} = join(",", reverse(@er));

    ($gkv{$a[1]}{s}) = $gkv{$a[1]}{er} =~ /^(\d+)/; 
}
close IN;


# output reversed VDJ info

foreach my $g (sort { $gkv{$a}{s} <=> $gkv{$b}{s} } keys %gkv) {
    print "$ref\t$g\t$gkv{$g}{f}\t$gkv{$g}{o}\t$gkv{$g}{er}\n";
}
