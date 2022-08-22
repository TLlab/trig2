#!/usr/bin/perl -w
# usage : CombineTRADCDR.pl

use strict;


# load TRAD gene names

my $ref;
my %tratrad;
my %trads;

open IN, "<$ARGV[0]" || die "open $ARGV[0]: $!\n";

while (<IN>) {
    my @a = split "\t"; chomp $a[-1];
    $ref = $a[0];
    if ($a[1] =~ /TRA(.+)D/) {
	$tratrad{"TRA$1"} = $a[1];
    } else {
	$tratrad{$a[1]} = $a[1];
    }
    ($trads{$a[1]}) = $a[4] =~ /^(\d+)/; 
}
close IN;


# load TRA and TRD cdr info

my %tradcdr;

open IN, "<$ARGV[1]" || die "open $ARGV[1]: $!\n";
while (<IN>) {
    my @a = split "\t";

    my $trad = $tratrad{$a[1]};
    $tradcdr{$trad} = join("\t", @a[2..$#a]);
			  
}
close IN;

open IN, "<$ARGV[2]" || die "open $ARGV[2]: $!\n";
while (<IN>) {
    my @a = split "\t";

    $tradcdr{$a[1]} = join("\t", @a[2..$#a]);
}
close IN;


# output combined cdr

foreach my $trad (sort { $trads{$a} <=> $trads{$b} } keys %tradcdr) {
    print "$ref\t$trad\t$tradcdr{$trad}";
}
