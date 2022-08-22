#!/usr/bin/perl -w
# usage : ExtractVDJCInfo.pl gb species_gene

use strict;


# load genbank file and extract VDJC information

my %gq;

open IN, "<$ARGV[0]" || die "open $ARGV[0]: $!\n";

while (<IN>) {

    # at the line of annotation 
    if (/_segment/ || /_region/) {
	my @a = split /\s+/; chomp $a[-1];

	# obtain coordinate
	my $coord = $a[-1];
	while ($coord =~ /\(/ && $coord !~ /\)/) {
	    $_ = <IN>;
	    my @b = split /\s+/; chomp $b[-1];
	    $coord .= $b[-1];
	}

	# process coordinate
	my $o = $coord =~ /complement/ ? "-" : "+";
	if ($coord =~ /\(<?([\d\.,]+)\)/) {
	    $coord = $1;
	}

	# and gene name
	$_ = <IN>;
	my ($g) = /gene=\"(.+)\"/;

	# check pesudo gene till the end of record
	my $f = "F";
	while (<IN>) {
	    my $line = $_; chomp $line;
	    my @dq = $line =~ /\"/g;
	    while (@dq==1) {
		$_ = <IN>; chomp $_;
		$line .= $_;
		@dq = $line =~ /\"/g;
	    }
	    last if $line !~ /\//;
	    $f = "P" if $line =~ /pseudo/;
	}

	# output information while avoid redundant record
	if (!$gq{$g}) {
	    print "$ARGV[1]\t$g\t$f\t$o\t$coord\n";
	    $gq{$g} = 1;
	}
    }
}
close IN;
