#!/usr/bin/perl -w
# usage : ExonSeq.pl fa vdj

use strict;


# load reference sequence

my $rseq;

open IN, "<$ARGV[0]" || die "open $ARGV[0]: $!\n";

<IN>;
while (<IN>) {
    chomp;
    $rseq .= $_;
}
close IN;


# load VDJC exon info

open IN, "<$ARGV[1]" || die "open $ARGV[1]: $!\n";

while (<IN>) {
    my @a = split "\t"; chomp $a[-1];

    my $seg;
    my @exon = split ",", $a[4];

    if ($a[1] =~ /^...[VJ]/) {
	foreach my $ex (@exon) {
	    my ($s, $e) = $ex =~ /(\d+)\.\.(\d+)/;
	    $seg .= substr($rseq, $s - 1, $e - $s + 1);
	}
    } else {
	my ($s, $e) = $exon[0] =~ /(\d+)\.\.(\d+)/;
	$seg .= substr($rseq, $s - 1, $e - $s + 1);
	$a[4] = "$s..$e";
    }

    $seg = ReverseComplement($seg) if $a[3] eq "-";
    
    print ">$a[0]:$a[4]:$a[3]:$a[1]\n$seg\n";
}
close IN;



######################################################################


sub ReverseComplement {
    my $s = reverse($_[0]);
    $s =~ tr/acgtACGT/tgcaTGCA/;
    return $s;
}
