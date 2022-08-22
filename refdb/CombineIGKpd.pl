#!/usr/bin/perl -w
# usage : CombineIGKpd.pl

use strict;


# load sequences of IGKp and IGKd and combine

my %rseq;
LoadFasta("hsa_igkp.fa", \%rseq);
LoadFasta("hsa_igkd.fa", \%rseq);

my $igkseq = $rseq{"hsa_igkp"} . ('n' x 1000) . $rseq{"hsa_igkd"};

open OUT, ">hsa_igk.fa" || die "open hsa_igk.fa: $!\n";

my @igkseq = unpack("(A70)*", $igkseq);

print OUT ">hsa_igk\n";
print OUT join("\n", @igkseq) . "\n";

close OUT;


# load VDJ of IGKp and IGKd

open OUT, ">hsa_igk.vdj" || die "open hsa_igk.vdj: $!\n";

my $igkpl = length($rseq{"hsa_igkp"});

open IN, "<hsa_igkp.vdj" || die "open hsa_igkp.vdj: $!\n";
while (<IN>) {
    $_ =~ s/igkp/igk/;
    print OUT $_;
}

open IN, "<hsa_igkd.vdj" || die "open hsa_igkd.vdj: $!\n";

while (<IN>) {
    my @a = split "\t"; chomp $a[-1];

    chop($a[0]);
    my @er = split ",", $a[4];
    foreach my $er (@er) {
	my ($s, $e) = $er =~ /(\d+)\.\.(\d+)/;

	$s += $igkpl + 1000;
	$e += $igkpl + 1000;
	$er = "$s..$e";
    }
    $a[4] = join(",", @er);

    print OUT join("\t", @a) . "\n";
}



######################################################################


sub LoadFasta {
    my ($fn, $pidseq) = @_;
    
    open IN, "<$fn" || die "open $fn: $!\n";

    my $id;
    while (<IN>) {
        if (/^>(\S+)/) {
            $id = $1;
        } else {
            chomp $_;
            $pidseq->{$id} .= $_;
        }
    }
    close IN;
}
