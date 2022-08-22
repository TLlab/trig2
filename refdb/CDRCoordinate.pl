#!/usr/bin/perl -w
# usage : CDRCoordinate.pl s_g_v_nt.aln s_g_j_nt.aln s_g_j_aa.aln s_g_exon.fa

use strict;


# FR and CDR regions based on IMGT unique numbering
# http://www.imgt.org/IMGTScientificChart/Nomenclature/IMGT-FRCDRdefinition.html

# FR1  : 1-26 (1st Cys(C) at 23)
# CDR1 : 27-38 
# FR2  : 39-55 (Trp(W) at 41)
# CDR2 : 56-65
# FR3  : 66-104 (2nd Cys(C) at 104)
# CDR3 : 105-117
# FR4  : 118-129 (Phe(F) or Trp(W) at 118)

my @frs  = (1, 39, 66, 118);
my @fre  = (26, 55, 104, 129);
my @cdrs = (27, 56, 105);
my @cdre = (38, 65, 117);


# load IMGT aligned nucleotide sequences

my %gaseq;

open IN, "<$ARGV[0]" || die "open $ARGV[0]: $!\n";
while (<IN>) {
    my ($g) = />(.+)\/?/;
    #my @g = split " and ", $g;
    my $seq = <IN>; chomp $seq;
    #foreach my $g (@g) {
	$gaseq{$g} = $seq;
    #}
}
close IN;

open IN, "<$ARGV[1]" || die "open $ARGV[1]: $!\n";
while (<IN>) {
    my ($g) = />(\S+)/;
    $gaseq{$g} = <IN>; chomp $gaseq{$g};
}
close IN;


# get numbering on J genes

my %jfin;

open IN, "<$ARGV[2]" || die "open $ARGV[2]: $!\n";
if ($ARGV[2] !~ /igh/) {
    while (<IN>) {
	my $seq = <IN>;
	my $i = index($seq, "FG");
	$jfin{$i}++ if $i != -1;
    }
} else {
    while (<IN>) {
	my $seq = <IN>;
	my $i = index($seq, "WG");
	$jfin{$i}++ if $i != -1;
    }
}
close IN;

my @jfi = sort { $jfin{$b} <=> $jfin{$a} } keys %jfin;
my $jfi = 3 * $jfi[0];


# load exon sequence

open IN, "<$ARGV[3]" || die "open $ARGV[3]: $!\n";

while (<IN>) {
    my ($ref, $er, $o, $g) = />(\S+):(\S+):([+-]):(\S+)/;
    my $seq = <IN>; chomp $seq;

    $ref =~ s/igkp/igk/;
    $ref =~ s/igkd/igk/;
    
    # skip non-V/J genes
    next if $g !~ /[VJ]/;
    
    # if no aligned sequence
    if (!$gaseq{$g}) {
	print "$ref\t$g" . ("\t---" x 8) ."\n";
	next;
    }

    # for V gene
    if ($g =~ /V/) {
	
	# get exon range
	my @er = split ",", $er;
	my @es;
	my @ee;
	for (my $i = 0; $i < @er; $i++) {
	    ($es[$i], $ee[$i]) = $er[$i] =~ /(\d+)\.\.(\d+)/;
	}
	
	# CDR regions
	# rf : reading frame
	my $rf;
	my @cdrcoord;
	for (my $i = 0; $i < @cdrs; $i++) {
	    my $s = 3 * ($cdrs[$i] - 1);
	    my $l = 3 * ($cdre[$i] - $cdrs[$i] + 1);
	    my $gaseg = uc(substr($gaseq{$g}, $s, $l));
	    $gaseg =~ s/\.//g;
	    my $gasegl = length($gaseg);

	    $s = index($seq, $gaseg);
	    while ($s == -1) {
		chop $gaseg;
		$s = index($seq, $gaseg);
	    }
	    $rf = $s % 3 if !(defined $rf);

	    my $e;

	    # if gene on the plus strand
	    if ($o eq "+") {

		# if two exons
		if ($es[1]) {
		    $s -= ($ee[0] - $es[0] + 1);
		    $s += $es[1];

		# if only one exon
		} else {
		    $s += $es[0];
		}
		$e = $s + $gasegl - 1;
		
	    # if on the minus strand (always two exons)
	    } else {
		$s -= ($ee[1] - $es[1] + 1) if $es[1];
		$s = $ee[0] - $s;
		$e = $s - $gasegl + 1;
	    }
	    my $coord = $i == 2 ? $s : "$s..$e";
	    push(@cdrcoord, $coord);
	}

	# FR regions
	my @frcoord;
	for (my $i = 0; $i < $#frs; $i++) {
	    my $s = 3 * ($frs[$i] - 1);
	    my $l = 3 * ($fre[$i] - $frs[$i] + 1);
	    my $gaseg = uc(substr($gaseq{$g}, $s, $l));
	    $gaseg =~ s/\.//g;
	    my $gasegl = length($gaseg);
	    
	    $s = index($seq, $gaseg);
	    while ($s == -1) {
		chop $gaseg;
		$s = index($seq, $gaseg);
	    }

	    my $e;
	    if ($o eq "+") {
		if ($es[1]) {
		    $s -= ($ee[0] - $es[0] + 1);
		    $s += $es[1];
		} else {
		    $s += $es[0];
		}
		$e = $s + $gasegl - 1;
		
	    } else {
		$s -= ($ee[1] - $es[1] + 1) if $es[1];
		$s = $ee[0] - $s;
		$e = $s - $gasegl + 1;
	    }
	    push(@frcoord, "$s..$e");
	}
	print "$ref\t$g\t$rf\t" . join("\t", @cdrcoord) ."\t" . join("\t", @frcoord) . "\t---\n";

    # for J gene
    } elsif ($g =~ /J/) {
	my ($es, $ee) = $er =~ /(\d+)\.\.(\d+)/;

	# CDR3 end
	my $gaseg = uc(substr($gaseq{$g}, $jfi));
	$gaseg =~ s/\.//g;
	while ((length($gaseg) % 3) != 0) {
	    chop $gaseg;
	}
	my $gasegl = length($gaseg);
	
	my $s = index($seq, $gaseg);
	while ($s == -1) {
	    chop $gaseg;
	    $s = index($seq, $gaseg);
	}
	my $rf = $s % 3;
	my $cdr3e = $es + $s - 1;

	# FR4
	my $fr4s = $cdr3e + 1;
	my $fr4e = $fr4s + $gasegl - 1;

	print "$ref\t$g\t$rf\t---\t---\t$cdr3e\t---\t---\t---\t$fr4s..$fr4e\n";
    }
}
close IN;
