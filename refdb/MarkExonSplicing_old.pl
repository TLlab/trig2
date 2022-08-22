#!/usr/bin/perl -w
# usage : MarkExonSplicing.pl vdj ref.fa

use strict;


# set IUPAC code

my %cn = (
    'A' => 'A', 'C' => 'C', 'G' => 'G', 'T' => 'T', 'U' => 'U',
    'R' => '[AG]', 'Y' => '[CT]', 'N' => '[ACGT]', 'W' => '[AT]',
    'S' => '[CG]', 'M' => '[AC]', 'K' => '[GT]', 'B' => '[CGT]',
    'H' => '[ACT]', 'D' => '[AGT]', 'V' => '[ACG]'
    );

my ($rid) = $ARGV[0] =~ /(.+)\.vdj/;
$rid =~ s/_r$//;
    
 
# load VDJC exon information

my @es;
my @ee;

open IN, "<$ARGV[0]" || die "open $ARGV[0]: $!\n";

while (<IN>) {
    my @a = split "\t"; chomp $a[-1];

    my @exon = split ",", $a[4];
    foreach my $ex (@exon) {
	my ($s, $e) = $ex =~ /(\d+)\.\.(\d+)/;
	push(@es, $s);
	push(@ee, $e);
    }
}
close IN;


# load reference sequence

my $rseq;

open IN, "<$ARGV[1]" || die "open $ARGV[1]: $!\n";

<IN>;
while (<IN>) {
    chomp;
    $rseq .= $_;
}
close IN;

my $rsl = length($rseq);


# mark exon, i.e., turn intron bases into lower case

my $rseqex;

my $lastee = 0;     # last exon end
for (my $i = 0; $i < @es; $i++) {
    $rseqex .= lc(substr($rseq, $lastee, $es[$i] - $lastee - 1));
    my $l = $ee[$i] - $es[$i] + 1;
    $rseqex .= substr($rseq, $es[$i] - 1, $l);
    $lastee = $ee[$i];
}
$rseqex .= lc(substr($rseq, $lastee, length($rseq) - $lastee));


# get potential splicing sites

my @donor = MotifLocus($rseq, "GGT");
@donor = map { $_ - 1 } @donor;

my @acceptor = MotifLocus($rseq, "YAG");
@acceptor = map { $_ + 3 } @acceptor;

my @splice = Unique(@donor, @acceptor);

@splice = sort { $a <=> $b } @splice;

my @ss = ($splice[0]);
my @se = ();
my $se = $splice[0] + 1;

for (my $i = 1; $i < @splice; $i++) {
    if ($splice[$i] > $se) {
	push(@se, $se);
	push(@ss, $splice[$i]);
    }
    $se = $splice[$i] + 1;
}
push(@se, $se);

while ($ss[-1] > $rsl || $se[-1] > $rsl) {
    if ($ss[-1] > $rsl) {
	pop(@ss);
	pop(@se);
    } else {
	$se[-1] = $rsl;
    }
}


# mark lable potential splicing site

my $rseqexsp = "";

$lastee = 0;
for (my $i = 0; $i < @ss; $i++) {
    $rseqexsp .= substr($rseqex, $lastee, $ss[$i] - $lastee - 1);
    $rseqexsp .= uc(substr($rseqex, $ss[$i] - 1, $se[$i] - $ss[$i] + 1));
    $lastee = $se[$i];
}
$rseqexsp .= substr($rseqex, $lastee, length($rseqex) - $lastee);


# output marked sequence

my @rseqexsp = unpack("(A70)*", $rseqexsp); 

print ">$rid\n" . join("\n", @rseqexsp) . "\n";



######################################################################


sub Unique {
    my %seen;
    return grep(!$seen{$_}++, @_);
}


sub MotifLocus {
    my ($seq, $motif) = @_;
    my @locus = ();
    
    # reformat sequence and motif
    $seq = uc($seq);
    my @m = split "", uc($motif);
    @m = map($cn{$_}, @m);
    $motif = join "", @m;

    # find motif pattern
    my @pat = $seq =~ /($motif)/g;

    # return if no pattern
    return @locus if !@pat;

    # get motif locus
    my $last = 0;
    for (my $i = 0; $i < @pat; $i++) {
        my $p = index($seq, $pat[$i], $last);
        push(@locus, $p+1);
        $last = $p + 1;
    }
    
    return @locus;
}
