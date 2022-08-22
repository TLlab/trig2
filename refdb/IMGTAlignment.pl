#!/usr/bin/perl -w
# usage : IMGTAlignment.pl species group

use strict;
use LWP::Simple qw(get);
use HTML::TableExtract;


# set parameters

my %scn = ("hsa" => "human", "mmu" => "mouse");
my %sfn = ("hsa" => "Homo%20sapiens", "mmu" => "Mus%20musculus");

my $species = $ARGV[0];
my $group = uc($ARGV[1]);


##### get alignments of all V genes

my $vg = $group . "V";
my $vurl = "http://www.imgt.org/IMGTrepertoire/index.php?section=LocusGenes&repertoire=genetable&species=$scn{$species}&group=$vg";
my $vhtml = get($vurl);
sleep(3);

my @headers = ("IMGT gene name");
my $te = HTML::TableExtract->new( headers => \@headers ); 
$te->parse($vhtml);

my @vag;
foreach my $ts ($te->tables) {
    foreach my $row ($ts->rows) {
	push(@vag, grep($_, @$row));
    }
}

open OUTvaa, ">$ARGV[0]\_$ARGV[1]\_v_aa.aln" || die "open $ARGV[0]\_$ARGV[1]\_v_aa.aln: $!\n";
open OUTvnt, ">$ARGV[0]\_$ARGV[1]\_v_nt.aln" || die "open $ARGV[0]\_$ARGV[1]\_v_nt.aln: $!\n";

foreach my $g (@vag) {
    my $url = "http://www.imgt.org/genedb/fasta?outputFormat=fastaAGaps&selectGenes=$g&selectSpecies=$sfn{$species}";
    my $html = get($url);
    sleep(3);

    my $aaseq = ParseAlignedSeq($html);
    print OUTvaa ">$g\n$aaseq\n" if $aaseq;
    
    $url = "http://www.imgt.org/genedb/fasta?outputFormat=fastaNGaps&selectGenes=$g&selectSpecies=$sfn{$species}";
    $html = get($url);
    sleep(3);

    my $ntseq = ParseAlignedSeq($html);
    print OUTvnt ">$g\n$ntseq\n" if $ntseq;
    
    #print "done: $g\n";
}

close OUTvaa;
close OUTvnt;


##### get alignments of all J genes

my @ajg;
my @ajaa;
my @ajnt;

my $jg = $group . "J";
my $jurl = "http://www.imgt.org/IMGTrepertoire/Proteins/alleles/index.php?species=$sfn{$species}&group=$jg&gene=$jg-overview";
my $jhtml = get($jurl);

my @line = split "\n", $jhtml;
for (my $i = 0; $i < @line; $i++) {
    my $j;
    
    # for TR*J
    if ($line[$i] =~ /<span class="fgxg">/) {
	$line[$i] = ">" . $line[$i] . "<";
	my (@aa) = $line[$i] =~ />([A-Z\s\*]+)</g;
	$_ =~ s/ //g for @aa;
	$i++;

	($j) = $line[$i] =~ />\s*($group.+)\*01<\/span><\/span>/;
	my ($nt) = $line[$i] =~ /gDNA\s+([a-z\s]+)\s*$/;
	my @nt = split " ", $nt;
	if (length($nt[0]) < 3) {
	    $nt[0] = ('.' x (3 - length($nt[0]))) . $nt[0];
	    $aa[0] = '.' . $aa[0];
	}
	if (length($nt[-1]) < 3) {
	    pop(@nt);
	}

	# take care of TRGJ alignments (not ending at the same A.A.)
	if ($j =~ /^TRGJ.$/) {
	    push(@aa, (".", "."));
	    push(@nt, ("...", "..."));
	}
	
	push(@ajg, $j);
	push(@ajaa, join("", @aa));
	push(@ajnt, join("", @nt));

    # for IG*J
    } elsif ($line[$i] =~ />($group.+)\*01<\/FONT>/) {
	$j = $1;
	my ($nt) = $line[$i] =~ /<\/FONT>(.+)\s*$/;
	$nt =~ s/\&nbsp;/ /g;
	$nt =~ s/^\s+//;
	my @nt = split " ", $nt;

	$line[$i-1] =~ s/^<FONT SIZE=-1>//;
	$line[$i-1] =~ s/\&nbsp;/ /g;
	$line[$i-1] = ">" . $line[$i-1] . "<"; 
	my (@aa) = $line[$i-1] =~ />([A-Z\s]+)</g;
	$_ =~ s/ //g for @aa;

	if (length($nt[0]) < 3) {
	    $nt[0] = ('.' x (3 - length($nt[0]))) . $nt[0];
	    $aa[0] = '.' . $aa[0];
	}

	push(@ajg, $j);
	push(@ajaa, join("", @aa));
	push(@ajnt, join("", @nt));
    }

    #print "done: $j\n" if $j;
}

my @jaal = sort { $b <=> $a } map(length($_), @ajaa);
my $maxjaal = $jaal[0];
foreach my $aa (@ajaa) {
    $aa = ('.' x ($maxjaal - length($aa))) . $aa;
}

my @jntl = sort { $b <=> $a } map(length($_), @ajnt);
my $maxjntl = $jntl[0];
foreach my $nt (@ajnt) {
    $nt = ('.' x ($maxjntl - length($nt))) . $nt;
}

open OUTjaa, ">$ARGV[0]\_$ARGV[1]\_j_aa.aln" || die "open $ARGV[0]\_$ARGV[1]\_j_aa.aln: $!\n";
open OUTjnt, ">$ARGV[0]\_$ARGV[1]\_j_nt.aln" || die "open $ARGV[0]\_$ARGV[1]\_j_nt.aln: $!\n";

for (my $i = 0; $i < @ajg; $i++) {
    print OUTjaa ">$ajg[$i]\n$ajaa[$i]\n";
    print OUTjnt ">$ajg[$i]\n$ajnt[$i]\n";
}

close OUTjaa;
close OUTjnt;



######################################################################


sub ParseAlignedSeq {
    my $html = $_[0];
    my $aseq = "";

    return $aseq if $html =~ /Number of results = 0/;

    $html =~ s/\015/\n/g;
    my @line = split "\n", $html;
    for (my $i = 0; $i < @line; $i++) {
        if ($line[$i] =~ /^>/) {
            $i++;
            while ($line[$i] && $line[$i] ne "" && $line[$i] !~ /^>/) {
                $aseq .= $line[$i];
                $i++;
            }
            last;
        }
    }
    
    return $aseq;
}
