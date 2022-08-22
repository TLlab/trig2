#!/usr/bin/perl -w
# usage : IMGTAlignment.pl species group
# note  : species: human, gene: tra, trb, trd, trg, igh, igl, igk

use strict;
use LWP::Simple qw(get);
use HTML::TableExtract;

# set parameters

my %scn = ("hsa" => "human", "mmu" => "mouse");
my %sfn = ("hsa" => "Homo%20sapiens", "mmu" => "Mus%20musculus");

my $species = $sfn{$ARGV[0]};
my $group = uc($ARGV[1]);


##### get alignments of all V genes
=pod
my $vg = $group . "V";
my @vag;
my @vaaseq;
my @vntseq;

GetAlignedSeq($ARGV[0], $vg, \@vag, \@vaaseq, \@vntseq);

open OUTvaa, ">$ARGV[0]\_$ARGV[1]\_v_aa.aln1" || die "open $ARGV[0]\_$ARGV[1]\_v_aa.aln1: $!\n";
open OUTvnt, ">$ARGV[0]\_$ARGV[1]\_v_nt.aln1" || die "open $ARGV[0]\_$ARGV[1]\_v_nt.aln1: $!\n";

for (my $i = 0; $i < @vag; $i++) {
    next if !$vaaseq[$i];
    print OUTvaa ">$vag[$i]\n$vaaseq[$i]\n";
    print OUTvnt ">$vag[$i]\n$vntseq[$i]\n";
}

close OUTvaa;
close OUTvnt;
=cut

##### get alignments of all J genes

my $jg = $group . "J";
my @jag;
my @jaaseq;
my @jntseq;

GetAlignedSeq($ARGV[0], $jg, \@jag, \@jaaseq, \@jntseq);

open OUTjaa, ">$ARGV[0]\_$ARGV[1]\_j_aa.aln1" || die "open $ARGV[0]\_$ARGV[1]\_j_aa.aln1: $!\n";
open OUTjnt, ">$ARGV[0]\_$ARGV[1]\_j_nt.aln1" || die "open $ARGV[0]\_$ARGV[1]\_j_nt.aln1: $!\n";

for (my $i = 0; $i < @jag; $i++) {
    next if !$jaaseq[$i];
    print OUTjaa ">$jag[$i]\n$jaaseq[$i]\n";
    print OUTjnt ">$jag[$i]\n$jntseq[$i]\n";
}

close OUTjaa;
close OUTjnt;



######################################################################


sub GetAlignedSeq {
    my ($species, $gene, $pag, $paaseq, $pntseq) = @_;
    my @ag;
    my @aaseq;
    my @ntseq;

    my $url = "http://www.imgt.org/IMGTrepertoire/index.php?section=LocusGenes&repertoire=genetable&species=$scn{$species}&group=$gene";
    my $html = get($url);
    sleep(3);

    my @headers = $gene =~ /V/ ? ("IMGT gene name") : ("IMGT allele name");
    my $te = HTML::TableExtract->new( headers => \@headers ); 
    $te->parse($html);
    
    foreach my $ts ($te->tables) {
	foreach my $row ($ts->rows) {
	    push(@ag, grep($_, @$row));
	}
    }
    
    foreach my $g (@ag) {
	my $url = "http://www.imgt.org/genedb/fasta?outputFormat=fastaAGaps&selectGenes=$g&selectSpecies=$sfn{$species}";
	my $html = get($url);
	sleep(3);
	
	my $aaseq = ParseAlignedSeq($html);
	push(@aaseq, $aaseq);

	$url = "http://www.imgt.org/genedb/fasta?outputFormat=fastaNGaps&selectGenes=$g&selectSpecies=$sfn{$species}";
	$html = get($url);
	sleep(3);
	
	my $ntseq = ParseAlignedSeq($html);
	push(@ntseq, $ntseq);

	print "done: $g\t$aaseq\t$ntseq\n";
    }

    @{$pag} = @ag;
    @{$paaseq} = @aaseq;
    @{$pntseq} = @ntseq;
}


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
