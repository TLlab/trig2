#!/usr/bin/perl -w

use strict;
use File::Basename qw(dirname);
use Cwd qw(abs_path);
use Getopt::Long qw(GetOptions);

use lib dirname(abs_path($0));
use VDJDelta qw(SetVDJInfo LoadOneDelta GetOptimalSet OutputAln);
use GroupAln qw(SetScoringFunction FilterAlignment AdjustOverlap SetAnnotation OutputGroupAln);
use Fastx qw(LoadFasta LoadThisFasta LoadOneFasta);


##### set parameters

my $species  = "hsa";
my $gene     = "trb";
my $minmatch = 15;
my $adjolq   = 1;
my $frac     = 0.5;

GetOptions(
    "species=s"  => \$species,
    "gene=s"     => \$gene,
    "minmatch=i" => \$minmatch,
    "adjolq=i"   => \$adjolq,
    "frac=f"   => \$frac,
    );

my $sg = "$species\_$gene";


##### set VDJ and CDR3 info

SetVDJInfo("$sg.vdj");
SetScoringFunction();
GroupAln::setparas($minmatch);


##### process nucmer alignments

open my $afh, "<$ARGV[0]" || die "open $ARGV[0]: $!\n";

# load refernece sequence and set query file handle
$_ = <$afh>; chomp; <$afh>;
my ($ref, $qry) = /^(\S+) (\S+)/;

my %rseq;
LoadFasta($ref, \%rseq);

open my $qfh, "<$qry" || die "open $qry: $!\n";

# load nucmer alignments
while( <$afh> ) {
    my ($r, $q, $rl, $ql) = /^>(\S+) (\S+) (\d+) (\d+)/;

    # load all alignments of a query
    my @aln = LoadOneDelta($afh, $r, $q);

    # get optimal alignments
    my @oaln = GetOptimalSet(@aln);
    $_->annotate() for @oaln;

    # filter alignments by annotation
    my @galn = FilterAlignment(@oaln);

    # skip if no alignment remained after filtering
    next if !@galn;

    # load query sequence for adjusting overlap and output
    my $qseq = LoadThisFasta($qfh, $q, 1);

    # adjust overlap
    if ($adjolq) {
	AdjustOverlap(\%rseq, \$qseq, \@galn);
    }
    
    # skip if no alignment remained after overlap adjustment
    if (!@galn) {
	print "$q\t$ql\t---\n";
	next;
    }

    # set annotation for whole read
    my ($reg, $mref, $ori, $vdj, $vdji, $ual, $al) = SetAnnotation(\@galn, $ql);
    my $faln = OutputGroupAln(@galn);

    # treat as unaligned if alignment fraction is small (i.e., <$frac)
    $reg = "---" if $al < $ql * $frac;

    # output
    print "$q\t$ql\t$reg\t$mref\t$ori\t$vdj\t$faln\t$vdji\t$ual\t$al\n";
}
close $afh;

# output the remaining unaligned reads
my ($q, $qseq) = LoadOneFasta($qfh);

while ($q) {
    my $ql = length($qseq);
    print "$q\t$ql\t---\n";
    ($q, $qseq) = LoadOneFasta($qfh);
}
close $qfh;
