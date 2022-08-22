#!/usr/bin/perl -w
# usage : ExtractCDR3.pl vdjdelta fq

use strict;
use File::Basename qw(dirname);
use Cwd qw(abs_path);
use Getopt::Long qw(GetOptions);

# use TRIg's modules
use lib dirname(abs_path($0));
use Fastx qw(LoadThisFastq);


#################### set parameters ####################

my $species = "hsa";
my $gene = "trb";

GetOptions(
    "species=s" => \$species,
    "gene=s"    => \$gene,
    );

my $sg = "$species\_$gene";

my %codonaa = qw( 
    TTT F TCT S TAT Y TGT C TTC F TCC S TAC Y TGC C TTA L TCA S TAA * TGA * TTG L TCG S TAG * TGG W 
    CTT L CCT P CAT H CGT R CTC L CCC P CAC H CGC R CTA L CCA P CAA Q CGA R CTG L CCG P CAG Q CGG R 
    ATT I ACT T AAT N AGT S ATC I ACC T AAC N AGC S ATA I ACA T AAA K AGA R ATG M ACG T AAG K AGG R 
    GTT V GCT A GAT D GGT G GTC V GCC A GAC D GGC G GTA V GCA A GAA E GGA G GTG V GCG A GAG E GGG G 
    );


#################### set CDR3 position ####################

my %vjcdr3p;

open IN, "<$sg.cdr" || die "open $sg.cdr: $!\n";

while (<IN>) {
    my @a = split "\t"; chomp $a[-1];
    
    # skip if no CDR3 coordinate (e.g., for pseudogene)                                                                                                                                        
    next if $a[5] eq "---";
    
    $vjcdr3p{$a[1]} = $a[5];
}
close IN;


#################### load VDJ delta and extract CDR3 ####################

# set read file handle
open my $qfh, "<$ARGV[1]" || die "open $ARGV[1]: $!\n";

open IN, "<$ARGV[0]" || die "open $ARGV[0]: $!\n";

while( <IN> ) {
    my @a = split "\t"; chomp $a[-1];

    # if no alignment or non-regular
    if ($a[2] eq "---" || $a[2] == 0) {
	print "$a[0]\t$a[2]\t---\n";
	next;
    }

    # get the VJ alignments
    my $q = $a[0];
    my $o = $a[4];
    my ($v, $d, $j) = split ":", $a[5];
    my @galn = split " ", $a[6];
    my ($vi, $di, $ji) = split ",", $a[7];
    my @valn = split /\|/, $galn[$vi];
    my @jaln = split /\|/, $galn[$ji];

    # extract CDR3
    my ($qseq, $qqual) = LoadThisFastq($qfh, $q);
    my @cdr3c = ();
    my @cdr3q = ();
    my @cdr3aa = ();
    ExtractCDR3(\@valn, \@jaln, $qseq, $qqual, \@cdr3c, \@cdr3q, \@cdr3aa);

    # output CDR3
    print "$q\t$a[2]\t";
    print join("|", @cdr3c) . "\t";
    print join("|", @cdr3q) . "\t";
    print join("|", @cdr3aa) . "\n";
}
close IN;



#################### subroutines ####################


sub ExtractCDR3 {
    my ($pvaln, $pjaln, $qseq, $qqual, $pcdr3c, $pcdr3q, $pcdr3aa) = @_;

    # for all VJ pairs
    foreach my $valn (@{$pvaln}) {
	my ($ref, $ve, $vrr, $vqr, $vmmg, $vgap) = split ":", $valn;
	my ($v) = $ve =~ /(.+)\_e/;
	my ($vrs, $vre) = split "-", $vrr;
	my ($vqs, $vqe) = split "-", $vqr;
	
	foreach my $jaln (@{$pjaln}) {
	    my ($ref, $je, $jrr, $jqr, $jmmg, $jgap) = split ":", $jaln;
	    my ($j) = $je =~ /(.+)\_e/;
	    my ($jrs, $jre) = split "-", $jrr;
	    my ($jqs, $jqe) = split "-", $jqr;
	    
	    my $cdr3seq = "---";
            my $cdr3qual = "---";
	    my $cdr3aa = "---";
            
            # check if the CDR3 positions on the reference are available (i.e., not pseudogene)
            # and the positions are covered by the V and J alignments (take care of V30 on the minus strand)
            my $cdr3q = 0;
            if ($vjcdr3p{$v} && $vjcdr3p{$j}) {
                if ($vrs <= $vjcdr3p{$v} && $vjcdr3p{$v} <= $vre &&
		    $jrs <= $vjcdr3p{$j} && $vjcdr3p{$j} <= $jre) {
                    $cdr3q = 1;
                }
            }
            
            # if good, first get the CDR3 starting and ending positions on the query
            if ($cdr3q == 1) {
                my $vqp = AlignmentRpQp($vrs, $vre, $vqs, $vqe, $vgap, $vjcdr3p{$v});
                my $jqp = AlignmentRpQp($jrs, $jre, $jqs, $jqe, $jgap, $vjcdr3p{$j});
                
                # then get the CDR3 segment on the plus strand
                if ($jqs < $jqe) {
                    $cdr3seq = substr($qseq, $vqp - 1, $jqp - $vqp + 1);
                    $cdr3qual = substr($qqual, $vqp - 1, $jqp - $vqp + 1);
                } else {
                    $cdr3seq = substr($qseq, $jqp - 1, $vqp - $jqp + 1);
                    $cdr3qual = substr($qqual, $jqp - 1, $vqp - $jqp + 1);
                    $cdr3seq = ReverseComplement($cdr3seq);
                    $cdr3qual = reverse($cdr3qual);
                }

		# translate CDR3
		$cdr3aa = Translate($cdr3seq);
            }
            push(@{$pcdr3c}, "$v:$cdr3seq:$j");
            push(@{$pcdr3q}, $cdr3qual);
	    push(@{$pcdr3aa}, $cdr3aa);
        }
    }
}


sub AlignmentRpQp {
    my ($rs, $re, $q1, $q2, $gap, $rp) = @_;
    my $qp;

    # set alignment info
    my $o = $q1 < $q2 ? 1 : -1;
    my @gap = split ",", $gap; pop(@gap);

    # bs : block size (idea from blat)
    # rg : reference gap
    my @bs = ();
    my @rg = ();
    foreach my $g (@gap) {
        if ($g > 0) {
            push(@bs, $g - 1);
            push(@rg, 1);
        } else {
            push(@bs, -$g - 1);
            push(@rg, 0);
        }
    }
    push(@bs, $re - $rs + 1 - Total(@bs, @rg));

    # calculate query position
    my $rl = $rp - $rs;
    my $ri = -1;
    my $qi = -1;
    my $bi = 0;
    while ($ri < $rl) {
        my $d = $rl - $ri;
        if ($d <= $bs[$bi]) {
            $ri += $d;
            $qi += $d;
            last;
        } else {
            $ri += $bs[$bi];
            $qi += $bs[$bi];
            if ($rg[$bi] == 1) {
                $ri++;
            } else {
                $qi++;
            }
        }
        $bi++;
    }
    $qp = $q1 + $qi * $o;

    return $qp;
}


sub ReverseComplement {
    my $s = reverse($_[0]);
    $s =~ tr/acgtACGT/tgcaTGCA/;
    return $s;
}


sub Total {
    my $s = 0;
    foreach my $e (@_) { $s += $e; }
    return $s;
}


sub Translate {
    my $dna = $_[0];
    
    my @codon = unpack('(A3)*', $dna);
    my @aa = map { $codonaa{$_} ? $codonaa{$_} : "_" } @codon;
    my $aa = join "", @aa;
    
    return $aa;
}
