package VDJDelta;

use strict;
use Exporter qw(import);

our @EXPORT_OK = qw(SetScoringFunction SetVDJInfo LoadOneDelta GetOptimalSet OutputAln);


############################## set scoring function ##############################

# default values based on nucmer
# msc  : match score
# mmsc : mismatch score
# gsc  : gap score

my $msc  = 3;
my $mmsc = -7;
my $gsc  = -7;


############################## set VDJ and CDR3 information ##############################

my %rvdje;
my %vdjekv;


sub SetVDJInfo {
    my $fn = $_[0];

    # load and set VDJ info
    open IN, "<$fn" || die "open $fn: $!\n";

    while (<IN>) {
        my @a = split "\t"; chomp $a[-1];

        # er : exon range
        # e0 : only one exon
        # e1, e2, ... : more than one exon
        # ns : numbering start
        my @er = split ",", $a[4];
        my $ns = @er > 1 ? 1 : 0;
        
        # set VDJ exon info
        for (my $i = 0; $i < @er; $i++) {
            my $ei = $a[3] eq "+" ? $i + $ns : $#er + $ns - $i;
            my $vdje = "$a[1]\_$ei";
            push(@{$rvdje{$a[0]}}, $vdje);
            
            my ($s, $e) = $er[$i] =~ /(\d+)\.\.(\d+)/;
            $vdjekv{$vdje}{s} = $s;
            $vdjekv{$vdje}{e} = $e;
            $vdjekv{$vdje}{o} = $a[3];
            $vdjekv{$vdje}{vdj} = $a[1];
        }
    }
    close IN;
    
    $vdjekv{"TRBI_0"}{vdj} = "TRBI";
}


############################## define class and methods ##############################


sub new {
    my ($class, $aln) = @_;

    # ref     : reference
    # r/q s/e : reference/query start/end
    # ro      : orientation on reference
    # oqs/e   : oriented query start/end (i.e., along query)
    # qal     : query aligned length
    # mmg     : number of mismatches and gaps
    # gap     : gap positions
    # id      : alignment identity
    # sc      : alignment score
    # vdje    : VDJ exon
    # vdj     : VDJ
    # ge      : gene exon (e.g., V2, J0, C1)
    # go      : orientation on gene
    # rseg    : reference segment
    # qseg    : query segment
    my ($ref, $rs, $re, $qs, $qe, $ro, $oqs, $oqe, $qal, $mmg, $gap,
	$id, $sc, $vdje, $vdj, $ge, $go, $rseg, $qseg);

    # initialize if alignment info is provided    
    if ($aln) {

	# if in the delta format
	if ($aln =~ /\n/) {
	    
	    # ai : alingment info
	    ($ref, my $ai, my @gap) = split "\n", $aln;
	    ($rs, $re, $qs, $qe, $mmg) = split " ", $ai;
	    ($ro, $oqs, $oqe) = $qs < $qe ? ("+", $qs, $qe) : ("-", $qe, $qs);
	    $go = $ro;
	    $qal = $oqe - $oqs + 1;
	    $gap = join ",", @gap;
	    
	    # tg : total gaps
	    # qg : query gaps
	    my $tg = $#gap;
	    my $qg = scalar(grep($_ > 0, @gap));
	    my $mm = $mmg - $tg;
	    my $m = $qal + $qg - $mmg;
	    $id = $m / ($qal + $qg);
	    $sc = $m * $msc + $mm * $mmsc + $tg * $gsc;

	# if in the vdjdelta format
	} elsif ($aln =~ /:/) {

	    ($ref, $vdje, my $rse, my $qse, $mmg, $gap) = split ":", $aln;
	    ($rs, $re) = split "-", $rse;
	    ($qs, $qe) = split "-", $qse;
	    ($ro, $oqs, $oqe) = $qs < $qe ? ("+", $qs, $qe) : ("-", $qe, $qs);
	    $qal = $oqe - $oqs + 1;

	    if ($vdje ne "---") {
		if ($vdje eq "TRBI_0") {
		    $vdj = "TRBI";
		} else {
		    my @vdje = split "~", $vdje;
		    my @vdj = map { $vdjekv{$_}{vdj} } @vdje;
		    $vdj = join "~", @vdj;
		    $vdjekv{$vdje}{o} = $vdjekv{$vdje[0]}{o} if !$vdjekv{$vdje}{o};
		}
		my $g = $vdje eq "TRBI" ? "I" : substr($vdje, 3, 1);
		my $e = $vdje =~ /_(\d)/ ? $1 : 0;
		$ge = "$g$e";
	    }
	    $go = ($vdjekv{$vdje}{o} && $vdjekv{$vdje}{o} eq "-") ? ($ro eq "+" ? "-" : "+") : $ro;
	}
    }
    
    # define attributes
    my $self = {"ref" => $ref, "rs" => $rs, "re" => $re, "qs" => $qs, "qe" => $qe,
		    "ro" => $ro, "oqs" => $oqs, "oqe" => $oqe, "qal" => $qal, 
		    "mmg" => $mmg, "gap" => $gap, "id" => $id, "sc" => $sc,
		    "vdje" => $vdje, "vdj" => $vdj, "ge" => $ge, "go" => $go,
		    "rseg" => $rseg, "qseg" => $qseg};
    
    bless $self, $class;
}


sub out {
    my $self = $_[0];

    my $out = "$self->{ref}:$self->{vdje}";
    $out .= ":$self->{rs}-$self->{re}";
    $out .= ":$self->{qs}-$self->{qe}";
    $out .= ":$self->{mmg}:$self->{gap}";

    return $out;
}


sub LoadOneDelta {
    my ($afh, $r, $q) = @_;

    # load all alignments of a query
    my @aln;
    my $delta;
    while (<$afh>) {
        if (/^>(\S+) (\S+)/) {
	    if ($2 eq $q) {
		$r = $1;
		next;
	    } else {
		seek($afh, -length($_), 1);
		last;
	    }
	} else {
            $delta .= $_;
            if ($_ eq "0\n") {
		my $aln = VDJDelta->new("$r\n$delta");
                push(@aln, $aln);
		$delta = "";
            }
        }
    }

    return @aln;
}


sub GetOptimalSet {
    my @aln = @_;
    my @oaln;
    
    @aln = sort { $b->{sc} <=> $a->{sc} || $b->{id} <=> $a->{id} } @aln;

    # iterate until no alignment remains
    while (@aln) {

        # get current best alignment
        my @cba = shift(@aln);
        while (@aln && $cba[-1]->{sc} == $aln[0]->{sc} && $cba[-1]->{id} == $aln[0]->{id}) {
            push(@cba, shift(@aln));
        }
        push(@oaln, @cba);
        
        # filter alignment that overlaps the current best alignment
        foreach my $a (@cba) {
	    @aln = grep { $a->overlapq($_) == 0 } @aln;
        }
    }

    return @oaln;
}


sub overlapq {
    my ($self, $a) = @_;

    my $ol = $self->{oqs} < $a->{oqs} ? $self->{oqe} - $a->{oqs} + 1 : $a->{oqe} - $self->{oqs} + 1;
    my $q = $ol > $a->{qal} * 0.5 ? 1 : 0;

    return $q;
}


sub annotate {
    my $self = $_[0];

    $self->{vdje} = RangeVDJE($self->{ref}, $self->{rs}, $self->{re});
        
    # take care of alignment that spans more than one exon
    if (!$vdjekv{$self->{vdje}}{vdj}) {
	my @vdje = split '~', $self->{vdje};
	
	$vdjekv{$self->{vdje}}{s} = $vdjekv{$vdje[0]}{s};
	$vdjekv{$self->{vdje}}{e} = $vdjekv{$vdje[-1]}{e};
	$vdjekv{$self->{vdje}}{o} = $vdjekv{$vdje[0]}{o};
	$_ =~ s/(_\d)// for @vdje;
	$vdjekv{$self->{vdje}}{vdj} = join '~', @vdje;
    }
    $self->{vdj} = $vdjekv{$self->{vdje}}{vdj};

    # set gene exon (e.g., V2, J0)
    my $g = $self->{vdje} eq "TRBI" ? "I" : substr($self->{vdje}, 3, 1);
    $g = "C" if $self->{vdje} =~ /IGH[MEGAD]_/;
    my $e = $self->{vdje} =~ /_(\d)/ ? $1 : 0;
    $self->{ge} = "$g$e";
    $self->{go} = $self->{ro} eq $vdjekv{$self->{vdje}}{o} ? "+" : "-" if $self->{vdj} ne "TRBI";
}


sub setrefqryseg {
    my ($self, $prseq, $pqseq) = @_;

    $self->{rseg} = substr($prseq->{$self->{ref}}, $self->{rs} - 1, $self->{re} - $self->{rs} + 1);
    $self->{qseg} = substr(${$pqseq}, $self->{oqs} - 1, $self->{oqe} - $self->{oqs} + 1);
}


sub RangeVDJE {
    my ($r, $s, $e) = @_;

    my @vdje = @{$rvdje{$r}};
    my $si = 0;
    my $ei = scalar(@vdje) - 1;
    my @annot = ();

    # if the range start is after start of the last vdj
    if ($s > $vdjekv{$vdje[$ei]}{s}) {
        my $annot = $s > $vdjekv{$vdje[$ei]}{e} ? "TRBI" : $vdje[$ei];
        return $annot;
    }

    # locate the range start using dichotomy
    my $mi = int(($si + $ei) / 2);
    while ($mi != $si) {
        if ($s < $vdjekv{$vdje[$mi]}{s}) {
            $ei = $mi;
        } else {
            $si = $mi;
        }
        $mi = int(($si + $ei) / 2);
    }
    
    # get vdj exon from the starting location
    if ($s <= $vdjekv{$vdje[$si]}{e}) {
        push(@annot, $vdje[$si]);
    }
    while ($vdje[++$si]) {
        last if $e < $vdjekv{$vdje[$si]}{s};
        push (@annot, $vdje[$si]);
    }

    # an intergenic region if no VDJ exon in the region
    @annot = ("TRBI") if !@annot;

    return join "~", @annot;
}


sub reset_aln {
    my ($self, $delta) = @_;

    my ($rse, $qse, $mmg, $gap) = split ":", $delta;
    my ($rs, $re) = split "-", $rse;
    my ($qs, $qe) = split "-", $qse;

    $self->{rs} = $rs;
    $self->{re} = $re;
    $self->{qs} = $qs;
    $self->{qe} = $qe;
    $self->{mmg} = $mmg;
    $self->{gap} = $gap;
    $self->{qal} = abs($qe-$qs) + 1;
}


sub OutputAln {
    my @aln = @_;

    my @out;
    foreach my $a (@aln) {
        push(@out, $a->out());
    }
    my $out = join " ", @out;

    return $out;
}


1;
