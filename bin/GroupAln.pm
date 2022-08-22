package GroupAln;

use strict;
use Exporter qw(import);

our @EXPORT_OK = qw(SetScoringFunction FilterAlignment AdjustOverlap SetAnnotation OutputGroupAln);


############################## set scoring function ##############################

# default values based on nucmer
# msc  : match score
# mmsc : mismatch score
# gsc  : gap score
# exb  : exon bonus

my $msc  = 3;
my $mmsc = -7;
my $gsc  = -7;
my $exb  = 1;
my $minmatch = 15;

sub setparas {
    $minmatch = $_[0]; 
}
    
my %bbsc;

sub SetScoringFunction {
    my @base = ("A", "C", "G", "T");

    foreach my $b1 (@base) {
	foreach my $b2 (@base) {
	    if ($b1 eq $b2) {
		$bbsc{$b1}{$b2} = $msc + $exb;
		$bbsc{lc($b1)}{$b2} = $msc;
	    } else {
		$bbsc{$b1}{$b2} = $mmsc;
		$bbsc{lc($b1)}{$b2} = $mmsc;
	    }
	}
    }

    foreach my $b (@base) {
	$bbsc{$b}{'-'} = $gsc;
	$bbsc{lc($b)}{'-'} = $gsc;
	$bbsc{'-'}{$b} = $gsc;
	$bbsc{'-'}{'-'} = 0;
    }
}


############################## define class and methods ##############################


sub new {
    my ($class, @aln) = @_;

    # n : number of alignments
    my $n = scalar(@aln);

    # refer to VDJDelta.pm for the following definitions
    my $ref  = UniqueJoin("ref", @aln);
    my $vdje = UniqueJoin("vdje", @aln);
    my $vdj  = UniqueJoin("vdj", @aln);
    my $ge   = UniqueJoin("ge", @aln);
    my $ro   = UniqueJoin("ro", @aln);
    my $go   = UniqueJoin("go", @aln);
    my $rs   = UniqueJoin("rs", @aln);
    my $re   = UniqueJoin("re", @aln);
    my $qs   = UniqueJoin("qs", @aln);
    my $qe   = UniqueJoin("qe", @aln);
    my $oqs  = UniqueJoin("oqs", @aln);
    my $oqe  = UniqueJoin("oqe", @aln);
    my @qal  = sort { $a <=> $b } map($_->{qal}, @aln);
    my $qal  = $qal[0];

    # cxx : consistency in xx 
    # cbq : for broken alignment
    # cmq : for missing alignment
    my $cref  = $ref =~ /\|/ ? 0 : 1;
    my $cge   = $ge =~ /\|/ ? 0 : 1;
    my $cro   = $ro =~ /\|/ ? 0 : 1;
    my $cgo   = $go =~ /\|/ ? 0 : 1;
    my $coqse = $oqs =~ /\|/ || $oqe =~ /\|/? 0 : 1;
    my $cbq   = $cref * $cro * $coqse;
    my $cmq   = $cref * $cgo * $cge * $coqse;
    
    # Majority of inconsistent group alignments are INTs, which can have different query starts/ends 
    # and/or orientation. Some are different V alignments with slightly different query starts/ends 
    # and/or orientation. Few have V and J alignments.
=pod
    my $cq = $cge * $cro * $cgo * $coqse;
    if ($cq == 0) {
	print "inconsistent group alignemnt:";
	foreach my $a (@aln) {
	    print " " . $a->out();
	}
	print "\n";
    }
=cut

    # set attributes
    my $self = {"aln" => \@aln, "n" => $n, "ref" => $ref, "vdje" => $vdje, "vdj" => $vdj, "ge" => $ge, "ro" => $ro, "go" => $go,
		    "rs" => $rs, "re" => $re, "qs" => $qs, "qe" => $qe, "oqs" => $oqs, "oqe" => $oqe, "qal" => $qal, 
		    "ref" => $ref, "cge" => $cge, "cro" => $cro, "cgo" => $cgo, "coqse" => $coqse, "cbq" => $cbq, "cmq" => $cmq};
    
    bless $self, $class;
}


sub UniqueJoin {
    my ($key, @aln) = @_;

    my @uak = Unique(map($_->{$key}, @aln));

    return join("|", @uak);
}


sub Unique {
    my %seen = ();

    return grep(!$seen{$_}++, @_);
}


sub out {
    my $self = $_[0];

    my @aln = @{$self->{aln}};
    $_ = $_->out() for @aln;

    return join "|", @aln;
}


sub FilterAlignment {
    my @oaln = @_;
    my @galn;
    
    # filter short (<30 bp) intergenic alignments
    @oaln = grep { !($_->{vdje} eq "INT" && $_->{qal} < 30) } @oaln;
    
    # exit if no alignment remains
    return @galn if !@oaln;

    # check number of references
    my @ref = Unique(map($_->{ref}, @oaln));
    
    # group overlapping alignments based on query region
    my @sa = sort { $a->{oqs} <=> $b->{oqs} } @oaln;
    my @ga = ($sa[0]);
    my $g = 0;
    for (my $i = 1; $i < @sa; $i++) {

	# if in the same group
        if ($sa[$i-1]->overlapq($sa[$i]) == 1) {
            push(@ga, $sa[$i]);

	# if in a different group
        } else {

	    # filter intergenic alignment in a group if there exist 
	    # both an exonic and non-intergenic alignment
	    @ga = FilterINTwithExon(@ga);
	    push(@galn, GroupAln->new(@ga));
            @ga = ($sa[$i]);
        }
    }
    @ga = FilterINTwithExon(@ga);
    push(@galn, GroupAln->new(@ga));

    # check consistent reference for all groups of alignments
    my $srq = 1;
    $srq *= $_->{cref} for @galn;
    
    # if not single reference, filter based on total alignment length on the reference
    if ($srq == 0) {
	my %rqal;
	foreach my $ga (@galn) {
	    $rqal{$_->{ref}} += $_->{qal} for @{$ga->{aln}};
	}
	foreach my $ga (@galn) {
	    next if $ga->{cref} == 1;
	    
	    my @a = sort { $rqal{$b->{ref}} <=> $rqal{$a->{ref}} } @{$ga->{aln}};
	    @a = grep { $rqal{$_->{ref}} == $rqal{$a[0]->{ref}} } @a;
	    $ga = GroupAln->new(@a) if @a < $ga->{n};
	}
    }

    # filter ambiguous V alignments (e.g., eqaully well alignments to TRBV12-3 and TRBV12-4) 
    if (@galn > 1) {

	# find V and ambiguous V alignments
	my %vin;
	for (my $i = 0; $i < @galn; $i++) {
	    next if $galn[$i]->{vdj} !~ /V/;
	    $vin{$i} = $galn[$i]->{n};
	}
	my @vi = keys %vin;
	my @avi = grep { $vin{$_} > 1 } @vi;
	
	# filter ambiguous V alignment based on the unique one
	if (@avi) {
	    
	    # for each ambiguous V alignment
	    foreach my $i (@avi) {
		
		# skip if no chance to reduce the ambiguity
		my @uvi = grep { $vin{$_} < $vin{$i} } @vi;
		next if !@uvi;
		
		# for each possible reduction
		foreach my $uvi (@uvi) {
		    my @ga;
		    
		    # get the reduced alignment
		    if ($vin{$uvi} == 1) {
			@ga = grep { $_->{vdj} eq $galn[$uvi]->{vdj} } @{$galn[$i]->{aln}};
		    } else {
			my %vq;
			$vq{$_->{vdj}} = 1 for @{$galn[$uvi]->{aln}};
			@ga = grep { $vq{$_->{vdj}} } @{$galn[$i]->{aln}};
		    }
		    
		    # if non-null reduced alignment, reset the alignment 
		    if (0 < @ga && @ga < $vin{$i}) {
			$galn[$i] = GroupAln->new(@ga);
			$vin{$i} = 1;
			push(@vi, $i);
		    }
		}
	    }
	}
    }

    # calculate total alignment length of each V, D, J gene for identifying potential CDR3 alignment
    # i.e., to avoid mistaking a broken alignment (thus the summed VDJ length is long) as CDR3. 
    my %vdjqal;
    foreach my $ga (@galn) {
	$vdjqal{$ga->{vdj}} += $ga->{qal} if $ga->{vdj} ne "INT";
    }

    # identify potentail CDR3 segment, i.e., flanked by V2 and J of the same reference
    my $pcdr3i = -1;
    for (my $i = 1; $i < $#galn; $i++) {
	
        # skip if the alignment is good, i.e., D (should be there),
	# J (for concatenations of two Js), and long (60 bp is about the maximal of CDR3) 
        next if $galn[$i]->{ge} eq "D0";
        next if $galn[$i]->{ge} eq "J0" && $galn[$i]->{qal} >= 30;
        next if $galn[$i]->{qal} >= 60;

	# if flanked by V2 and J0 and shorter than the total alignment length of the V and J
       	if ($galn[$i-1]->{ge} eq "V2" && $galn[$i+1]->{ge} eq "J0" && $galn[$i-1]->{ref} eq $galn[$i+1]->{ref}) {

	    # short intergenic alignment is likely CDR3
	    if ($galn[$i]->{vdj} eq "INT") {
		$pcdr3i = $i if $galn[$i]->{qal} <= $vdjqal{$galn[$i-1]->{vdj}};

	    # shorter V alignment is also likely CDR3 
	    } elsif($vdjqal{$galn[$i]->{vdj}} < $vdjqal{$galn[$i-1]->{vdj}}) {
		$pcdr3i = $i; 
	    }
	} elsif ($galn[$i-1]->{ge} eq "J0" && $galn[$i+1]->{ge} eq "V2" && $galn[$i-1]->{ref} eq $galn[$i+1]->{ref}) {

	    if ($galn[$i]->{vdj} eq "INT") {
		$pcdr3i = $i if $galn[$i]->{qal} <= $vdjqal{$galn[$i+1]->{vdj}};

	    } elsif($vdjqal{$galn[$i]->{vdj}} < $vdjqal{$galn[$i+1]->{vdj}}) {
		$pcdr3i = $i;
	    }
	}
    }

    splice(@galn, $pcdr3i, 1) if $pcdr3i != -1;

    # return filtered group alignment
    return @galn;
}


sub FilterINTwithExon {
    my @ga = @_;
    
    my @non = grep { $_->{vdj} ne "INT" } @ga;
    my @int = grep { $_->{vdj} eq "INT" } @ga;
    @ga = @non if @non && @int;

    return @ga;
}


sub AdjustOverlap {
    my ($prseq, $pqseq, $pgaln) = @_;

    # check overlap for each pair of neighboring alignments
    my %iolq;
    my %iol;
    my @galn = @{$pgaln};
    for (my $i = 1; $i < @galn; $i++) {

	# skip if not on the same reference or not consistent reference
	next if $galn[$i-1]->{ref} ne $galn[$i]->{ref} || $galn[$i-1]->{cref} * $galn[$i]->{cref} == 0;
	
	# skip if no consistent query aligned positions
	next if $galn[$i-1]->{coqse} * $galn[$i]->{coqse} == 0;
	
	# skip if no overlap
	my $ol = $galn[$i-1]->{oqe} - $galn[$i]->{oqs} + 1;
	next if $ol <= 0;

	$iolq{$i-1} = 1;
	$iolq{$i} = 1;
	$iol{$i} = $ol;
    }

    # end if no overlap
    return if !%iol;

    # set aligned segments if there is an overlap
    foreach my $i (keys %iolq) {
	foreach my $aln (@{$galn[$i]->{aln}}) {
	    $aln->setrefqryseg($prseq, $pqseq);
	}
    }
	
    # adjust alignments
    for (my $i = 1; $i < @galn; $i++) {
	next if !$iol{$i};

	my $ol = $iol{$i};
	my %a1a2kv;
	my @a1 = @{$galn[$i-1]->{aln}};
	my @a2 = @{$galn[$i]->{aln}};
	foreach my $a1 (@a1) {
	    foreach my $a2 (@a2) {
		($a1a2kv{$a1}{$a2}{aa1}, $a1a2kv{$a1}{$a2}{aa2}) = AdjustAlignments($a1, $a2, $ol);
	    }
	}

	# skip if the adjusted alignments are not consistent
        my $cq = 1;
        foreach my $a1 (@a1) {
            next if @a2 == 1;
            my @aa1 = map($a1a2kv{$a1}{$_}{aa1}, @a2);
            @aa1 = Unique(@aa1);
            if (@aa1 > 1) {
                $cq = 0;
                last;
            }
        }
        if ($cq == 1) {
            foreach my $a2 (@a2) {
                next if @a1 == 1;
                my @aa2 = map($a1a2kv{$_}{$a2}{aa2}, @a1);
                @aa2 = Unique(@aa2);
                if (@aa2 > 1) {
                    $cq = 0;
                    last;
                }
            }
        }
        next if $cq == 0;

	# reset the adjusted alignment
	foreach my $a1 (@{$galn[$i-1]->{aln}}) {
	    $a1->reset_aln($a1a2kv{$a1}{$a2[0]}{aa1});
	}
	foreach my $a2 (@{$galn[$i]->{aln}}) {
	    $a2->reset_aln($a1a2kv{$a1[0]}{$a2}{aa2});
	}
	$galn[$i-1]->reset_aln(@{$galn[$i-1]->{aln}});
	$galn[$i]->reset_aln(@{$galn[$i]->{aln}});
    }

    # filter short adjusted alignment
    my $k = 0;
    for (my $i = 0; $i < @galn; $i++) {
	if ($galn[$i]->{qal} < $minmatch || ($galn[$i]->{vdje} eq "INT" && $galn[$i]->{qal} < 30)) {
	    splice(@galn, $i, 1);
	    $i--;
	    $k++;
	} else {
	    $pgaln->[$i] = $galn[$i];
	}
    }
    foreach (1..$k) {
	pop(@{$pgaln});
    }
}


sub AdjustAlignments {
    my ($a1, $a2, $ol) = @_;

    # get aligned segments in the overlap and make the two aligments consistent
    my ($ra1ol, $qa1ol) = OverlappingAlignment($a1, "t", $ol);
    my ($ra2ol, $qa2ol) = OverlappingAlignment($a2, "h", $ol);
    ($ra1ol, $qa1ol, $ra2ol, $qa2ol) = MakeConsistent($ra1ol, $qa1ol, $ra2ol, $qa2ol);

    # cut overlapping alignment
    my ($cra1, $cqa1, $cra2, $cqa2) = CutOverlappingAlignment($ra1ol, $qa1ol, $ra2ol);
    
    # modify alignment info
    my $aa1 = ModifyAlignment($a1, 't', $cra1, $cqa1);
    my $aa2 = ModifyAlignment($a2, 'h', $cra2, $cqa2);
    
    return ($aa1, $aa2);
}


sub OverlappingAlignment {
    my ($aln, $qht, $ol) = @_;
    my $raseg;
    my $qaseg;

    # get aligned segment oriented along the query
    my $rseg = $aln->{ro} eq "+" ? $aln->{rseg} : ReverseComplement($aln->{rseg});
    my $qseg = $aln->{qseg};
    my @gap = split ",", $aln->{gap}; pop( @gap );
    if ( ($qht eq "h" && $aln->{ro} eq "-") || ($qht eq "t" && $aln->{ro} eq "+") ) {
	@gap = ReverseGap($aln->{qal}, @gap);
    }

    # get the overlapping aligned segments
    if ($qht eq "h") {
	($raseg, $qaseg) = LeftEndAlignment($rseg, $qseg, $ol, @gap);
    } else {
	($raseg, $qaseg) = RightEndAlignment($rseg, $qseg, $ol, @gap);
    }
	
    return ($raseg, $qaseg);
}


sub ReverseComplement {
    my $s = reverse($_[0]);
    $s =~ tr/acgtACGT/tgcaTGCA/;
    return $s;
}


sub ReverseGap {
    my ($qal, @gap) = @_;

    return @gap if !@gap;

    # bs: block size (including the gap)
    # sign: - for ref and + for query
    my @bs = map(abs($_), @gap);
    my @sign = map { $_ > 0 ? 1 : -1 } @gap;
    my $qg = grep { $_ > 0 } @gap;
    my $qtbs = $qal - Total(@bs) + $qg + 1;

    push(@bs, $qtbs);
    @bs = reverse(@bs); pop(@bs);
    @sign = reverse(@sign);
    for (my $i = 0; $i < @bs; $i++) {
	$gap[$i] = $sign[$i] * $bs[$i]; 
    }

    return @gap;
}


sub Total {
    my $s = 0;
    $s += $_ for @_;
    return $s;
}


sub LeftEndAlignment {
    my ($rseg, $qseg, $ol, @gap) = @_;
    my $rhseg = "";
    my $qhseg = "";
    
    # if no gap, simply get the overlapped segments
    if (!@gap) {
	$rhseg = substr($rseg, 0, $ol);
	$qhseg = substr($qseg, 0, $ol);
	return ($rhseg, $qhseg);
    }
    
    # if there is a gap
    my $ri = 0;
    my $qi = 0;
    my $gi = 0;
    while ($qi < $ol) {

	# if not the last block yet
	if ($gi < @gap) {
	    my $bs = abs($gap[$gi]);
	
	    # if the gap is on the query
	    if ($gap[$gi] > 0) {
		
		# if the overlap ends within this block
		if ($ol <= ($qi + $bs - 1)) {
		    $rhseg .= substr($rseg, $ri, $ol - $qi);
		    $qhseg .= substr($qseg, $qi, $ol - $qi);
		    $ri += $ol - $qi;
		    $qi += $ol - $qi;
		    
		# if the overlap goes further
		} else {
		    $rhseg .= substr($rseg, $ri, $bs);
		    $qhseg .= substr($qseg, $qi, $bs - 1) . '-';
		    $ri += $bs;
		    $qi += $bs - 1;
		    $gi++;
		}
		
	    # if the gap is on the reference
	    } else {
		
		# if the overlap ends within this block
		if ($ol <= ($qi + $bs)) {
		    
		    # if the overlap extends to the end of this block
		    if ($ol == ($qi + $bs)) {
			$rhseg .= substr($rseg, $ri, $bs - 1) . '-';
			$ri += $bs - 1;
		    } else {
			$rhseg .= substr($rseg, $ri, $ol - $qi);
			$ri += $ol - $qi;
		    }
		    $qhseg .= substr($qseg, $qi, $ol - $qi);
		    $qi += $ol - $qi;
		    
	        # if the overlap goes further
		} else {
		    $rhseg .= substr($rseg, $ri, $bs - 1) . '-';
		    $qhseg .= substr($qseg, $qi, $bs);
		    $ri += $bs - 1;
		    $qi += $bs;
		    $gi++;
		}
	    }

	# if the last block    
	} else {
	    $rhseg .= substr($rseg, $ri, $ol - $qi);
	    $qhseg .= substr($qseg, $qi, $ol - $qi);
	    $ri += $ol - $qi;
	    $qi += $ol - $qi;
	}
    }

    return ($rhseg, $qhseg);
}


sub RightEndAlignment {
    my ($rseg, $qseg, $ol, @gap) = @_;
    my $rtseg = "";
    my $qtseg = "";
    
    # if no gap, simply get the overlapped segments
    if (!@gap) {
	$rtseg = substr($rseg, -$ol);
	$qtseg = substr($qseg, -$ol);
	return ($rtseg, $qtseg);
    }

    # if there is a gap
    my $ri = 0;
    my $qi = 0;
    my $gi = 0;
    while (-$ol < $qi) {
	
	# if not the last block yet
	if ($gi < @gap) {
	    my $bs = abs($gap[$gi]);
	
	    # if the gap is on the query
	    if ($gap[$gi] > 0) {
		
		# if the overlap ends within this block
		if (($qi - $bs + 1) < -$ol) {
		    my $d = $qi + $ol;
		    $ri -= $d;
		    $qi -= $d;
		    $rtseg = substr($rseg, $ri, $d) . $rtseg;
		    $qtseg = substr($qseg, $qi, $d) . $qtseg;
		    
		# if the overlap goes further
		} else {
		    $ri -= $bs;
		    $qi -= $bs - 1;
		    $rtseg = substr($rseg, $ri, $bs) . $rtseg;
		    $qtseg = '-' . substr($qseg, $qi, $bs - 1) . $qtseg;
		    $gi++;
		}

	    # if the gap is on the reference
	    } else {
		
		# if the overlap ends within this block
		if (($qi - $bs) <= -$ol) {
		    
		    # if the overlap extends to the end of this block
		    if (($qi - $bs) == -$ol) {
			$ri -= $bs - 1;
			$rtseg = '-' . substr($rseg, $ri, $bs - 1) . $rtseg;
		    } else {
			my $d = $qi + $ol;
			$ri -= $d;
			$rtseg = substr($rseg, $ri, $d) . $rtseg;
		    }
		    my $d = $qi + $ol;
		    $qi -= $d;
		    $qtseg = substr($qseg, $qi, $d) . $qtseg;
		    
		# if the overlap goes further
		} else {
		    $ri -= $bs - 1;
		    $qi -= $bs;
		    $rtseg = '-' . substr($rseg, $ri, $bs - 1) . $rtseg;
		    $qtseg = substr($qseg, $qi, $bs) . $qtseg;
		    $gi++;
		}
	    }

	# if the last block
	} else {
	    my $d = $qi + $ol;
	    $ri -= $d;
	    $qi -= $d;
	    $rtseg = substr($rseg, $ri, $d) . $rtseg;
	    $qtseg = substr($qseg, $qi, $d) . $qtseg;
        }
    }

    return ($rtseg, $qtseg);
}


sub MakeConsistent {
    my ($ra1ol, $qa1ol, $ra2ol, $qa2ol) = @_;

    if( $qa1ol ne $qa2ol ) {
        my %gpl1 = GapPositionLength($qa1ol);
        my %gpl2 = GapPositionLength($qa2ol);
        my %agp1 = GapComplement(\%gpl2, \%gpl1);
        my %agp2 = GapComplement(\%gpl1, \%gpl2);
        ($ra1ol, $qa1ol) = AddAlignmentGap($ra1ol, $qa1ol, \%agp1) if %agp1;
        ($ra2ol, $qa2ol) = AddAlignmentGap($ra2ol, $qa2ol, \%agp2) if %agp2;
    }
    
    return ($ra1ol, $qa1ol, $ra2ol, $qa2ol);
}


sub GapPositionLength {
    my $aseg = $_[0];
    my %gpl = ();
    
    my @gap = $aseg =~ /(-+)/g;
    my $p = -1;
    foreach my $g (@gap) {
        $p = index($aseg, $g, $p+1);
        $gpl{$p} = length($g);
    }

    return %gpl;
}


sub GapComplement {
    my %gpl1 = %{$_[0]};
    my %gpl2 = %{$_[1]};
    my %agpl;
 
    # get gaps only in aln1 or longer in aln1
    foreach my $p1 (keys %gpl1) {
        if (!$gpl2{$p1}) {
            $agpl{$p1} = $gpl1{$p1};
        } else {
            my $d = $gpl1{$p1} - $gpl2{$p1};
            $agpl{$p1} = $d if $d > 0;
        }
    }

    return %agpl;
}


sub AddAlignmentGap {
    my ($raseg, $qaseg, $pagpl) = @_;
    my %agpl = %{$pagpl};

    my @agp = sort { $b <=> $a } keys %agpl;
    foreach my $gp (@agp) {
        my $g = '-' x $agpl{$gp};
        $raseg = substr($raseg, 0, $gp) . $g. substr($raseg, $gp);
        $qaseg = substr($qaseg, 0, $gp) . $g. substr($qaseg, $gp);
    }
    
    return ($raseg, $qaseg);
}


sub CutOverlappingAlignment {
    my ($ra1, $qa, $ra2) = @_;
    
    my @ra1 = split "", $ra1;
    my @qa  = split "", $qa;
    my @ra2 = split "", $ra2;

    # get score difference at each cut point
    my @scd = (0);
    for (my $i = 0; $i < @qa; $i++) {
        my $d = $scd[-1] + $bbsc{$ra1[$i]}{$qa[$i]} - $bbsc{$ra2[$i]}{$qa[$i]};
        push(@scd, $d);
    }
    
    # get maximal cut point
    my $max = 0;
    my $maxi = 0;
    for (my $i = 1; $i < @scd; $i++) {
        if ($scd[$i] > $max) { 
            $max = $scd[$i];
            $maxi = $i;
        }
    }
    
    # alignment that is cut
    my $cra1 = join "", @ra1[$maxi..$#ra1];
    my $cra2 = join "", @ra2[0..($maxi-1)];
    my $cqa1 = join "", @qa[$maxi..$#ra1];
    my $cqa2 = join "", @qa[0..($maxi-1)];

    return ($cra1, $cqa1, $cra2, $cqa2);
}


sub ModifyAlignment {
    my ($aln, $qht, $cra, $cqa) = @_;
    my $maln;
    
    my @gap = split ",", $aln->{gap}; pop( @gap );
    my @ai = ($aln->{rs}, $aln->{re}, $aln->{qs}, $aln->{qe}, $aln->{mmg});

    my ($cmm, $crg, $cqg, $chb, $ctb) = AlignmentMismatchGap($cra, $cqa);
    my $crb = length($cra) - $crg;
    my $cqb = length($cqa) - $cqg;
    my $cmmg = $cmm + $crg + $cqg;
    if ($aln->{ro} eq '+') {
        if ($qht eq 't') {
            $ai[1] -= $crb;
            $ai[3] -= $cqb;
            foreach my $g (1..($crg+$cqg)) {
                pop( @gap );
            }
        } else {
            $ai[0] += $crb;
            $ai[2] += $cqb;
            foreach my $g (1..($crg+$cqg)) {
                shift( @gap );
            }
            if (@gap) {
                my $l = abs( $gap[0] ) - $ctb;
                $gap[0] = $gap[0] > 0 ? $l : -$l;
            }
        }
    } else {
        if ($qht eq 't') {
            $ai[0] += $crb;
            $ai[2] -= $cqb;
            foreach my $g (1..($crg+$cqg)) {
                shift( @gap );
            }
            if (@gap) {
                my $l = abs( $gap[0] ) - $chb;
                $gap[0] = $gap[0] > 0 ? $l : -$l;
            }
        } else {
            $ai[1] -= $crb;
            $ai[3] += $cqb;
            foreach my $g (1..($crg+$cqg)) {
                pop( @gap );
            }
        }
    }
    $ai[4] -= $cmmg;
    $ai[5] -= $cmmg;
    push(@gap, 0);
    
    $maln = "$ai[0]-$ai[1]:$ai[2]-$ai[3]:$ai[4]:";
    $maln .= join(",", @gap);

    return $maln;
}


sub AlignmentMismatchGap {
    my ($ra, $qa) = @_;
    my @ra = split "", uc($ra);
    my @qa = split "", $qa;

    my $mm = 0;
    my @rg = ();
    my @qg = ();
    for (my $i = 0; $i < @ra; $i++) {
        next if $ra[$i] eq $qa[$i];
        if ($ra[$i] eq '-') {
            push(@rg, $i);
        } elsif ($qa[$i] eq '-') {
            push(@qg, $i);
        } else {
            $mm++;
        }
    }

    my $al = length($ra);
    my $rg = scalar(@rg);
    my $qg = scalar(@qg);
    my @g = sort { $a <=> $b } (@rg, @qg);
    my $hb = @g ? $g[0] : $al;
    my $tb = @g ? $al - $g[-1] - 1 : $al;

    return ($mm, $rg, $qg, $hb, $tb);
}


sub reset_aln {
    my ($self, @aln) = @_;
    
    $self->{rs}  = UniqueJoin("rs", @aln);
    $self->{re}  = UniqueJoin("re", @aln);
    $self->{qs}  = UniqueJoin("qs", @aln);
    $self->{qe}  = UniqueJoin("qe", @aln);
    $self->{oqs} = UniqueJoin("oqs", @aln);
    $self->{oqe} = UniqueJoin("oqe", @aln);
    $self->{qal} = UniqueJoin("qal", @aln);
}


sub SetAnnotation {
    my ($pgaln, $ql) = @_;
    my @galn = @{$pgaln};

    my $reg = -1;
    my $mref = "---";
    my $ori = "o";
    my ($vi, $di, $ji) = (-1, -1, -1);
    my ($v, $d, $j) = ("---", "---", "---");

    my @ref;
    my @ge;
    my @go;
    my @vdj;
    push(@ref, $_->{ref}) for @galn;
    push(@ge, $_->{ge}) for @galn;
    push(@go, $_->{go}) for @galn;
    push(@vdj, $_->{vdj}) for @galn;
    my @uref = Unique(@ref);
    my $ge = join "", @ge;
    my $rge = join "", reverse(@ge);
    
    # if possibly regular
    if ($ge =~ /V[02](D0)?J0/ || $rge =~ /V[02](D0)?J0/) {
	
	# determine orientation
	if ($ge =~ /V[02](D0)?J0/) {
	    $ori = "+";
	    
	} else {
	    $ori = "-";
	    $ge = $rge;
	    @ref = reverse(@ref);
	    @ge = reverse(@ge);
	    @go = reverse(@go);
	    @vdj = reverse(@vdj); 
	}
       
	# determine the VDJ index
	my ($vdj) = $ge =~ /(V[02](D0)?J0)/;
	my $i = index($ge, $vdj);
	$vi = $i / 2;
	$v = $vdj[$vi];
	if ($vdj =~ /D0/) {
	    $di = $vi + 1; 
	    $ji = $vi + 2;
	    $d = $vdj[$di];
	    $j = $vdj[$ji];
	} else {
	    $ji = $vi + 1;
	    $j = $vdj[$ji];
	}

	# confirm that the V and J genes are on the same reference
	if ($ref[$vi] eq $ref[$ji]) {
	    
	    # confirm consistency of orientation
	    my @o = Unique(@go[$vi..$ji]);
	    if (@o == 1 && $o[0] eq $ori) {
		$reg = 1;
		$mref = $ref[$vi];
		
		# if regular, check if strictly regular, i.e., follow the VDJ(C) order, 
		# all alignments on the same gene strand, with only one V and J, and
		# on the same reference
		if ($reg == 1 && $ge =~ /^(V1)?V[02](D0)?J0(C1)?$/) {
		    my @o = Unique(@go);
		    my @v = Unique(@vdj[0..$vi]);
		    $reg = 2 if @o == 1 && @v == 1 && $j !~ /~/ && @uref == 1;
		}
	    }
	}

	# reset if not regular
	if ($reg == -1) {
	    if($ori eq "-") {
		@ref = reverse(@ref);
		@ge = reverse(@ge);
		@go = reverse(@go);
		@vdj = reverse(@vdj);
		$ge = join "", @ge;
		$ori = "o";
	    }
	    $vi = -1; $di = -1; $ji = -1;
	    $v = "---"; $d = "---"; $j = "---";
	}
    }

    # if not regular
    if ($reg == -1) {
	$reg = 0;

	# set orientation if not set yet (i.e., not even reg=1)
	if ($ori eq "o") {

	    # if there is a JC, set orientation based on the JC
	    if ($ge =~ /(J0C1)/) {
		my $jc = $1;
		my $i = index($ge, $jc);
		$ji = $i / 2;
		
		# confirm the consistency of the J and C orientation and reference
		if ($go[$ji] eq "+" && $go[$ji+1] eq "+" && $ref[$ji] eq $ref[$ji+1]) {
		    $ori = "+";
		    $j = $vdj[$ji];
		    $mref = $ref[$ji];
		} else {
		    $ji = -1;
		}

	    } elsif ($rge =~ /(J0C1)/) {
		my $jc = $1;
		my $i = index($rge, $jc);
		$ji = $#ge - $i/2;
		
		if ($go[$ji] eq "-" && $go[$ji-1] eq "-" && $ref[$ji] eq $ref[$ji-1]) {
		    $ori = "-";
		    $j = $vdj[$ji];
		    $mref = $ref[$ji];
		} else {
		    $ji = -1;
		}
	    }

	    # if the orientation is still not set
	    if ($ori eq "o") {

		# find the reference with a longer alignment length
		my %rqal;
		$rqal{$_->{ref}} += $_->{qal} for @galn;
		my @r = sort { $rqal{$b} <=> $rqal{$a} } keys %rqal;
		$mref = $r[0];
		my @rgaln = grep { $_->{ref} eq $mref } @galn;

		my @rge = map($_->{ge}, @rgaln);
		my @rgo = map($_->{go}, @rgaln);
		my %rgei;
		for (my $i = 0; $i < @rge; $i++) {
		    push(@{$rgei{$rge[$i]}}, $i);
		}

		# set it based on J, V2, V0, V1, D0, INT, C in order
		my @g = ("J0", "V2", "V0", "V1", "D0", "I0", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9");
		my @ro = map { $rgei{$_} ? $rgo[$rgei{$_}[0]] : "" } @g;
		@ro = grep($_, @ro);
		if (@ro) {
		    $ori = $ro[0];
		} else {
		    $ori = $rgaln[0]->{aln}[0]->{go};
		}
	    }
	}

	#if (!$ori) {
	#    my $out = OutputGroupAln(@galn);
	#    print "error: $out\n";
	#}

	# set the VDJ index
	if ($ori eq "-") {
	    @ref = reverse(@ref);
	    @ge = reverse(@ge);
	    @vdj = reverse(@vdj);
	}
	my %gei;
	for (my $i = 0; $i < @ge; $i++) {
	    push(@{$gei{$ge[$i]}}, $i);
	}
	if ($gei{"V2"} && $ref[$gei{"V2"}[0]] eq $mref) {
	    $vi = $gei{"V2"}[0];
	    $v = $vdj[$vi];
	}
	if ($gei{"V0"} && $ref[$gei{"V0"}[0]] eq $mref) {
	    $vi = $gei{"V0"}[0];
	    $v = $vdj[$vi];
	}
	if ($gei{"D0"} && $ref[$gei{"D0"}[0]] eq $mref) {
	    $di = $gei{"D0"}[0];
	    $d = $vdj[$di];
	}
	if ($gei{"J0"} && $ref[$gei{"J0"}[0]] eq $mref && $ji == -1) {
	    $ji = $gei{"J0"}[0];
	    $j = $vdj[$ji];
	}
    }

    # re-orient the alignments
    if ($ori eq "-") {
	@galn = reverse(@galn);
	@{$pgaln} = @galn;
    }

    # calculate alignment length
    my %oqsoqe;
    if ($reg > 0) {
	if ($ori eq "+") {
	    $oqsoqe{$galn[$vi]->{aln}[0]->{oqs}} = $galn[$ji]->{aln}[0]->{oqe};
	} else {
	    $oqsoqe{$galn[$ji]->{aln}[0]->{oqs}} = $galn[$vi]->{aln}[0]->{oqe};
	}
	for (my $i = 0; $i < @galn; $i++) {
	    next if $i == $vi || $i == $ji;
	    $oqsoqe{$galn[$i]->{aln}[0]->{oqs}} = $galn[$i]->{aln}[0]->{oqe};
	}
    } else {
	for (my $i = 0; $i < @galn; $i++) {
	    $oqsoqe{$galn[$i]->{aln}[0]->{oqs}} = $galn[$i]->{aln}[0]->{oqe};
	}
    }
    my @oqs = sort { $a <=> $b } keys %oqsoqe;

    my $ual = ($oqs[0] - 1) . "," . ($ql - $oqsoqe{$oqs[-1]});
    
    my $al = 0;
    my $oqs = $oqs[0];
    my $oqe = $oqsoqe{$oqs};
    for (my $i = 1; $i < @oqs; $i++) {
	if ($oqs[$i] > $oqe) {
	    $al += ($oqe - $oqs + 1);
	    $oqs = $oqs[$i];
	    $oqe = $oqsoqe{$oqs};
	} else {
	    $oqe = $oqsoqe{$oqs[$i]};
	}
    }
    $al += ($oqe - $oqs + 1);

    return ($reg, $mref, $ori, "$v:$d:$j", "$vi,$di,$ji", $ual, $al);
}


sub OutputGroupAln {
    my @galn = @_;

    my @out;
    foreach my $a (@galn) {
	push(@out, $a->out());
    }
    my $out = join " ", @out;

    return $out;
}


1;
