#!/usr/bin/perl -w
# usage : CombinePEVDJDelta.pl r1.vdjdelta r2.vdjdelta

use strict;


#################### combine paired-end results ####################

open IN1, "<$ARGV[0]" || die "open $ARGV[0]: $!\n";
open IN2, "<$ARGV[1]" || die "open $ARGV[1]: $!\n";

while (<IN1>) {
    my $a1 = $_; chomp $a1;
    my @a1 = split "\t", $a1;
    my $a2 = <IN2>; chomp $a2;
    my @a2 = split "\t", $a2;

    # exit if the query IDs are different
    die "error: inconsistent paired reads!\n" if $a1[0] ne $a2[0];

    # if one of the paired reads is not aligned
    my @a = ();
    if ($a1[2] eq "---" || $a2[2] eq "---") {
	if ($a2[2] eq "---") {
            @a = ($a1[0], "1", @a1[1..$#a1], "2", @a2[1..$#a2]);
	} else {
            @a = ($a2[0], "2", @a2[1..$#a2], "1", @a1[1..$#a1]);
	}
        print join("\t", @a) . "\n";
	next;
    }
    
    # combine VDJ annotation based on the better alignment
    my ($v1, $d1, $j1) = split ":", $a1[5];
    my ($v2, $d2, $j2) = split ":", $a2[5]; 
    my @aln1 = split " ", $a1[6];
    my @aln2 = split " ", $a2[6];
    my ($v1i, $d1i, $j1i) = split ",", $a1[7];
    my ($v2i, $d2i, $j2i) = split ",", $a2[7];

    # if alignment of read 1 is better
    if ($a1[8] >= $a2[8]) {
	
	# if strictly regular
	if ($a1[2] == 2) {
	    
	    # if orientation of read 2 is consistent
	    if ($a1[4] ne $a2[4]) {
		my ($cv, $caln) = ConsensusV($v1, $v2, $aln1[$v1i]);
		my $cd = ConsensusD($d1, $d2);
		my $cj = ConsensusJ($j1, $j2);
		
		# if the Vs, Ds, or Js are inconsistent, lower the regularity
		if (!$cv || !$cd || !$cj) {
		    $a1[2] = 1;
		    
		# if the Vs and Js are consistent, deambiguity if possible
		} elsif ($cv ne $v1) {
		    $a1[5] = "$cv:$cd:$cj";
		    $aln1[$v1i] = $caln;
		    $a1[6] = join " ", @aln1;
		}
	    }
	    
	# if non-regular
	} elsif ($a1[2] == 0) {

	    # set possible VJ recombination
	    my $cd = ConsensusD($d1, $d2);
	    if ($v1i == -1 && $v2i != -1 && $j1i != -1 && $j2i == -1) {
		$a1[5] = "$v2:$cd:$j1";
	    } elsif ($v1i != -1 && $v2i == -1 && $j1i == -1 && $j2i != -1) {
		$a1[5] = "$v1:$cd:$j2";
	    }
	}

    # if alignment of read 2 is better
    } else {
	
	# if strictly regular
	if ($a2[2] == 2) {
	    
	    # if orientation of read 2 is consistent
	    if ($a2[4] ne $a1[4]) {
		my ($cv, $caln) = ConsensusV($v2, $v1, $aln2[$v2i]);
		my $cd = ConsensusD($d2, $d1);
		my $cj = ConsensusJ($j2, $j1);
		
		# if the Vs or Js are inconsistent, lower the regularity
		if (!$cv || !$cd || !$cj) {
		    $a2[2] = 1;
		    
		# if the Vs and Js are consistent, deambiguity if possible
		} elsif ($cv ne $v2) {
		    $a2[5] = "$cv:$cd:$cj";
		    $aln2[$v2i] = $caln;
		    $a2[6] = join " ", @aln2;
		}
	    }
	    
	# if non-regular
	} elsif ($a2[2] == 0) {

	    # set possible VJ recombination
	    my $cd = ConsensusD($d1, $d2);
	    if ($v2i == -1 && $v1i != -1 && $j2i != -1 && $j1i == -1) {
		$a2[5] = "$v1:$cd:$j2";
	    } elsif ($v2i != -1 && $v1i == -1 && $j2i == -1 && $j1i != -1) {
		$a2[5] = "$v2:$cd:$j1";
	    }
	}
    }

    # resolve read1,2 ambiguity 
    if ($a1[6] =~ /\|/ && $a2[6] !~ /\|/) {
        my %g2c;
        foreach my $ta (@aln2) {
            if ($ta =~ /([VDJC][\d-]+)/) {
                $g2c{$1}++;
            }
        }
        if (%g2c) {
            foreach my $ta (@aln1) {
                next if $ta !~ /\|/;
                my @ta = split /\|/, $ta;
                my @g;
                foreach my $t (@ta) {
                    if ($t =~ /([VDJC][\d-]+)/) {
                        push(@g, $1);
                    }
                }
                if (@g) {
                    foreach my $e (@g) { $g2c{$e} = 0 if !$g2c{$e}; }
                    @g = sort { $g2c{$b} <=> $g2c{$a} } @g;
                    if ($g2c{$g[0]} > $g2c{$g[1]}) {
                        @ta = grep(/$g[0]/, @ta);
                        $ta = $ta[0];
                    } elsif ($ta =~ /C/ && $g2c{"J2"}) {
                        @ta = grep(/C2/, @ta);
                        $ta = $ta[0];
                    }
                }
            }
            $a1[6] = join(" ", @aln1);
        }
    } elsif ($a1[6] !~ /\|/ && $a2[6] =~ /\|/) {
        my %g1c;
        foreach my $ta (@aln1) {
            if ($ta =~ /([VDJC][\d-]+)/) {
                $g1c{$1}++;
            }
        }
        if (%g1c) {
            foreach my $ta (@aln2) {
                next if $ta !~ /\|/;
                my @ta = split /\|/, $ta;
                my @g;
                foreach my $t (@ta) {
                    if ($t =~ /([VDJC][\d-]+)/) {
                        push(@g, $1);
                    }
                }
                if (@g) {
                    foreach my $e (@g) { $g1c{$e} = 0 if !$g1c{$e}; }
                    @g = sort { $g1c{$b} <=> $g1c{$a} } @g;
                    if ($g1c{$g[0]} > $g1c{$g[1]}) {
                        @ta = grep(/$g[0]/, @ta);
                        $ta = $ta[0];
                    } elsif ($ta =~ /C/ && $g1c{"J2"}) {
                        @ta = grep(/C2/, @ta);
                        $ta = $ta[0];
                    }
                }
            }
            $a2[6] = join(" ", @aln2);
        }
    }  
    
    # get combined results
    # if alignment of read 1 is better
    if ($a1[8] >= $a2[8]) {
        @a = ($a1[0], "1", @a1[1..$#a1], "2", @a2[1..$#a2]);
    } else {
        @a = ($a2[0], "2", @a2[1..$#a2], "1", @a1[1..$#a1]);
    }
        
    # output
    print join("\t", @a) . "\n";
}



############################################################


sub ConsensusV {
    my ($v1, $v2, $aln1) = @_;

    if ($v2 eq "---") {
	return ($v1, $aln1);

    } else {
	my @v1 = split /\|/, $v1;
	my @aln1 = split /\|/, $aln1;
	my @v2 = split /\|/, $v2;
	my %v2q;
	foreach my $e (@v2) { $v2q{$e} = 1; }
	my @cv;
	my @caln;
	for (my $i = 0; $i < @v1; $i++) {
	    if ($v2q{$v1[$i]}) {
		push(@cv, $v1[$i]);
		push(@caln, $aln1[$i]);
	    }
	}
	my $cv = join("|", @cv);
	my $caln = join("|", @caln);
	return ($cv, $caln);
    }
}


sub ConsensusD {
    my ($d1, $d2) = @_;

    if ($d1 eq "---" && $d2 eq "---") {
	return "---";

    } elsif ($d2 eq "---") {
	return $d1;

    } elsif ($d1 eq "---") {
	return $d2;

    } elsif ($d1 eq $d2) {
	return $d1;

    } else {
	return "";
    }
}


sub ConsensusJ {
    my ($j1, $j2) = @_;

    if ($j2 eq "---") {
	return $j1;

    } elsif ($j1 eq $j2) {
	return $j1;

    } else {
	return "";
    }
}
