#!/usr/bin/perl -w
# usage : LabelRecombination.pl -s hsa -g trb [vdjdelta]
# update: *using for new delta format in Trig2
#       : *add CorrectAmbiguousRegion subroutine
#       : *correct delta information and defined category order
#       : *corrcet C1 -> C2 if align to J2
#       : *unified field

use strict;
use File::Basename qw(dirname);
use Cwd qw(abs_path);
use Getopt::Long qw(GetOptions);

# set parameters
my $s = "hsa";
my $g = "trb";

GetOptions(
	"s=s" => \$s,
	"g=s" => \$g,
);

my $sg = "$s\_$g";
my $dir = dirname(dirname(abs_path($0)));


# defined category order
my %Cat_order = ("CH" => 1, "AS" => 2, "AR" => 3, "CR" => 4, "PCR" => 5, 
				"PR" => 6, "PPR" => 7, "NR" => 8, "NS" => 9, "UC" => 10);

# load VDJ info

my %gs;
open IN, "<$dir/gene/$sg.vdj" || die "open $dir/gene/$sg.vdj: $!\n";

while (<IN>) {
    my @a = split "\t"; chomp $a[-1];

    ($gs{$a[1]}) = $a[4] =~ /^(\d+)/;
}
close IN;


# load RSS info

my %lrss;
open IN, "<$dir/gene/$sg\_rss.txt" || die "open $dir/gene/$sg\_rss.txt: $!\n";

while (<IN>) {
    my @a = split "\t"; chomp $a[-1];

    next if $a[0] eq "RSS12";
    
    my $l = $a[1];
    for (my $i = $l - 30; $i <= $l + 3; $i++) {
        $lrss{$i} = $a[0];
    }
}

close IN;


# load vdj delta

open IN, "<$ARGV[0]" || die "open $ARGV[0]: $!\n";

while (<IN>) {
	my @a = split "\t"; chomp $a[-1];
	next if @a == 4 || @a == 7;

	if (@a == 10) {
		my ($cat, $g, $j, $da) = Category($a[3], $a[7]);
		$a[7] = $da;
		my $record = join "\t", @a;
		print "$cat\t$g\t$j\t$record\n";

	} elsif (@a == 13) {
		my ($cat, $g, $j, $da) = Category($a[3], $a[7]);
		$a[7] = $da;
		my $record = join "\t", @a;
		print "$cat\t$g\t$j\t$record\n";

	} elsif (@a == 19) {
		my ($_cat1, $_g1, $_j1, $_da1) = Category($a[3], $a[7]);
		my ($_cat2, $_g2, $_j2, $_da2) = Category($a[12], $a[16]);

		# sort category
		my ($cat) = sort{ $Cat_order{$a} <=> $Cat_order{$b} } ($_cat1, $_cat2);

		my $g = $cat eq $_cat1 ? $_g1 : $_g2;
		my $j = $cat eq $_cat1 ? $_j1 : $_j2;
		$a[7]  = $_da1;
		$a[16] = $_da2;
		my $record = join "\t", @a;
		print "$cat\t$g\t$j\t$record\n";

	} else {
		print STDERR "Warning: undifined filed!\n";
	}
}


#####

sub Category {
	my $fg = shift;
	my $da = shift;

	#### new version delta : hsa_trb:TRBI_0:548436-548538:59-161:2:0 hsa_trb:TRBC1_1:554193-554215:162-184:0:0
	my @vdje = split " ", $da;
	my @vdjeqr = map { (split ":")[3] } @vdje;
	my @o = map { &ori($_) } @vdjeqr;
	my @uo = Unique(@o);

	# remove "hsa_trb:" to get the same array index to conform old vdjdelta format
	my @delta = split " ", $da =~ s/hsa_trb://gr;

	# corrected VDJ order
	if (@delta > 2 and @uo == 1) {
		my $s1 = (split ":", $delta[0])[1];
		my $s2 = (split ":", $delta[1])[1];
		@delta = reverse @delta if (split "-", $s1)[0] > (split "-", $s2)[0];
	}

	# process delta alignments
	
	# J1 C1|C2 -> J1 C1
	# J2 C1|C2 -> J2 C2
	@delta = CorrectAmbiguousRegion(\@delta, $uo[0]) if @uo == 1;

	# C1 J2 -> J2 C2 (also change reference position)
	@delta = CorrectC1C2Ambiguous(\@delta, $uo[0]) if @uo == 1;
	@delta = MergeBrokenAlignment(\@delta, $uo[0]) if @uo == 1;

	foreach my $d (@delta) {
		if ($d =~ /V/ && $d =~ /_1/ && $d !~ /~/) {
			$d =~ s/V/v/;
		}
	}

	# check if V is good (one V gene without broken pieces)
	my %vc;
	my %vec;
	foreach my $d (@delta) {
		my @b = split ":", $d;

		if ($b[0] =~ /V/) {
			my ($vdj, $e) = $b[0] =~ /TRB(.+)\_([12])/ ? ($1, $2) : ("V8-1", "1");
			$vc{$vdj}++;
			$vec{"$vdj$e"}++;
		}

	}
	my $vn = scalar(keys %vc);
	my @tmp = sort { $b <=> $a } (values %vec);
	my $vemaxc = $tmp[0];
	my $vgood = ($vn == 1 && $vemaxc <= 1) ? 1 : 0;

	# get gene count
	my @g = map(substr($_, 3, 1), @delta);
	my $g = join("", @g);
	my %gc;
	foreach my $g (@g) {
		$gc{$g}++;
	}
	foreach my $g ("v", "V", "D", "J", "C","I") {
		$gc{$g} = 0 if !$gc{$g};
	}

	my $cat;

	if ($fg == 2) {
		$cat = "CR";

	} elsif (@uo > 1 || $gc{"C"} > 1 || IncreasingQ(@delta) == 0) {
		$cat = "CH";

	} else {
		if ($g =~ /^v?VD?JC?$/) {
			$cat = "CR";

		} elsif ($g =~ /^DJC?$/ || $g eq "D") {
			if ($g ne "D") {
				$cat = "PR";
			} else {
				$cat = "PPR";
			}

		} elsif ($g =~ /^JC?$/) {
			$cat = ProcessJC(\@delta);

		} elsif ($g =~ /^v?VD?$/ && $vgood == 1) {
			$cat = "PCR";

		} elsif ($g eq "C") {
			$cat = ProcessC(\@delta);

		} elsif ($g eq "I" || $g eq "v") {
			$cat = "UC";

		} else {
			if ($g =~ /^IJC?$/) {
				$cat = ProcessIJC(\@delta, $uo[0]);

			} elsif ($g =~ /^JJC?$/) {
				#$cat = ProcessIJC(\@delta, $uo[0]);
				$cat = "AS";

			} else {
				$cat = "AS";
			}
		}
	}

	my $j = ($cat eq "UC" || $cat eq "CH") ? 0 : 1;

	$da = join(" ", @delta);
	
	return ($cat, $g, $j, $da);
}

sub ori {
	my ($s, $e) = shift =~ /(\d+)-(\d+)/;
	return $s < $e ? "+" : "-";
}


sub Unique {
    my %seen;
    return grep(!$seen{$_}++, @_);
}


sub ModifyC {
    my @d = @_;

    my @s;
    foreach my $d (@d) {
	my @b = split ":", $d;
	my ($s) = $b[1] =~ /^(\d+)/;
	push(@s, $s);
    }

    for (my $i = 1; $i < @d; $i++) {
	if (($d[$i-1] =~ /[DJ]2/ && $d[$i] =~ /C1/) ||
	    ($d[$i-1] =~ /I/ && $d[$i] =~ /C1/ && $s[$i-1] > $s[$i])) {
	    my @a = split ":", $d[$i];
	    $a[0] =~ s/C1/C2/;
	    my ($s, $e) = $a[1] =~ /(\d+)-(\d+)/;
	    $s += 9346;
	    $e += 9346;
	    $a[1] = "$s-$e";
	    $d[$i] = join(":", @a);
	}
    }

    return @d;
}

sub CorrectAmbiguousRegion{
	my @delta = @{$_[0]};
	my $o = $_[1];
	for (my $i = 1; $i < @delta; $i++) {
		my @ambiguous = split '\|', $delta[$i];
		if (@ambiguous > 1) {
			my @ld = split ":", $delta[$i-1];
			if ($ld[0] =~ /TRBJ1/) {
				$delta[$i] = $ambiguous[0] if $ambiguous[0] =~ /TRBC1/;
			} elsif ($ld[0] =~ /TRBJ2/) {
				$delta[$i] = $ambiguous[1] if $ambiguous[1] =~ /TRBC2/;
			}
		}
	}
	return @delta;
}

sub CorrectC1C2Ambiguous{
	my @delta = @{$_[0]};
	my $o = $_[1];
	for (my $i = 1; $i < @delta; $i++) {
		my @ld = split ":", $delta[$i-1];
		my @rd = split ":", $delta[$i];
		if ($ld[0] =~ /TRBJ2/ && $rd[0] eq "TRBC1_1") {
			$rd[0] = "TRBC2_1";
			$rd[1] = join "-", map{ $_ + 9346 } split "-", $rd[1];
			$delta[$i] = join ":", @rd;
			last;
		}
	}
	return @delta;
}

sub MergeBrokenAlignment {
    my @delta = @{$_[0]};
    my $o = $_[1];
    
    for (my $i = 1; $i < @delta; $i++) {
	my @ld = split ":", $delta[$i-1];
	my ($lts, $lte) = $ld[1] =~ /(\d+)-(\d+)/;
	my ($lqs, $lqe) = $ld[2] =~ /(\d+)-(\d+)/;
	my @rd = split ":", $delta[$i];
	my ($rts, $rte) = $rd[1] =~ /(\d+)-(\d+)/;
	my ($rqs, $rqe) = $rd[2] =~ /(\d+)-(\d+)/;

	my $td = $rts - $lte - 1;
	next if $td < 0;
	my $qd = $o eq "+" ? $rqs - $lqe - 1 : $lqe - $rqs - 1; 
	my $bd = $td > $qd ? $td : $qd;
	next if $bd == 0;
	my $dd = abs($td - $qd);
	my $ddp = $dd / $bd;

	if ($ddp <= 0.5 || $bd <= 5) {
	    my $lg = $ld[0];
	    my $rg = $rd[0];

	    my $ng;
	    if ($lg eq $rg) {
		$ng = $lg;
	    } elsif ($lg =~ /I/ && $rg !~ /I/) {
		$ng = $rg;
	    } elsif ($lg !~ /I/ && $rg =~ /I/) {
		$ng = $lg;
	    } else {
		$ng = "$lg~$rg";
	    }

	    my $mm = $ld[3] + $rd[3];
	    my $gap = "$ld[4],$rd[4]";
	    my $nd = "$ng:$lts-$rte:$lqs-$rqe:$mm:$gap";
	    $delta[$i-1] = $nd;
	    splice(@delta, $i, 1);
	    $i--;
	}
    }

    return @delta;
}


sub IncreasingQ {
    my @delta = @_;

    my @start;
    foreach my $d (@delta) {
	my @b = split ":", $d;
	my ($s, $e) = $b[1] =~ /^(\d+)-(\d+)/;
	push(@start, $s);
    }

    for (my $i = 0; $i < @delta; $i++) {
	if ($delta[$i] =~ /V30/) {
	    splice(@delta, $i, 1);
	    splice(@start, $i, 1);
	    $i++;
	}
    }
    
    my $q = 1;
    for (my $i = 1; $i < @delta; $i++) {
	if ($start[$i] < $start[$i-1]) {
	    $q = 0;
	    last;
	}
    }
    
    return $q;
}    


sub ProcessJC {
    my @delta = @{$_[0]};

    my @d = split ":", $delta[0];
	my ($ts, $te) = $d[1] =~ /(\d+)-(\d+)/;

    my $g = $d[0]; $g =~ s/_\d//;
    if ($g =~ /~/) {
	my @tmp = split "~", $g;
	$g = $tmp[-1]; $g =~ s/_\d//;
    }
    my $gul = $gs{$g} - $ts;

	if ($gul >= 15 || $d[0] =~ /J2-2P/) {
	return "NR";
    } else {
	return "UC";
    }
}


sub ProcessC {
    my @delta = @{$_[0]};

    my @d = split ":", $delta[0];
	#my ($ts) = $d[1] =~ /(\d+)-\d+/;
	my ($ts, $te) = $d[1] =~ /(\d+)-(\d+)/;

	my $g = $d[0];
	#$g =~ s/_e\d//;
	$g =~ s/_\d//;
    my $gul = $gs{$g} - $ts;

    if ($gul >= 15) {
	return "NS";
    } else {
	return "UC";
    }
}


sub ProcessIJC {
    my @delta = @{$_[0]};
    my $o = $_[1];

    my @j = split ":", $delta[1];
	my ($qjs, $qje) = $j[2] =~ /(\d+)-(\d+)/;
	#my ($tjs) = $j[1] =~ /(\d+)-\d+/;
	my ($tjs, $tje) = $j[1] =~ /(\d+)-(\d+)/;

    my @i = split ":", $delta[0];
	my ($qis, $qie) = $i[2] =~ /(\d+)-(\d+)/;
	#my ($tie) = $i[1] =~ /\d+-(\d+)/;
	my ($tis, $tie) = $i[1] =~ /(\d+)-(\d+)/;

    if ((548682<=$tie && $tie<=548711) || (558174<=$tie && $tie<=558203)) {
	return "PPR";

    } elsif($lrss{$tie}) {
        return "AR";
        
    } else {
	return "AS";
    }
}
