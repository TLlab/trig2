#!/usr/bin/perl -w
# usage : PatchAlignment.pl vdjdelta read.fa

use strict;
use File::Basename qw(dirname);
use Cwd qw(abs_path);
use Getopt::Long qw(GetOptions);

use lib dirname(abs_path($0));
use VDJDelta qw(SetVDJInfo);
use GroupAln qw(SetAnnotation OutputGroupAln);
use Fastx qw(LoadFasta LoadThisFasta);


##### set parameters

my $species  = "hsa";
my $gene     = "trb";
my $minmatch = 15;
my $minid    = 0.8;
my $maxe     = 0.01;

GetOptions(
    "species=s"  => \$species,
    "gene=s"     => \$gene,
    "minmatch=i" => \$minmatch,
    );

my $sg = "$species\_$gene";
my @gene = split "_", $gene;

SetVDJInfo("$sg.vdj");


##### load VDJ delta

my @bmq;
my %qkv;

open IN, "<$ARGV[0]" || die "open $ARGV[0]: $!\n";

while (<IN>) {
    my @a = split "\t"; chomp $a[-1];

    # skip unaligned read
    next if $a[2] eq "---";

    my ($q, $ql, $reg, $sg, $ori, $vdj, $vdjdelta) = @a;
    
    # load alignment info
    my @galn = split " ", $vdjdelta;
    foreach my $ga (@galn) {
	my @aln = split /\|/, $ga;
	foreach my $a (@aln) {
	    $a = VDJDelta->new($a);
	}
	$ga = GroupAln->new(@aln);
    }

    # find potential broken alignments
    my $bflag = 0;
    my $bref = "";
    my ($bsi, $bei);
    my ($brs, $bre);
    my ($bqs, $bqe);
    
    # for each pair of neighboring alignments
    for (my $i = 1; $i < @galn; $i++) {
	
	# check if the grouped alignments are consistent
	my $cbq = $galn[$i-1]->{cbq} * $galn[$i-1]->{cbq};
	$cbq *= ($galn[$i-1]->{ref} eq $galn[$i]->{ref} ? 1 : 0);
	$cbq *= ($galn[$i-1]->{n} == $galn[$i]->{n} ? 1 : 0);
	
	# calculate alignment distances (based on the first of a grouped alignment)
	my $ro = $galn[$i]->{ro};
	my $rad = $galn[$i]->{aln}[0]->{rs} - $galn[$i-1]->{aln}[0]->{re} - 1;
	my $qad = $ro eq "+" ? $galn[$i]->{aln}[0]->{qs} - $galn[$i-1]->{aln}[0]->{qe} - 1 : $galn[$i-1]->{aln}[0]->{qe} - $galn[$i]->{aln}[0]->{qs} - 1;
	my $add = abs($rad - $qad);
        my $adm = $rad > $qad ? $rad : $qad;

	# if broken alignments
	if ($add <= $adm*0.5 && $galn[$i-1]->{ro} eq $ro && $cbq == 1) {
	    
            # if the first broken alignment
            if ($bflag == 0) {
		$bref = $galn[$i-1]->{ref};
                ($bsi, $bei) = ($i - 1, $i);
                ($brs, $bre) = ($galn[$i-1]->{rs}, $galn[$i]->{re});
                ($bqs, $bqe) = ($galn[$i-1]->{qs}, $galn[$i]->{qe});
                $bflag = 1;

            # if the following broken alignment
            } else {
                $bei = $i;
                $bre = $galn[$i]->{re};
                $bqe = $galn[$i]->{qe};
            }

	# if not broken alignments
	} else {
	    
            # record the previous broken alignments if exist
            if ($bflag == 1) {
		push(@{$qkv{$q}{bref}}, $bref);
                push(@{$qkv{$q}{bsei}}, "$bsi,$bei");
                push(@{$qkv{$q}{brse}}, "$brs,$bre");
                push(@{$qkv{$q}{bqse}}, $ro eq "+" ? "$bqs,$bqe" : "$bqe,$bqs");
		$bflag = 0;
            }
	    
	    # check potential missing J alignment
	    my $cmq = $galn[$i-1]->{cmq} * $galn[$i-1]->{cmq};
            if ($qad > $minmatch && $cmq == 1) {

		# skip if not possible J
		next if $galn[$i-1]->{ge} !~ /[VD][02]/ || $galn[$i]->{ge} ne "C1";

		# skip if orentations of the VD and C are not consistent
		next if $galn[$i-1]->{go} ne $galn[$i]->{go};

		my $qs = $ro eq "+" ? ( ($galn[$i-1]->{ro} eq "+") ? $galn[$i-1]->{qe} + 1 : $galn[$i-1]->{qs} + 1 ) : $galn[$i]->{qs} + 1;
		my $qe = $ro eq "+" ? $galn[$i]->{qs} - 1 : ( ($galn[$i-1]->{ro} eq "-") ? $galn[$i-1]->{qe} - 1 : $galn[$i-1]->{qs} - 1 );
		push(@{$qkv{$q}{mi}}, $i);
		push(@{$qkv{$q}{mref}}, $galn[$i-1]->{ref});
		push(@{$qkv{$q}{mqse}}, "$qs,$qe");
		push(@{$qkv{$q}{mge}}, "J0");
	    }
	}
    }

    # record the last broken alignments if exist
    if ($bflag == 1) {
	push(@{$qkv{$q}{bref}}, $bref);
        push(@{$qkv{$q}{bsei}}, "$bsi,$bei");
        push(@{$qkv{$q}{brse}}, "$brs,$bre");
        push(@{$qkv{$q}{bqse}}, $galn[-1]->{ro} eq "+" ? "$bqs,$bqe" : "$bqe,$bqs");
    }
    
    # check potential missing C1 alignment at the right
    if (($galn[-1]->{ge} eq "J0" || $galn[-1]->{ge} eq "I0") && $galn[-1]->{cmq} == 1) {
        my $qs = $galn[-1]->{ro} eq "+" ? $galn[-1]->{qe} + 1 : 1;
        my $qe = $galn[-1]->{ro} eq "+" ? $ql : $galn[-1]->{qe} - 1;
        my $qad = $qe - $qs + 1;
        
        if ($qad >= $minmatch) {
            push(@{$qkv{$q}{mi}}, $#galn + 1);
	    push(@{$qkv{$q}{mref}}, $galn[-1]->{ref});
            push(@{$qkv{$q}{mqse}}, "$qs,$qe");
            push(@{$qkv{$q}{mge}}, "C1");
        }
    }

    if ($qkv{$q}) {
	push(@bmq, $q);
	$qkv{$q}{ql} = $ql;
	$qkv{$q}{record} = $_;
	$qkv{$q}{ori} = $ori;
	@{$qkv{$q}{galn}} = @galn;
    }
}
close IN;


##### get sequence segments and do alignment

# set file handles of reference and query

my %rseq;
LoadFasta("$sg.fa", \%rseq);

open my $qfh, "<$ARGV[1]" || die "open $ARGV[1]: $!\n";

my %rgmsegfh;
foreach my $g (@gene) {
    my $sg = "$species\_$g";
    open $rgmsegfh{$sg}{j}, ">$ARGV[0].$sg.jqseg.fa" || die "open $ARGV[0].$sg.jqseg.fa: $!\n";
    open $rgmsegfh{$sg}{c}, ">$ARGV[0].$sg.cqseg.fa" || die "open $ARGV[0].$sg.cqseg.fa: $!\n";
}

foreach my $q (@bmq) {
    my $qseq = LoadThisFasta($qfh, $q, 0);

    # if there is a broken alignment
    if ($qkv{$q}{bref}) {
	my @bref = @{$qkv{$q}{bref}};
        my @bsei = @{$qkv{$q}{bsei}};
        my @brse = @{$qkv{$q}{brse}};
        my @bqse = @{$qkv{$q}{bqse}};

        # for each broken alignment
        for (my $i = 0; $i < @bref; $i++) {
	    my $bref = $bref[$i];

            # prepare reference and query segment
            open OUTBR, ">$ARGV[0].$bref.brseg.fa" || die "open $ARGV[0].$bref.brseg.fa: $!\n";
            open OUTBQ, ">$ARGV[0].$bref.bqseg.fa" || die "open $ARGV[0].$bref.bqseg.fa: $!\n";
	    
	    my ($rs, $re) = split ",", $brse[$i];
            if ($rs !~ /\|/) {
                my $rseg = substr($rseq{$bref}, $rs - 1, $re - $rs + 1);
                print OUTBR ">$bref[$i]:$rs-$re\n$rseg\n";
            
            } else {
                my @rs = split /\|/, $rs;
                my @re = split /\|/, $re;
                for (my $i = 0; $i < @rs; $i++) {
		    my $rseg = substr($rseq{$bref}, $rs[$i] - 1, $re[$i] - $rs[$i] + 1);
                    print OUTBR ">$bref:$rs[$i]-$re[$i]\n$rseg\n";
                }
	    }
            my ($qs, $qe) = split ",", $bqse[$i];
	    my $qseg = substr($qseq, $qs - 1, $qe - $qs + 1);
	    print OUTBQ ">$q:$qs-$qe\n$qseg\n";

            # do global alignment using usearch
            my $cmd = "usearch -usearch_global $ARGV[0].$bref.bqseg.fa -db $ARGV[0].$bref.brseg.fa"; 
            $cmd .= " -strand both -id $minid -evalue $maxe -maxaccepts 4 -userout $ARGV[0].$bref.bqseg.usearch";
            $cmd .= " -userfields query+target+qstrand+qlo+qhi+tlo+thi+mism+gaps+qrowdots+trowdots+alnlen+evalue+bits";
            $cmd .= " 2> /dev/null";
            `$cmd`;
            
            # get filled broken alignment
	    FillBrokenAlignment("$ARGV[0].$bref.bqseg.usearch");
	}
    }

    # if there is a missing alignment
    if ($qkv{$q}{mi}) {
        my @mi = @{$qkv{$q}{mi}};
	my @mref = @{$qkv{$q}{mref}};
        my @mqse = @{$qkv{$q}{mqse}};
        my @mge = @{$qkv{$q}{mge}};

        for (my $i = 0; $i < @mi; $i++) {
	    my ($qs, $qe) = split ",", $mqse[$i];
	    my $qseg = substr($qseq, $qs - 1, $qe - $qs + 1);

	    if ($mge[$i] eq "J0") {
		print { $rgmsegfh{$mref[$i]}{j} } ">$q:$qs-$qe\n$qseg\n";
	    } elsif ($mge[$i] eq "C1") {
		print {$rgmsegfh{$mref[$i]}{c}} ">$q:$qs-$qe\n$qseg\n";
	    }
	}
    }
}
foreach my $g (@gene) {
    my $sg = "$species\_$g";
    close $rgmsegfh{$sg}{j};
    close $rgmsegfh{$sg}{c};
}


# for each gene
foreach my $g (@gene) {
    my $sg = "$species\_$g";
    
    # align potential missing J0 alignment
    if (-s "$ARGV[0].$sg.jqseg.fa") {

	# do local alignment using usearch
	my $cmd = "usearch -ublast $ARGV[0].$sg.jqseg.fa -db $sg\_j.fa"; 
	$cmd .= " -strand both -wordlength 4 -lopen 1.0 -lext 0.5 -evalue $maxe -userout $ARGV[0].$sg.jqseg.ublast";
	$cmd .= " -userfields query+target+qstrand+qlo+qhi+tlo+thi+mism+gaps+qrowdots+trowdots+alnlen+evalue+bits";
	$cmd .= " 2> /dev/null";
	`$cmd`;
	
	GetMissingAlignment("$ARGV[0].$sg.jqseg.ublast");
    }

    # align potential missing C1 alignment
    if (-s "$ARGV[0].$sg.cqseg.fa") {
	
	# do local alignment using usearch
	my $cmd = "usearch -ublast $ARGV[0].$sg.cqseg.fa -db $sg\_c.fa"; 
	$cmd .= " -strand both -wordlength 4 -lopen 1.0 -lext 0.5 -evalue $maxe -userout $ARGV[0].$sg.cqseg.ublast";
	$cmd .= " -userfields query+target+qstrand+qlo+qhi+tlo+thi+mism+gaps+qrowdots+trowdots+alnlen+evalue+bits";
	$cmd .= " 2> /dev/null";
	`$cmd`;
	
	GetMissingAlignment("$ARGV[0].$sg.cqseg.ublast");
    }
}

# remove intermediate files
`rm $ARGV[0]*seg*`;


##### output patched alignments

my %qreg;

open IN, "<$ARGV[0]" || die "open $ARGV[0]: $!\n";
open OUT, ">$ARGV[0].1" || die "open $ARGV[0].1: $!\n";

while( <IN> ) {
    my @a = split "\t"; chomp $a[-1];
    my $q = $a[0];

    # output if not patched
    if (!$qkv{$q}) {
	print OUT $_;
        next;
    }
	
    # if there is a patched broken or missing alignment
    my %sil;
    my %sia;
    
    # set the starting index and length
    if ($qkv{$q}{faln}) {
        my @bsei = @{$qkv{$q}{bsei}};
        my @faln = @{$qkv{$q}{faln}};
        for (my $i = 0; $i < @bsei; $i++) {
            my ($si, $ei) = split ",", $bsei[$i];
            $sil{$si} = $ei - $si + 1 if $faln[$i];
            $sia{$si} = $faln[$i];
        }
    } elsif ($qkv{$q}{maln}) {
        my @mi = @{$qkv{$q}{mi}};
        my @maln = @{$qkv{$q}{maln}};
        for (my $i = 0; $i < @mi; $i++) {
            $sil{$mi[$i]} = 0 if $maln[$i];
            $sia{$mi[$i]} = $maln[$i];
        }
    }
    
    # modify alignment
    if (%sia) {
        my @galn = @{$qkv{$q}{galn}};
        my @si = sort { $b <=> $a } keys %sil;
        foreach my $si (@si) {
            splice(@galn, $si, $sil{$si}, $sia{$si});
        }
        
        # re-annotate modified alignments
        @galn = reverse(@galn) if $qkv{$q}{ori} eq "-";
        my ($reg, $ref, $ori, $vdj, $vdji, $ual, $al) = SetAnnotation(\@galn, $qkv{$q}{ql});
        my @a = split "\t", $qkv{$q}{record};
        $qreg{$q} = $reg if $reg != $a[2];
        $a[2] = $reg;
        $a[3] = $ref;
        $a[5] = $vdj;
        $a[6] = OutputGroupAln(@galn);
        #$a[7] = $vdji;
        #$a[8] = $ual;
        #$a[9] = $al;
        print OUT join("\t", @a) . "\n";
        
    } else {
        print OUT $_;
    }
}
close IN;
close OUT;

# modify CDR3 info

my $fn = $ARGV[0]; $fn =~ s/\.vdjdelta/\.cdr3/;

open IN, "<$fn" || die "open $fn: $!\n";
open OUT, ">$fn.1" || die "open $fn.1: $!\n";

while (<IN>) {
    my @a = split "\t";

    if (!$qreg{$a[0]}) {
        print OUT $_;
        next;
    }

    $a[1] = $qreg{$a[0]};
    $_ = join("\t", @a);
    print OUT $_;
}
close IN;
close OUT;



######################################################################


sub FillBrokenAlignment {
    my $fn = $_[0];

    open IN, "<$fn" || die "open $fn: $!\n";
    while( <IN> ) {
	my @a = split "\t"; chomp $a[-1];
	my ($q) = $a[0] =~ /^(\S+):\d+-\d+$/;
	my $bbits = $a[-1];
	
	# get best filled alignment
	my @faln;
	push(@faln, FormatAlignment(@a));
	while (<IN>) {
	    my @a = split "\t"; chomp $a[-1];
	    last if $a[-1] < $bbits;
	    push(@faln, FormatAlignment(@a));
	}

	my $gfa = @faln ? GroupAln->new(@faln) : "";
	push(@{$qkv{$q}{faln}}, $gfa);
    }
    close IN;
}


sub FormatAlignment {
    my @a = @_;
        
    # obtain alignment info
    my ($qs) = $a[0] =~ /:(\d+)-\d+$/;
    my ($ref, $rs) = $a[1] =~ /^(\S+):(\d+)/;
    my $aqs = $a[2] eq "+" ? $qs + $a[3] - 1 : $qs + $a[4] - 1;
    my $aqe = $a[2] eq "+" ? $qs + $a[4] - 1 : $qs + $a[3] - 1;
    my $ars = $rs + $a[5] - 1;
    my $are = $rs + $a[6] - 1;
    my $mmg = $a[7] + $a[8];
    my @gap = (0);
    if ($a[8] > 0) {
	@gap = GapInfo($a[10], $a[9]);
    }
    my $gap = join ",", @gap;
    my $fa = "$ref:---:$ars-$are:$aqs-$aqe:$mmg:$gap";

    # format the filled alignment
    $fa = VDJDelta->new($fa);
    $fa->annotate();
    
    return $fa;
}


sub GapInfo {
    my ($ta, $qa) = @_;
    my @gap = ();
    
    # obtain gap positions
    my @tgp = GapPosition($ta);
    my @qgp = GapPosition($qa);
    my @gp = sort { $a <=> $b } (@tgp, @qgp);
    my %gptq;
    foreach my $gp (@tgp) {
        $gptq{$gp} = -1;
    }
    foreach my $gp (@qgp) {
        $gptq{$gp} = 1;
    }

    # convert into nucmer delta format
    push(@gap, $gptq{$gp[0]} * $gp[0]);
    for (my $i = 1; $i < @gp; $i++) {
        my $l = $gp[$i] - $gp[$i-1];
        push(@gap, $l * $gptq{$gp[$i]});
    }
    push(@gap, 0);

    return @gap;
}


sub GapPosition {
    my $a = $_[0];
    my @gp = ();

    my $i = index($a, '-');
    while ($i != -1) {
        push(@gp, $i + 1);
        $i = index($a, '-', $i + 1);
    }

    return @gp;
}


sub GetMissingAlignment {
    my $fn = $_[0];

    open IN, "<$fn" || die "open $fn: $!\n";
    while (<IN>) {
	my @a = split "\t"; chomp $a[-1];
	my ($q) = $a[0] =~ /(\S+):/;
	my $bbits = $a[-1];
	
	# get best missing alignment
	my @maln;
	push(@maln, FormatAlignment(@a));
	while (<IN>) {
	    my @a = split "\t"; chomp $a[-1];
	    if ($a[0] !~ /^$q:/) {
		seek(IN, -length($_), 1);
		last;
	    }
	    next if $a[-1] < $bbits;
	    push(@maln, FormatAlignment(@a));
	}

	# discard too short alignments
	@maln = grep($_->{qal} >= 15, @maln);

	my $gfa = @maln ? GroupAln->new(@maln) : "";
	push(@{$qkv{$q}{maln}}, $gfa);
    }
    close IN;
}
