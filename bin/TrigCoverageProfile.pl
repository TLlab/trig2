#!/usr/bin/perl -w
# usage : TrigCoverageProfile.pl vdjdelta

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


# load reference length

my $rl = 0;

open IN, "<$dir/gene/$sg.fa" || die "open $dir/gene/$sg.fa: $!\n";

<IN>;
while (<IN>) {
    chomp;
    $rl += length($_);
}
close IN;


# load VDJ exon info
# pa : position annotation

my %pa;

open IN, "<$dir/gene/$sg.vdj" || die "open $dir/gene/$sg.vdj: $!\n";

my $le = 0;
my $lg = "START";

while (<IN>) {
    my @a = split "\t"; chomp $a[-1];
    
    # er : exon range
    # e0 : only one exon
    # e1, e2, ... : more than one exon
    # i1, i2, ... : intron within a V or C gene
    # ns : numbering start
    my @er = split ",", $a[4];
    my $ns = @er > 1 ? 1 : 0;
    
    # set VDJ exon info
    for (my $i = 0; $i < @er; $i++) {
	my $ei = $a[3] eq "+" ? $i + $ns : $#er + $ns - $i;
	my $vdje = "$a[1]\_e$ei";
	my ($s, $e) = $er[$i] =~ /(\d+)\.\.(\d+)/;
	for (my $i = $s; $i <= $e; $i++) {
	    $pa{$i} = $vdje;
	}
    }
    
    # set VDJ intron info
    if (@er > 1) {
	for (my $i = 1; $i < @er; $i++) {
	    my $ii = $a[3] eq "+" ? $i + $ns - 1 : $#er + $ns - $i - 1;
	    my $vdji = "$a[1]\_i$ii";
	    my ($s) = $er[$i-1] =~ /\.\.(\d+)/; $s++;
	    my ($e) = $er[$i] =~ /(\d+)\.\./; $e--;
	    for (my $i = $s; $i <= $e; $i++) {
		$pa{$i} = $vdji;
	    }
	}
    }

    # set intergenic info
    my ($ts) = $er[0] =~ /(\d+)\.\./;
    my $int = "INT_$lg:$a[1]";
    for (my $i = $le + 1; $i < $ts; $i++) {
	$pa{$i} = $int;
    }
    ($le) = $er[-1] =~ /\.\.(\d+)/;
    $lg = $a[1];
}
close IN;

for (my $i = $le + 1; $i <= $rl; $i++) {
    $pa{$i} = "INT_$lg:END";
}


# load TRIg annotation
# pi : position increment

my %pi;

open IN, "<$ARGV[0]" || die "open $ARGV[0]: $!\n";

while (<IN>) {
    my @a = split "\t"; chomp $a[-1];
    
    next if @a == 3 || @a == 7;

    if (@a == 9) {
        SetCoverage($a[6]);
    } elsif (@a == 13) {
        SetCoverage($a[7]);
    } elsif (@a == 19) {
        SetCoverage($a[7]);
        SetCoverage($a[16]);
    }
}
close IN;


# calculate coverage profile

my @pc = map(0, (1..$rl));

my @p = sort { $a <=> $b } keys %pi;

$pc[$p[0]-1] = $pi{$p[0]};
    
for (my $i = 1; $i < @p; $i++) {
    for (my $k = $p[$i-1] + 1; $k < $p[$i]; $k++) {
	$pc[$k-1] = $pc[$p[$i-1]-1];
    }
    $pc[$p[$i]-1] = $pc[$p[$i-1]-1] + $pi{$p[$i]};
}


# output coverage profile

for (my $i = 1; $i <= $rl; $i++) {
    my $pc = sprintf("%.1f", $pc[$i-1]);
    print "$i\t$pc\t$pa{$i}\n";
}



######################################################################


sub SetCoverage {
    my @aln = split " ", $_[0];
    
    foreach my $aln (@aln) {
        my @aa = split /\|/, $aln;
        my $n = 1 / scalar(@aa);

        foreach my $aa (@aa) {
            my @b = split ":", $aa;
            my ($s, $e) = $b[2] =~ /(\d+)-(\d+)/;
            $pi{$s} += $n;
            $pi{$e+1} -= $n;
        }
    }
}
