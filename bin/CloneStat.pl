#!/usr/bin/perl -w
# usage : CloneStat.pl clone.txt [prefix]
# to do : statistics clone (vpc, jpc, vjpc, nlpc, cnpc, aanpc, vjnpc)

use strict;
use File::Basename qw(dirname);
use Cwd qw(abs_path);
use Getopt::Long qw(GetOptions);


##### set parameters

my $s = "hsa";
my $g = "trb";

GetOptions(
    "s=s" => \$s,
    "g=s" => \$g,
    );

my $sg = "$s\_$g";
my $dir = dirname(dirname(abs_path($0)));

# set prefix
my $bs = `basename -s .txt $ARGV[0]`; chomp $bs;
$bs  = $ARGV[1] if defined $ARGV[1];


#####  load all VJ genes

my @v;
my @j;

open IN, "<$dir/gene/$sg.vdj" || die "open $dir/gene/$sg.vdj: $!\n";

while (<IN>) {
    my @F = split /\s+/, $_;
    my $vdjc = substr($F[1], 3, 1);
    if ($vdjc eq "V") {
	push(@v, $F[1]);
    } elsif ($vdjc eq "J" && $F[1] !~ /P/) {
	push(@j, $F[1]);
    }
}


##### load clone data

my ($count, %vn, %jn, %vjn, %aan, %nln);

open IN, "<$ARGV[0]" || die "$!\n";
<IN>;
while (<IN>) {
    chomp;
    my @F = split "\t", $_;
    next if $F[3] =~ /[\*|\_]/;  # stop condon, not triple
    $count += $F[0];
    $vn{$F[4]} += $F[0];
    $jn{$F[6]} += $F[0];
    $vjn{"$F[4]:$F[6]"} += $F[0];
    $aan{$F[3]} += $F[0];
    $nln{length($F[2])} += $F[0];
}


##### output clone stat

open OVNPC, ">$bs.vnpc" || die "$!\n";
open OJNPC, ">$bs.jnpc" || die "$!\n"; 
open OVJPC, ">$bs.vjpc" || die "$!\n"; 
open ONLPC, ">$bs.nlpc" || die "$!\n";
open OCNPC, ">$bs.cnpc" || die "$!\n";
open OVJNPC, ">$bs.vjnpc" || die "$!\n";
open OAANPC, ">$bs.aanpc" || die "$!\n";

if (!$count) {
    close OVNPC;
    close OJNPC;
    close OVJPC;
    close ONLPC;
    close OCNPC;
    close OVJNPC;
    close OAANPC;
    exit 0;
}

# vnpc
for (@v) {
    $vn{$_} = 0 unless $vn{$_};
    my $vpc = sprintf("%.4f", $vn{$_}/$count*100);
    print OVNPC "$_\t$vn{$_}\t$vpc\n";
}

# jnpc
for (@j) {
    $jn{$_} = 0 unless $jn{$_};
    my $jpc = sprintf("%.4f", $jn{$_}/$count*100);
    print OJNPC "$_\t$jn{$_}\t$jpc\n";
}

# nlpc
for (sort {$a <=> $b} keys %nln) {
    my $nlpc = sprintf( "%.4f", $nln{$_}/$count*100 );
    print ONLPC "$_\t$nln{$_}\t$nlpc\n";
}

# cnpc
seek(IN, 0, 0);
<IN>;    # skip title
while (<IN>){
    chomp;
    my @F = split "\t", $_;
    next if $F[3] =~ /[\*|\_]/;  # stop condon, not triple
    print OCNPC "$F[4]:$F[2]:$F[6]\t$F[0]\t" . sprintf("%.4f", $F[0]/$count*100) . "\n";
}
close IN;

# vjpc/vjnpc
print OVJPC "V\\J\t" . join("\t", @j) . "\n";
for my $v (@v){
    print OVJPC "$v\t";
    for my $j (@j){
	$vjn{"$v:$j"} = 0 unless $vjn{"$v:$j"};
	my $vjpc = sprintf("%.4f", $vjn{"$v:$j"}/$count*100);
	print OVJPC "$vjpc\t";
	print OVJNPC "$v:$j\t" . $vjn{"$v:$j"} . "\t$vjpc\n";
    }
    seek(OVJPC, -1, 1); print OVJPC "\n"; # replace \t -> \n
}

# aanpc
for (sort { $aan{$b} <=> $aan{$a} } keys %aan){
    my $aapc = sprintf("%.4f", $aan{$_}/$count*100);
    print OAANPC "$_\t$aan{$_}\t$aapc\n";
}

close OVNPC;
close OJNPC;
close OVJPC;
close ONLPC;
close OCNPC;
close OVJNPC;
close OAANPC;
