#!/usr/bin/perl -w
# usage : CollectCloneStat.pl folder

use strict;


##### set parameter #####

my $prefix = "all";


##### load trig folders and get sample names #####

my @f;

if (-e $ARGV[0]) {
    open IN, "<$ARGV[0]" || die "open $ARGV[0]: $!\n";
    @f = <IN>;
    chomp @f;
    close IN;
    
} else {
    @f = @ARGV;
}

my @s;

foreach my $f (@f) {
    my $s = $f;
    if ($f =~ /trig_(\S+)/) {
	$s = $1;
    }
    push (@s, $s);
}


##### load statistics #####

# V gene usage

my @v;
my %vsp;

open IN, "<$f[0]/clone.vnpc" || die "open $f[0]/clone.vnpc: $!\n";
while (<IN>) {
    my @a = split "\t"; chomp $a[-1];
    push (@v, $a[0]);
    $vsp{$a[0]}{$s[0]} = $a[2];
}
close IN;

for (my $i = 1; $i < @f; $i++) { 
    open IN, "<$f[$i]/clone.vnpc" || die "open $f[$i]/clone.vnpc: $!\n";
    while (<IN>) {
	my @a = split "\t"; chomp $a[-1];
	$vsp{$a[0]}{$s[$i]} = $a[2];
    }
    close IN;
}

# J gene usage

my @j;
my %jsp;

open IN, "<$f[0]/clone.jnpc" || die "open $f[0]/clone.jnpc: $!\n";
while (<IN>) {
    my @a = split "\t"; chomp $a[-1];
    push (@j, $a[0]);
    $jsp{$a[0]}{$s[0]} = $a[2];
}
close IN;

for (my $i = 1; $i < @f; $i++) { 
    open IN, "<$f[$i]/clone.jnpc" || die "open $f[$i]/clone.jnpc: $!\n";
    while (<IN>) {
	my @a = split "\t"; chomp $a[-1];
	$jsp{$a[0]}{$s[$i]} = $a[2];
    }
    close IN;
}

# VJ frequency

my %vjtn;
my %vjsn;
my %vjsp;

for (my $i = 0; $i < @f; $i++) { 
    open IN, "<$f[$i]/clone.vjnpc" || die "open $f[$i]/clone.vjnpc: $!\n";
    while (<IN>) {
	my @a = split "\t"; chomp $a[-1];
	$vjtn{$a[0]} += $a[1];
	$vjsn{$a[0]}{$s[$i]} = $a[1];
	$vjsp{$a[0]}{$s[$i]} = $a[2];
    }
    close IN;
}

my @vj;
foreach my $v (@v) {
    foreach my $j (@j) {
	push (@vj, "$v:$j");
    }
}

@vj = sort { $vjtn{$b} <=> $vjtn{$a} } @vj;

# clone frequency

my %ctn;
my %csn;
my %csp;

for (my $i = 0; $i < @f; $i++) { 
    open IN, "<$f[$i]/clone.cnpc" || die "open $f[$i]/clone.cnpc: $!\n";
    while (<IN>) {
	my @a = split "\t"; chomp $a[-1];
	$ctn{$a[0]} += $a[1];
	$csn{$a[0]}{$s[$i]} = $a[1];
	$csp{$a[0]}{$s[$i]} = $a[2];
    }
    close IN;
}

my @c = sort { $ctn{$b} <=> $ctn{$a} } keys %ctn;

foreach my $c (@c) {
    foreach my $s (@s) {
	$csn{$c}{$s} = 0 if !defined( $csn{$c}{$s} );
	$csp{$c}{$s} = 0 if !defined( $csp{$c}{$s} );
    }
}

# AA frequency

my %aatn;
my %aasn;
my %aasp;

for (my $i = 0; $i < @f; $i++) { 
    open IN, "<$f[$i]/clone.aanpc" || die "open $f[$i]/clone.aanpc: $!\n";
    while (<IN>) {
	my @a = split "\t"; chomp $a[-1];
	$aatn{$a[0]} += $a[1];
	$aasn{$a[0]}{$s[$i]} = $a[1];
	$aasp{$a[0]}{$s[$i]} = $a[2];
    }
    close IN;
}

my @aa = sort { $aatn{$b} <=> $aatn{$a} } keys %aatn;

foreach my $aa (@aa) {
    foreach my $s (@s) {
	$aasn{$aa}{$s} = 0 if !defined( $aasn{$aa}{$s} );
	$aasp{$aa}{$s} = 0 if !defined( $aasp{$aa}{$s} );
    }
}

# CDR3 length

my %cdr3nlsp;

for (my $i = 0; $i < @f; $i++) { 
    open IN, "<$f[$i]/clone.cdr3nlpc" || die "open $f[$i]/clone.cdr3nlpc: $!\n";
    while (<IN>) {
	my @a = split "\t"; chomp $a[-1];
	$cdr3nlsp{$a[0]}{$s[$i]} = $a[2];
    }
    close IN;
}

my @cdr3nl = sort { $a <=> $b } keys %cdr3nlsp;

foreach my $nl (@cdr3nl) {
    foreach my $s (@s) {
	$cdr3nlsp{$nl}{$s} = 0 if !defined($cdr3nlsp{$nl}{$s});
    }
}

# read statistics

my %sinn;
my %sregn;
my %spcdr3n;

for (my $i = 0; $i < @f; $i++) { 
    my $tmp = `tail -1 $f[$i]/clone.div`;
    my @a = split /\s+/, $tmp;
    $sinn{$s[$i]} = $a[0];
    $sregn{$s[$i]} = $a[1];
    $spcdr3n{$s[$i]} = $a[2];
}


##### output collected data #####

# V gene usage

open OUT, ">$prefix.vp" || die "open $prefix.vp: $!\n";
print OUT "V\t" . join("\t", @s) ."\n";
foreach my $v (@v) {
    print OUT "$v\t" . join("\t", map($vsp{$v}{$_}, @s)) ."\n";
}
close OUT;

# J gene usage

open OUT, ">$prefix.jp" || die "open $prefix.jp: $!\n";
print OUT "J\t" . join("\t", @s) ."\n";
foreach my $j (@j) {
    print OUT "$j\t" . join("\t", map($jsp{$j}{$_}, @s)) ."\n";
}
close OUT;

# VJ number and frequency

open OUT, ">$prefix.vjn" || die "open $prefix.vjn: $!\n";
print OUT "VJ\t" . join("\t", @s) ."\n";
foreach my $vj (@vj) {
    print OUT "$vj\t" . join("\t", map($vjsn{$vj}{$_}, @s)) ."\n";
}
close OUT;

open OUT, ">$prefix.vjp" || die "open $prefix.vjp: $!\n";
print OUT "VJ\t" . join("\t", @s) ."\n";
foreach my $vj (@vj) {
    print OUT "$vj\t" . join("\t", map($vjsp{$vj}{$_}, @s)) ."\n";
}
close OUT;

# clone number and frequency

open OUT, ">$prefix.cn" || die "open $prefix.cn: $!\n";
print OUT "Clone\t" . join("\t", @s) ."\n";
foreach my $c (@c) {
    print OUT "$c\t" . join("\t", map($csn{$c}{$_}, @s)) ."\n";
}
close OUT;

open OUT, ">$prefix.cp" || die "open $prefix.cp: $!\n";
print OUT "Clone\t" . join("\t", @s) ."\n";
foreach my $c (@c) {
    print OUT "$c\t" . join("\t", map($csp{$c}{$_}, @s)) ."\n";
}
close OUT;

# AA number and frequency

open OUT, ">$prefix.aan" || die "open $prefix.aan: $!\n";
print OUT "AA\t" . join("\t", @s) ."\n";
foreach my $aa (@aa) {
    print OUT "$aa\t" . join("\t", map($aasn{$aa}{$_}, @s)) ."\n";
}
close OUT;

open OUT, ">$prefix.aap" || die "open $prefix.aap: $!\n";
print OUT "AA\t" . join("\t", @s) ."\n";
foreach my $aa (@aa) {
    print OUT "$aa\t" . join("\t", map($aasp{$aa}{$_}, @s)) ."\n";
}
close OUT;

# CDR3 length frequency

open OUT, ">$prefix.cdr3nlp" || die "open $prefix.cdr3nlp: $!\n";
print OUT "CDR3L\t" . join("\t", @s) ."\n";
foreach my $nl (@cdr3nl) {
    print OUT "$nl\t" . join("\t", map($cdr3nlsp{$nl}{$_}, @s)) ."\n";
}
close OUT;

# read statistics

open OUT, ">$prefix.read_stat" || die "open $prefix.read_stat: $!\n";
print OUT "Sample\tInput\tRegular\t%\tProductive_CDR3\t%\n";
foreach my $s (@s) {
    my $regp = sprintf "%.2f", $sregn{$s}/$sinn{$s}*100;
    my $pdr3p = sprintf "%.2f", $spdr3n{$s}/$sregn{$s}*100;
    print OUT "$s\t$sinn{$s}\t$sregn{$s}\t$regp\t$spcdr3n{$s}\t$pcdr3p\n";
}
close OUT;
