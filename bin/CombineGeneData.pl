#!/usr/bin/perl -w
# usage : CombineGeneData.pl g1.prefix g2.prefix ...

use strict;


# set file handles

my @g;
my %gvdjdeltafh;
my %gcdr3fh;

foreach my $fp (@ARGV) {
    my ($g) = $fp =~ /^([^\.]+)/;
    push(@g, $g);

    open $gvdjdeltafh{$g}, "<$fp.vdjdelta" || die "open $fp.vdjdelta: $!\n";
    open $gcdr3fh{$g}, "<$fp.cdr3" || die "open $fp.cdr3: $!\n";
}

my ($prefix) = $ARGV[0] =~ /$g[0]\.(.+)/;

open OUTvdjdelta, ">$prefix.vdjdelta" || die "open $prefix.vdjdelta: $!\n";
open OUTcdr3, ">$prefix.cdr3" || die "open $prefix.cdr3: $!\n";


# load VDJ delta and CDR3

while (readline($gvdjdeltafh{$g[0]})) {
    my @a = split "\t"; chomp $a[-1];
    
    my %gkv;
    $gkv{$g[0]}{vdjdelta} = $a[2] eq "---" ? join("\t", @a[0..2]) : join("\t", @a[0..9]);
    $gkv{$g[0]}{al} = $a[2] eq "---" ? 0 : $a[9];

    foreach my $g (@g[1..$#g]) {
	$_ = readline($gvdjdeltafh{$g});
	my @a = split "\t"; chomp $a[-1];
	
	$gkv{$g}{vdjdelta} = $a[2] eq "---" ? join("\t", @a[0..2]) : join("\t", @a[0..9]);
	$gkv{$g}{al} = $a[2] eq "---" ? 0 : $a[9];
    }

    foreach my $g (@g) {
	$gkv{$g}{cdr3} = readline($gcdr3fh{$g});
    }

    # find the gene that gives the longest alignment
    my @sg = sort { $gkv{$b}{al} <=> $gkv{$a}{al} } @g;
    
    # output data
    print OUTvdjdelta $gkv{$sg[0]}{vdjdelta} . "\n";
    print OUTcdr3 $gkv{$sg[0]}{cdr3};
}

foreach my $g (@g) {
    close $gvdjdeltafh{$g};
    close $gcdr3fh{$g};
}

close OUTvdjdelta;
close OUTcdr3;
