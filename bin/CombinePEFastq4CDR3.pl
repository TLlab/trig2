#!/usr/bin/perl -w
# usage : CombinePEFastq4CDR3.pl um_read.vdjdelta um1_read.fq um2_read.fq

use strict;
use File::Basename qw(dirname);
use Cwd qw(abs_path);
use lib dirname(abs_path $0);
use Fastx qw(LoadOneFastq);


# load combined VDJ delta and the read 1 and 2

open IN, "<$ARGV[0]" || die "open $ARGV[0]: $!\n";
open my $fh1, "<$ARGV[1]" || die "open $ARGV[1]: $!\n";
open my $fh2, "<$ARGV[2]" || die "open $ARGV[2]: $!\n";

# for each paired-end
while (<IN>) {
    my @a = split "\t"; chomp $a[-1];
    my ($q1, $seq1, $qual1) = LoadOneFastq($fh1);
    my ($q2, $seq2, $qual2) = LoadOneFastq($fh2);
    
    # print the major read
    if ($a[2] eq "---" || $a[10] == 2) {
	print "\@$q1\n$seq1\n+\n$qual1\n";
    } else {
	print "\@$q2\n$seq2\n+\n$qual2\n";
    }
}
close IN;
close $fh1;
close $fh2;
