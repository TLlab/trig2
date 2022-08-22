#!/usr/bin/perl -w
# usage : CombinePECDR3.pl read1.cdr3 read2.cdr3

use strict;


# load CDR3 info of read 1 and read 2

open IN1, "<$ARGV[0]" || die "open $ARGV[0]: $!\n";
open IN2, "<$ARGV[1]" || die "open $ARGV[1]: $!\n";

while (<IN1>) {
    my $cdr31 = $_;
    my $cdr32 = <IN2>;

    if ($cdr31 eq $cdr32) {
        print $cdr31;
    } else {
        my @a1 = split "\t", $cdr31; chomp $a1[-1];
        my @a2 = split "\t", $cdr32; chomp $a2[-1];

        if (!$a2[2]) {
            print $cdr31;
            print $cdr32;
            exit 0;
        }
        
        if ($a1[2] eq "---" && $a2[2] eq "---") {
            if ($a1[1] ne "---") {
                print $cdr31;
            } else {
                print $cdr32;
            }
        } elsif ($a1[2] ne "---" && $a2[2] ne "---") {
            print $cdr31;
        } else {
            if ($a1[2] ne "---") {
                print $cdr31;
            } else {
                print $cdr32;
            }
        }
    }
}
