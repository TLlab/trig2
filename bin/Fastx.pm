#!/usr/bin/perl -w

package Fastx;

use strict;
use Exporter qw(import);

our @EXPORT_OK = qw(Fastq2Fasta Fasta2Fastq SplitFastq SplitFasta LoadFasta LoadThisFasta LoadOneFasta LoadThisFastq LoadOneFastq);


#################### define subroutines ####################


sub Fastq2Fasta {
    my ($ifn, $ofn) = @_;
    
    open IN, "<$ifn" || die "open $ifn: $!\n";
    open OUT, ">$ofn" || die "open $ofn: $!\n";

    while (<IN>) {
	$_ = ">" . substr( $_, 1 );
	print OUT $_;
	
	$_ = <IN>;
	print OUT $_;
	
	<IN>;
	<IN>;
    }
    close IN;
    close OUT;
}


sub Fasta2Fastq {
    my ($ifn, $ofn) = @_;
    
    open IN, "<$ifn" || die "open $ifn: $!\n";
    open OUT, ">$ofn" || die "open $ofn: $!\n";

    while (<IN>) {
	my $id = substr($_, 1); chomp $id;
	my $seq;
	while (<IN>) {
	    if (/^>/) {
		seek(IN, -length($_), 1);
		last;
	    }
	    chomp;
	    $seq .= $_;
	}
	my $qual = 'I' x length($seq);
	print OUT "\@$id\n$seq\n+\n$qual\n";
    }
    close IN;
    close OUT;
}


sub SplitFastq {
    my ($fn, $n) = @_;

    # set prefix
    my $bn = `basename $fn`;          
    my ($pre) = $bn =~ /(.+)\.f(ast)*q/;

    # link file if no need to split
    if( $n == 1 ) {
	`ln -s $fn $pre.1.fq`;
	return;
    }

    # get number of sequences per split file
    my $nos = `wc -l $fn`;
    ($nos) = $nos =~ /^(\d+)/; $nos /= 4;

    # when too few reads
    if ($nos < $n) {
	`cp $fn $pre.1.fq`;
	foreach my $i (2..$n) {
	    `echo -n > $pre.$i.fq`;
	}
	return;
    }
    
    $nos = int($nos/$n);
    
    # load sequences and output
    my $k = 1;
    open OUT, ">$pre.$k.fq" || die "open $pre.$k.fq: $!\n";
    
    my $i = 0;
    open IN, "<$fn" || die "open $fn: $!\n";
    while (<IN>) {
	$i++;
	if ($i > $nos) {
	    if ($k == $n) {
		print OUT $_;
		last;
	    }
	    close OUT;
	    $k++;
	    open OUT, ">$pre.$k.fq" || die "open $pre.$k.fq: $!\n";
	    $i = 1;
	}
	print OUT $_;
	$_ = <IN>; print OUT $_;
	$_ = <IN>; print OUT $_;
	$_ = <IN>; print OUT $_;
    }

    # output the remaining reads
    while (<IN>) {
	print OUT $_;
    }
    close IN;
    close OUT;
}


sub SplitFasta {
    my ($fn, $n) = @_;

    # set prefix
    my $bn = `basename $fn`;          
    my ($pre) = $bn =~ /(.+)\.f(ast)*a/;

    # link file if no need to split
    if( $n == 1 ) {
	`ln -s $fn $pre.1.fa`;
	return;
    }

    # count number of sequences
    my $nos = `grep -c \\> $fn`;
    ($nos) = $nos =~ /(\d+)/;

    # when too few reads
    if ($nos < $n) {
	`cp $fn $pre.1.fa`;
	foreach my $i (2..$n) {
	    `echo -n > $pre.$i.fa`;
	}
	return;
    }

    $nos = int($nos/$n);

    # load sequences and output
    my $k = 1;
    open OUT, ">$pre.$k.fa" || die "open $pre.$k.fa: $!\n";

    my $i = 0;
    open IN, "<$fn" || die "open $fn: $!\n";
    while( <IN> ) {
	$i++ if /^>/;
	if( $i > $nos ) {
	    if ($k == $n) {
		print OUT $_;
		last;
	    }
	    close OUT;
	    $k++;
	    open OUT, ">$pre.$k.fa" || die "open $pre.$k.fa: $!\n";
	    $i = 1;
	}
	print OUT $_;
    }

    # output the remaining reads
    while (<IN>) {
	print OUT $_;
    }
    close IN;
    close OUT;
}


sub LoadFasta {
    my ($fn, $pidseq) = @_;
    
    open IN, "<$fn" || die "open $fn: $!\n";

    my $id;
    while (<IN>) {
	if (/^>(\S+)/) {
	    $id = $1;
	} else {
	    chomp $_;
	    $pidseq->{$id} .= $_;
	}
    }
    close IN;
}


sub LoadThisFasta {
    my ($fh, $q, $pq) = @_;
    my $seq = "";
    
    while (<$fh>) {
	if (/^>(\S+)/) {
	    my $id = $1;

	    # load sequence
	    while (<$fh>) {
		if (/^>/) {
		    seek($fh, -length($_), 1);
		    last;
		} else {
		    chomp;
		    $seq .= $_;
		}
	    }

	    # if the desired query
	    if ($id eq $q) {
		return $seq;

	    # output if not the desired
	    } else {
		if ($pq == 1) {
		    my $l = length($seq);
		    print "$id\t$l\t---\n";
		}
		$seq = "";
	    }
	}
    }
}


sub LoadOneFasta {
    my ($fh) = $_[0];
    my ($id, $seq);
    
    while (<$fh>) {
	if (/^>(\S+)/) {
	    $id = $1;

	    # load sequence
	    while (<$fh>) {
		if (/^>/) {
		    seek($fh, -length($_), 1);
		    last;
		} else {
		    chomp;
		    $seq .= $_;
		}
	    }
	    last;
	}
    }

    return ($id, $seq);
}


sub LoadThisFastq {
    my ($fh, $id) = @_;
    my ($seq, $qual);

    while (<$fh>) {
        if (/^\@$id\s+/) {
            $seq = <$fh>; chomp $seq;
            <$fh>;
            $qual = <$fh>; chomp $qual;
            last;
        } else {
            <$fh>;
            <$fh>;
            <$fh>;
        }
    }

    return ($seq, $qual);
}


sub LoadOneFastq {
    my ($fh) = $_[0];
    my ($id, $seq, $qual);

    $_ = <$fh>;
    ($id) = /^\@(\S+)/;
    $seq = <$fh>; chomp $seq;
    <$fh>;
    $qual = <$fh>; chomp $qual;

    return ($id, $seq, $qual);
}


sub LoadRefSeq {
    my $fn = $_[0];
    my $seq;

    open IN, "<$fn" || die "open $fn: $!\n";
    <IN>;
    while (<IN>) {
	chomp;
	$seq .= $_;
    }
    close IN;

    return $seq;
}
