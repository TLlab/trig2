#!/usr/bin/perl -w

use strict;


# load html and get amino acid and codon

my @aa;
my @codon;

open IN, "<$ARGV[0]" || die "open $ARGV[0]: $!\n";

while (<IN>) {

    # at the number line
    if (/<span class="strand">\s*\d+\s*<\/span>/) {

	# get amino acid
	$_ = <IN>;
	/(<span class=.+<\/span>)/;
	my @a = split /<\/span>/, $1;
	foreach my $e (@a) {
	    if ($e =~ /([A-Z\s])$/) {
		my $aa = $1 eq " " ? "." : $1;
		push(@aa, $aa);
	    }
	}

	# get codon
	$_ = <IN>; chomp;
	/(<span class=.+<\/span>)/;
	@a = split /> </, $1;
	foreach my $e (@a) {
	    my @n = $e =~ />(.)<\/span/g;
	    push(@codon, join("", @n));
	}
    }
}

print scalar(@aa) . "\n";
print join("", @aa) . "\n";
print scalar(@codon) . "\n";
print join("", @codon) . "\n";
