#!/usr/bin/perl -w
# usage : CorrectError.pl read.cdr3 > output
# to do : filter low qulity, clustering, error corrcet and convet to vdjtools format

use 5.010;
use strict;
use POSIX;
use List::MoreUtils qw(indexes);

##### find core clone and defer low quality reads for further processing

# clone : clone count
# lqseq : sequence contains low quality
# defer : deferred read count
my %clone;
my %lqseq;
my %defer;

open IN, "<$ARGV[0]" || die "$!\n";

while (<IN>) {
	my @F = split "\t"; chomp $F[-1];
	
	# skip if no alignment or non-regular
	next if $F[1] eq "---" || $F[1] != 2;
	
	next if grep(/---/, @F[2,3]);         # skip ---
	my $qua = (split '\|', $F[3])[0];     # use first qua
	my @seq =  split '\|', $F[2];
	my $wt = 1/@seq;

	# position of low quality (q<10) bases
	my @lqb = indexes { $_ - 33 < 10 } unpack("C*", $qua);

	# core clone if no low quality base
	unless (@lqb) {
		$clone{$_} += $wt foreach @seq;
	}

	# if < 70% of the bases are of low quality
	elsif (@lqb / length($qua) < 0.7) {
		$lqseq{$F[0]} = {
			WT  => $wt,
			SEQ => \@seq,
			LQB => \@lqb,
		};
	}
}
close IN;

foreach my $id (keys %lqseq){
	foreach (@{ $lqseq{$id}->{SEQ} }) {

		# accept it if the clone exists
		if (exists $clone{$_}) {
			$clone{$_} += $lqseq{$id}->{WT};
		}

		# otherwise, mask the low quality bases (i.e., make them Ns)
		else {
			my @F = split ":";
			substr($F[1], $_, 1) = 'N' foreach @{ $lqseq{$id}->{LQB} };
			$defer{"$F[0]:$F[1]:$F[2]"} += $lqseq{$id}->{WT};
		}
	}
}

##### rescue deferred reads using usearch

# prepare core clone sequences and deferred reads
open DB, ">core.fa" || die "$!\n";
open QY, ">defer.fa" || die "$!\n";

foreach my $c (sort { $clone{$b} <=> $clone{$a} } keys %clone) {
	state $id++;
	my @F = split ":", $c;
	print DB ">core$id\t$clone{$c}\t$F[0]:$F[2]\n$F[1]\n";
}

foreach my $d (sort { $defer{$b} <=> $defer{$a} } keys %defer) {
	state $id++;
	my @F = split ":", $d;
	print QY ">defer$id\t$defer{$d}\t$F[0]:$F[2]\n$F[1]\n";
}

close DB;
close QY;

# align deferred reads to core clone sequencs
`usearch -usearch_global defer.fa -db core.fa -strand plus -id 1.0 -userout hits.ur -userfields query+target+ql+tl+tseq 2>&1`;

# assign the masked deferred reads to the perfectly aligned core clone
# if they are of the same VJ pair and length
open IN, "<hits.ur" || die "$!\n"; 
while (<IN>) {
	chomp;
	my @F = split /[\t\:]/;
	$F[10] = uc $F[10];
	$clone{"$F[6]:$F[10]:$F[7]"} += $F[1] if (@F[2,3,8] eq @F[6,7,9]);
}
close IN;
unlink("core.fa", "defer.fa", "hits.ur");

##### remove errors via clustering clones 

# vjrid : clone IDs of a VJ pair
my %ridkv;
my %vjrid;

foreach my $c (sort { $clone{$b} <=> $clone{$a} } keys %clone) {
	state $id = 0;
	$id++;
	my @F = split ":", $c;    # v:seq:j
	$ridkv{$id} = {
		VJ  => "$F[0]:$F[2]",
		SIZ => $clone{$c},
		LEN => length($F[1]),
		SEQ => $F[1],
	};
	push(@{ $vjrid{"$F[0]:$F[2]"} }, $id);
}

# for each VJ pair
for my $vj (keys %vjrid){
	my @child;

	# cluster until no clone remains
	while (@{ $vjrid{$vj} }) {
		my $id = shift(@{ $vjrid{$vj} });

		# stop cluster if the clone size is <100, i.e., cannot have child 
		last if $ridkv{$id}->{SIZ} < 100;

		# identify potential child clone
		my @qrid = grep($ridkv{$_}->{SIZ} < ($ridkv{$id}->{SIZ}/100) &&
			abs($ridkv{$_}->{LEN} - $ridkv{$id}->{LEN}) <= 1 &&
			!$ridkv{$_}->{PAR}, @{ $vjrid{$vj} });

		# skip if no potential child clone
		next unless @qrid;

		# prepare parent and child sequences
		open DB, ">par.fa" || die "$!\n";
		open QY, ">chi.fa" || die "$!\n";
		print DB ">$id\n$ridkv{$id}->{SEQ}\n";
		print QY ">$_\n$ridkv{$_}->{SEQ}\n" for @qrid;
		close DB;
		close QY;

		# align child to parent sequences
		`usearch -usearch_global chi.fa -db par.fa -strand plus -id 0.9 -gapopen 2.0 -userout hits.ur -userfields query+target+id+mism+gaps 2>&1`;

		# assign child to the parent if the edit distance is less than 1
		open IN, "hits.ur" || die "$!\n";
		while (<IN>) {
			my @F = split "\t"; chomp $F[-1];
			my $ed = $F[3]+$F[4];
			$ridkv{$F[0]}->{PAR} = $F[1], push(@child, $F[0]) if $ed <= 1;
		}
		close IN;
	}

	# skip if no child
	next unless @child;

	# increase size of the parent clone by the child size
	for my $chi (reverse @child){
		my $par = $ridkv{$chi}->{PAR};
		$ridkv{$par}->{SIZ} += $ridkv{$chi}->{SIZ};
	}
}
unlink ("par.fa","chi.fa","hits.ur");


##### put clone info in vdjtools format

my $total_count = 0;
my $TRBD1 = "ACAGGG";
my $TRBD2 = "ACTAGC";
my %aacode = ( 
	TTT => "F", TTC => "F", TTA => "L", TTG => "L",
	TCT => "S", TCC => "S", TCA => "S", TCG => "S",
	TAT => "Y", TAC => "Y", TAA => "*", TAG => "*",
	TGT => "C", TGC => "C", TGA => "*", TGG => "W",
	CTT => "L", CTC => "L", CTA => "L", CTG => "L",
	CCT => "P", CCC => "P", CCA => "P", CCG => "P",
	CAT => "H", CAC => "H", CAA => "Q", CAG => "Q",
	CGT => "R", CGC => "R", CGA => "R", CGG => "R",
	ATT => "I", ATC => "I", ATA => "I", ATG => "M",
	ACT => "T", ACC => "T", ACA => "T", ACG => "T",
	AAT => "N", AAC => "N", AAA => "K", AAG => "K",
	AGT => "S", AGC => "S", AGA => "R", AGG => "R",
	GTT => "V", GTC => "V", GTA => "V", GTG => "V",
	GCT => "A", GCC => "A", GCA => "A", GCG => "A",
	GAT => "D", GAC => "D", GAA => "E", GAG => "E",
	GGT => "G", GGC => "G", GGA => "G", GGG => "G",
);

my @rid = grep(!$ridkv{$_}->{PAR}, (1..keys %clone));

# round up
for (@rid) {
	$ridkv{$_}->{SIZ} = ceil($ridkv{$_}->{SIZ});
	$total_count += $ridkv{$_}->{SIZ};
}

print "count\tfreq\tcdr3nt\tcdr3aa\tv\td\tj\n"; # title

foreach (@rid){
	my $count = $ridkv{$_}->{SIZ};
	my $freq = $count / $total_count;
	my $seq = $ridkv{$_}->{SEQ};
	my @codons = unpack('(A3)*', $seq);
	my $cdr3aa = join '',  map { exists $aacode{$_} ? $aacode{$_} : "*" } @codons;
	my @F = split ':', $ridkv{$_}->{VJ};
	my $D = "."; $D = "TRBD1" if $seq =~ $TRBD1; $D = "TRBD2" if $seq =~ $TRBD2;
	print "$count\t$freq\t$seq\t$cdr3aa\t$F[0]\t$D\t$F[1]\n";
}
