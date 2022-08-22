#!/usr/bin/perl -w
# usage : IMGTAlignment.pl species group

use strict;
use LWP::Simple qw(get);


# set parameters

my %sfn = ("hsa" => "Homo%20sapiens", "mmu" => "Mus%20musculus");

my $species = $sfn{$ARGV[0]};
my $group = uc($ARGV[1]);


##### get alignments of all V genes

open OUTvaa, ">$ARGV[0]\_$ARGV[1]\_v_aa.aln" || die "open $ARGV[0]\_$ARGV[1]\_v_aa.aln: $!\n";
open OUTvnt, ">$ARGV[0]\_$ARGV[1]\_v_nt.aln" || die "open $ARGV[0]\_$ARGV[1]\_v_nt.aln: $!\n";

# get html and select the line with gene names
my $vg = $group . "V";
my $vurl = "http://www.imgt.org/IMGTrepertoire/Proteins/alleles/list_alleles.php?species=$species&group=$vg";
my $vhtml = get($vurl);
my @h = split "\n", $vhtml;
@h = grep { /<li>/ && /$vg/ } @h; 

# get all gene names
my @avg;
my ($list) = $h[0] =~ /(<li>.+)<\/a><\/li>/;
my @a = split /<\/a>/, $list;
foreach my $e (@a) {
    my ($g) = $e =~ /([^>]+)$/;
    push(@avg, $g) if $g !~ /OR/ || $g =~ /TRBV21/;     # TRBV21/OR9-2 is used for TRBV21
}
sleep(3);

# get alignments for each gene
foreach my $g (@avg) {
    my $cg = $g;
    $cg =~ s/ /\%20/g;
    my $url = "http://www.imgt.org/IMGTrepertoire/Proteins/alleles/index.php?species=$species&group=$vg&gene=$cg";
    my $html = get($url);
    my ($aa, $nt) = ParseVAlignment($html);

    # take care special genes
    if ($g =~ /TRBV21/) {
	$g = "TRBV21-1";
	$nt =~ s/aaaagacat/aaagcacat/;
	$nt =~ s/gccaacagcaaagc/gccagcagcaaagc/;
    } elsif ($g eq "IGKV3D-11") { 
	$nt =~ s/cagagtgtt/cagggtgtt/;
    }

    print "done: $g\n";
    print OUTvaa ">$g\n$aa\n";
    print OUTvnt ">$g\n$nt\n";
    sleep(3);
}

close OUTvaa;
close OUTvnt;


##### get alignments of all J genes

my @ajg;
my @ajaa;
my @ajnt;

my $jg = $group . "J";
my $jurl = "http://www.imgt.org/IMGTrepertoire/Proteins/alleles/index.php?species=$species&group=$jg&gene=$jg-overview";
my $jhtml = get($jurl);

my @line = split "\n", $jhtml;
for (my $i = 0; $i < @line; $i++) {
    
    if ($line[$i] =~ /<span class="fgxg">/) {
	$line[$i] = ">" . $line[$i] . "<";
	my (@aa) = $line[$i] =~ />([A-Z\s\*]+)</g;
	$_ =~ s/ //g for @aa;
	$i++;

	my ($j) = $line[$i] =~ />\s*($group.+)\*01<\/span><\/span>/;
	my ($nt) = $line[$i] =~ /gDNA\s+([a-z\s]+)\s*$/;
	my @nt = split " ", $nt;
	if (length($nt[0]) < 3) {
	    $nt[0] = ('.' x (3 - length($nt[0]))) . $nt[0];
	    $aa[0] = '.' . $aa[0];
	}
	if (length($nt[-1]) < 3) {
	    pop(@nt);
	}

	# take care of TRGJ alignments (not ending at the same A.A.)
	if ($j =~ /^TRGJ.$/) {
	    push(@aa, (".", "."));
	    push(@nt, ("...", "..."));
	}
	
	push(@ajg, $j);
	push(@ajaa, join("", @aa));
	push(@ajnt, join("", @nt));

    # for IG*J
    } elsif ($line[$i] =~ />($group.+)\*01<\/FONT>/) {
	my $j = $1;
	my ($nt) = $line[$i] =~ /<\/FONT>(.+)\s*$/;
	$nt =~ s/\&nbsp;/ /g;
	$nt =~ s/^\s+//;
	my @nt = split " ", $nt;

	$line[$i-1] =~ s/^<FONT SIZE=-1>//;
	$line[$i-1] =~ s/\&nbsp;/ /g;
	$line[$i-1] = ">" . $line[$i-1] . "<"; 
	my (@aa) = $line[$i-1] =~ />([A-Z\s]+)</g;
	$_ =~ s/ //g for @aa;

	if (length($nt[0]) < 3) {
	    $nt[0] = ('.' x (3 - length($nt[0]))) . $nt[0];
	    $aa[0] = '.' . $aa[0];
	}

	push(@ajg, $j);
	push(@ajaa, join("", @aa));
	push(@ajnt, join("", @nt));
    }
}

my @jaal = sort { $b <=> $a } map(length($_), @ajaa);
my $maxjaal = $jaal[0];
foreach my $aa (@ajaa) {
    $aa = ('.' x ($maxjaal - length($aa))) . $aa;
}

my @jntl = sort { $b <=> $a } map(length($_), @ajnt);
my $maxjntl = $jntl[0];
foreach my $nt (@ajnt) {
    $nt = ('.' x ($maxjntl - length($nt))) . $nt;
}

open OUTjaa, ">$ARGV[0]\_$ARGV[1]\_j_aa.aln" || die "open $ARGV[0]\_$ARGV[1]\_j_aa.aln: $!\n";
open OUTjnt, ">$ARGV[0]\_$ARGV[1]\_j_nt.aln" || die "open $ARGV[0]\_$ARGV[1]\_j_nt.aln: $!\n";

for (my $i = 0; $i < @ajg; $i++) {
    print OUTjaa ">$ajg[$i]\n$ajaa[$i]\n";
    print OUTjnt ">$ajg[$i]\n$ajnt[$i]\n";
}

close OUTjaa;
close OUTjnt;



######################################################################


sub ParseVAlignment {
    my $html = $_[0];
    my @aa;
    my @codon;
    
    my @h = split "\n", $html;
    for(my $i = 0; $i < @h; $i++) {
	
	# at the number line
	if ($h[$i] =~ /<span class="strand">\s*\d+\s*<\/span>/) {
	    
	    # get amino acids in the next line
	    $i++;
	    $h[$i] =~ /(<span class=.+<\/span>)/;
	    my @a = split /<\/span>/, $1;
	    foreach my $e (@a) {
		if ($e =~ /([A-Z\s])$/) {
		    my $aa = $1 eq " " ? "." : $1;
		    push(@aa, $aa);
		} elsif ($e =~ /\*<\/a>$/) {
		    push(@aa, "*");
		}
	    }
	    
	    # then get codons in the next line
	    $i++;
	    $h[$i] =~ /(<span class=.+<\/span>)/;
	    @a = split /> </, $1;
	    foreach my $e (@a) {
		my @n = $e =~ />(.)<\/span/g;
		push(@codon, join("", @n));
	    }
	}
    }

    my $aa = join("", @aa);
    my $nt = join("", @codon);

    return ($aa, $nt);
}
