#!/usr/bin/perl -w
# usage : trig.pl [option] read.fa/q [read2.fq]
# e.g.  : trig.pl -s hsa -g trb -o test read_1.fq read_2.fq

use strict;
use File::Basename qw(dirname basename);
use Cwd qw(abs_path);
use Getopt::Long qw(GetOptions);

use lib dirname(abs_path($0));
use Fastx qw(Fastq2Fasta Fasta2Fastq SplitFastq SplitFasta);

# to do : try on different genes, e.g., igh
# to do : add utilities (calculate diversity, D50, hydrophobicity, etc)
# to do : table and graphics (use VDJtools)


############################## set parameters ##############################

# minmatch : minimal match for nucmer alignment
# adjolq   : adjust overalp Q

my $species  = "hsa";
my $gene     = "trb";
my $outdir   = "";
my $minmatch = 15;
my $mergeq   = 1;
my $adjolq   = 1;
my $patchq   = 0;
my $frac     = 0.5;
my $thread   = 1;
my $help;

GetOptions(
    "species=s"  => \$species,
    "gene=s"     => \$gene,
    "outdir=s"   => \$outdir,
    "minmatch=i" => \$minmatch,
    "mergeq=i"   => \$mergeq,
    "adjolq=i"   => \$adjolq,
    "patchq=i"   => \$patchq,
    "frac=f"     => \$frac,
    "thread=i"   => \$thread,
    "help"       => \$help,
    );

$help = 1 if !@ARGV;

Usage(), exit 0 if $help;


##### confirm input and output

# single-end reads can be either fasta or fastq
die "input file error!\n" if !-e $ARGV[0] || $ARGV[0] !~ /[qa]$/;

# paired-end fastq
my $peq = 1 if $ARGV[1] ? 1 : 0;
if ($peq) {
    die "input file error!\n" if $ARGV[0] !~ /q$/ || !-e $ARGV[1] || $ARGV[1] !~ /q$/;
}

# set output directory
if (!$outdir) {
    my $fn = basename $ARGV[0];

    my $fp;
    if (!$peq) {
	($fp) = $fn =~ /(.+)\.f/;
    } else {
	($fp) = $fn =~ /(.+)\_(R)?1\.f/;
    }
    $fp = "output" if !$fp;

    $outdir = "trig_$fp";
}

die "output folder $outdir exists!\n" if -d $outdir;

# get file info (handle linked files)
my $file = -l $ARGV[0] ? readlink($ARGV[0]) : abs_path($ARGV[0]);
my $fn = basename($file);
my $fd = dirname(abs_path($file));
my $f2n;
my $f2d;

if ($peq) {
    my $file2 = -l $ARGV[1] ? readlink($ARGV[1]) : abs_path($ARGV[1]);
    $f2n = basename($file2);
    $f2d = dirname(abs_path($file2));
}

# confirm species and gene
Usage() if $species ne "hsa" && $species ne "mmu";

my @gene = split ",", $gene;
if ($gene eq "all") {
    @gene = qw(trad trb trg igh igl igk);
} elsif($gene eq "tr") {
    @gene = qw(trad trb trg);
} elsif($gene eq "ig") {
    @gene = qw(igh igl igk);
} else {
    foreach my $g (@gene) {
	$g = "trad" if $g eq "tra" || $g eq "trd";
	Usage() if $g !~ /^tr[bg]$/ && $g ne "trad" && $g !~ /^ig[hlk]$/;
    }
}
$gene = join("_", @gene);

my $sg = "$species\_$gene";

# get directory of TRIg
my $trigdir = dirname(dirname(abs_path($0)));


############################## prepare data ##############################

# create output directory and log
`mkdir $outdir`;
chdir($outdir);

my $date = `date`;
open  LOG,">log" || die "open log: $!\n";
print LOG "input file name   : @ARGV\n";
print LOG "species           : $species\n";
print LOG "gene              : $gene\n";
print LOG "nucmer minmatch   : $minmatch\n";
print LOG "adjust overlap    : $adjolq\n";
print LOG "threads           : $thread\n";
print LOG "program start     : $date";


##### link and preprocess data

# link reference data
if (@gene == 1) {
    `ln -s $trigdir/gene/$sg.$_` for ("fa", "vdj", "cdr");
    if ($patchq) {
	`ln -s $trigdir/gene/$sg\_$_.fa` for ("j", "c");
    }
} else {
    foreach my $ext ("fa", "vdj", "cdr") {
	my $command = "cat " . join(" ", map("$trigdir/gene/$species\_$_.$ext", @gene)) . " > $sg.$ext";
	`$command`;
    }
    if ($patchq) {
	foreach my $g (@gene) {
	    `ln -s $trigdir/gene/$species\_$g\_$_.fa` for ("j", "c");
	}
    }
}

# default data type
my $ext = "fq";

# if single-end
if (!$peq) {

    # link input sequence
    $ext = "fa" if $fn =~ /a$/;
    `ln -s $fd/$fn read.$ext`;

    # convert into fasta for nucmer and split
    if ($ext eq "fq") {
	Fastq2Fasta("read.fq", "read.fa");
    } else {
	Fasta2Fastq("read.fa", "read.fq");
    }
    SplitFasta("read.fa", $thread);
    SplitFastq("read.fq", $thread);
    

# if paired-end
} else {

    # link paired-end input
    `ln -s $fd/$fn read_1.fq`;
    `ln -s $f2d/$f2n read_2.fq`;

    # merge paired-end reads
    if ($mergeq) {
        my $cmd = "usearch -fastq_mergepairs read_1.fq -reverse read_2.fq";
        $cmd .= " -threads $thread -fastq_minovlen 20 -fastq_minmergelen 30 -fastq_pctid 90 -fastq_maxdiffs 10 -fastq_trunctail 0";
        #$cmd .= " -threads $thread -fastq_minovlen 10 -fastq_pctid 75 -fastq_maxdiffs 75 -fastq_trunctail 0";
        $cmd .= " -fastqout read.fq -fastqout_notmerged_fwd um1_read.fq -fastqout_notmerged_rev um2_read.fq"; 
        $cmd .= " 2> usearch_fastqmergepairs.log";
        `$cmd`;
    } else {
        `ln -s read_1.fq um1_read.fq`;
        `ln -s read_2.fq um2_read.fq`;
    }

    # convert fastq into fasta then split file
    Fastq2Fasta("read.fq", "read.fa") if $mergeq;
    Fastq2Fasta("um1_read.fq", "um1_read.fa");
    Fastq2Fasta("um2_read.fq", "um2_read.fa");
    SplitFastq("read.fq", $thread) if $mergeq;
    SplitFastq("um1_read.fq", $thread);
    SplitFastq("um2_read.fq", $thread);
    SplitFasta("read.fa", $thread) if $mergeq;
    SplitFasta("um1_read.fa", $thread);
    SplitFasta("um2_read.fa", $thread);
}


############################## TRIg pipeline ##############################

# run TRIg in parallel

my @child = ();

for (my $i = 1; $i <= $thread; $i++) {
    my $pid = fork();
    
    if ($pid) {
        push(@child, $pid);
	
    } elsif ($pid == 0) {
	
	# analyze single or merged reads
        if ($peq == 0 || $mergeq) {

            # run nucmer while taking care of null file
            if (-s "read.$i.fa") {
                `nucmer --maxmatch -l $minmatch -c $minmatch -b $minmatch -p initial.$i $sg.fa read.$i.fa 2> /dev/null`;
            } else {
                my $command = "echo \"$sg.fa read.$i.fa\" > initial.$i.delta";
                `$command`;
                $command = "echo \"NUCMER\" >> initial.$i.delta";
                `$command`;
            }
            
            # run TRIg kernel
            if ($patchq) {
                `ProcessAlignment -s $species -g $gene -m $minmatch -a $adjolq -f $frac -o processed.$i initial.$i.delta`;
                `PatchAlignment.pl -s $species -g $gene -m $minmatch processed.$i.vdjdelta read.$i.fa`;
                `mv processed.$i.vdjdelta.1 read.$i.vdjdelta`;
                `mv processed.$i.cdr3.1 read.$i.cdr3`;
            } else {
                `ProcessAlignment -s $species -g $gene -m $minmatch -a $adjolq -f $frac -o read.$i initial.$i.delta`;
            }
        }
	
	# analyze un-merged paired-end reads
	if ($peq) {
	    if (-s "um1_read.$i.fa") {
		`nucmer --maxmatch -l $minmatch -c $minmatch -b $minmatch -p um1_initial.$i $sg.fa um1_read.$i.fa 2> /dev/null`;
		`nucmer --maxmatch -l $minmatch -c $minmatch -b $minmatch -p um2_initial.$i $sg.fa um2_read.$i.fa 2> /dev/null`;
	    } else {
		my $command = "echo \"$sg.fa um1_read.$i.fa\" > um1_initial.$i.delta"; `$command`;
		$command = "echo \"NUCMER\" >> um1_initial.$i.delta"; `$command`;
		$command = "echo \"$sg.fa um2_read.$i.fa\" > um2_initial.$i.delta"; `$command`;
		$command = "echo \"NUCMER\" >> um2_initial.$i.delta"; `$command`;
	    }
	    if($patchq) {
		`ProcessAlignment -s $species -g $gene -m $minmatch -a $adjolq -f $frac -o um1_processed.$i um1_initial.$i.delta`;
		`PatchAlignment.pl -s $species -g $gene -m $minmatch um1_processed.$i.vdjdelta um1_read.$i.fa`;
                `mv um1_processed.$i.vdjdelta.1 um1_read.$i.vdjdelta`;
                `mv um1_processed.$i.cdr3.1 um1_read.$i.cdr3`;
                `ProcessAlignment -s $species -g $gene -m $minmatch -a $adjolq -f $frac -o um2_processed.$i um2_initial.$i.delta`;
		`PatchAlignment.pl -s $species -g $gene -m $minmatch um2_processed.$i.vdjdelta um2_read.$i.fa`;
                `mv um2_processed.$i.vdjdelta.1 um2_read.$i.vdjdelta`;
                `mv um2_processed.$i.cdr3.1 um2_read.$i.cdr3`;
	    } else {
		`ProcessAlignment -s $species -g $gene -m $minmatch -a $adjolq -f $frac -o um1_read.$i um1_initial.$i.delta`;
                `ProcessAlignment -s $species -g $gene -m $minmatch -a $adjolq -f $frac -o um2_read.$i um2_initial.$i.delta`;
	    }
	    `CombinePEVDJDelta.pl um1_read.$i.vdjdelta um2_read.$i.vdjdelta > um_read.$i.vdjdelta`;
            `CombinePECDR3.pl um1_read.$i.cdr3 um2_read.$i.cdr3 > um_read.$i.cdr3`;
	    #`CombinePEFastq4CDR3.pl um_read.$i.vdjdelta um1_read.$i.fq um2_read.$i.fq > um_read.$i.fq`;
	    #`ExtractCDR3.pl -s $species -g $gene um_read.$i.vdjdelta um_read.$i.fq > um_read.$i.cdr3`;
	}
	exit 0;

    } else {
        print "fork: $!\n";
    }
}

foreach (@child) {
    waitpid($_, 0);
}


##### combine results and remove intermediate files

my $command;

if ($peq == 0 || $mergeq) {
    $command = "cat " . join(" ", map("read.$_.vdjdelta", (1..$thread))) . " > read.vdjdelta";
    `$command`;
    `rm read.*.vdjdelta`;

    $command = "cat " . join(" ", map("read.$_.cdr3", (1..$thread))) . " > read.cdr3";
    `$command`;
    `rm read.*.cdr3`;
}

if ($peq) {
    $command = "cat " . join(" ", map("um_read.$_.vdjdelta", (1..$thread))) . " > unmerged.vdjdelta";
    `$command`;
    `rm um*_read.*.vdjdelta`;
    
    $command = "cat " . join(" ", map("um_read.$_.cdr3", (1..$thread))) . " > unmerged.cdr3";
    `$command`;
    `rm um*_read.*.cdr3`;
    
    `mv read.vdjdelta merged.vdjdelta`;
    `mv read.cdr3 merged.cdr3`;
    if ($mergeq) {
        `cat merged.vdjdelta unmerged.vdjdelta > read.vdjdelta`;
        `cat merged.cdr3 unmerged.cdr3 > read.cdr3`;
    } else {
        `ln -s unmerged.vdjdelta read.vdjdelta`;
        `ln -s unmerged.cdr3 read.cdr3`;
    }
}

`rm *read.*.fa`;
`rm *read.*.fq`;
`rm *initial.*.delta`;
`rm *processed.*.vdjdelta*` if $patchq == 1;

$command = "CorrectCDR3Error.pl read.cdr3 > clone.txt";
`$command`;

$command = "CloneStat.pl clone.txt";
`$command`;

# finish log
$date = `date`;
print LOG "program end       : $date";

chdir("..");



################################################################################


sub Usage {
    print "usage  : trig.pl [options] read.fq/a\n";
    print "e.g.   : trig.pl -s hsa -g trb -o test read.fq\n\n";
    print "option : -species  <str>    species name                 [hsa*, mmu] (*default)\n";
    print "         -gene     <str>    immune receptor gene         [tra, trb*, trd, trg, igh, igl, igk]\n";
    print "                            can also be combined, e.g.,  [tra,trb]\n";
    print "         -outdir   <str>    output directory name        [output*]\n";
    print "         -minmatch <int>    minimal match of nucmer      [15*]\n";
    print "         -adjolq   <int>    adjust overlap Q             [0, 1*]\n";
    print "         -patchq   <int>    patch missing alignment Q    [0*, 1]\n";
    print "         -frac     <float>  alignment length fraction    [0.5]\n";
    print "         -thread   <int>    number of processors         [1*]\n";
    print "\n";
    exit 0;
}
