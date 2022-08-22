#include <iostream>
#include <fstream>
#include <getopt.h>
#include <vector>
#include <unordered_map>
#include "delta.hpp"
#include "fastx_read.hpp"
#include "extractCDR3.hpp"

using namespace std;

//====================================================Options===
string    OPT_Delta_file;
string    OPT_Species    = "hsa";
string    OPT_Gene       = "trb";
int       OPT_Minmatch   = 15;
int       OPT_Adjolq     = 1;
float     OPT_Frac       = 0.5;
string    OPT_Output     = "read";

//===================================================Function===
void ParseArgs(int argc, char ** argv);
void help();

// initialization static object in DeltaFilter_t
unordered_map<string, string> DeltaFilter_t::refseq_m;
vector<VDJInfo_t> DeltaFilter_t::VDJInfo_m;

// initialization static object in ExtractCDR3_t
map<std::string, int> ExtractCDR3_t::cdr3p_m;

//=======================================================Main===
int main(int argc, char **argv) {

	// Command line parsing
	ParseArgs(argc, argv);
	const string vdjpath = OPT_Species + "_"  + OPT_Gene + ".vdj";
	const string deltapath = OPT_Delta_file;
	const string cdr3path = OPT_Species + "_"  + OPT_Gene + ".cdr";
	const string vdjdeltaout = OPT_Output + ".vdjdelta";
	const string cdr3out = OPT_Output + ".cdr3";
	
	// output file stream
	ofstream OUT_V;
	ofstream OUT_C;
	OUT_V.open(vdjdeltaout);
	OUT_C.open(cdr3out);

	// load MUMmer delta file
	DeltaReader_t dr;
	dr.open(deltapath);
	
	// delta header
	const string refpath = dr.getReferencePath();
	const string qrypath_fa = dr.getQueryPath();
	const string qrypath_fq = qrypath_fa.substr(0, qrypath_fa.length()-2) + "fq";

	// load reference, query sequence and vdj info to static object in deltafilter_t;
	DeltaFilter_t::refseq_m = FASTA_t::getfasta(refpath);
	VDJReader_t vr;
	DeltaFilter_t::VDJInfo_m = vr.getallVDJInfo(vdjpath);

	// load cdr3 info
	CDRReader_t cr;
	ExtractCDR3_t::cdr3p_m = cr.getCDR3(cdr3path);

	// load fastq file using FastqReader_t
	FastqReader_t fr;
	fr.open(qrypath_fq);

	// with gap information
	while(dr.readNext(true)) {

		DeltaFilter_t df(dr.getRecord());

		// fasta without delta information
		while(fr.readNext()) {
			const string uid = fr.getUID();
			const string seq = fr.getSEQ();
			if (dr.getRecord().idQ == uid) {
				df.qryseq_m = seq;
				break;
			} else {
				OUT_V << uid << "\t" << seq.length() << "\t" << "---" << endl;
				OUT_C << uid << "\t" << "---" << "\t" << "---" << endl;
			}
		}

		df.getOptimalSet();
                df.printResult();
		df.annotateVDJ();
		df.groupAlignment();
		df.filterAlignment();
                df.setRecombCode();
		if (OPT_Adjolq && df.rc_m != "CH")
			df.adjustOverlap();
		df.annotateQuery();

                //df.printResult(OUT_V);
                //if (df.alf_m >= OPT_Frac) {
                if (df.al_m >= 30) {
			df.printResult(OUT_V);
		} else {
			OUT_V << dr.getRecord().idQ << "\t" << dr.getRecord().lenQ << "\t" << "---" << endl;
		}
                
                // CDR3
		if (df.getREG() == 1 || df.getREG() == 2) {
			ExtractCDR3_t excdr(df.getREC(), df.getREG(), df.getORI(), df.getVDJ(), df.getVi(), df.getJi());
                        excdr.inputFastq(fr.getSEQ(), fr.getQUA());
			excdr.extractCDR3();
			excdr.printResult(OUT_C);
		} else {
			OUT_C << fr.getUID() << "\t" << 0 << "\t" << "---" << endl;
		}
	}

	// else fasta without delta information
        while(fr.readNext()) {
		OUT_V << fr.getUID() << "\t" << fr.getSEQ().length() << "\t" << "---" << endl;
		OUT_C << fr.getUID() << "\t---\t---" << endl;
	}

	OUT_V.close();
	OUT_C.close();
	return 0;
}

//==================================================ParseArgs===
void ParseArgs(int argc, char ** argv) {
	int opt, errflg = 0;
	const char *optstring = "s:g:m:a:f:o:";
	const struct option int_opts[] = {
		{"species",  1, NULL, 's'},
		{"gene",     1, NULL, 'g'},
		{"minmatch", 1, NULL, 'm'},
		{"adjolq",   1, NULL, 'a'},
		{"frac",     1, NULL, 'f'},
		{"output",   1, NULL, 'o'},
	};

	while((opt = getopt_long(argc, argv, optstring, int_opts, NULL)) != -1) {
		switch(opt) {
			case (int)'s':
				OPT_Species = optarg;
				break;
			case (int)'g':
				OPT_Gene = optarg;
				break;
			case (int)'m':
				OPT_Minmatch = atoi(optarg);
				break;
			case (int)'a':
				OPT_Adjolq = atoi(optarg);
				break;
			case (int)'f':
				OPT_Frac = atof(optarg);
			case (int)'o':
				OPT_Output = optarg;
				break;
			default:
				errflg++;
		}
	}

	if (errflg > 0 || optind != argc -1) help();
	OPT_Delta_file = argv[optind++];
}

//=======================================================Help===
void help() {
	cout << "usage  : ProcAlgn [option] initial.delta\n\n" <<
		"option : -s | --species  species name              [hsa*, mmu] (*default)\n" <<
		"         -g | --gene     immune receptor gene      [tra, trb*, trd, trg, igh, igl, igk]\n" <<
		"         -m | --minmatch minimal match of nucmer   [15*]\n" <<
		"         -a | --adjolq   adjust overlap Q          [0, 1*]\n" <<
		"         -f | --frac     alignment length fraction [0.5*]\n" <<
		"         -o | --output   output filenames prefix   [read*] (ext: .vdjdelta .cdr3)\n\n";
	exit(0);
}
