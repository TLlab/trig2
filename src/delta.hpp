#ifndef DELTA_H
#define DELTA_H

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <unordered_map>

#define MSC  3
#define MMSC -7
#define GSC  -7
#define EXB  2
#define SPC  1


//=======================================================

// information of a VDJ exon
struct VDJInfo_t {
	std::string species_gene;
	std::string species;
	std::string gene;
	std::string vdj;
	std::string vdj_exon;
	char strand;
	int exon_start;
	int exon_end;
	
	VDJInfo_t() {
		clear();
	}
	
	void clear() {
		species_gene.erase();
		species.erase();
		gene.erase();
		vdj.erase();
		vdj_exon.erase();
		strand = '\0';
		exon_start = exon_end = 0;
	}
};

//====================================================VDJReader_t=====

// for reading VDJ information file
class VDJReader_t
{
private:
	std::string vdj_path_m;          // vdj file path
	std::ifstream vdj_stream_m;      // vdj file input stream
	std::vector<VDJInfo_t> info_m;   // vdj info of all exons of a gene
	bool is_open_m;                  // vdj stream is open
	
	void CheckStream() {
		if(!vdj_stream_m.good()) {
			std::cerr << "\033[31mERROR:\033[0m Could not parse vdj file, "
				  << vdj_path_m << std::endl;
			exit(1);
		}
	}

public:
	VDJReader_t() {
		is_open_m = false;
	}
	~VDJReader_t() {
		vdj_stream_m.close();
		is_open_m = false;
	}

	void open(const std::string &vdj_path);
	bool readNext();

	const std::vector<VDJInfo_t> &getperVDJInfo() const {
		return info_m;
	}

	std::vector<VDJInfo_t> getallVDJInfo(const std::string &vdj_path) {
		open(vdj_path);
		std::vector<VDJInfo_t> VDJInfo;

		while(readNext()) {
			std::vector<VDJInfo_t> eachVDJ = getperVDJInfo();
			VDJInfo.insert(VDJInfo.end(), eachVDJ.begin(), eachVDJ.end());
		}
		return VDJInfo;
	}
};

//=======================================================

// information of an alignment (modified from MUMMER's source)
struct DeltaAlignment_t {
	int sR;     // start coordinate on the reference
	int eR;     // end coordinate on the reference
	int sQ;     // start coordinate on the query
	int eQ;     // end coordinate on the query
	int mmgp;   // number of mismatches and gaps in the alignment
	int simc;   // number of similarity scores < 1 in the alignment
	int stpc;   // number of stop codons in the alignment

	char ro;    // orientation along reference
        char go;    // orientation along gene
	int osQ;    // orientated start
	int oeQ;    // orientated end
	int alQ;    // query alignment len
  	int gcQ;    // count of gaps in the aligned query
	int sc;     // score
	float id;   // identity of alignment
        bool tbf;   // to be filtered

        std::string rseg;   // reference segment
        std::string rlfk;   // reference left flanking bases
        std::string rrfk;   // reference right flanking bases
        std::string qseg;   // query segment
	std::string vdj;    // annotated by VDJInfo_t
	std::string vdje;   // annotated by VDJInfo_t
	std::string ge;     // annotated by VDJInfo_t
        
	std::vector<int> deltas;            // gap positions in delta alignment
	std::vector<DeltaAlignment_t> gm;   // group member
	
	DeltaAlignment_t() {
		clear();
	}

	void clear() {
		sR = eR = sQ = eQ = 0;
		mmgp = simc = stpc = 0;
		ro = 'o';
                go = 'o';
		osQ = oeQ = alQ = gcQ = 0;
		id = sc = 0;
                tbf = false;
                rseg.erase();
                rlfk.erase();
                rrfk.erase();
                qseg.erase();
		vdj.erase();
		deltas.clear();
		gm.clear();
	}

	int overlapQ(const DeltaAlignment_t &);

	std::string concise_form();

        std::vector<std::string> getEndAlignment(const int &, bool);
        void cutEndAlignment(const int &);
};

//=======================================================

// all alignments and info of a query
struct DeltaRecord_t {
	std::string idR;   // reference contig ID
	std::string idQ;   // query contig ID
	int lenR;          // length of the reference contig
	int lenQ;          // length of the query contig

	std::vector<DeltaAlignment_t> aligns;   // alignments of query segments to reference

	DeltaRecord_t() {
		clear();
	}
	
	void clear() {
		idR.erase();
		idQ.erase();
		lenR = lenQ = 0;
		aligns.clear();
	}
};

//====================================================DeltaReader_t===

// for reading a delta file (incorporated from MUMMER)
class DeltaReader_t
{
private:
	std::string delta_path_m;     // delta input file
	std::ifstream delta_stream_m; // delta file input stream
	std::string reference_path_m; // reference file
	std::string query_path_m;     // query file
	std::string data_type_m;      // type of data
	DeltaRecord_t record_m;       // current delta information record
	bool is_record_m;             // valid record
	bool is_open_m;               // delta stream is open

	bool readNextRecord (const bool read_deltas);
	void readNextAlignment (DeltaAlignment_t & align, const bool read_deltas);

	void CheckStream() {
		if(!delta_stream_m.good()) {
			std::cerr << "\033[31mERROR:\033[0m Could not parse delta file, "
				<< delta_path_m << std::endl;
			exit(1);
		}
	}

public:
	DeltaReader_t() {
		is_open_m = false;
		is_record_m = false;
	}
	~DeltaReader_t() {
		close();
	}

	void open(const std::string &delta_path);
	
	void close() {
		delta_path_m.erase();
		delta_stream_m.close();
		reference_path_m.erase();
		query_path_m.erase();
		data_type_m.erase();
		record_m.clear();
		is_open_m = false;
		is_record_m = false;
	}

	bool readNext (bool getdeltas = true) {
		return readNextRecord (getdeltas);
	}

	const std::string &getReferencePath() const {
		return reference_path_m;
	}

	const std::string &getQueryPath() const {
		return query_path_m;
	}

	DeltaRecord_t &getRecord() {
		return record_m;
	}
};

//====================================================DeltaFilter_t===

// for processing alignments of a query, i.e., kernel of TRIg
class DeltaFilter_t
{
private:
	DeltaRecord_t rec_m;

	int reg_m;                  // regularity of a query
        char ori_m;                 // orientation of a query
	std::string CombineVDJ_m;   // combined VDJ annotation
	std::vector<int> lndisi;    // index of longest non-decreasing sub-alignments

	int vi, di, ji;   // index of V, D, and J gene in all alignments of a query

	void LNDIS();     // get longest non-decreasing or non-increasing sub-alignments

	// find maximal score position in the overlapping alignments
        int maxScorePosition(std::vector<std::string> &adjseq1, std::vector<std::string> &adjseq2,
                              std::string &rfk1, std::string &rfk2);

        // update vdj_index (whenever sorting rec_m.aligns) (to be discarded)
        void update_VDJ_index();

public:
        int al_m;                                                      // alignment length
	float alf_m;                                                   // aligned length fraction
        std::string rc_m;                                              // recombination code
	static std::unordered_map<std::string, std::string> refseq_m;  // load reference sequence
	std::string qryseq_m;                                          // load query sequence by one
	static std::vector<VDJInfo_t> VDJInfo_m;                       // load vdj information

	DeltaFilter_t(DeltaRecord_t &rec) {
                clear();
		rec_m = rec;
	}
	~DeltaFilter_t() {
		clear();
	}

	void clear() {
		rec_m.clear();
		reg_m = -1;
		ori_m = 'o';
                al_m = 0;
		alf_m = 0;
                rc_m.erase();
		vi = di = ji = 0;
		CombineVDJ_m.erase();
		lndisi = {};
	}

	// TRIg functions
	void getOptimalSet();
	void annotateVDJ();
	void groupAlignment();
	void filterAlignment();
	void adjustOverlap();
        void setRecombCode();
	void annotateQuery();
	void printResult();
	void printResult(std::ofstream &cout);

	const DeltaRecord_t &getREC() const {
		return rec_m;
	}
	const int &getREG() const {
		return reg_m;
	}
	const char &getORI() const {
		return ori_m;
	}
	const std::string &getVDJ() const {
		return CombineVDJ_m;
	}
	const int &getVi() const {
		return vi;
	}
	const int &getDi() const {
		return di;
	}
	const int &getJi() const {
		return ji;
	}
};

#endif /* delta.h */
