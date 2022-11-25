#ifndef EXTRACTCDR3_H
#define EXTRACTCDR3_H

#include "delta.hpp"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <map>
#include <vector>
#include <numeric>

//============================================CDRReader_t

class CDRReader_t
{
private:
	std::string cdr_path_m;
	std::ifstream cdr_stream_m;
	std::string vdj_m;
	int cdr3p_m;
	std::map<std::string, int> allcdr3p_m;

	bool is_open_m;
	void CheckStream() {
		if(!cdr_stream_m.good()) {
			std::cerr << "\033[31mERROR:\033[0m Could not parse cdr file, "
				<< cdr_path_m << std::endl;
			exit(1);
		}
	}

public:
	CDRReader_t() {
		is_open_m = false;
	}
	~CDRReader_t() {
		cdr_stream_m.close();
		is_open_m = false;
	}
	void open(const std::string &cdr_path);
	bool readNext();
	std::map<std::string, int> getCDR3(const std::string &cdr_path);

	const std::string &getVDJ() const {
		return vdj_m;
	}
	const int &getCDR3p() const {
		return cdr3p_m;
	}
};

//=============================================ExtractCDR3_t

class ExtractCDR3_t
{
private:
	std::string idQ_m;
	int reg_m;
	char ori_m;
	std::string v_m, d_m, j_m;

	DeltaAlignment_t valn_m;
	DeltaAlignment_t jaln_m;

	std::string fqseq_m;
	std::string fqqua_m;

	std::string cdr3seq_m;
	std::string cdr3qua_m;
	std::string cdr3aa_m;

	std::vector<std::string> cdr3c_m;
	std::vector<std::string> cdr3q_m;
	std::vector<std::string> cdr3a_m;

	int AlignmentRpQp(DeltaAlignment_t align);

public:
	static std::map<std::string, int> cdr3p_m;
	ExtractCDR3_t(const DeltaRecord_t &rec, const int &reg, const char &ori, const std::string &vdj, const int &vi, const int &ji) {
		idQ_m = rec.idQ;
		reg_m = reg;
		ori_m = ori;

		std::stringstream ss(vdj);
		getline(ss, v_m, ':');
		getline(ss, d_m, ':');
		getline(ss, j_m, ':');

		valn_m = rec.aligns[vi];
		valn_m.gm.push_back(valn_m);

		jaln_m = rec.aligns[ji];
		jaln_m.gm.push_back(jaln_m);
		cdr3seq_m = cdr3qua_m = cdr3aa_m = "---";
	}
	~ExtractCDR3_t() {
		clear();	
	}

	void clear() {
		idQ_m.erase();
		reg_m = 0;
		ori_m = 0;
		v_m.erase(); d_m.erase(); j_m.erase();
		valn_m.clear();
		jaln_m.clear();

		fqseq_m.erase();
		fqqua_m.erase();

		cdr3seq_m.erase();
		cdr3qua_m.erase();
		cdr3aa_m.erase();

		cdr3c_m.clear();
		cdr3q_m.clear();
		cdr3a_m.clear();
	}

	void inputFastq(std::string seq, std::string qua) {
		fqseq_m = seq;
		fqqua_m = qua;
	}
	void extractCDR3();
	void printResult(std::ofstream &cout);
	friend std::ostream& operator<< (std::ostream &out, const ExtractCDR3_t &et);
};

#endif /* extractCDR3.h */
