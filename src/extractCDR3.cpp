#include "delta.hpp"
#include "fastx_read.hpp"
#include "extractCDR3.hpp"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <unordered_map>
#include <algorithm>

//================================Translate
const std::unordered_map<std::string, std::string> aacode = {                                                                                                                   
    {"TTT", "F"}, {"TTC", "F"}, {"TTA", "L"}, {"TTG", "L"},
    {"TCT", "S"}, {"TCC", "S"}, {"TCA", "S"}, {"TCG", "S"},
    {"TAT", "Y"}, {"TAC", "Y"}, {"TAA", "*"}, {"TAG", "*"},
    {"TGT", "C"}, {"TGC", "C"}, {"TGA", "*"}, {"TGG", "W"},
    {"CTT", "L"}, {"CTC", "L"}, {"CTA", "L"}, {"CTG", "L"},
    {"CCT", "P"}, {"CCC", "P"}, {"CCA", "P"}, {"CCG", "P"},
    {"CAT", "H"}, {"CAC", "H"}, {"CAA", "Q"}, {"CAG", "Q"},
    {"CGT", "R"}, {"CGC", "R"}, {"CGA", "R"}, {"CGG", "R"},
    {"ATT", "I"}, {"ATC", "I"}, {"ATA", "I"}, {"ATG", "M"},
    {"ACT", "T"}, {"ACC", "T"}, {"ACA", "T"}, {"ACG", "T"},
    {"AAT", "N"}, {"AAC", "N"}, {"AAA", "K"}, {"AAG", "K"},
    {"AGT", "S"}, {"AGC", "S"}, {"AGA", "R"}, {"AGG", "R"},
    {"GTT", "V"}, {"GTC", "V"}, {"GTA", "V"}, {"GTG", "V"},
    {"GCT", "A"}, {"GCC", "A"}, {"GCA", "A"}, {"GCG", "A"},
    {"GAT", "D"}, {"GAC", "D"}, {"GAA", "E"}, {"GAG", "E"},
    {"GGT", "G"}, {"GGC", "G"}, {"GGA", "G"}, {"GGG", "G"},
};

std::string Translate(const std::string seq, int str = 0) {
    int i;
    std::string aa; 
    for (i = str; i < seq.length()/3*3; i += 3) {
            if (aacode.find(seq.substr(i, 3)) == aacode.end()) {
                    aa += "*";
            } else {
                    aa += aacode.at(seq.substr(i, 3));
            }
    }

    // not a multiple of three
    if (i < seq.size())
        aa += "_";

    return aa; 
}

//================================CDRReader_t

std::map<std::string, int> CDRReader_t::getCDR3(const std::string &cdr_path) {
	open(cdr_path);
	while(readNext()) {
		allcdr3p_m[getVDJ()] = getCDR3p();
	}
	return allcdr3p_m; 
}

void CDRReader_t::open(const std::string &cdr_path) {
	cdr_path_m = cdr_path;
	cdr_stream_m.open(cdr_path_m);
	CheckStream();
	is_open_m = true;
}

bool CDRReader_t::readNext(){
	if (cdr_stream_m.peek() == EOF) {
		return false;
	}

	// read cdr3 information (skip psudo gene)
	std::string speces_gene, vdj, n, cdr1, cdr2, cdr3;
	do {
		cdr_stream_m >> speces_gene;   // e.g. hsa_trb
		cdr_stream_m >> vdj;           // e.g. TRBV1
		cdr_stream_m >> n;             // 0 1 2 ---
		cdr_stream_m >> cdr1;
		cdr_stream_m >> cdr2;
		cdr_stream_m >> cdr3;

		// flush string after cdr3
		while(cdr_stream_m.get() != '\n');
	} while (cdr3 == "---");

	vdj_m = vdj;
	cdr3p_m = stoi(cdr3);

	return true;
}

void ExtractCDR3_t::extractCDR3() {
	
	// for all vj pairs
	for (auto v = valn_m.gm.begin(); v != valn_m.gm.end(); v++) {
		for (auto j = jaln_m.gm.begin(); j != jaln_m.gm.end(); j++) {

			cdr3seq_m = "---";
			cdr3qua_m = "---";
			cdr3aa_m = "---";

			// check if the CDR3 positions on the reference are available (i.e., not pseudogene)
			// and the positions are covered by the V and J alignments (take care of V30 on the minus strand)
			bool cdr3q = false;
			if (cdr3p_m.find(v->vdj) != cdr3p_m.end() && cdr3p_m.find(j->vdj) != cdr3p_m.end()) {
				if (v->sR <= cdr3p_m.at(v->vdj) && cdr3p_m.at(v->vdj) <= v->eR &&
						j->sR <= cdr3p_m.at(j->vdj) && cdr3p_m.at(j->vdj) <= j->eR) {
					cdr3q = true;
				}
			}

			// if good, first get the CDR3 starting and ending positions on the query
			if (cdr3q) {
				int vqp = AlignmentRpQp(*v);
				int jqp = AlignmentRpQp(*j);

				// then get the CDR3 segment on the plus strand
				if (j->sQ < j->eQ) {
					cdr3seq_m = fqseq_m.substr(vqp -1, jqp - vqp +1);
					cdr3qua_m = fqqua_m.substr(vqp -1, jqp - vqp +1);
				} else {
					cdr3seq_m = fqseq_m.substr(jqp -1, vqp - jqp +1);
					cdr3qua_m = fqqua_m.substr(jqp -1, vqp - jqp +1);
					cdr3seq_m = FASTA_t::revcom(cdr3seq_m);
					std::reverse(cdr3qua_m.begin(), cdr3qua_m.end());
				}

				// translate cdr3
				cdr3aa_m = Translate(cdr3seq_m, 0);
			}

			// push to array
			cdr3c_m.push_back(v->vdj + ":" + cdr3seq_m + ":" + j->vdj);
			cdr3q_m.push_back(cdr3qua_m);
			cdr3a_m.push_back(cdr3aa_m);

		}
	}
}

int ExtractCDR3_t::AlignmentRpQp(DeltaAlignment_t align) {

	// ro = align.ro (+/-)
	int o = align.ro == '+' ? 1 : -1;
	// gap = align.deltas = delta (vector<int>)
	align.deltas.pop_back();

	// bs : block size (idea from blat)
	// rg : reference gap
	std::vector<int> bs;
	std::vector<int> rg;

	for (const int &d : align.deltas) {
		if (d > 0) {
			bs.push_back(d-1);
			rg.push_back(1);
		} else {
			bs.push_back(-d-1);
			rg.push_back(0);
		}
	}
	bs.push_back(align.eR - align.sR + 1 
			- (std::accumulate(bs.begin(), bs.end(), 0) + std::accumulate(rg.begin(), rg.end(), 0)));

	// calculate query position
	int rl = cdr3p_m.at(align.vdj) - align.sR;
	int ri = -1;
	int qi = -1;
	int bi = 0;
	while(ri < rl) {
		int d = rl - ri;
		if (d <= bs[bi]) {
			ri += d;
			qi += d;
			break;
		} else {
			ri += bs[bi];
			qi += bs[bi];
			if (rg[bi] == 1) {
				ri++;
			} else {
				qi++;
			}
		}
		bi++;
	}
	return align.sQ + qi * o;
}

void ExtractCDR3_t::printResult() {
	std::cout << idQ_m << "\t" << reg_m << "\t";
	std::cout << cdr3c_m[0];
	for(auto c = cdr3c_m.begin()+1; c != cdr3c_m.end(); c++)
		std::cout << "|" << *c;
	std::cout << "\t";
	std::cout << cdr3q_m[0];
	for(auto q = cdr3q_m.begin()+1; q != cdr3q_m.end(); q++)
		std::cout << "|" << *q;
	std::cout << "\t";
	std::cout << cdr3a_m[0];
	for(auto a = cdr3a_m.begin()+1; a != cdr3a_m.end(); a++)
		std::cout << "|" << *a;
	std::cout << std::endl;
}

void ExtractCDR3_t::printResult(std::ofstream &cout) {
	cout << idQ_m << "\t" << reg_m << "\t";
	cout << cdr3c_m[0];
	
	for(auto c = cdr3c_m.begin()+1; c != cdr3c_m.end(); c++)
		cout << "|" << *c;
	cout << "\t";
	cout << cdr3q_m[0];
	for(auto q = cdr3q_m.begin()+1; q != cdr3q_m.end(); q++)
		cout << "|" << *q;
	cout << "\t";
	cout << cdr3a_m[0];
	for(auto a = cdr3a_m.begin()+1; a != cdr3a_m.end(); a++)
		cout << "|" << *a;
	
	cout << std::endl;
}
