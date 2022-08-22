#include "delta.hpp"
#include "fastx_read.hpp"

#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <cctype>
#include <vector>
#include <map>
#include <unordered_map>

// IGH constant region MEGAD issue to be solved
// IGK distal orientation issue to be solved
std::map<char, int> GGO = { {'V', 0}, {'J', 1}, {'D', 2}, {'C', 3}, {'I', 4} };
std::map<std::string, int> GEO = { {"V0", 0}, {"V1", 1}, {"V2", 2}, {"D0", 3}, {"J0", 4},
				   {"C1", 5}, {"C2", 6}, {"C3", 7}, {"C4", 8}};

//====================================================VDJReader_t=====

void VDJReader_t::open(const std::string &vdj_path) {
	vdj_path_m = vdj_path;
	vdj_stream_m.open(vdj_path_m);
	CheckStream();
	is_open_m = true;
}

bool VDJReader_t::readNext() {
	if (vdj_stream_m.peek() == EOF)
		return false;
	
	// make way for the new vdj info
	info_m.clear();

	// read the vdj information
	VDJInfo_t info;
	std::string exon_range;
	std::vector<int> each_sted;   // each range of start and end

	vdj_stream_m >> info.species_gene;   // e.g., hsa_trb, mmu_igh
	vdj_stream_m >> info.vdj;            // e.g., V3-1, D1, J2-3
	vdj_stream_m.ignore(2);              // i.e., F, P
	vdj_stream_m >> info.strand;         // i.e., +, -
	vdj_stream_m >> exon_range;          // e.g., 1..41,159..450

	info.species = info.species_gene.substr(0, 3);
	info.gene = info.species_gene.substr(4, 3);
	std::transform(info.gene.begin(), info.gene.end(), info.gene.begin(), ::toupper);
	
	// flush the remaining whitespace
	while(vdj_stream_m.get() != '\n');

	// split string by [.,]
	int num_of_exon = 1;
	std::string temp;
	for(int i = 0; i < (int)exon_range.size(); i++) {
		if(exon_range[i] == '.') {
			each_sted.push_back(stoi(temp));
			temp = "";
			i++;   // skip next dot
		} else if (exon_range[i] == ',') {
			each_sted.push_back(stoi(temp));
			temp = "";
			num_of_exon++;
		} else {
			temp += exon_range[i];
		}
	}
	each_sted.push_back(stoi(temp));

	// push each vdj info to vector<VDJInfo_t>
	if (num_of_exon == 1) {
		info.vdj_exon = info.vdj + "_0";
		info.exon_start = each_sted[0];
		info.exon_end = each_sted[1];
		info_m.push_back(info);
	} else {
		if (info.strand == '+') {
			for (int i = 1; i <= num_of_exon; i++) {
				info.vdj_exon = info.vdj + "_" + std::to_string(i);
				info.exon_start = each_sted[i*2-2];
				info.exon_end = each_sted[i*2-1];
				info_m.push_back(info);
			}
		} else { // strand == '-'
			for (int i = 1; i <= num_of_exon; i++) {
				info.vdj_exon = info.vdj + "_" + std::to_string(num_of_exon-i+1);
				info.exon_start = each_sted[i*2-2];
				info.exon_end = each_sted[i*2-1];
				info_m.push_back(info);
			}
		}
	}
	return true;
}

//====================================================DeltaAlignment_t===

int DeltaAlignment_t::overlapQ(const DeltaAlignment_t &b) {
	int ol = this->osQ < b.osQ ? this->oeQ - b.osQ + 1 : b.oeQ - this->osQ + 1;
	int q = ol > b.alQ * 0.5 ? 1 : 0;
	return q;
}

std::string DeltaAlignment_t::concise_form() {
	std::string cf = "hsa_trb:" + this->vdje + ":";
	cf += std::to_string(this->sR) + "-" + std::to_string(this->eR) + ":";
	cf += std::to_string(this->sQ) + "-" + std::to_string(this->eQ) + ":";
	cf += std::to_string(this->mmgp) + ":";
	cf += std::to_string(this->deltas[0]);
	for (auto i = this->deltas.begin()+1; i != this->deltas.end(); i++) {
		cf += "," + std::to_string(*i);
	}
	return cf;
}

std::vector<std::string> DeltaAlignment_t::getEndAlignment(const int &n, bool rev) {
        std::vector<std::string> aseg;

        // starting loci
        int ri = 0;
        int qi = 0;
        int gi = 0;
        int gs = 0;
        if (n < 0) {
                int qs = this->alQ + n;
                while (qi < qs) {
                        int g = this->deltas[gi];
                        int d = qs - qi;
                        if (g > 0) {
                                if (d < g-1) {
                                        ri += d;
                                        qi += d;
                                        gs = d;
                                } else {
                                        ri += g;
                                        qi += g-1;
                                        gi += 1;
                                }
                        } else if (g == 0) {
                                ri += d;
                                qi += d;
                                gs = d;
                        } else {
                                g = -g;
                                if (d < g) {
                                        ri += d;
                                        qi += d;
                                        gs = d;
                                } else {
                                        ri += g-1;
                                        qi += g;
                                        gi += 1;
                                }
                        }
                }
        }

        // get aligned segment on reference and query
        std::string raseg;
        std::string qaseg;
        
        // if no gap
        if (this->deltas[gi] == 0) {
                raseg = this->rseg.substr(ri, n);
                qaseg = this->qseg.substr(qi, n);
                if (rev == true) {
                        raseg = FASTA_t::revcom(raseg);
                        qaseg = FASTA_t::revcom(qaseg);
                }
                aseg.push_back(raseg);
                aseg.push_back(qaseg);
                return aseg;
        }

        // if there is a gap
        int d = abs(n);   // distance to desired length
        while (d > 0) {
                int g = this->deltas[gi];
                if (g > 0) {
                        g -= gs;
                        if (d < g) {
                                raseg += this->rseg.substr(ri, d);
                                qaseg += this->qseg.substr(qi, d);
                                d = 0;
                        } else {
                                raseg += this->rseg.substr(ri, g);
                                qaseg += this->qseg.substr(qi, g-1) + '_';
                                d -= g-1;
                                ri += g;
                                qi += g-1;
                                gi += 1;
                                gs = 0;
                        }
                } else if (g == 0) {
                        raseg += this->rseg.substr(ri, d);
                        qaseg += this->qseg.substr(qi, d);
                        d = 0;
                } else {
                        g = -g - gs;
                        if (d < g) {
                                raseg += this->rseg.substr(ri, d);
                                qaseg += this->qseg.substr(qi, d);
                                d = 0;
                        } else {
                                raseg += this->rseg.substr(ri, g-1) + '_';
                                qaseg += this->qseg.substr(qi, g);
                                d -= g;
                                ri += g-1;
                                qi += g;
                                gi += 1;
                                gs = 0;
                        }
                }
        }

        if (rev == true) {
                raseg = FASTA_t::revcom(raseg);
                qaseg = FASTA_t::revcom(qaseg);
        }
        aseg.push_back(raseg);
        aseg.push_back(qaseg);
        return aseg;
}

void DeltaAlignment_t::cutEndAlignment(const int &n) {
        if (n == 0)
                return;
        
        // get end alignment, count mismatch and gap, and update
        int cutoffmmgp = 0;
        std::vector<std::string> cutoffseq = getEndAlignment(n, false);
        for (size_t i = 0; i < cutoffseq[0].size(); i++) {
                if ( (char)std::toupper(cutoffseq[0][i]) != (char)cutoffseq[1][i] ) cutoffmmgp++;
        }
        this->mmgp -= cutoffmmgp;

        // move to the desired cut point
        int d = n > 0 ? n : this->alQ + n;   // distance to the desired cut point
        int ri = 0;
        int qi = 0;
        int gi = 0;
        int gs = 0;
        std::vector<int> newdeltas;

        while (d > 0) {
                int g = this->deltas[gi];
                if (g > 0) {
                        if (d < g-1) {
                                ri += d;
                                qi += d;
                                gs = d;
                                d = 0;
                        } else {
                                ri += g;
                                qi += g-1;
                                gi += 1;
                                newdeltas.push_back(g);
                                d -= g-1;
                        }
                } else if (g == 0) {
                        ri += d;
                        qi += d;
                        d = 0;
                } else {
                        g = -g;
                        if (d < g) {
                                ri += d;
                                qi += d;
                                gs = d;
                                d = 0;
                        } else {
                                ri += g-1;
                                qi += g;
                                gi += 1;
                                newdeltas.push_back(-g);
                                d -= g;
                        }
                }
        }

        // cut from head
        if (n > 0) {
                this->sR   += ri;
                this->alQ  -= n;
                this->rseg  = this->rseg.substr(ri);
                this->qseg  = this->qseg.substr(qi);
                if (this->deltas[gi] == 0) {
                        newdeltas = {0};
                } else {
                        newdeltas = {};
                        int g = this->deltas[gi];
                        g = g > 0 ? g - gs : g + gs;
                        newdeltas.push_back(g);
                        for (int i = gi+1; i < this->deltas.size(); i++) {
                                newdeltas.push_back(this->deltas[i]);
                        }
                }
                this->deltas = newdeltas;
                if (this->ro == '+') {
                        this->sQ  += qi;
                        this->osQ += qi;
                } else {
                        this->sQ  -= qi;
                        this->oeQ -= qi;
                }
                
        // cut from tail
        } else {
                int rl = (this->eR - this->sR + 1) - ri;
                int ql = this->alQ - qi;
                this->eR   -= rl;
                this->alQ  += n;
                this->rseg  = this->rseg.substr(0, ri);
                this->qseg  = this->qseg.substr(0, qi);
                newdeltas.push_back(0);
                this->deltas = newdeltas;
                if (this->ro == '+') {
                        this->eQ  -= ql;
                        this->oeQ -= ql;
                } else {
                        this->eQ  += ql;
                        this->osQ += ql;
                }
        }

        // Note that simc, stpc, gcQ, sc, and id are not modified.
        // But these will not be used afterward.
}

//====================================================DeltaReader_t===

void DeltaReader_t::open(const std::string &delta_path) {
	delta_path_m = delta_path;
	
	// open delta file
	delta_stream_m.open(delta_path_m);
	CheckStream();

	// read file header
	delta_stream_m >> reference_path_m;
	delta_stream_m >> query_path_m;
	delta_stream_m >> data_type_m;
	is_open_m = true;

	while(delta_stream_m.peek() != '>')
		if (delta_stream_m.get() == EOF)
			break;
}

bool DeltaReader_t::readNextRecord(const bool read_deltas) {
	
	// EOF or or any other abnormality
	if(delta_stream_m.peek() != '>')
		return false;

	// make way for the new record
	record_m.clear();
	is_record_m = true;

	// read the record header
	delta_stream_m.get();             // remove '>'
	delta_stream_m >> record_m.idR;
	delta_stream_m >> record_m.idQ;
	delta_stream_m >> record_m.lenR;
	delta_stream_m >> record_m.lenQ;

	// flush the remaining whitespace
	while(delta_stream_m.get () != '\n');

	// for each alignment...
	DeltaAlignment_t align;
	while(delta_stream_m.peek () != '>' && delta_stream_m.peek () != EOF) {
		readNextAlignment(align, read_deltas);
		record_m.aligns.push_back(align);
	}

	return true;
}

void DeltaReader_t::readNextAlignment(DeltaAlignment_t &align, const bool read_deltas) {
	int delta;      // indel pos
	int gapT = 0;   // total gaps
	int gapQ = 0;   // query gaps

	// make way for the new alignment
	align.clear();

	// read the alignment header
	delta_stream_m >> align.sR;
	delta_stream_m >> align.eR;
	delta_stream_m >> align.sQ;
	delta_stream_m >> align.eQ;
	delta_stream_m >> align.mmgp;
	delta_stream_m >> align.simc;
	delta_stream_m >> align.stpc;

	align.ro = align.sQ < align.eQ ? '+' : '-';
        align.go = align.ro;
        align.osQ = align.sQ < align.eQ ? align.sQ : align.eQ;
	align.oeQ = align.sQ < align.eQ ? align.eQ : align.sQ;
	align.alQ = align.oeQ - align.osQ + 1;

	// get gap info
	do {
		delta_stream_m >> delta;
		if (delta > 0) {
			gapQ++;
			gapT++;
		}
		if (delta < 0)
			gapT++;
		if (read_deltas)
			align.deltas.push_back(delta);
	} while(delta != 0);

	int mismatch = align.mmgp - gapT;
	int match    = align.alQ + gapQ - align.mmgp;
	align.id     = match / (float)(align.alQ + gapQ);
	align.sc     = match * MSC + mismatch * MMSC + gapT * GSC;

	// flush the remaining whitespace
	while(delta_stream_m.get () != '\n');
}

//====================================================DeltaFilter_t===

// comparing function for sorting alignments by score and identity
struct aligns_SC_Cmp_t {
	bool operator() (const DeltaAlignment_t &i, const DeltaAlignment_t &j) const {
		if (i.sc > j.sc)
			return true;
		else if (i.sc < j.sc)
			return false;
		else if (i.id > j.id)
			return true;
		else
			return false;
	}
};

void DeltaFilter_t::getOptimalSet() {
	std::sort(rec_m.aligns.begin(), rec_m.aligns.end(), aligns_SC_Cmp_t());
	std::vector<DeltaAlignment_t> oaln;

	while(!rec_m.aligns.empty()) {

		// get current best alignment
		std::vector<DeltaAlignment_t> cba = {rec_m.aligns[0]};
		rec_m.aligns.erase(rec_m.aligns.begin());

		while(!rec_m.aligns.empty() && rec_m.aligns[0].sc == cba.back().sc && rec_m.aligns[0].id == cba.back().id) {
			cba.push_back(rec_m.aligns[0]);
			rec_m.aligns.erase(rec_m.aligns.begin());
		}
		oaln.insert(oaln.end(), cba.begin(), cba.end());

		// filter alignment that overlaps the current best alignment
		for (auto a : cba) {
			std::vector<DeltaAlignment_t> aln_tmp;
			for (auto b : rec_m.aligns) {
				if (a.overlapQ(b) == 0)
					aln_tmp.push_back(b); 
			}
			rec_m.aligns.assign(aln_tmp.begin(), aln_tmp.end());
		}
	}
	rec_m.aligns.assign(oaln.begin(), oaln.end());
}

void DeltaFilter_t::annotateVDJ() {
	for(auto i = rec_m.aligns.begin(); i != rec_m.aligns.end(); i++) {

		const int rs = i->sR;          // alignment range start
		const int re = i->eR;          // alignment range end
		int si = 0;                    // start of info
		int ei = VDJInfo_m.size()-1;   // end of info
		std::string annot;
		std::string annot_exon;

		// range start is after start of the last vdj
		if (rs > VDJInfo_m[ei].exon_start) {
			i->vdj = rs > VDJInfo_m[ei].exon_end ? VDJInfo_m[ei].gene + "I" : VDJInfo_m[ei].vdj;
			i->vdje = rs > VDJInfo_m[ei].exon_end ? VDJInfo_m[ei].gene + "I_0" : VDJInfo_m[ei].vdj_exon;
			i->ge = "I0";
			continue;
		}

		// locate the range start
		int mid = (si+ei) / 2;
		while (mid != si) {
			rs < VDJInfo_m[mid].exon_start ? ei = mid : si = mid;
			mid = (si+ei) / 2;
		}

		// get vdj exon from the starting location
		if (rs <= VDJInfo_m[si].exon_end) {
			annot = VDJInfo_m[si].vdj;
			annot_exon = VDJInfo_m[si].vdj_exon;
		}
		while(++si) {
			if (si > VDJInfo_m.size()-1 || re < VDJInfo_m[si].exon_start)
				break;
			if (!annot_exon.empty()) {
				annot += "~";
				annot_exon += "~";
			}
			annot += VDJInfo_m[si].vdj;
			annot_exon += VDJInfo_m[si].vdj_exon;
		}
		if (annot_exon.empty()) {
			i->vdj = VDJInfo_m[si].gene + "I";
			i->vdje = VDJInfo_m[si].gene + "I_0";
			i->ge = "I0";
			continue;
		}
		i->vdj = annot;
		i->vdje = annot_exon;
		i->ge = std::string(1, annot[3]) + annot_exon[annot_exon.length()-1];
                if (annot == "TRBV30") {
                        i->go = i->ro == '+' ? '-' : '+';
                }
	}
}

void DeltaFilter_t::groupAlignment() {
	std::vector<DeltaAlignment_t> galn;

        // sort along the query
        std::sort(rec_m.aligns.begin(), rec_m.aligns.end(),
		  [](DeltaAlignment_t a, DeltaAlignment_t b){ return a.osQ < b.osQ; });

	// filter short (<30 bp) intergenic alignment
	for (int i = 0; i < rec_m.aligns.size(); i++) {
		if (rec_m.aligns[i].ge == "I0" && rec_m.aligns[i].alQ < 30) {
			rec_m.aligns.erase(rec_m.aligns.begin()+i);
			i--;
		}
	}
	
	while(!rec_m.aligns.empty()) {

		// initialize grouped alignment
		std::vector<DeltaAlignment_t> ga = {rec_m.aligns[0]};
		rec_m.aligns.erase(rec_m.aligns.begin());

		// find overlapping alignment and put it in a group
		while(!rec_m.aligns.empty() && ga[0].overlapQ(rec_m.aligns[0]) == 1) {
			ga.push_back(rec_m.aligns[0]);
			rec_m.aligns.erase(rec_m.aligns.begin());
		}

		// find representative alignment of a group
		if (ga.size() > 1) {
			std::sort(ga.begin(), ga.end(),
				  [](DeltaAlignment_t a, DeltaAlignment_t b){ return GGO[a.vdj[3]] < GGO[b.vdj[3]]; });
		}
		DeltaAlignment_t ra = ga[0];
		ra.gm.assign(ga.begin()+1, ga.end());
		galn.push_back(ra);
	}

	rec_m.aligns.assign(galn.begin(), galn.end());
}

void DeltaFilter_t::filterAlignment() {

	// exit if no alignment remains
	if (rec_m.aligns.size() == 0)
		return;

	// count alignments of each V, D, or J
	std::map<std::string, int> VDJ_c;
	for (auto ra : rec_m.aligns) {
                VDJ_c[ra.vdj]++;
		for (auto m : ra.gm)
			VDJ_c[m.vdj]++;
	}

	// filter ambiguous alignment in a group
	for (auto i = rec_m.aligns.begin(); i != rec_m.aligns.end(); i++) {

		// skip if no ambiguity
		if (i->gm.size() == 0)
			continue;

                // if V or C ambiguity
                if (i->vdj[3] == 'V' || i->vdj[3] == 'C') {

                        // keep only V or C alignments in the group
                        std::vector<DeltaAlignment_t> va = {*i};
                        for (auto m : i->gm) {
                                if (m.vdj[3] == 'V' || m.vdj[3] == 'C') 
                                        va.push_back(m);
                        }
                        if (va.size() == 1) {
                                *i = va[0];
                                i->gm.assign(va.begin()+1, va.end());
                                continue;
                        }

                        // find V or C alignment with more counts
                        int maxn = 0;
                        for (auto m : va) {
                                if (VDJ_c[m.vdj] > maxn)
                                        maxn = VDJ_c[m.vdj];
                        }
                        std::vector<DeltaAlignment_t> tva = {};
                        for (auto a : va) {
                                if (VDJ_c[a.vdj] == maxn) 
                                        tva.push_back(a);
                        }
                        
                        // reset the grouped alignment if the ambiguity can be resolved
                        if (tva.size() > 0 && tva.size() < i->gm.size()+1) {
                                *i = tva[0];
                                i->gm.assign(tva.begin()+1, tva.end());
                        }

                } else if ((i->vdj[2] == 'B' || i->vdj[2] == 'G') && i->vdj[3] == 'C') {   // if TRB/G C ambiguity   

                        // filter the C1 if J2-C1|C2
                        std::vector<DeltaAlignment_t> tca = {};
                        if ((i->ro == '+' && i > rec_m.aligns.begin() && (i-1)->vdj[3] == 'J' && (i-1)->vdj[4] == '2') ||
                            (i->ro == '-' && i < rec_m.aligns.end()-1 && (i+1)->vdj[3] == 'J' && (i+1)->vdj[4] == '2')) {
                                if (i->vdj[4] == '2')
                                        tca.push_back(*i);
                                for (auto a : i->gm) {
                                        if (a.vdj[4] == '2') 
                                                tca.push_back(a);
                                }
                        }
                        if (tca.size() == 1) {
                                *i = tca[0];
                                i->gm.assign(tca.begin()+1, tca.end());
                        }
                        
                } else if (i->vdj[2] == 'L' && i->vdj[3] == 'C') {   // if IGL C ambiguity
                        // 12367
                }
        }
        
	// filter alignment of potential CDR3 segment
	int pcdr3i = -1;
	for (auto i = rec_m.aligns.begin()+1; i < rec_m.aligns.end()-1; i++) {
		
        // skip if the alignment is good
		if (i->ge == "D0")
			continue;
		if (i->alQ >= 60)
			continue;
		if (i->alQ >= 30 && i->id > 0.9)
			continue;

		// if flanked by V2 and J0 and not broken V or J
                if ( ((i-1)->ge == "V2" && (i+1)->ge == "J0" &&
                      i->vdje != (i-1)->vdje && i->vdje != (i+1)->vdje) ||
                     ((i-1)->ge == "J0" && (i+1)->ge == "V2" &&
                      i->vdje != (i+1)->vdje && i->vdje != (i-1)->vdje) ) { 
                        pcdr3i = i - rec_m.aligns.begin();
                }
	}

	// reset alignment if a potential CDR3 is found
	if (pcdr3i != -1) {
                rec_m.aligns.erase(rec_m.aligns.begin() + pcdr3i);
	}
}

int DeltaFilter_t::maxScorePosition(std::vector<std::string> &adjseq1, std::vector<std::string> &adjseq2,
                                     std::string &rfk1, std::string &rfk2) {
        std::vector<int> sps1;
        std::vector<int> sps2;
        bool spq = (rfk1.size() == 2 && rfk2.size() == 2);
        if (spq) {
                // get splicing bonus 1
                std::string rseg1 = adjseq1[0] + rfk1;
                std::transform(rseg1.begin(), rseg1.end(), rseg1.begin(), ::toupper);
                for (size_t i = 0; i < rseg1.size()-1; i++) {
                        if (rseg1[i] == 'G' && rseg1[i+1] == 'T') {
                                sps1.push_back(SPC);
                        } else {
                                sps1.push_back(0);
                        }
                }
                // get splicing bonus 2
                std::string rseg2 = rfk2 + adjseq2[0];
                std::transform(rseg2.begin(), rseg2.end(), rseg2.begin(), ::toupper);
                for (size_t i = 2; i < rseg2.size()+1; i++) {
                        if (rseg2[i-2] == 'A' && rseg2[i-1] == 'G') {
                                sps2.push_back(SPC);
                        } else {
                                sps2.push_back(0);
                        }
                }
        }

        /*
        std::cout << "sps1:";
        for (int i: sps1) {
                std::cout << i << " ";
        }
        std::cout << std::endl;
        std::cout << "sps2:";
        for (int i: sps2) {
                std::cout << i << " ";
        }
        std::cout << std::endl;
        */
        
        // calculate score of alignment 1
        std::vector<int> qps1;
        int cs1 = 0;
        int ri1 = -1;
        if (spq) {
                qps1.push_back(sps1[ri1+1]);
        } else {
                qps1.push_back(0);
        }
        for (size_t i = 0; i < adjseq1[0].size(); i++) {
                if (adjseq1[0][i] != '_')
                        ri1 += 1;
                if (adjseq1[1][i] != '_') {
                        if (adjseq1[0][i] != '_') {
                                if (adjseq1[0][i] == adjseq1[1][i]) {
                                        cs1 += MSC + EXB;
                                } else if (adjseq1[0][i]-32 == adjseq1[1][i]) {
                                        cs1 += MSC;
                                } else {
                                        cs1 += MMSC;
                                }
                        } else {
                                cs1 += GSC;
                        }
                        if (spq) {
                                qps1.push_back(cs1 + sps1[ri1+1]);
                        } else {
                                qps1.push_back(cs1);
                        }
                } else {
                        cs1 += GSC;
                }
        }

        // calculate score of alignment 2
        std::vector<int> qps2;
        int cs2 = 0;
        int ri2 = -1;
        if (spq) {
                qps2.push_back(-sps2[ri2+1]);
        } else {
                qps2.push_back(0);
        }
        for (size_t i = 0; i < adjseq2[0].size(); i++) {
                if (adjseq2[0][i] != '_')
                        ri2 += 1;
                if (adjseq2[1][i] != '_') {
                        if (adjseq2[0][i] != '_') {
                                if (adjseq2[0][i] == adjseq2[1][i]) {
                                        cs2 += MSC + EXB;
                                } else if (adjseq2[0][i]-32 == adjseq2[1][i]) {
                                        cs2 += MSC;
                                } else {
                                        cs2 += MMSC;
                                }
                        } else {
                                cs2 += GSC;
                        }
                        if (spq) {
                                qps2.push_back(cs2 - sps2[ri2+1]);
                        } else {
                                qps2.push_back(cs2);
                        }
                } else {
                        cs2 += GSC;
                }
        }

        // combine two scores
        std::vector<int> qps;
        for (size_t i = 0; i < qps1.size(); i++) {
                qps.push_back(qps1[i] - qps2[i]);
        }

        /*
        std::cout << "qps1:";
        for (int i: qps1) {
                std::cout << i << " ";
        }
        std::cout << std::endl;
        std::cout << "qps2:";
        for (int i: qps2) {
                std::cout << i << " ";
        }
        std::cout << std::endl;
        std::cout << "qps:";
        for (int i: qps) {
                std::cout << i << " ";
        }
        std::cout << std::endl;
        */
        
        // find the max score, if same score choose the left one
        return std::distance(qps.begin(), std::max_element(qps.begin(), qps.end()));
}

void DeltaFilter_t::adjustOverlap() {

	// if no alignment remains
	if (rec_m.aligns.size() == 0)
		return;
	
	// sort by orientated query coordinate
	//std::sort(rec_m.aligns.begin(), rec_m.aligns.end(),
        //          [](DeltaAlignment_t a, DeltaAlignment_t b){ return a.osQ < b.osQ; });
	
	// update vi di ji
	//update_VDJ_index();

	// check overlap for each pair of neighboring alignments
	for (auto i = rec_m.aligns.begin()+1; i < rec_m.aligns.end(); i++) {
		
		// skip if more than one group member
		if ((i-1)->gm.size() > 0 || i->gm.size() > 0)
			continue;

		// skip if no overlap
		int ol = (i-1)->oeQ - i->osQ +1;
		if (ol <= 0)
                        continue;

                // label (i-1) to be removed if it is enclosed in i after being cut
                if (((i-1)->osQ + (i-1)->oeQ) >= 2 * i->osQ) {
                        (i-1)->tbf = true;
                        continue;
                }

		// skip if the aligned segment is at the beginning of the reference 
		if ((i-1)->sR == 1 || i->sR == 1)
                        continue;
		
                // load reference and query segments
                if ((i-1)->rseg.length()==0) {
                        (i-1)->rseg = FASTA_t::subseq(refseq_m[rec_m.idR], (i-1)->sR, (i-1)->eR);
                        (i-1)->rlfk = FASTA_t::subseq(refseq_m[rec_m.idR], (i-1)->sR-2, (i-1)->sR-1);
                        (i-1)->rrfk = FASTA_t::subseq(refseq_m[rec_m.idR], (i-1)->eR+1, (i-1)->eR+2);
                        (i-1)->qseg = FASTA_t::subseq(qryseq_m, (i-1)->osQ, (i-1)->oeQ);
                        if ((i-1)->ro == '-')
                                (i-1)->qseg = FASTA_t::revcom((i-1)->qseg);
                }
                i->rseg = FASTA_t::subseq(refseq_m[rec_m.idR], i->sR, i->eR);
                i->rlfk = FASTA_t::subseq(refseq_m[rec_m.idR], i->sR-2, i->sR-1);
                i->rrfk = FASTA_t::subseq(refseq_m[rec_m.idR], i->eR+1, i->eR+2);
                i->qseg = FASTA_t::subseq(qryseq_m, i->osQ, i->oeQ);
                if (i->ro == '-')
                        i->qseg = FASTA_t::revcom(i->qseg);

		// check adjustment of overlapping alignments
                /*
		std::cout << "***before_adj*****\n";
		std::cout << (i-1)->osQ << "-" << (i-1)->oeQ << " "
			  << i->osQ << "-" << i->oeQ << std::endl;
		std::cout << (i-1)->sQ << "-" << (i-1)->eQ << " "
			  << i->sQ << "-" << i->eQ << std::endl;
		std::cout << (i-1)->sR << "-" << (i-1)->eR << " "
			  << i->sR << "-" << i->eR << std::endl << std::endl;
		for (int d : i->deltas) {
			std::cout << d << " " << std::endl;
		}
		std::cout << "*****************" << std::endl;
                */

                // get aligned segments
                std::vector<std::string> adjseq1;
		std::vector<std::string> adjseq2;

                // along reference
                if ((i-1)->ro == i->ro) {
                        if (i->ro == '+') {
                                adjseq1 = (i-1)->getEndAlignment(-ol, false);
                                adjseq2 = i->getEndAlignment(ol, false);
                        } else {
                                adjseq1 = i->getEndAlignment(-ol, false);
                                adjseq2 = (i-1)->getEndAlignment(ol, false);
                        }
                        
                        int max_pos = maxScorePosition(adjseq1, adjseq2, (i-1)->rrfk, i->rlfk);
                        int n1 = -(ol-max_pos);
                        int n2 = max_pos;

                        if (i->ro == '+') {
                                (i-1)->cutEndAlignment(n1);
                                i->cutEndAlignment(n2);
                        } else {
                                i->cutEndAlignment(n1);
                                (i-1)->cutEndAlignment(n2);
                        }
                        
                // along query
                } else {
                        adjseq1 = (i-1)->ro == '+' ? (i-1)->getEndAlignment(-ol, false) :
                                (i-1)->getEndAlignment(ol, true);
                        adjseq2 = i->ro == '+' ? i->getEndAlignment(ol, false) :
                                i->getEndAlignment(-ol, true);

                        std::string rfk1 = "";
                        std::string rfk2 = "";
                        int max_pos = maxScorePosition(adjseq1, adjseq2, rfk1, rfk2);

                        int n1 = (i-1)->ro == '+' ? -(ol-max_pos) : (ol-max_pos);
                        (i-1)->cutEndAlignment(n1);
                        int n2 = i->ro == '+' ? max_pos : -max_pos;
                        i->cutEndAlignment(n2);
                }

                /*
		std::cout << "------\n";
		std::cout << adjseq1[0] << std::endl;
		std::cout << adjseq1[1] << std::endl;
		std::cout << adjseq2[0] << std::endl;
		std::cout << adjseq2[1] << std::endl;
		std::cout << "------\n";

		std::cout << "after adjust: " << std::endl;
		std::cout << (i-1)->osQ << "-" << (i-1)->oeQ << " " << (i-1)->alQ << " "
			  << i->osQ << "-" << i->oeQ << " " << i->alQ << std::endl;
		std::cout << (i-1)->sQ << "-" << (i-1)->eQ << " "
			  << i->sQ << "-" << i->eQ << std::endl;
		std::cout << (i-1)->sR << "-" << (i-1)->eR << " "
			  << i->sR << "-" << i->eR << std::endl << std::endl;
                for (int d : i->deltas) {
			std::cout << d << " ";
		}
		std::cout << std::endl;
                */
        }

	// filter short intergenic alignment again after adjusting overlap
	for (int i = 0; i < rec_m.aligns.size(); i++) {
		if ((rec_m.aligns[i].ge == "I0" && rec_m.aligns[i].alQ < 30) || rec_m.aligns[i].tbf == true) {
			rec_m.aligns.erase(rec_m.aligns.begin()+i);
			i--;
		}
	}

        // sort by the overall orientation
        /*
        if (ori_m == '-') {
		std::reverse(rec_m.aligns.begin(), rec_m.aligns.end());
                update_VDJ_index();
	}
        */
}

void DeltaFilter_t::setRecombCode() {

	// exit if no alignment remains
	if (rec_m.aligns.size() == 0)
		return;

	// set recombination code
        int cn = 0;    // number of C segments
        std::unordered_map<char, int> oc = { {'+', 0}, {'-', 0} };
        
	for (int i = 0; i < rec_m.aligns.size(); i++) {
                if (rec_m.aligns[i].vdj[3] == 'C')
                        cn++;
                if (rec_m.aligns[i].vdj != "TRBV30") {
                        oc[rec_m.aligns[i].ro] += 1;
                } else {
                        char o = rec_m.aligns[i].ro == '+' ? '-' : '+';
                        oc[o] += 1;
                }
        }
        int soq = 1;   // same orientation Q
        if (oc['+'] > 0 && oc['-'] > 0)
                soq = 0;
        
        // set recombination code
        rc_m = "NCH";
        if (cn > 1 || soq == 0) {
                rc_m = "CH";
                
        } else {   
                int siq = 1;   // check increasing Q
                
                // remove V30 from the alignment
                std::vector<DeltaAlignment_t> nov30aln;
                for (auto a : rec_m.aligns) {
                        if (a.vdj != "TRBV30") {
                                nov30aln.push_back(a);
                        }
                }
                if (oc['-'] > 0) {
                        std::reverse(nov30aln.begin(), nov30aln.end());
                }

                // allow C1 to be replaced by C2
                for (int i = 1; i < nov30aln.size(); i++) {
                        if (nov30aln[i-1].sR > nov30aln[i].sR) {
                                if (nov30aln[i].vdj == "TRBC1") {
                                        if (nov30aln[i-1].sR > (nov30aln[i].sR + 9346)) {
                                                siq = 0;
                                        }
                                }
                        }
                }
                
                if (siq == 0)
                        rc_m = "CH"; 
        }
}

void DeltaFilter_t::annotateQuery() {

	// exit if no alignment remains
	if (rec_m.aligns.size() == 0)
		return;
                
        // get longest non-decreasing or non-increasing sub-alignments
	LNDIS();
        
	// calculate fraction of aligned length
        
        al_m = 0;
        int qs = rec_m.aligns[0].osQ;
        int qe = rec_m.aligns[0].oeQ;
        for (int i = 1; i < rec_m.aligns.size(); i++) {
                if (rec_m.aligns[i].osQ > qe) {
                        al_m += (qe - qs + 1);
                        qs = rec_m.aligns[i].osQ;
                        qe = rec_m.aligns[i].oeQ;
                } else {
                        qe = rec_m.aligns[i].oeQ;
                }
        }
        al_m += (qe - qs + 1);

        //al_m = rec_m.aligns.back().oeQ - rec_m.aligns[0].osQ + 1;
	alf_m = ((float) al_m) / rec_m.lenQ;

	// orient alignments in the order of VDJ
	if (ori_m == '-') {
		std::reverse(rec_m.aligns.begin(), rec_m.aligns.end());
	}

	// get index of V, D, and J alignments
	vi = -1;
	di = -1;
	ji = -1;
        //int cn = 0;    // number of C segments
        //int soq = 1;   // same orientation Q
        
	for (int i = 0; i < rec_m.aligns.size(); i++) {
		if (rec_m.aligns[i].ge == "V0" || rec_m.aligns[i].ge == "V2") {
			vi = i;
                } else if (rec_m.aligns[i].ge == "D0" && (vi != -1 || di == -1)) {
			di = i;
		} else if (rec_m.aligns[i].ge == "J0" && (vi != -1 || di != -1 || ji == -1)) {
			ji = i;
		}
                /*
                if (rec_m.aligns[i].vdj[3] == 'C')
                        cn++;
                if (rec_m.aligns[i].ro != ori_m) {
                        if (rec_m.aligns[i].vdj != "TRBV30")
                                soq = 0;
                }
                */
        }

        /*
        // check increasing Q while allowing C1 to be replaced by C2
        int siq = 1;
        std::vector<DeltaAlignment_t> nov30aln;
        for (auto a : rec_m.aligns) {
                if (a.vdj != "TRBV30") {
                        nov30aln.push_back(a);
                }
        }
        for (int i = 1; i < nov30aln.size(); i++) {
                if (nov30aln[i-1].sR > nov30aln[i].sR) {
                        if (nov30aln[i].vdj[3] == 'C' && nov30aln[i].vdj[4] == '1') {
                                if (nov30aln[i-1].sR > (nov30aln[i].sR + 9346)) {
                                        siq = 0;
                                }
                        }
                }
        }
        
        // set recombination code
        rc_m = "NCH";
        if (cn > 1 || soq == 0 || siq == 0)
                rc_m = "CH";
        */    

	// set regularity
	std::string combv = "---";
	std::string combd = "---";
	std::string combj = "---";
	if (vi != -1) {
		combv = rec_m.aligns[vi].vdj;
		for (auto m = rec_m.aligns[vi].gm.begin(); m != rec_m.aligns[vi].gm.end(); m++)
			combv += "|" + m->vdj;
	}
	if (di != -1) {
		combd = rec_m.aligns[di].vdj;
		for (auto m = rec_m.aligns[di].gm.begin(); m != rec_m.aligns[di].gm.end(); m++)
			combd += "|" + m->vdj;
	}
	if (ji != -1) {
		combj = rec_m.aligns[ji].vdj;
		for (auto m = rec_m.aligns[ji].gm.begin(); m != rec_m.aligns[ji].gm.end(); m++)
			combj += "|" + m->vdj;
	}
	if (vi != -1 && ji != -1) {
		if ((ji - vi) == 1 || ((ji-vi) == 2 && di == (vi+ji)/2)) {
                        reg_m = 1;
                        int vq = 1;
                        if (vi > 1) {
                                vq = 0;
                        } else if (vi == 1) {
                                if (rec_m.aligns[0].ge != "V1" || rec_m.aligns[0].vdj != rec_m.aligns[1].vdj)
                                        vq = 0;
                        }
                        int cq = 1;
                        if (ji < rec_m.aligns.size()-2) {
                                cq = 0;
                        } else if (ji == rec_m.aligns.size()-2) {
                                if (rec_m.aligns[ji+1].ge != "C1")
                                        cq = 0;
                        }
                        if (vq == 1 && cq == 1)
                                reg_m = 2;
                } else {
                        reg_m = 0;
                }
	} else {
		reg_m = 0;
	}

        if (alf_m < 0.5) {
                reg_m = -1;
        }
        
	CombineVDJ_m = combv + ":" + combd + ":" + combj;
}

void DeltaFilter_t::LNDIS() {
	int lndsl[rec_m.aligns.size()]  = {0};
	int lnisl[rec_m.aligns.size()]  = {0};
	int ndfrom[rec_m.aligns.size()] = {0};
	int nifrom[rec_m.aligns.size()] = {0};

	for (int i = 0; i < rec_m.aligns.size(); i++) {
		lndsl[i] = rec_m.aligns[i].alQ;
		lnisl[i] = rec_m.aligns[i].alQ;
		ndfrom[i] = -1;
		nifrom[i] = -1;
	}

	/*
	std::cout << "alQ: ";
	for (int i = 0; i < rec_m.aligns.size(); i++) {
		std::cout << " " << rec_m.aligns[i].alQ << " " << rec_m.aligns[i].ge;
	}
	std::cout << std::endl;
	*/

	for (int i = 1; i < rec_m.aligns.size(); i++) {
		for (int j = 0; j < i; j++) {
			if (rec_m.aligns[j].ge != "I0" && rec_m.aligns[i].ge != "I0" &&
			    GEO[rec_m.aligns[j].ge] <= GEO[rec_m.aligns[i].ge] &&
			    (lndsl[j] + rec_m.aligns[i].alQ) > lndsl[i]) {
				lndsl[i] = lndsl[j] + rec_m.aligns[i].alQ;
				ndfrom[i] = j;
			}
			if (rec_m.aligns[j].ge != "I0" && rec_m.aligns[i].ge != "I0" &&
			    GEO[rec_m.aligns[j].ge] >= GEO[rec_m.aligns[i].ge] &&
			    (lnisl[j] + rec_m.aligns[i].alQ) > lnisl[i]) {
				lnisl[i] = lnisl[j] + rec_m.aligns[i].alQ;
				nifrom[i] = j;
			}
		}
	}

	/*
	std::cout << "lndsl: ";
	for (int i = 0; i < rec_m.aligns.size(); i++) {
		std::cout << " " << lndsl[i];
	}
	std::cout << std::endl;

	std::cout << "lnisl: ";
	for (int i = 0; i < rec_m.aligns.size(); i++) {
		std::cout << " " << lnisl[i];
	}
	std::cout << std::endl;
	*/

	int ndmaxl = 0;
	int nimaxl = 0;
	int ndmaxi = -1;
	int nimaxi = -1;
	for (int i = 0; i < rec_m.aligns.size(); i++) {
		if (lndsl[i] >= ndmaxl) {
			ndmaxl = lndsl[i];
			ndmaxi = i;
		}
		if (lnisl[i] >= nimaxl) {
			nimaxl = lnisl[i];
			nimaxi = i;
		}
	}

	if (ndmaxl >= nimaxl) {
		int ndi = ndmaxi;
		lndisi.insert(lndisi.begin(), ndi);
		while (ndfrom[ndi] != -1) {
			ndi = ndfrom[ndi];
			lndisi.insert(lndisi.begin(), ndi);
		}
		if (ndmaxl > nimaxl && lndisi.size() > 1) {
			ori_m = '+';
		} else {
			ori_m = rec_m.aligns[ndmaxi].ro;
		}
	} else {
		int nii = nimaxi;
		lndisi.insert(lndisi.begin(), nii);
		while (nifrom[nii] != -1) {
			nii = nifrom[nii];
			lndisi.insert(lndisi.begin(), nii);
		}
		ori_m = '-';
	}
}

void DeltaFilter_t::update_VDJ_index() {
	for (int i = 0; i < rec_m.aligns.size(); i++) {
		if (rec_m.aligns[i].ge == "V0" || rec_m.aligns[i].ge == "V2") {
			vi = i;
		} else if (rec_m.aligns[i].ge == "D0" && (vi != -1 || di == -1)) {
			di = i;
		} else if (rec_m.aligns[i].ge == "J0" && (vi != -1 || di != -1 || ji == -1)) {
			ji = i;
		}
	}
}

//=============================================

void DeltaFilter_t::printResult() {

	// if no alignment remains
	if (rec_m.aligns.size() == 0)
		return;
	
	std::cout << rec_m.idQ << "\t" << rec_m.lenQ;
	std::cout << "\t" << reg_m << "\t" << rec_m.idR << "\t" << ori_m;
	std::cout << "\t" << CombineVDJ_m;

	std::cout << "\t" << rec_m.aligns[0].concise_form();
	for (auto m = rec_m.aligns[0].gm.begin(); m != rec_m.aligns[0].gm.end(); m++) {
		std::cout << "|" << m->concise_form();
	}
	for (auto i = rec_m.aligns.begin()+1; i != rec_m.aligns.end(); i++) {
		std::cout << " ";
		std::cout << i->concise_form();
		for (auto m = i->gm.begin(); m != i->gm.end(); m++) {
			std::cout << "|" << m->concise_form();
		}
	}
	std::cout << std::endl;
}

void DeltaFilter_t::printResult(std::ofstream &cout) {

	// if no alignment remains
	if (rec_m.aligns.size() == 0)
		return;
	
	cout << rec_m.idQ << "\t" << rec_m.lenQ;
	cout << "\t" << reg_m << "\t" << rec_m.idR << "\t" << ori_m;
	cout << "\t" << CombineVDJ_m;

	cout << "\t" << rec_m.aligns[0].concise_form();

	for (auto m = rec_m.aligns[0].gm.begin(); m != rec_m.aligns[0].gm.end(); m++) {
		cout << "|" << m->concise_form();
	}
	for (auto i = rec_m.aligns.begin()+1; i != rec_m.aligns.end(); i++) {
		cout << " ";
		cout << i->concise_form();
		for (auto m = i->gm.begin(); m != i->gm.end(); m++) {
			cout << "|" << m->concise_form();
		}
	}
        cout << "\t" << vi << "," << di << "," << ji << "\t" << al_m;
	cout << std::endl;
}
