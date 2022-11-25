#include "vdjreader.hpp"

#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <algorithm>


void VDJReader_t::open(const std::string &vdj_path) {
    vdj_path_m = vdj_path;
    vdj_stream_m.open(vdj_path_m);
    CheckStream();
    is_open_m = true;
}

bool VDJReader_t::readNext() {
	if (vdj_stream_m.peek() == EOF)
		return false;

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
		info_m[info.species_gene].push_back(info);
	} else {
		if (info.strand == '+') {
			for (int i = 1; i <= num_of_exon; i++) {
				info.vdj_exon = info.vdj + "_" + std::to_string(i);
				info.exon_start = each_sted[i*2-2];
				info.exon_end = each_sted[i*2-1];
				info_m[info.species_gene].push_back(info);
			}
		} else { // strand == '-'
			for (int i = 1; i <= num_of_exon; i++) {
				info.vdj_exon = info.vdj + "_" + std::to_string(num_of_exon-i+1);
				info.exon_start = each_sted[i*2-2];
				info.exon_end = each_sted[i*2-1];
				info_m[info.species_gene].push_back(info);
			}
		}
	}
	return true;
}

