#ifndef _VDJREADER_H_
#define _VDJREADER_H_

#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <algorithm>

/* 
  usage: 
  VDJReader_t vr;
  std::unordered_map<std::string, std::vector<VDJInfo_t>> vdj = vr.getallVDJInfo();
*/

// information of a VDJ gene detail                                                                                                                                             
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

// for reading VDJ information file
class VDJReader_t
{
private:
    std::string vdj_path_m;          // vdj file path
    std::ifstream vdj_stream_m;      // vdj file input stream
	
	std::unordered_map< std::string, std::vector<VDJInfo_t> > info_m;   // vdj info of all exons of a gene
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

	std::unordered_map< std::string, std::vector<VDJInfo_t> > getallVDJInfo(const std::string &vdj_path) {
		open(vdj_path);
		while(readNext()){};
		return info_m;
	}

};

#endif /* vdjreader.h */
