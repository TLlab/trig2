#ifndef FASTX_READ_HPP
#define FASTX_READ_HPP

#include <iostream>
#include <string>
#include <fstream>
#include <unordered_map>

namespace FASTA_t
{
	// load fasta file and convert to unordered map
	std::unordered_map<std::string, std::string> getfasta(const std::string &fasta_path);

	// subsequence (position start from 1)
	std::string subseq(const std::string &seq, const int str, const int end);

	// reverse complement
	std::string revcom(std::string seq);
}

// load fasta file using "while"
class FastaReader_t
{
private:
	std::string fasta_path_m;      // path
	std::ifstream fasta_stream_m;  // stream
	bool is_open_m;                // stream is open

	std::string uid_m;             // fasta id(>)
	std::string seq_m;             // fasta sequnece

	void CheckStream() {
		if (!fasta_stream_m.good()) {
			std::cerr << "\033[31mERROR:\033[0m Could not parse fasta file, "
				<< fasta_path_m << std::endl;
			exit(1);
		}
	}

public:
	FastaReader_t() {
		is_open_m = false;
	}
	~FastaReader_t() {
		fasta_stream_m.close();
		is_open_m = false;
	}

	void open(const std::string &fasta_path);
	bool readNext();

	const std::string &getUID() const {
		return uid_m;
	}
	const std::string &getSEQ() const {
		return seq_m;
	}
};

class FastqReader_t
{
private:
	std::string fastq_path_m;
	std::ifstream fastq_stream_m;
	bool is_open_m;

	std::string uid_m;
    std::string seq_m;
	std::string qid_m;
	std::string qua_m;

	void CheckStream() {
		if (!fastq_stream_m.good()) {
			std::cerr << "\033[31mERROR:\033[0m Could not parse fastq file, "
				<< fastq_path_m << std::endl;
			exit(1);
		}
	}

public:
	FastqReader_t() {
		is_open_m = false;
	}
	~FastqReader_t() {
		fastq_stream_m.close();
		is_open_m = false;
	}

	void open(const std::string &fastq_path);
	bool readNext();

	const std::string &getUID() const {
		return uid_m;
	}
	const std::string &getSEQ() const {
		return seq_m;
	}
	const std::string &getQID() const {
		return qid_m;
	};
	const std::string &getQUA() const {
		return qua_m;
	};
};

#endif /* fastx_read.hpp  */
