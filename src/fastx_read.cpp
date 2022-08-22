#include "fastx_read.hpp"

#include <iostream>
#include <string>
#include <fstream>
#include <unordered_map>
#include <algorithm>


std::unordered_map<std::string, std::string> FASTA_t::getfasta(const std::string &fasta_path) {
	std::unordered_map<std::string, std::string> fasta;
	std::ifstream fasta_stream;
	fasta_stream.open(fasta_path);

	// check stream
	if (!fasta_stream.good()) {
		std::cerr << "\033[31mERROR:\033[0m Could not parse fasta file, " 
			<< fasta_path << std::endl;
		exit(1);
	}

	// load fasta
	std::string uid, seq;
	while (fasta_stream.peek() != EOF) {
		fasta_stream.get();
		fasta_stream >> uid;
		fasta_stream >> seq;
		fasta[uid] = seq;
		while(fasta_stream.get() != '\n');
		while(fasta_stream.peek() != EOF && fasta_stream.peek() != '>') {
			fasta_stream >> seq;
			fasta[uid] += seq;
			while(fasta_stream.get() != '\n');
		}
	}
	fasta_stream.close();
	return fasta;
}

// subsequence (position start from 1)
std::string FASTA_t::subseq(const std::string &seq, const int str, const int end) {
	int pos = str - 1;
	int len = end - str + 1;
	return seq.substr(pos, len);
}

std::string FASTA_t::revcom(std::string seq) {
	reverse(seq.begin(), seq.end());
	for (size_t i = 0; i < seq.length(); i++) {
		switch (seq[i]) {
			case 'A' :
				seq[i] = 'T';
				break;
			case 'C' :
				seq[i] = 'G';
				break;
			case 'G' :
				seq[i] = 'C';
				break;
			case 'T' :
				seq[i] = 'A';
				break;
			case 'a' :
				seq[i] = 't';
				break;
			case 'c' :
				seq[i] = 'g';
				break;
			case 'g' :
				seq[i] = 'c';
				break;
			case 't' :
				seq[i] = 'a';
				break;
		}
	}
	return seq;
}

void FastaReader_t::open(const std::string &fasta_path) {
	fasta_path_m = fasta_path;
	fasta_stream_m.open(fasta_path_m);
	CheckStream();
	is_open_m = true;
}

bool FastaReader_t::readNext() {
	if (fasta_stream_m.peek() == EOF)
		return false;

	// make way for the new ID and sequence
	uid_m.clear();
	seq_m.clear();

	// load one fasta
	std::string seq;
	fasta_stream_m.get();
	fasta_stream_m >> uid_m;
	fasta_stream_m >> seq_m;
	while(fasta_stream_m.get() != '\n');
	while(fasta_stream_m.peek() != EOF && fasta_stream_m.peek() != '>') {
		fasta_stream_m >> seq;
		seq_m += seq;
		while(fasta_stream_m.get() != '\n');
	}
	return true;
}

void FastqReader_t::open(const std::string &fastq_path) {
	fastq_path_m = fastq_path;
	fastq_stream_m.open(fastq_path_m);
	CheckStream();
	is_open_m = true;
}

bool FastqReader_t::readNext() {
    if (fastq_stream_m.peek() == EOF)
        return false;

    // make way for the new ID and sequence
    uid_m.clear();
    seq_m.clear();
    qid_m.clear();
    qua_m.clear();

    // load one fastq
    fastq_stream_m.get();
    fastq_stream_m >> uid_m;
    while(fastq_stream_m.get() != '\n');
    fastq_stream_m >> seq_m;
    while(fastq_stream_m.get() != '\n');
    fastq_stream_m >> qid_m;
    while(fastq_stream_m.get() != '\n');
    fastq_stream_m >> qua_m;
    while(fastq_stream_m.get() != '\n');

    return true;
}
