//============================================================================
//                                  I B E X
// File        : ibex_Cov.cpp
// Author      : Gilles Chabert
// Copyright   : IMT Atlantique (France)
// License     : See the LICENSE file
// Created     : Nov 07, 2018
// Last update : Dec 24, 2018
//============================================================================

#include "ibex_Cov.h"

#include <sstream>
#include <cassert>
#include <string.h>

using namespace std;

namespace ibex {

const unsigned int Cov::subformat_level = 0;
const unsigned int Cov::subformat_number = 0;
const size_t Cov::SIGNATURE_LENGTH = 20;
const char* Cov::SIGNATURE = "IBEX COVERING FILE ";
const uint32_t Cov::FORMAT_VERSION = 1;
const string Cov::separator = "+-------------------+-----------------------------------------------------------\n";
const string Cov::space     = "|                   |";

Cov::Cov(size_t n) : n(n) {

}

Cov::Cov(const char* filename) : n(0 /* tmp */) {
	stack<unsigned int> format_seq;
	ifstream* f = Cov::read(filename, *this, format_seq);
	f->close();
	delete f;
}

void Cov::save(const char* filename) const {
	stack<unsigned int> format_seq;
	ofstream* of=Cov::write(filename, *this, format_seq);
	of->close();
	delete of;
}

Cov::~Cov() {

}

ifstream* Cov::read(const char* filename, Cov& cov, stack<unsigned int>& format_seq) {

	ifstream* f = new ifstream();

	f->open(filename, ios::in | ios::binary);

	if (f->fail()) ibex_error("[Cov]: cannot open input file.\n");

	int input_format_version = read_signature(*f);

	if (input_format_version!=FORMAT_VERSION)
		ibex_error("[Cov] unsupported format version");

	read_format_seq(*f, format_seq);

	if (format_seq.top()!=subformat_number) return f;
	else format_seq.pop();

	size_t _n = read_pos_int(*f);

	(size_t&) cov.n = _n;

	return f;
}

ofstream* Cov::write(const char* filename, const Cov& cov, std::stack<unsigned int>& format_seq) {

	ofstream* f = new ofstream();

	f->open(filename, ios::out | ios::binary);

	if (f->fail())
		ibex_error("[Cov]: cannot create output file.\n");

	write_signature(*f);

	format_seq.push(subformat_number);

	write_format_seq(*f, format_seq);

	write_int(*f, cov.n);

	return f;
}

int Cov::read_signature(ifstream& f) {
	char* sig=new char[SIGNATURE_LENGTH];
	f.read(sig, SIGNATURE_LENGTH*sizeof(char));
	if (f.eof()) ibex_error("[manifold]: unexpected end of file.");
	if (strcmp(sig,SIGNATURE)!=0)
	ibex_error("[Cov]: not an Ibex \"cover\" file.");
	uint32_t format_version=read_pos_int(f); // manifold version
	if (format_version>FORMAT_VERSION) {
		cout << "format version = " << format_version << endl;
		ibex_error("[Cov]: wrong format version");
		return -1;
	} else
		return format_version;
}

void Cov::read_format_seq(std::ifstream& f, stack<unsigned int>& format_seq) {
	size_t format_level = read_pos_int(f);

	list<unsigned int> tmp;
	for (size_t i=0; i<=format_level; i++)
		tmp.push_front(read_pos_int(f));

	for (list<unsigned int>::const_iterator it=tmp.begin(); it!=tmp.end(); ++it)
		format_seq.push(*it);
}

unsigned int Cov::read_pos_int(ifstream& f) {
	uint32_t x;
	f.read((char*) &x, sizeof(x)); //f >> x;
	if (f.eof()) ibex_error("[Cov]: unexpected end of file.");
	if (x<0) ibex_error("[Cov]: negative integer unexpected");
	return x;
}

double Cov::read_double(ifstream& f) {
	double x;
	f.read((char*) &x, sizeof(x)); //f >> x;
	if (f.eof()) ibex_error("[Cov]: unexpected end of file.");
	return x;
}

void Cov::write_signature(ofstream& f) {
	f.write(SIGNATURE, SIGNATURE_LENGTH*sizeof(char));
	write_int(f, FORMAT_VERSION);
}

void Cov::write_format_seq(std::ofstream& f, std::stack<unsigned int>& format_seq) {
	write_int(f, format_seq.size()-1); // level = size-1 (starts from 0)

	while (!format_seq.empty()) {
		write_int(f, format_seq.top());
		format_seq.pop();
	}
}

void Cov::write_int(ofstream& f, uint32_t x) {
	f.write((char*) &x, sizeof(uint32_t));
}

void Cov::write_double(ofstream& f, double x) {
	f.write((char*) &x, sizeof(x));
}

void Cov::format(stringstream& ss, const string& title, stack<unsigned int>& format_seq) {

	format_seq.push(subformat_number);

	ss
	<< "\n"
	<< "+-------------------------------------------------------------------------------\n"
	<< "|                          " << title << " file format v" << FORMAT_VERSION << "\n"
	"|\n"
	//			"The " << title << " text format (obtained with --txt) is described below.\n"
	//			"The " << title << " binary format (.cov) is exactly the same except that:\n"
	//			"  - all separating characters (space, line return) are removed except\n"
	//			"    those inside the signature (line 1 in text format)\n"
	"| Note:\n"
	"| - integer values are unsigned 32 bits integer (uint32_t)\n"
	"| - real values are 64 bits double\n"
	"| - all indices start from 0.\n"
	<< separator
	<< space << " - the signature: the null-terminated sequence of " << SIGNATURE_LENGTH <<  "\n"
	<< space <<	"                  characters \"" << SIGNATURE << "\"\n"
	<< "|        Cov        |" <<
	"                  (mind the space at the end)\n"
	<< space << "                  followed by the format version number: " << FORMAT_VERSION << "\n"
	<< space << " - 1 integer:     the subformat level L (=" << format_seq.size()-1 << " in the case of\n"
	<< space << "                  " << title << ")\n"
	<< space << " - L integers:    the subformat identifying sequence:\n"
	<< space << "                      ";
	while (!format_seq.empty()) {
		ss << format_seq.top();
		format_seq.pop();
		if (!format_seq.empty()) ss << ' ';
	}
	ss
	<< "\n"
	<< space << "                  (in the case of " << title << ")\n"
	<< space << " - 1 integer:     the dimension n of boxes (# of variables)\n"
	<< separator;
}

string Cov::format() {
	stringstream ss;
	stack<unsigned int> format_seq;
	format(ss, "COV", format_seq);
	return ss.str();
}

} /* namespace ibex */