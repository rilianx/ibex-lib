//============================================================================
//                                  I B E X
// File        : ibex_Json.cpp
// Author      : Ignacio Araya
// License     : See the LICENSE file
// Created     : 2026
//============================================================================

#include "ibex_Json.h"

#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sstream>

using namespace std;

namespace ibex {

/*================================ JsonValue ================================*/

JsonValue::JsonValue() : type(NUL), boolean(false), number(0.0) { }

namespace {

/* Recursive-descent parser over a null-terminated buffer. */
class Parser {
public:
	Parser(const string& s) : s(s), i(0) { }

	JsonValue parse_document() {
		JsonValue v = parse_value();
		skip_ws();
		if (i != s.size()) fail("trailing characters");
		return v;
	}

private:
	const string& s;
	size_t i;

	void fail(const string& msg) {
		ostringstream o;
		o << "JSON error at offset " << i << ": " << msg;
		throw JsonError(o.str());
	}

	void skip_ws() {
		while (i<s.size() && (s[i]==' ' || s[i]=='\t' || s[i]=='\n' || s[i]=='\r')) i++;
	}

	bool literal(const char* lit) {
		size_t n = strlen(lit);
		if (s.compare(i, n, lit)==0) { i += n; return true; }
		return false;
	}

	JsonValue parse_value() {
		skip_ws();
		if (i>=s.size()) fail("unexpected end of input");

		switch (s[i]) {
		case '{' : return parse_object();
		case '[' : return parse_array();
		case '"' : {
			JsonValue v; v.type=JsonValue::STR; v.str=parse_string(); return v;
		}
		case 't' : {
			if (!literal("true")) fail("invalid literal");
			JsonValue v; v.type=JsonValue::BOOL; v.boolean=true; return v;
		}
		case 'f' : {
			if (!literal("false")) fail("invalid literal");
			JsonValue v; v.type=JsonValue::BOOL; v.boolean=false; return v;
		}
		case 'n' : {
			if (!literal("null")) fail("invalid literal");
			return JsonValue();
		}
		case 'N' : {
			if (!literal("NaN")) fail("invalid literal");
			JsonValue v; v.type=JsonValue::NUM; v.number=NAN; return v;
		}
		case 'I' : {
			if (!literal("Infinity")) fail("invalid literal");
			JsonValue v; v.type=JsonValue::NUM; v.number=POS_INFINITY; return v;
		}
		default :
			return parse_number();
		}
	}

	JsonValue parse_number() {
		// "-Infinity" is the only negative literal that is not a plain number
		if (s[i]=='-' && s.compare(i,9,"-Infinity")==0) {
			i += 9;
			JsonValue v; v.type=JsonValue::NUM; v.number=NEG_INFINITY; return v;
		}
		size_t start = i;
		if (i<s.size() && (s[i]=='-' || s[i]=='+')) i++;
		while (i<s.size() && (isdigit((unsigned char) s[i]) || s[i]=='.' ||
				s[i]=='e' || s[i]=='E' ||
				((s[i]=='-' || s[i]=='+') && (s[i-1]=='e' || s[i-1]=='E')))) i++;
		if (i==start) fail("invalid value");

		string tok = s.substr(start, i-start);
		char* end = NULL;
		double d = strtod(tok.c_str(), &end);
		if (end==NULL || *end!='\0') fail("invalid number '" + tok + "'");

		JsonValue v; v.type=JsonValue::NUM; v.number=d; return v;
	}

	string parse_string() {
		if (s[i]!='"') fail("expected '\"'");
		i++;
		string out;
		while (true) {
			if (i>=s.size()) fail("unterminated string");
			char c = s[i++];
			if (c=='"') break;
			if (c!='\\') { out.push_back(c); continue; }
			if (i>=s.size()) fail("unterminated escape");
			char e = s[i++];
			switch (e) {
			case '"' : out.push_back('"');  break;
			case '\\': out.push_back('\\'); break;
			case '/' : out.push_back('/');  break;
			case 'b' : out.push_back('\b'); break;
			case 'f' : out.push_back('\f'); break;
			case 'n' : out.push_back('\n'); break;
			case 'r' : out.push_back('\r'); break;
			case 't' : out.push_back('\t'); break;
			case 'u' : {
				if (i+4>s.size()) fail("truncated \\u escape");
				unsigned int cp = 0;
				for (int k=0; k<4; k++) {
					char h = s[i++];
					cp <<= 4;
					if (h>='0' && h<='9') cp |= (unsigned int) (h-'0');
					else if (h>='a' && h<='f') cp |= (unsigned int) (h-'a'+10);
					else if (h>='A' && h<='F') cp |= (unsigned int) (h-'A'+10);
					else fail("invalid \\u escape");
				}
				// UTF-8 encoding of the BMP code point (surrogates not decoded)
				if (cp<0x80) out.push_back((char) cp);
				else if (cp<0x800) {
					out.push_back((char) (0xC0 | (cp>>6)));
					out.push_back((char) (0x80 | (cp & 0x3F)));
				} else {
					out.push_back((char) (0xE0 | (cp>>12)));
					out.push_back((char) (0x80 | ((cp>>6) & 0x3F)));
					out.push_back((char) (0x80 | (cp & 0x3F)));
				}
				break;
			}
			default: fail("invalid escape");
			}
		}
		return out;
	}

	JsonValue parse_array() {
		JsonValue v; v.type=JsonValue::ARR;
		i++; // '['
		skip_ws();
		if (i<s.size() && s[i]==']') { i++; return v; }
		while (true) {
			v.array.push_back(parse_value());
			skip_ws();
			if (i>=s.size()) fail("unterminated array");
			if (s[i]==',') { i++; continue; }
			if (s[i]==']') { i++; break; }
			fail("expected ',' or ']'");
		}
		return v;
	}

	JsonValue parse_object() {
		JsonValue v; v.type=JsonValue::OBJ;
		i++; // '{'
		skip_ws();
		if (i<s.size() && s[i]=='}') { i++; return v; }
		while (true) {
			skip_ws();
			string k = parse_string();
			skip_ws();
			if (i>=s.size() || s[i]!=':') fail("expected ':'");
			i++;
			v.object.push_back(make_pair(k, parse_value()));
			skip_ws();
			if (i>=s.size()) fail("unterminated object");
			if (s[i]==',') { i++; continue; }
			if (s[i]=='}') { i++; break; }
			fail("expected ',' or '}'");
		}
		return v;
	}
};

} // anonymous namespace

JsonValue JsonValue::parse(const string& text) {
	Parser p(text);
	return p.parse_document();
}

const JsonValue* JsonValue::get(const string& key) const {
	if (type!=OBJ) return NULL;
	for (size_t k=0; k<object.size(); k++)
		if (object[k].first==key) return &object[k].second;
	return NULL;
}

bool JsonValue::is_null() const { return type==NUL; }

double JsonValue::as_double() const {
	if (type!=NUM) throw JsonError("expected a number");
	return number;
}

int JsonValue::as_int() const { return (int) as_double(); }

long JsonValue::as_long() const { return (long) as_double(); }

bool JsonValue::as_bool() const {
	if (type==BOOL) return boolean;
	if (type==NUM) return number!=0.0;
	throw JsonError("expected a boolean");
}

const string& JsonValue::as_string() const {
	if (type!=STR) throw JsonError("expected a string");
	return str;
}

double JsonValue::get_double(const std::string& key, double def) const {
	const JsonValue* v = get(key);
	return (v==NULL || v->is_null()) ? def : v->as_double();
}

int JsonValue::get_int(const std::string& key, int def) const {
	const JsonValue* v = get(key);
	return (v==NULL || v->is_null()) ? def : v->as_int();
}

long JsonValue::get_long(const std::string& key, long def) const {
	const JsonValue* v = get(key);
	return (v==NULL || v->is_null()) ? def : v->as_long();
}

bool JsonValue::get_bool(const std::string& key, bool def) const {
	const JsonValue* v = get(key);
	return (v==NULL || v->is_null()) ? def : v->as_bool();
}

string JsonValue::get_string(const std::string& key, const std::string& def) const {
	const JsonValue* v = get(key);
	return (v==NULL || v->is_null()) ? def : v->as_string();
}

IntervalVector JsonValue::as_box(int expected_size) const {
	if (type!=ARR) throw JsonError("expected a box (array)");
	int n = (int) array.size();
	if (expected_size>=0 && n!=expected_size) {
		ostringstream o;
		o << "expected a box of size " << expected_size << ", got " << n;
		throw JsonError(o.str());
	}
	IntervalVector box(n>0 ? n : 1);
	if (n==0) throw JsonError("empty box");

	for (int j=0; j<n; j++) {
		const JsonValue& e = array[j];
		if (e.type==NUM) {
			box[j] = Interval(e.number);
		} else if (e.type==ARR) {
			if (e.array.size()!=2) throw JsonError("an interval must be [lb,ub]");
			double lb = e.array[0].as_double();
			double ub = e.array[1].as_double();
			if (lb>ub) throw JsonError("an interval must have lb<=ub");
			box[j] = Interval(lb, ub);
		} else
			throw JsonError("an interval must be [lb,ub] or a real");
	}
	return box;
}

/*================================= JsonOut =================================*/

JsonOut::JsonOut(std::ostream& os) : os(os), expect_value(false) { }

void JsonOut::write_double(std::ostream& os, double d) {
	if (std::isnan(d)) { os << "NaN"; return; }
	if (d==POS_INFINITY) { os << "Infinity"; return; }
	if (d==NEG_INFINITY) { os << "-Infinity"; return; }
	char buf[40];
	// 17 significant digits round-trip through an IEEE-754 double
	snprintf(buf, sizeof(buf), "%.17g", d);
	os << buf;
}

void JsonOut::write_string(std::ostream& os, const std::string& s) {
	os << '"';
	for (size_t k=0; k<s.size(); k++) {
		unsigned char c = (unsigned char) s[k];
		switch (c) {
		case '"' : os << "\\\""; break;
		case '\\': os << "\\\\"; break;
		case '\b': os << "\\b";  break;
		case '\f': os << "\\f";  break;
		case '\n': os << "\\n";  break;
		case '\r': os << "\\r";  break;
		case '\t': os << "\\t";  break;
		default:
			if (c<0x20) {
				char buf[8];
				snprintf(buf, sizeof(buf), "\\u%04x", c);
				os << buf;
			} else
				os << (char) c;
		}
	}
	os << '"';
}

void JsonOut::sep() {
	if (expect_value) { expect_value=false; return; } // value of a key: no comma
	if (!count.empty()) {
		if (count.back()>0) os << ',';
		count.back()++;
	}
}

JsonOut& JsonOut::obj()     { sep(); os << '{'; count.push_back(0); return *this; }
JsonOut& JsonOut::end_obj() { os << '}'; count.pop_back(); return *this; }
JsonOut& JsonOut::arr()     { sep(); os << '['; count.push_back(0); return *this; }
JsonOut& JsonOut::end_arr() { os << ']'; count.pop_back(); return *this; }

JsonOut& JsonOut::key(const char* k) {
	sep();
	write_string(os, k);
	os << ':';
	expect_value = true;
	return *this;
}

JsonOut& JsonOut::val(double d)  { sep(); write_double(os, d); return *this; }
JsonOut& JsonOut::val(int v)     { sep(); os << v; return *this; }
JsonOut& JsonOut::val(long v)    { sep(); os << v; return *this; }
JsonOut& JsonOut::val(size_t v)  { sep(); os << v; return *this; }
JsonOut& JsonOut::val(bool b)    { sep(); os << (b ? "true" : "false"); return *this; }
JsonOut& JsonOut::val(const char* s)        { sep(); write_string(os, s); return *this; }
JsonOut& JsonOut::val(const std::string& s) { sep(); write_string(os, s); return *this; }
JsonOut& JsonOut::null()         { sep(); os << "null"; return *this; }

JsonOut& JsonOut::val(const Interval& x) {
	if (x.is_empty()) return null();
	arr(); val(x.lb()); val(x.ub()); return end_arr();
}

JsonOut& JsonOut::val(const IntervalVector& box) {
	if (box.is_empty()) return null();
	arr();
	for (int j=0; j<box.size(); j++) val(box[j]);
	return end_arr();
}

JsonOut& JsonOut::val(const Vector& v) {
	arr();
	for (int j=0; j<v.size(); j++) val(v[j]);
	return end_arr();
}

JsonOut& JsonOut::kv(const char* k, double d)  { return key(k).val(d); }
JsonOut& JsonOut::kv(const char* k, int v)     { return key(k).val(v); }
JsonOut& JsonOut::kv(const char* k, long v)    { return key(k).val(v); }
JsonOut& JsonOut::kv(const char* k, size_t v)  { return key(k).val(v); }
JsonOut& JsonOut::kv(const char* k, bool b)    { return key(k).val(b); }
JsonOut& JsonOut::kv(const char* k, const char* s)        { return key(k).val(s); }
JsonOut& JsonOut::kv(const char* k, const std::string& s) { return key(k).val(s); }
JsonOut& JsonOut::kv(const char* k, const Interval& x)       { return key(k).val(x); }
JsonOut& JsonOut::kv(const char* k, const IntervalVector& b) { return key(k).val(b); }
JsonOut& JsonOut::kv(const char* k, const Vector& v)         { return key(k).val(v); }
JsonOut& JsonOut::knull(const char* k) { return key(k).null(); }

} // end namespace ibex
