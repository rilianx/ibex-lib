//============================================================================
//                                  I B E X
// File        : ibex_Json.h
// Author      : Ignacio Araya
// License     : See the LICENSE file
// Created     : 2026
//
// Minimal dependency-free JSON reader/writer, used by the machine-learning
// data-collection server (see #ibex::MLNodeServer).
//
// Non-standard extension (on purpose): the literals "Infinity", "-Infinity"
// and "NaN" are both accepted and emitted, because they are what Python's
// "json" module produces and consumes by default, and interval bounds are
// routinely infinite.
//============================================================================

#ifndef __IBEX_JSON_H__
#define __IBEX_JSON_H__

#include "ibex_IntervalVector.h"

#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace ibex {

/**
 * \ingroup ml
 * \brief Exception raised on malformed JSON or on a type mismatch.
 */
class JsonError : public std::runtime_error {
public:
	explicit JsonError(const std::string& msg) : std::runtime_error(msg) { }
};

/**
 * \ingroup ml
 * \brief A parsed JSON value.
 */
class JsonValue {
public:
	typedef enum { NUL, BOOL, NUM, STR, ARR, OBJ } Type;

	JsonValue();

	/**
	 * \brief Parse a complete JSON document.
	 *
	 * \throw JsonError if the text is malformed or has trailing garbage.
	 */
	static JsonValue parse(const std::string& text);

	/**
	 * \brief Look up a key.
	 *
	 * \return NULL if this is not an object or the key is absent.
	 */
	const JsonValue* get(const std::string& key) const;

	bool is_null() const;

	/* Strict accessors (throw JsonError on a type mismatch). */
	double as_double() const;
	int as_int() const;
	long as_long() const;
	bool as_bool() const;
	const std::string& as_string() const;

	/* Lenient accessors: an absent key or a null value yields "def". */
	double get_double(const std::string& key, double def) const;
	int get_int(const std::string& key, int def) const;
	long get_long(const std::string& key, long def) const;
	bool get_bool(const std::string& key, bool def) const;
	std::string get_string(const std::string& key, const std::string& def) const;

	/**
	 * \brief Read a box given as [[lb,ub],...] (or as a flat list of reals,
	 *        each one being read as a degenerate interval).
	 *
	 * \throw JsonError if the value is not a box of the expected size
	 *        (any size is accepted if expected_size is -1).
	 */
	IntervalVector as_box(int expected_size=-1) const;

	Type type;
	bool boolean;
	double number;
	std::string str;   // named "str", not "string", so that "std::string" stays usable in member definitions
	std::vector<JsonValue> array;
	std::vector<std::pair<std::string,JsonValue> > object;
};

/**
 * \ingroup ml
 * \brief Streaming JSON writer.
 *
 * Commas and the object/array nesting are handled automatically, so that large
 * feature dumps can be written without ever materializing a value tree.
 */
class JsonOut {
public:
	explicit JsonOut(std::ostream& os);

	JsonOut& obj();
	JsonOut& end_obj();
	JsonOut& arr();
	JsonOut& end_arr();

	/** \brief Write a key (only valid inside an object). */
	JsonOut& key(const char* k);

	JsonOut& val(double d);
	JsonOut& val(int v);
	JsonOut& val(long v);
	JsonOut& val(size_t v);
	JsonOut& val(bool b);
	JsonOut& val(const char* s);
	JsonOut& val(const std::string& s);
	JsonOut& null();

	/** \brief An interval, written as the pair [lb,ub]. */
	JsonOut& val(const Interval& x);

	/** \brief A box, written as [[lb,ub],...]. An empty box is written as null. */
	JsonOut& val(const IntervalVector& box);

	/** \brief A vector, written as [x0,x1,...]. */
	JsonOut& val(const Vector& v);

	/* "key: value" shortcuts. */
	JsonOut& kv(const char* k, double d);
	JsonOut& kv(const char* k, int v);
	JsonOut& kv(const char* k, long v);
	JsonOut& kv(const char* k, size_t v);
	JsonOut& kv(const char* k, bool b);
	JsonOut& kv(const char* k, const char* s);
	JsonOut& kv(const char* k, const std::string& s);
	JsonOut& kv(const char* k, const Interval& x);
	JsonOut& kv(const char* k, const IntervalVector& box);
	JsonOut& kv(const char* k, const Vector& v);
	JsonOut& knull(const char* k);

	/** \brief Write a double the way this writer does (Infinity/NaN aware). */
	static void write_double(std::ostream& os, double d);

	/** \brief Write a JSON string literal (with escaping). */
	static void write_string(std::ostream& os, const std::string& s);

protected:
	/* Emit a separator if an item was already written at the current level. */
	void sep();

	std::ostream& os;
	std::vector<int> count;  // number of items written, per nesting level
	bool expect_value;       // true right after key(): no separator, no counting
};

} // end namespace ibex

#endif // __IBEX_JSON_H__
