//============================================================================
//                                  I B E X                                   
// File        : ibex_Timer.h
// Author      : Jordan Ninin
// Copyright   : IMT Atlantique (France)
// License     : See the LICENSE file
// Created     : May 13, 2012
// Last Update : May 13, 2012
//============================================================================

#ifndef __IBEX_TIMER_H__
#define __IBEX_TIMER_H__

#include "ibex_Exception.h"

#ifndef _WIN32
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>

#else
/* A time value that is accurate to the nearest
   microsecond but also has a range of years.  */
struct mytimeval
  {
	long int tv_sec;		/* Seconds.  */
	long int tv_usec;	/* Microseconds.  */
  };
#endif // _MSC_VER

namespace ibex {

/** \ingroup tools
 *
 * \brief Timeout exception.
 *
 */
class TimeOutException : public Exception { };

class Timer;
/**
 * \brief Time limit of the running search, checked from inside contractors.
 *
 * The optimizer checks its time limit between nodes only; a contraction that
 * does not converge (fixpoint of the polytope hull, ACID, propagation) can then
 * run forever. When set, check_search_deadline() throws TimeOutException once
 * search_timer exceeds search_timeout. It is checked every 256 calls (the
 * clock is a system call).
 */
extern Timer* search_timer;
extern double search_timeout;
void check_search_deadline_slow();
inline void check_search_deadline() {
	static unsigned int n = 0;
	if (search_timer && ((++n) & 255) == 0) check_search_deadline_slow();
}


/** \ingroup tools
 *
 * \brief Timer.
 */
class Timer {
 public:
	Timer();

	void start();
	void stop();
	void restart();
	double get_time();
	void check(double timeout);
private:
	  double start_time;
	  bool active;
};

/** \ingroup tools
 *
 * \brief Timer.
 */
class StaticTimer {
 public:
  static void start();

 private:
  friend class Timer;

  typedef double Time;

  static Time get_localtime();

  static Time local_time;

#ifdef _WIN32
  static struct mytimeval tp;
  static Time real_lapse;
  static Time real_time;

#else
  static struct rusage res;
  static Time virtual_ulapse;
  static Time virtual_slapse;
  static Time virtual_utime;
  static Time virtual_stime;
  static long resident_memory;
#endif
};

} // end namespace ibex
#endif // __IBEX_TIMER_H__
