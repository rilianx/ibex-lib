//============================================================================
//                                  I B E X                                   
// File        : ibex_CellSet.h
// Author      : Ignacio Araya
// Copyright   : Ecole des Mines de Nantes (France)
// License     : See the LICENSE file
// Created     : Oct 28, 2014
// Last Update : Oct 28, 2014
//============================================================================

#ifndef __IBEX_CELL_SET_H__
#define __IBEX_CELL_SET_H__

#include "ibex_CellBuffer.h"
#include <set>

namespace ibex {

/** \ingroup strategy
 *
 * \brief Cell Set.
 *
 * For depth-first search. \see #CellBuffer
 */
class CellStack : public CellBuffer {
 public:
  /** Flush the buffer.
   * All the remaining cells will be *deleted* */
  void flush();

  /** Return the size of the buffer. */
  int size() const;

  /** Return true if the buffer is empty. */
  bool empty() const;

  /** push a new cell on the stack. */
  void push(Cell* cell);

  /** Pop a cell from the stack and return it.*/
  Cell* pop();

  /** Return the next box (but does not pop it).*/
  Cell* top() const;

 private:
  /* Stack of cells */
  std::set<Cell*> cstack;
};

} // end namespace ibex
#endif // __IBEX_CELLSET_H__
