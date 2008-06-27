/***************************************************************************
 *            newton_solver.h
 *
 *  Copyright  2006-8  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */
 
/*! \file newton_solver.h
 *  \brief Newton and interval Newton methods.
 */

#ifndef ARIADNE_NEWTON_SOLVER_H
#define ARIADNE_NEWTON_SOLVER_H

#include <exception>
#include <stdexcept>
#include <string>

#include "numeric/traits.h"
#include "function/declarations.h"
#include "geometry/declarations.h"

#include "solver_interface.h"
#include "solver_base.h"

namespace Ariadne {
  
      
    /*! \ingroup Solvers
     *  \brief Interval Newton solver.
     */
    template<class R>
    class IntervalNewtonSolver : public SolverBase<R>
    {
      typedef typename traits<R>::interval_type I;
     public:
      /*! \brief Constructor. */
      IntervalNewtonSolver(R max_error, uint max_steps) : SolverBase<R>(max_error,max_steps) { }
      
      /*! \brief Solve \f$f(x)=0\f$, using the interval Newton method. */
      Point<I>
      solve(const FunctionInterface<R>& f, 
            const Point<I>& pt); 
    };            
    
    
} // namespace Ariadne

#endif /* ARIADNE_NEWTON_SOLVER_H */