/***************************************************************************
 *            vector_field_evolver_interface.h
 *
 *  Copyright  2006-7  Alberto Casagrande, Pieter Collins
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
 
/*! \file vector_field_evolver_interface.h
 *  \brief Interface for computing the evolution of sets.
 */

#ifndef ARIADNE_VECTOR_FIELD_EVOLVER_INTERFACE_H
#define ARIADNE_VECTOR_FIELD_EVOLVER_INTERFACE_H

#include <boost/smart_ptr.hpp>

#include "base/types.h"
#include "base/declarations.h"
#include "geometry/declarations.h"
#include "system/declarations.h"
#include "evaluation/declarations.h"


namespace Ariadne {
  namespace Evaluation {

    /*! \brief A class for computing the evolution of a discrete-time autonomous system.
     *  \ingroup Applicators
     */
    template<class R>
    class VectorFieldEvolverInterface {
     public:
      /*! \brief Virtual destructor. */
      virtual ~VectorFieldEvolverInterface() { }

      //@{
      //! \name Evaluation of vector_fields on abstract sets

      /*! \brief Compute an approximation to the set obtained by iterating \a steps times \a vector_field starting in \a initial_set. */
      virtual
      Geometry::SetInterface<R>*
      lower_evolve(const System::VectorField<R>& vector_field, 
                   const Geometry::SetInterface<R>& initial_set,
                   const Numeric::Rational& steps) const = 0;
    
      /*! \brief Compute an approximation to the reachable set of \a vector_field starting in \a initial_set iterating at most \a steps times. */
      virtual
      Geometry::SetInterface<R>*
      lower_reach(const System::VectorField<R>& vector_field, 
                  const Geometry::SetInterface<R>& initial_set,
                  const Numeric::Rational& steps) const = 0;
    
      /*! \brief Compute an approximation to the set obtained by iterating \a steps times \a vector_field starting in \a initial_set. */
      virtual
      Geometry::SetInterface<R>*
      upper_evolve(const System::VectorField<R>& vector_field, 
                   const Geometry::SetInterface<R>& initial_set,
                   const Numeric::Rational& steps) const = 0;
    
      /*! \brief Compute an approximation to the reachable set of \a vector_field starting in \a initial_set iterating at most \a steps times. */
      virtual
      Geometry::SetInterface<R>*
      upper_reach(const System::VectorField<R>& vector_field, 
                  const Geometry::SetInterface<R>& initial_set,
                  const Numeric::Rational& steps) const = 0;
    
      /*! \brief Compute an outer-approximation to the chain-reachable set of \a vector_field starting in \a initial_set. */
      virtual
      Geometry::SetInterface<R>*
      chainreach(const System::VectorField<R>& vector_field, 
                 const Geometry::SetInterface<R>& initial_set) const = 0;
    
   
      /*! \brief Compute an outer-approximation to the viability kernel of \a vector_field within \a bounding_set. */
      virtual
      Geometry::SetInterface<R>* 
      viable(const System::VectorField<R>& vector_field, 
             const Geometry::SetInterface<R>& bounding_set) const = 0;
    
      /*! \brief Attempt to verify that the reachable set of \a vector_field starting in \a initial_set remains in \a safe_set. */
      virtual
      tribool
      verify(const System::VectorField<R>& vector_field, 
             const Geometry::SetInterface<R>& initial_set, 
             const Geometry::SetInterface<R>& safe_set) const = 0;
      //@}


      //@{
      //! \name Methods for computing discretizations of vector_fields on grids

    };



  }
}

 


#endif /* ARIADNE_VECTOR_FIELD_EVOLVER_INTERFACE_H */
