/***************************************************************************
 *            integrator_interface.h
 *
 *  Copyright  2006-7  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, pieter.collins@cwi.nl
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
 
/*! \file integrator_interface.h
 *  \brief Class for computing the image of a basic set under a map.
 */

#ifndef ARIADNE_INTEGRATOR_INTERFACE_H
#define ARIADNE_INTEGRATOR_INTERFACE_H

#include <boost/shared_ptr.hpp>

#include "base/types.h"
#include "base/declarations.h"
#include "linear_algebra/declarations.h"
#include "geometry/declarations.h"
#include "system/declarations.h"

namespace Ariadne {
  namespace Evaluation {


    /*! \brief A class for computing the image of a basic set under a map. 
     *  \ingroup Integrators
     */
    template<class BS>
    class IntegratorInterface
    {
      typedef Numeric::Rational T;
      typedef typename BS::real_type R;
      typedef Numeric::Interval<R> I;
     public:
      //@{ 
      //! \name Constructors and cloning operations.
      /*! \brief Virtual destructor. */
      virtual ~IntegratorInterface() { }
      /*! \brief Make a dynamically-allocated copy. */
      virtual IntegratorInterface<BS>* clone() const = 0;
      //@}


      //@{ 
      //! \name Methods for applying a system to a basic set.

      /*! \brief Compute an integration time and a bounding box. */
      virtual 
      std::pair< Numeric::Rational, Geometry::Box<R> >
      flow_bounds(const System::VectorField<R>& f, 
                  const Geometry::Box<R>& s,
                  const Numeric::Rational& t) const = 0; 
      
      /*! \brief Compute the evolution of a basic set for times in the range [t1,t2]. */
      virtual 
      BS
      evolution_step(const System::VectorField<R>& f, 
                     const BS& s,
                     const Numeric::Interval<R>& t1, 
                     const Numeric::Interval<R>& t2, 
                     const Geometry::Box<R>& bb) const
      {
        ARIADNE_ASSERT(t2>=t1);
        BS es=this->integration_step(f,s,t1,bb); 
        if(t1==t2) { return es; }
        else { return this->reachability_step(f,es,Numeric::Interval<R>(t2-t1),bb); }
      }
          
      /*! \brief Compute the image of a basic set under a continuous function. */
      virtual 
      BS
      integration_step(const System::VectorField<R>& f, 
                       const BS& s,
                       const Numeric::Interval<R>& t, 
                       const Geometry::Box<R>& bb) const = 0; 
      
      /*! \brief Compute the image of a basic set under a continuous function. */
      virtual 
      BS
      reachability_step(const System::VectorField<R>& f, 
                        const BS& s,
                        const Numeric::Interval<R>& t, 
                        const Geometry::Box<R>& bb) const = 0;
      //@}

      //@{ 
      //! \name Input/output operators. */
      /*! \brief Write to an output stream. */
      virtual std::ostream& write(std::ostream&) const = 0;
      //@}
    };

    template<class BS> std::ostream& operator<<(std::ostream& os, const IntegratorInterface<BS>& i);




    /*! \brief A class for computing the flow of a vector field.
     *  \ingroup Integrators
     */
    template<class R>
    class FlowerInterface
    {
      typedef Numeric::Interval<R> I;
     public:
      /*! \brief Make a dynamically-allocated copy. */
      virtual FlowerInterface<R>* clone() const = 0;

      /*! \brief Compute the flow of a point. */
      virtual 
      Geometry::Point<I> 
      flow_step(const System::VectorField<R>& f, 
                const Geometry::Point<I>& s, 
                const Numeric::Interval<R>& t, 
                const Geometry::Box<R>& bb) const = 0;
    };


    /*! \brief A class for computing the flow of a differentiable vector field.
     *  \ingroup Integrators
     */
    template<class R>
    class DifferentiableFlowerInterface
      : public FlowerInterface<R>
    {
      typedef Numeric::Interval<R> I;
     public:
      /*! \brief Make a dynamically-allocated copy. */
      virtual DifferentiableFlowerInterface<R>* clone() const = 0;
      /*! \brief Compute the spacial jacobian over a flow step of time \a t starting at \a p assuming that the flow remains within \a bb. */
      virtual LinearAlgebra::Matrix<I> flow_step_jacobian(const System::VectorField<R>& vf,
                                                          const Geometry::Point<I>& p,
                                                          const Numeric::Interval<R>& t,
                                                          const Geometry::Box<R>& bb) const = 0;
    };

    template<class R> std::ostream& operator<<(std::ostream& os, const FlowerInterface<R>& i);


  }
}

#include "integrator_interface.inline.h"

#endif /* ARIADNE_INTEGRATOR_INTERFACE_H */
