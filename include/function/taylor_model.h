/***************************************************************************
 *            taylor_model.h
 *
 *  Copyright  2007 Pieter Collins
 *  pieter.collins@cwi.nl
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
 
/*! \file taylor_model.h
 *  \brief TaylorModels.
 */

#ifndef ARIADNE_TAYLOR_MODEL_H
#define ARIADNE_TAYLOR_MODEL_H

#include <iosfwd>
#include <string>
#include <sstream>

#include "base/types.h"
#include "base/array.h"

#include "linear_algebra/declarations.h"

#include "function/taylor_derivative.h"

namespace Ariadne {
  
   template<class R> class Point; 
   template<class R> class Box; 


  
    class MultiIndex;
    template<class R> class FunctionInterface;
    template<class R> class TaylorVariable;
    template<class R> class TaylorDerivative;

    template<class R> class TaylorModel;

    template<class R> TaylorModel<R> recentre(const TaylorModel<R>&, const Box<R>& bx, const Vector<R>& pt);
    template<class R> TaylorModel<R> truncate(const TaylorModel<R>&, const Box<R>&, size_type, size_type);

    template<class R> TaylorModel<R> compose(const TaylorModel<R>&, const TaylorModel<R>&);
    template<class R> TaylorModel<R> inverse(const TaylorModel<R>&, const Vector<R>&);
    template<class R> TaylorModel<R> implicit(const TaylorModel<R>&, const Vector<R>&);
    template<class R> TaylorModel<R> derivative(const TaylorModel<R>&, size_type);

 


    /*! \brief A taylor_model with multivalued output, using a den.
     *  \ingroup FunctionModel
     */
    template<class R>
    class TaylorModel {
      typedef Interval<R> I;
     public:
      /*! \brief The type of denotable real number used for the corners. */
      typedef R real_type;
     public:
      /*! \brief Default constructor constructs a Taylor model of order zero with no arguments and no result variables. */
      TaylorModel();
      /*! \brief The zero Taylor model in \a as variables with size \a rs image, order \a o and smoothness \a s, defined on the whole space with centre at the origin. */
      TaylorModel(size_type rs, size_type as, smoothness_type o, smoothness_type s);

      /*! \brief Construct from a domain, centre, and two derivative expansions, one for the centre and one over the entire domain. */
      TaylorModel(const Box<R>& domain, const Point<R>& centre, 
                  const TaylorDerivative<I>& centre_derivatives, const TaylorDerivative<I>& domain_derivatives);
      TaylorModel(const Box<R>& domain, const Vector<R>& centre, 
                  const TaylorDerivative<I>& centre_derivatives, const TaylorDerivative<I>& domain_derivatives);
      TaylorModel(const Vector<I>& domain, const Vector<R>& centre, 
                  const TaylorDerivative<I>& centre_derivatives, const TaylorDerivative<I>& domain_derivatives);

      /*! \brief Construct from a domain, centre, an order and a function. */
      TaylorModel(const Box<R>& domain, const Point<R>& centre,
                  smoothness_type order, smoothness_type smoothness,
                  const FunctionInterface<R>& function);
      TaylorModel(const Vector<I>& domain, const Vector<R>& centre,
                  smoothness_type order, smoothness_type smoothness,
                  const FunctionInterface<R>& function);

       
      /*! \brief Equality operator. */
      bool operator==(const TaylorModel<R>& p) const;
      /*! \brief Inequality operator. */
      bool operator!=(const TaylorModel<R>& p) const;

      // Data access
      /*! \brief The data used to define the centre of the Taylor model. */
      const TaylorDerivative<I>& centre_derivatives() const;
      /*! \brief The bounds on the derivative values over the domain of the Taylor model. */
      const TaylorDerivative<I>& domain_derivatives() const;

      // Data access
      /*! \brief The order of the Taylor model. */
      smoothness_type order() const;
      /*! \brief The smoothness of the function. */
      smoothness_type smoothness() const;
      /*! \brief The size of the argument. */
      size_type argument_size() const;
      /*! \brief The size of the result. */
      size_type result_size() const;

      /*! \brief The domain of validity of the Taylor model. */
      Box<R> domain() const;
      /*! \brief The centre of the derivative expansion. */
      Vector<R> centre() const;
      /*! \brief The range of values the Taylor model can take. */
      Box<R> range() const;
     
      /*! \brief Evaluate the Taylor model at the point \a x. */
      Vector<I> evaluate(const Vector<I>& x) const;
      Vector<I> evaluate(const Vector<R>& x) const;
      
      /*! \brief Compute the derivate of the map at a point. */
      Matrix<I> jacobian(const Vector<I>& s) const;
      Matrix<I> jacobian(const Vector<R>& s) const;

      /*! \brief Truncate to a model of lower order and/or smoothness, possibly on a different domain. */
      TaylorModel<R> truncate(const Box<R>& domain, const Vector<R>& centre, 
                              smoothness_type order, smoothness_type smoothness) const;

      /*! \brief The zero Taylor model with result size \a rs and argument size \a as. */
      static TaylorModel<R> zero(size_type rs, size_type as);
      /*! \brief The unit Taylor model with result size 1 and argument size \a as. */
      static TaylorModel<R> one(size_type as);
      /*! \brief The constant Taylor model with result size 1 and argument size \a as. */
      static TaylorModel<R> constant(size_type as, const R& c);

      /*! \brief Write to an output stream. */
      std::ostream& write(std::ostream& os) const;

#ifdef DOXYGEN
      /*! \brief Addition. */
      friend template<class R> TaylorModel<R> operator+(const TaylorModel<R1>&, const TaylorModel<R2>&);
      /*! \brief Subtraction. */
      friend template<class R> TaylorModel<R> operator-(const TaylorModel<R1>&, const TaylorModel<R2>&);
      /*! \brief Multiplication. At least one argument must be scalar-valued. */
      friend template<class R> TaylorModel<R> operator*(const TaylorModel<R1>&, const TaylorModel<R2>&);

      /*! \brief Multiplication by a scalar. */
      friend template<class R> TaylorModel<R> operator*(const R1&, const TaylorModel<R2>&);
      /*! \brief Multiplication by a scalar. */
      friend template<class R> TaylorModel<R> operator*(const TaylorModel<R1>&, const R2&);
      /*! \brief Division by a scalar. */
      friend template<class R> TaylorModel<R> operator/(const TaylorModel<R1>&, const R2&);

      /*! \brief Composition \f$p\circ q(x)=p(q(x))\f$. */
      friend template<class R> TaylorModel<R> compose(const TaylorModel<R1>&, const TaylorModel<R2>&);
      /*! \brief Power of a scalar Taylor model. */
      friend template<class R> TaylorModel<R> pow(const TaylorModel<R>& p, const unsigned int& n);
      /*! \brief Derivative with respect to variable \a k. */
      friend template<class R> TaylorModel<R> derivative(const TaylorModel<R>&, size_type k);
      /*! \brief Truncate within \a r to a Taylor model of order at most \a d, putting the error into terms of order \a s. */
      friend template<class R> TaylorModel<R> truncate(const TaylorModel<R>& p, const Rectangle<R>& bb, size_type d, size_type s);
#endif
     private:
      static void instantiate();
      array< array<I> > _powers(const Vector<I>&) const;
      void _compute_jacobian() const;
      void _set_argument_size(size_type n);
      size_type _compute_maximum_component_size() const;
     private:
      friend TaylorModel<R> recentre<>(const TaylorModel<R>&, const Box<R>& bx, const Vector<R>&);
      friend TaylorModel<R> inverse<>(const TaylorModel<R>&, const Vector<R>&);
      //template<class R> friend void add(TaylorModel<R>&,const TaylorModel<R>&,const TaylorModel<R>&);
      //template<class R> friend void sub(TaylorModel<R>&,const TaylorModel<R>&,const TaylorModel<R>&);
      //template<class R> friend void mul(TaylorModel<R>&,const TaylorModel<R>&,const TaylorModel<R>&);
      //template<class R> friend void div(TaylorModel<R>&,const TaylorModel<R>&,const TaylorModel<R>&);
      //template<class R> friend void compose(TaylorModel<R>&,const TaylorModel<R>&,const TaylorModel<R>&);
      //template<class R> friend void scale(TaylorModel<R>&,const R&);
     private:
      /* Domain of definition. */
      Box<R> _domain;
      /* The centre of the derivative expansion. */
      Vector<R> _centre;
      size_type _smoothness; 
      TaylorDerivative<I> _centre_derivatives;
      TaylorDerivative<I> _domain_derivatives;
   };
    

    template<class R> TaylorModel<R> add(const TaylorModel<R>&, const TaylorModel<R>&);
    template<class R> TaylorModel<R> sub(const TaylorModel<R>&, const TaylorModel<R>&);

    template<class R> TaylorModel<R> combine(const TaylorModel<R>&, const TaylorModel<R>&);
    template<class R> TaylorModel<R> join(const TaylorModel<R>&, const TaylorModel<R>&);


    template<class R> std::ostream& operator<<(std::ostream&, const TaylorModel<R>&);
  


} // namespace Ariadne


#include "taylor_model.inline.h"

#endif /* ARIADNE_TAYLOR_MODEL_H */
