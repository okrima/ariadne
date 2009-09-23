/***************************************************************************
 *            user_function.h
 *
 *  Copyright 2008-9  Pieter Collins
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

/*! \file user_function.h
 *  \brief Wrappers for user functions
 */

#ifndef ARIADNE_USER_FUNCTION_H
#define ARIADNE_USER_FUNCTION_H

#include <cstdarg>
#include <iosfwd>

#include "function_interface.h"

#include "macros.h"
#include "pointer.h"
#include "container.h"

#include "vector.h"
#include "matrix.h"
#include "polynomial.h"
#include "affine.h"
#include "taylor_model.h"
#include "differential.h"
#include "operators.h"
#include "expression.h"
#include "formula.h"

namespace Ariadne {


template<uint AS, uint PS=0u, uint SM=255u>
struct ScalarFunctionData
{
    const uint argument_size() const { return AS; }
    const uint parameter_size() const { return PS; }
    const uint smoothness() const { return SM; }
};


//! \brief A wrapper for converting templated C++ functions to %Ariadne functions.
//!
//! Given a C++ class T with a (static) template method
//!   <code>template<class R, class A, class P> compute(R& r, const A& a, const P& p);</code>
//! the type <code>UserExpression<T></code> is an Ariadne expression defined by \f$r=f(a)\f$.
//! The constructor for UserExpression<T> takes a Vector<Float> or Vector<Interval> argument which is used for \a p.
//!
//! The class T must also define meta-data <c>argument_size(), parameter_size()
//! and smoothness()</c>. These are most easily defined by inheriting from the
//! <tt>ExpressionData<AS,PS,SM=SMOOTH></tt> class.
//!
//! The constant \a SMOOTH is used for an arbitrarily-differentiable function.
template<class T> class ScalarUserFunction
    : public ScalarFunctionInterface
{
  public:
    ScalarUserFunction() : _p(this->parameter_size()) { }
    ScalarUserFunction(const Vector<Float>& p) : _p(p) { }
    ScalarUserFunction(const Vector<Interval>& p) : _p(p) { }

    const Vector<Interval> parameters() const { return _p; }

    virtual SizeType argument_size() const { return T::argument_size(); }
    virtual SizeType parameter_size() const { return T::parameter_size(); }
    virtual SmoothnessType smoothness() const { return T::smoothness(); }

    virtual Float evaluate(const Vector<Float>& x) const {
        Float r=0; T::compute(r,x,midpoint(_p)); return r; }
    virtual Interval evaluate(const Vector<Interval>& x) const {
        Interval r=0; T::compute(r,x,_p); return r; }

    virtual TaylorModel evaluate(const Vector<TaylorModel>& x) const {
        TaylorModel r(x[0].argument_size(),x[0].accuracy_ptr()); T::compute(r,x,_p); return r; }

    virtual Differential<Float> evaluate(const Vector< Differential<Float> >& x) const {
        Differential<Float> r(x[0].argument_size(),x[0].degree()); T::compute(r,x,midpoint(_p)); return r; }
    virtual Differential<Interval> evaluate(const Vector< Differential<Interval> >& x) const {
        Differential<Interval> r(x[0].argument_size(),x[0].degree()); T::compute(r,x,_p); return r; }

    virtual ScalarUserFunction<T>* derivative(uint j) { ARIADNE_NOT_IMPLEMENTED; }

    virtual Vector<Float> gradient(const Vector<Float>& x) const {
        return this->evaluate(Differential<Float>::variables(1u,x)).gradient(); }
    virtual Vector<Interval> gradient(const Vector<Interval>& x) const {
        return this->evaluate(Differential<Interval>::variables(1u,x)).gradient(); }

    virtual std::ostream& write(std::ostream& os) const  {
        return os << "ScalarUserFunction( argument_size="<<this->argument_size()<<" )"; }
  private:
    Vector<Interval> _p;
};



//! \brief A convenience class defining the meta-data of an %Ariadne function.
template<uint RS, uint AS, uint PS=0u, uint SM=255u>
class VectorFunctionData
{
  public:
    //!
    static const uint result_size() { return RS; }
    //!
    static const uint argument_size() { return AS; }
    //!
    static const uint parameter_size() { return PS; }
    //!
    static const uint smoothness() { return SM; }
};


//! \brief A wrapper for converting templated C++ functions to %Ariadne functions.
//!
//! Given a C++ class T with a (static) template method
//!   <code>template<class R, class A, class P> compute(R& r, const A& a, const P& p);</code>
//! the type <code>Function<T></code> is an Ariadne function defined by \f$r=f(a)\f$.
//! The constructor for Function<T> takes a Vector<Float> argument which is used for \a p.
//!
//! The class T must also define meta-data <c>result_size(), argument_size(), parameter_size()
//! and smoothness()</c> as static functions. These are most easily defined by inheriting from the
//! <tt>ScalarFunctionData<RS,AS,PS,SM=SMOOTH></tt> class.
//!
//! The constant \a SMOOTH is used for an arbitrarily-differentiable function.
template<class T> class VectorUserFunction
    : public VectorFunctionInterface
{
  public:
    VectorUserFunction() : _p(this->parameter_size()) { }
    VectorUserFunction(const Vector<Float>& p) : _p(p) { }
    VectorUserFunction(const Vector<Interval>& p) : _p(p) { }

    const Vector<Interval> parameters() const { return _p; }

    virtual VectorUserFunction<T>* clone() const { return new VectorUserFunction<T>(*this); }

    virtual SizeType result_size() const { return T::result_size(); }
    virtual SizeType argument_size() const { return T::argument_size(); }
    virtual SizeType parameter_size() const { return T::parameter_size(); }
    virtual SmoothnessType smoothness() const { return T::smoothness(); }

    virtual Vector<Float> evaluate(const Vector<Float>& x) const {
        Vector<Float> r(this->result_size(),0.0); T::compute(r,x,midpoint(_p)); return r; }
    virtual Vector<Interval> evaluate(const Vector<Interval>& x) const {
        Vector<Interval> r(this->result_size(),0.0); T::compute(r,x,_p); return r; }

    virtual Vector<TaylorModel> evaluate(const Vector<TaylorModel>& x) const {
        Vector<TaylorModel> r(this->result_size(),TaylorModel(x[0].argument_size(),x[0].accuracy_ptr()));
        T::compute(r,x,_p); return r; }

    virtual Vector< Differential<Float> > evaluate(const Vector< Differential<Float> >& x) const {
        Vector< Differential<Float> > r(this->result_size(),Differential<Float>(x[0].argument_size(),x[0].degree()));
        T::compute(r,x,midpoint(_p)); return r; }
    virtual Vector< Differential<Interval> > evaluate(const Vector< Differential<Interval> >& x) const {
        Vector< Differential<Interval> > r(this->result_size(),Differential<Interval>(x[0].argument_size(),x[0].degree()));
        T::compute(r,x,_p); return r; }

    virtual Matrix<Float> jacobian(const Vector<Float>& x) const {
        return Ariadne::jacobian(this->evaluate(Differential<Float>::variables(1u,x))); }
    virtual Matrix<Interval> jacobian(const Vector<Interval>& x) const {
        return Ariadne::jacobian(this->evaluate(Differential<Interval>::variables(1u,x))); }

    // TODO: Find a better way for writing functions which can handle transformations which may not have a
    // write() method or operator<<.
    virtual std::ostream& write(std::ostream& os) const  {
        return os << "VectorUserFunction( result_size="<<this->result_size()<<", argument_size="<<this->argument_size()<<" )"; }
  private:
    Vector<Interval> _p;
};


} // namespace Ariadne

#endif /* ARIADNE_USER_FUNCTION_H */
