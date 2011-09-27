/***************************************************************************
 *            expression_set.h
 *
 *  Copyright 2011  Pieter Collins
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

/*! \file expression_set.h
 *  \brief Sets defined using expressions over real variables.
 */

#ifndef ARIADNE_EXPRESSION_SET_H
#define ARIADNE_EXPRESSION_SET_H

#include <iostream>

#include "container.h"
#include "variables.h"
#include "expression.h"
#include "assignment.h"
#include "float.h"
#include "real.h"


namespace Ariadne {

class Real;
class RealInterval;

class RealBox;
class RealBoundedConstraintSet;

class Box;
class BoundedConstraintSet;

template<class X> class Variable;
typedef Variable<Real> RealVariable;
template<class X> class Expression;
typedef Expression<Real> RealExpression;
template<class X> class Space;
typedef Space<Real> RealSpace;

class RealVariableInterval;
class RealVariableBox;
class ExpressionConstraintSet;

struct RealVariableLowerInterval {
    Real _lower; RealVariable _variable;
    RealVariableLowerInterval(const Real& l, const RealVariable& v) : _lower(l), _variable(v) { }
    Real lower() const { return _lower; }
    const RealVariable& variable() const { return _variable; }
    operator Expression<Tribool>() const { return ( Real(this->_lower) <= RealExpression(this->_variable) ); };
};

struct RealVariableUpperInterval {
    RealVariable _variable; Real _upper;
    RealVariableUpperInterval(const RealVariable& v, const Real& u) : _variable(v), _upper(u)  { }
    const RealVariable& variable() const { return _variable; }
    Real upper() const { return _upper; }
    operator Expression<Tribool>() const { return ( RealExpression(this->_variable) <= Real(this->_upper) ); }
};


class RealVariableInterval {
  private:
    Real _lower;
    Variable<Real> _variable;
    Real _upper;
  public:
    RealVariableInterval(const Real& l, const Variable<Real>& v, const Real& u)
        : _lower(l), _variable(v), _upper(u) { ARIADNE_ASSERT_MSG(l<=u,"Interval("<<l<<","<<u<<") not provably nonempty"); }
    RealVariableInterval(const RealVariableLowerInterval& lv)
        : _lower(lv._lower), _variable(lv._variable), _upper(infty) { }
    RealVariableInterval(const RealVariableUpperInterval& vu)
        : _lower(-infty), _variable(vu._variable), _upper(vu._upper) { }
    Variable<Real> const& variable() const { return this->_variable; }
    const Interval approximate_interval() const;
    const RealInterval interval() const;
    const Real lower() const { return this->_lower; }
    const Real upper() const { return this->_upper; }
};

OutputStream& operator<<(OutputStream& os, const RealVariableInterval& eivl);

template<class T> template<class XL, class XU> inline RealVariableInterval Variable<T>::in(const XL& l, const XU& u) {
    //ARIADNE_FAIL_MESSAGE("Can't create interval in variable "<<*this<<" of type "<<name<T>()<<"\n");
    ARIADNE_FAIL_MESSAGE("Can't create interval in variable "<<*this<<"\n");
    return RealVariableInterval(l,*this,u);
}

template<> template<class XL, class XU> inline RealVariableInterval Variable<Real>::in(const XL& l, const XU& u) {
    return RealVariableInterval(l,*this,u);
}

inline RealVariableInterval operator<=(const RealVariableLowerInterval& lv, const Real& u) {
    return RealVariableInterval(lv.lower(),lv.variable(),u); }
inline RealVariableInterval operator>=(const Real& u, const RealVariableLowerInterval& lv) {
    return RealVariableInterval(lv.lower(),lv.variable(),u); }
inline RealVariableInterval operator<=(const Real& l, const RealVariableUpperInterval& vu) {
    return RealVariableInterval(l,vu.variable(),vu.upper()); }
inline RealVariableInterval operator>=(const RealVariableUpperInterval& vu, const Real& l) {
    return RealVariableInterval(l,vu.variable(),vu.upper()); }

inline RealVariableLowerInterval operator<=(const Real& l, const RealVariable& v) {
    return RealVariableLowerInterval(l,v); }
inline RealVariableLowerInterval operator>=(const RealVariable& v, const Real& l) {
    return RealVariableLowerInterval(l,v); }
inline RealVariableUpperInterval operator<=(const RealVariable& v, const Real& u) {
    return RealVariableUpperInterval(v,u); }
inline RealVariableUpperInterval operator>=(const Real& u, const RealVariable& v) {
    return RealVariableUpperInterval(v,u); }

inline RealVariableInterval operator==(const RealVariable& v, const Real& x) {
    return RealVariableInterval(x,v,x); }
inline RealVariableInterval operator==(const Real& x, const RealVariable& v) {
    return RealVariableInterval(x,v,x); }

inline RealVariableLowerInterval operator<=(double l, const RealVariable& v) { return Real(l)<=v; }
inline RealVariableLowerInterval operator>=(const RealVariable& v, double l) { return v>=Real(l); }
inline RealVariableUpperInterval operator<=(const RealVariable& v, double u) { return v<=Real(u); }
inline RealVariableUpperInterval operator>=(double u, const RealVariable& v) { return Real(u)>=v; }

inline RealVariableInterval operator==(const RealVariable& v, double x) { return v==Real(x); }
inline RealVariableInterval operator==(double x, const RealVariable& v) { return Real(x)==v; }

class RealVariableBox {
    Map<RealVariable,RealInterval> _bounds;
  public:
    RealVariableBox(const List<RealVariableInterval>& lst);
    RealVariableBox(const Map<RealVariable,RealInterval>& bnds) : _bounds(bnds) { }
    RealVariableBox(const RealSpace& spc, const RealBox& bx);
    Set<RealVariable> variables() const { return _bounds.keys(); }
    Map<RealVariable,RealInterval> bounds() const { return this->_bounds; }
    const RealInterval& operator[](const RealVariable& v) const { return this->_bounds[v]; }
    RealBox box(const List<RealVariable>& spc) const;
    friend RealBox euclidean_set(const RealVariableBox& ebx, const List<RealVariable>& spc);
    friend Box approximate_euclidean_set(const RealVariableBox& set, const RealSpace& space);
    friend OutputStream& operator<<(OutputStream& os, const RealVariableBox& ebx);
};



//! \ingroup GeometryModule ExactSetSubModule
//! \brief A set defined as the image of a box under a continuous function.
//! The set is described as \f$S=h(D) = \{ h(s) \mid s \in D\}\f$ where \f$D\f$ is the domain and \f$h\f$ the function.
class RealExpressionSet
{
    Map<RealVariable,RealInterval> _bounds;
    List<ContinuousPredicate> _constraints;
  public:
    RealExpressionSet(const List<RealVariableInterval>& domain);
    RealExpressionSet(const List<RealVariableInterval>& domain, const List<ContinuousPredicate>& constraints);
    RealExpressionSet(const RealVariableBox& box) : _bounds(box.bounds()) { }
    Map<RealVariable,RealInterval> bounds() const { return this->_bounds; }
    List<ContinuousPredicate> const& constraints() const { return this->_constraints; }
    friend RealBoundedConstraintSet euclidean_set(const RealExpressionSet& set, const RealSpace& space);
    friend BoundedConstraintSet approximate_euclidean_set(const RealExpressionSet& set, const RealSpace& space);
    friend std::ostream& operator<<(std::ostream& os, const RealExpressionSet& eset);
};


}



#endif