/***************************************************************************
 *            python/operators.h
 *
 *  Copyright  2005-7  Alberto Casagrande, Pieter Collins
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

/*! \file python/operators.h
 *  Commonly used inline methods for the Python interface.
 */
 
#ifndef ARIADNE_PYTHON_OPERATORS_H
#define ARIADNE_PYTHON_OPERATORS_H

#include <cstring>
#include <functional>

#include "numeric/integer.h"
#include "numeric/rational.h"
#include "python/float.h"

namespace Ariadne {
namespace Python {




  template<class Res, class Arg1, class Arg2, class Op>
  inline
  Res
  evaluate(const Arg1& a1, const Arg2& a2)
  {
    return Res(Op()(a1,a2));
  }


  template<class Res, class Arg> inline
  Res floor(const Arg& a) {
    return Res(floor(a));
  }

  template<class Res, class Arg> inline
  Res ceil(const Arg& a) {
    return Res(ceil(a));
  }

  template<class Res, class Arg1, class Arg2> inline
  Res min(const Arg1& a1, const Arg2& a2) {
    return std::min(Res(a1),Res(a2));
  }

  template<class Res, class Arg1, class Arg2> inline
  Res max(const Arg1& a1, const Arg2& a2) {
    return std::max(Res(a1),Res(a2));
  }

  template<class Res, class Arg> inline
  Res abs(const Arg& a) {
    return Res(abs(a));
  }

  template<class Res, class Arg> inline
  Res pos(const Arg& a) {
    return Res(+a);
  }

  template<class Res, class Arg> inline
  Res neg(const Arg& a) {
    return Res(-a);
  }

  template<class Res, class Arg1, class Arg2, class Tmp1, class Tmp2> inline
  Res add(const Arg1& a1, const Arg2& a2) {
    return Res(Tmp1(a1)+Tmp2(a2));
  }

  template<class Res, class Arg1, class Arg2> inline
  Res add(const Arg1& a1, const Arg2& a2) {
    return Res(a1+a2);
  }

  template<class Res, class Arg1, class Arg2, class Tmp1, class Tmp2> inline
  Res radd(const Arg1& a1, const Arg2& a2) {
    return Res(Tmp2(a2)+Tmp1(a1));
  }

  template<class Res, class Arg1, class Arg2> inline
  Res radd(const Arg1& a1, const Arg2& a2) {
    return Res(a2+a1);
  }


  template<class Res, class Arg1, class Arg2, class Tmp1, class Tmp2> inline
  Res sub(const Arg1& a1, const Arg2& a2) {
    return Res(Tmp1(a1)-Tmp2(a2));
  }

  template<class Res, class Arg1, class Arg2> inline
  Res sub(const Arg1& a1, const Arg2& a2) {
    return Res(a1-a2);
  }

  template<class Res, class Arg1, class Arg2, class Tmp1, class Tmp2> inline
  Res rsub(const Arg1& a1, const Arg2& a2) {
    return Res(Tmp2(a2)-Tmp1(a1));
  }

  template<class Res, class Arg1, class Arg2> inline
  Res rsub(const Arg1& a1, const Arg2& a2) {
    return Res(a2-a1);
  }


  template<class Res, class Arg1, class Arg2, class Tmp1, class Tmp2> inline
  Res mul(const Arg1& a1, const Arg2& a2) {
   return Res(Tmp1(a1)*Tmp2(a2));
  }

  template<class Res, class Arg1, class Arg2> inline
  Res mul(const Arg1& a1, const Arg2& a2) {
    return Res(a1*a2);
  }

  template<class Res, class Arg1, class Arg2, class Tmp1, class Tmp2> inline
  Res rmul(const Arg1& a1, const Arg2& a2) {
   return Res(Tmp2(a2)*Tmp1(a1));
  }

  template<class Res, class Arg1, class Arg2> inline
  Res rmul(const Arg1& a1, const Arg2& a2) {
    return Res(a2*a1);
  }


  template<class Res, class Arg1, class Arg2, class Tmp1, class Tmp2> inline
  Res div(const Arg1& a1, const Arg2& a2) {
    return Res(Tmp1(a1)/Tmp2(a2));
  }

  template<class Res, class Arg1, class Arg2> inline
  Res div(const Arg1& a1, const Arg2& a2) {
    return Res(a1/a2);
  }

  template<class Res, class Arg1, class Arg2, class Tmp1, class Tmp2> inline
  Res rdiv(const Arg1& a1, const Arg2& a2) {
    return Res(Tmp2(a2)/Tmp1(a1));
  }

  template<class Res, class Arg1, class Arg2> inline
  Res rdiv(const Arg1& a1, const Arg2& a2) {
    return Res(a2/a1);
  }

  template<class Res, class Arg1, class Arg2, class Tmp1, class Tmp2> inline
  Res pow(const Arg1& a1, const Arg2& a2) {
    return Res(pow(Tmp1(a1),Tmp2(a2)));
  }

  template<class Res, class Arg1, class Arg2> inline
  Res pow(const Arg1& a1, const Arg2& a2) {
    return Res(pow(a1,a2));
  }



  template<class Res,class Arg1,class Arg2> inline
  Res eq(const Arg1& a1, const Arg2& a2) {
    return a1==a2;
  }

  template<class Res,class Arg1,class Arg2> inline
  Res ne(const Arg1& a1, const Arg2& a2) {
    return a1!=a2;
  }

  template<class Res,class Arg1,class Arg2> inline
  Res lt(const Arg1& a1, const Arg2& a2) {
    return a1< a2;
  }

  template<class Res,class Arg1,class Arg2> inline
  Res le(const Arg1& a1, const Arg2& a2) {
    return a1<=a2;
  }

  template<class Res,class Arg1,class Arg2> inline
  Res gt(const Arg1& a1, const Arg2& a2) {
    return a1> a2;
  }

  template<class Res,class Arg1,class Arg2> inline
  Res ge(const Arg1& a1, const Arg2& a2) {
    return a1>=a2;
  }


  template<class Res,class Arg,class Tmp> inline
  Res sqrt(const Arg& a) {
    return sqrt(static_cast<Tmp>(a));
  }


  template<class Res,class Arg,class Tmp> inline
  Res exp(const Arg& a) {
    return exp(static_cast<Tmp>(a));
  }

  template<class Res,class Arg,class Tmp> inline
  Res log(const Arg& a) {
    return log(static_cast<Tmp>(a));
  }

  template<class Res,class Arg,class Tmp> inline
  Res sin(const Arg& a) {
    return sin(static_cast<Tmp>(a));
  }

  template<class Res,class Arg,class Tmp> inline
  Res cos(const Arg& a) {
    return cos(static_cast<Tmp>(a));
  }

  template<class Res,class Arg,class Tmp> inline
  Res tan(const Arg& a) {
    return tan(static_cast<Tmp>(a));
  }

  template<class Res,class Arg,class Tmp> inline
  Res asin(const Arg& a) {
    return asin(static_cast<Tmp>(a));
  }

  template<class Res,class Arg,class Tmp> inline
  Res acos(const Arg& a) {
    return acos(static_cast<Tmp>(a));
  }

  template<class Res,class Arg,class Tmp> inline
  Res atan(const Arg& a) {
    return atan(static_cast<Tmp>(a));
  }

}}

#endif /* ARIADNE_PYTHON_OPERATORS_H */
