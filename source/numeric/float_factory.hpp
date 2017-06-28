/***************************************************************************
 *            float_factory.hpp
 *
 *  Copyright 2008-17  Pieter Collins
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

/*! \file float_factory.hpp
 *  \brief Factories for creating floating-point numbers.
 */

#ifndef ARIADNE_FLOAT_FACTORY_HPP
#define ARIADNE_FLOAT_FACTORY_HPP

#include "float.decl.hpp"

namespace Ariadne {

template<class PR> class FloatFactory {
    PR _pr;
  public:
    typedef PR PrecisionType;
    typedef PR PropertiesType;
    FloatFactory(PR const& pr) : _pr(pr) { }
    PR precision() const { return this->_pr; }
    PR properties() const { return this->_pr; }
  public:
    FloatApproximation<PR> create(Number<ApproximateTag> const& y);
    FloatLowerBound<PR> create(Number<ValidatedLowerTag> const& y);
    FloatUpperBound<PR> create(Number<ValidatedUpperTag> const& y);
    FloatBounds<PR> create(Number<ValidatedTag> const& y);
    FloatBounds<PR> create(Number<EffectiveTag> const& y);
    FloatBounds<PR> create(Number<ExactTag> const& y);
    FloatBounds<PR> create(Real const& y);
    FloatBounds<PR> create(Rational const& y);
    FloatBounds<PR> create(Dyadic const& y);
    FloatBounds<PR> create(Integer const& y);
    FloatValue<PR> create(Dyadic const& y, ExactTag);
    FloatValue<PR> create(Integer const& y, ExactTag);
    template<class N, EnableIf<IsBuiltinSignedIntegral<N>> =dummy> FloatValue<PR> create(N const& y);
    template<class M, EnableIf<IsBuiltinUnsignedIntegral<M>> =dummy> PositiveFloatValue<PR> create(M const& y);
    template<class D, EnableIf<IsBuiltinFloatingPoint<D>> =dummy> FloatApproximation<PR> create(D const& y);
};

template<class Y, class PR> using ConcreteType = decltype(declval<FloatFactory<PR>>().create(declval<Y>()));
template<class Y, class PR> inline decltype(auto) make_float(Y const& y, PR pr) { return float_factory(pr).create(y); }

template<class PR> inline FloatFactory<PR> float_factory(PR pr) { return FloatFactory<PR>(pr); }
template<class PR> inline FloatFactory<PR> factory(FloatApproximation<PR> const& flt);
template<class PR> inline FloatFactory<PR> factory(FloatLowerBound<PR> const& flt);
template<class PR> inline FloatFactory<PR> factory(FloatUpperBound<PR> const& flt);
template<class PR> inline FloatFactory<PR> factory(FloatBounds<PR> const& flt);
template<class PR> inline FloatFactory<PR> factory(FloatBall<PR> const& flt);
template<class PR> inline FloatFactory<PR> factory(FloatValue<PR> const& flt);

/*
template<class PR> inline FloatFactory<PR> factory(FloatApproximation<PR> const& flt) { return FloatFactory<PR>(flt.precision()); }
template<class PR> inline FloatFactory<PR> factory(FloatLowerBound<PR> const& flt) { return FloatFactory<PR>(flt.precision()); }
template<class PR> inline FloatFactory<PR> factory(FloatUpperBound<PR> const& flt) { return FloatFactory<PR>(flt.precision()); }
template<class PR> inline FloatFactory<PR> factory(FloatBounds<PR> const& flt) { return FloatFactory<PR>(flt.precision()); }
template<class PR, class PRE> inline FloatFactory<PR> factory(FloatBall<PR,PRE> const& flt) { return FloatFactory<PR>(flt.precision()); }
template<class PR> inline FloatFactory<PR> factory(FloatValue<PR> const& flt) { return FloatFactory<PR>(flt.precision()); }

template<class PR> inline FloatApproximation<PR> FloatFactory<PR>::(Number<ApproximateTag> const& y) { return FloatApproximation<PR>(y,_pr); }
template<class PR> inline FloatLowerBound<PR> FloatFactory<PR>::(Number<ValidatedLowerTag> const& y) { return FloatLowerBound<PR>(y,_pr); }
template<class PR> inline FloatUpperBound<PR> FloatFactory<PR>::(Number<ValidatedUpperTag> const& y) { return FloatUpperBound<PR>(y,_pr); }

template<class PR> inline FloatBounds<PR> FloatFactory<PR>::(Number<ValidatedTag> const& y) { return FloatBounds<PR>(y,_pr); }
template<class PR> inline FloatBounds<PR> FloatFactory<PR>::(Number<EffectiveTag> const& y) { return FloatBounds<PR>(y,_pr); }
template<class PR> inline FloatBounds<PR> FloatFactory<PR>::(Number<ExactTag> const& y) { return FloatBounds<PR>(y,_pr); }
template<class PR> inline FloatBounds<PR> FloatFactory<PR>::(Real const& y) { return FloatBounds<PR>(y,_pr); }
template<class PR> inline FloatBounds<PR> FloatFactory<PR>::(Rational const& y) { return FloatBounds<PR>(y,_pr); }
template<class PR> inline FloatBounds<PR> FloatFactory<PR>::(Dyadic const& y) { return FloatBounds<PR>(y,_pr); }
template<class PR> inline FloatBounds<PR> FloatFactory<PR>::(Integer const& y) { return FloatBounds<PR>(y,_pr); }

template<class PR> inline FloatValue<PR> FloatFactory<PR>::(Dyadic const& y, ExactTag) { return FloatValue<PR>(y,_pr); }
template<class PR> inline FloatValue<PR> FloatFactory<PR>::(Integer const& y, ExactTag) { return FloatValue<PR>(y,_pr); }

template<class PR> inline template<class N, EnableIf<IsBuiltinSignedIntegral<N>> =dummy> FloatValue<PR> FloatFactory<PR>::(N const& y) { return FloatValue<PR>(y,_pr); }
template<class PR> inline template<class M, EnableIf<IsBuiltinUnsignedIntegral<M>> =dummy> PositiveFloatValue<PR> FloatFactory<PR>::(M const& y) { return PositiveFloatValue<PR>(y,_pr); }
template<class PR> inline template<class D, EnableIf<IsBuiltinFloatingPoint<D>> =dummy> FloatApproximation<PR> FloatFactory<PR>::(D const& y) { return FloatApproximation<PR>(RawFloat<PR>(y,_pr)); }
*/

} // namespace Ariadne

#endif
