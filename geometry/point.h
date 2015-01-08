/***************************************************************************
 *            point.h
 *
 *  Copyright 2008  Alberto Casagrande, Pieter Collins
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

/*! \file point.h
 *  \brief Points in Euclidean space.
 */

#ifndef ARIADNE_POINT_H
#define ARIADNE_POINT_H

#include "numeric/numeric.h"
#include "algebra/vector.h"

#include "output/graphics_interface.h"

namespace Ariadne {

template<class X> class Point;
typedef Point<ExactNumber> ExactPoint;
typedef Point<EffectiveNumber> EffectivePoint;
typedef Point<ValidatedNumber> ValidatedPoint;
typedef Point<ApproximateNumber> ApproximatePoint;

//! A point in Euclidean space.
template<class X>
class Point
    : public Vector<X>
    , public DrawableInterface
{
  public:
    typedef X RealType;
    //! Default constructor contructs the singleton point in zero dimensions.
    Point() : Vector<RealType>() { }
    //! The origin in \a n dimensions.
    explicit Point(Nat n) : Vector<RealType>(n) { }
    //! Construct from a string literal of the form "(x1,x2,...,xd)".
    explicit Point(const StringType& str);
    Point(const Vector<RealType>& v) : Vector<RealType>(v) { }
    template<class T, EnableIf<IsConvertible<T,X>> =dummy> Point(const Point<T>& pt) : Vector<RealType>(pt.vector()) { }
    //! Construct from an initializer list of floating-point values.
    template<class N, class T> Point(const N& n, const T& t) : Vector<RealType>(n,RealType(t)) { }
    //! Construct from an initializer list of floating-point values.
    explicit Point(InitializerList<double> lst);
    //! The origin in \a n dimensions.
    static Point origin(Nat n) { return Point(n,RealType(0)); }
    //! A dynamically-allocated copy.
    virtual Point<X>* clone() const;
    //! The dimension of the point.
    Nat dimension() const { return this->size(); }
    //! An explicit cast to a float vector. Useful to prevent ambiguous function overloads.
    const Vector<RealType>& vector() const { return *this; }

    Vector<RealType> centre() const { return *this; }

    //! Write to an output stream.
    virtual OutputStream& write(OutputStream& os) const {
        return os << static_cast<const Vector<RealType>&>(*this); }

    virtual Void draw(CanvasInterface& c, const Projection2d& p) const;
    virtual ExactBox bounding_box() const;
};

ExactPoint make_point(const StringType&);

} // namespace Ariadne

#endif // ARIADNE_POINT_H