/***************************************************************************
 *            python/export_point.cc
 *
 *  21 October 2005
 *  Copyright  2005-6  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, Pieter.Collins@cwi.nl
 ****************************************************************************/

/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is diself_ns::stributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */


#include "geometry/point.h"
#include "geometry/point_list.h"
#include "geometry/rectangle.h"
#include "linear_algebra/vector.h"

#include <boost/python.hpp>

#include "python/typedefs.h"
using namespace Ariadne;

#include <boost/python.hpp>
using namespace boost::python;

inline void rpoint_setitem_from_double(RPoint& p, uint i, double x) {
  p.set(i,Ariadne::convert_to<Real>(x));
}

inline RPoint rpoint_add_rVector(const RPoint& p, const RVector& v) {
  return p+v;
}

inline RVector rpoint_sub_rpoint(const RPoint& p, const RPoint& q) {
  return p-q;
}

inline RRectangle rectangle_expanded(const RPoint& p, const Real& r) {
  return (((RRectangle)p).expand_by(r));
}

inline RPoint point_list_get(const RPointList& pl, const size_type& n) {
  return pl[n];
}

void export_point() {
  class_<RPoint>("Point",init<int>())
    .def(init<RPoint>())
    .def(init<std::string>())
    .def("dimension", &RPoint::dimension)
    .def("__len__", &RPoint::dimension)
    .def("__getitem__", &RPoint::get)
    .def("__setitem__", &RPoint::set)
    .def("__setitem__", &rpoint_setitem_from_double)
    .def("__eq__", &RPoint::operator==)
    .def("__ne__", &RPoint::operator!=)
    .def("__add__", &rpoint_add_rVector)
    .def("__sub__", &rpoint_sub_rpoint)
    .def("rectangle_expanded_by", &rectangle_expanded)
    .def(self_ns::str(self))
  ;
}

void export_point_list() {
  class_<RPointList>("PointList",init<>())
    .def("size", &RPointList::size)
    .def("append", &RPointList::push_back)
    .def("__getitem__", &point_list_get)
  ;
}
