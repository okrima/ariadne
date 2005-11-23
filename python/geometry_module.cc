/***************************************************************************
 *            python/geometry_module.cc
 *
 *  21 October 2005
 *  Copyright  2005  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, Pieter.Collins@cwi.nl
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

#include "numerical_type.h"
#include "interval.h"
#include "state.h"
#include "rectangle.h"
#include "list_set.h"
#include "grid_set.h"
#include "binary_word.h"
#include "binary_tree.h"
#include "partition_tree_set.h"

#include <boost/python.hpp>

#include "real_typedef.h"
#include "python_utilities.h"

using Ariadne::BooleanArray;
using Ariadne::BinaryWord;
using Ariadne::BinaryTree;
using Ariadne::BinaryWord;
using Ariadne::SubdivisionSequence;

void export_state();
void export_interval();
void export_rectangle();
void export_polyhedron();
void export_list_set();
void export_grid();
void export_partition_tree();

typedef Ariadne::Interval<Real> RInterval;
typedef Ariadne::Geometry::Rectangle<Real> RRectangle;

typedef Ariadne::Geometry::ListSet<Real,Ariadne::Geometry::Rectangle> RRectangleListSet;

typedef Ariadne::Geometry::FiniteGrid<Real> RFiniteGrid;
typedef Ariadne::Geometry::GridRectangleListSet<Real> RGridRectangleListSet;
typedef Ariadne::Geometry::GridMaskSet<Real> RGridMaskSet;

typedef Ariadne::Geometry::PartitionScheme<Real> RPartitionScheme;
typedef Ariadne::Geometry::PartitionTree<Real> RPartitionTree;
typedef Ariadne::Geometry::PartitionTreeSet<Real> RPartitionTreeSet;

using Ariadne::Geometry::regular_intersection;
using Ariadne::Geometry::interiors_intersect;
using Ariadne::Geometry::disjoint;
using Ariadne::Geometry::inner_subset;
using Ariadne::Geometry::subset;

using boost::python::class_;
using boost::python::init;
using boost::python::self;
using boost::python::def;
using boost::python::self_ns::str;
using boost::python::return_value_policy;
using boost::python::copy_const_reference;
  
BOOST_PYTHON_MODULE(geometry)
{
  class_<BooleanArray>("BooleanArray")
    .def(init<int>())
    .def("__len__", &BooleanArray::size)
    .def("__getitem__", &get<BooleanArray>)
    .def("__setitem__", &set<BooleanArray>)
    .def(str(self))    // __str__
  ;

  class_<BinaryTree>("BinaryTree")
    .def(init<BooleanArray>())
    .def("__len__", &BinaryTree::size)
    .def(str(self))    // __str__
  ;

  class_<BinaryWord>("BinaryWord")
    .def(init<BinaryWord>())
    .def("__len__", &BinaryWord::size)
    .def("__getitem__", &BinaryWord::get)
    .def(str(self))    // __str__
  ;

  class_<SubdivisionSequence>("SubdivisionSequence")
    .def("__len__", &BinaryTree::size)
    .def("__getitem__", &SubdivisionSequence::get)
    .def(str(self))    // __str__
  ;

  export_state();
  export_interval();
  export_rectangle();
  export_polyhedron();
  export_list_set();
  export_grid();
  export_partition_tree();
}
