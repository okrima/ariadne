/***************************************************************************
 *            python/export_orbiter.cc
 *
 *  Copyright  2007  Alberto Casagrande, Pieter Collins
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

#include "python/float.h"

#include "geometry/box.h"
#include "geometry/box_list_set.h"
#include "geometry/rectangle.h"
#include "geometry/zonotope.h"
#include "geometry/orbit.h"
#include "system/map.h"
#include "evaluation/evolution_parameters.h"
#include "evaluation/standard_applicator.h"
#include "evaluation/standard_approximator.h"
#include "evaluation/map_orbiter.h"

using namespace Ariadne;

using namespace Ariadne::Geometry;
using namespace Ariadne::System;
using namespace Ariadne::Evaluation;
using namespace Ariadne::Python;

#include <boost/python.hpp>
using namespace boost::python;

template<class BS>
MapOrbiter<BS> make_orbiter(double mbsr) {
  typedef typename BS::real_type R;
  EvolutionParameters<R> parameters;
  parameters.set_maximum_basic_set_radius(mbsr);
  StandardApplicator<R> applicator;
  StandardApproximator<BS> approximator;
  return MapOrbiter<BS>(parameters,applicator,approximator);
}

template<class T1, class T2>
boost::python::tuple
make_tuple(const std::pair<T1,T2>& p)
{
  return make_tuple(p.first,p.second);
}

template<class BS, class R>
boost::python::tuple
lower_evolve(const MapOrbiter<BS>& orb, const Map<R>& f, const Box<R>& bx, const Integer& n) {
  return make_tuple(orb.lower_evolve(f,bx,n));
}


template<class R>
void export_orbiter() 
{
  typedef Integer N;
  typedef Zonotope<R,ExactTag> ZBS;
  typedef Zonotope<R,UniformErrorTag> EZBS;
  typedef Map<R> F;

  def("zonotope_orbiter",&make_orbiter<ZBS>);
  def("parallelotope_orbiter",&make_orbiter<EZBS>);

  class_< MapOrbiter<ZBS> > orbiter_class("ZonotopeOrbiter",no_init);
  orbiter_class.def("lower_evolve",&lower_evolve<ZBS,R>);
  orbiter_class.def("lower_reach",&MapOrbiter<ZBS>::lower_reach);
  orbiter_class.def("upper_evolve",&MapOrbiter<ZBS>::upper_evolve);
  orbiter_class.def("upper_reach",&MapOrbiter<ZBS>::upper_reach);
  orbiter_class.def("orbit",(Orbit<N,ZBS>*(MapOrbiter<ZBS>::*)(const F&,const Box<R>&,const N&)const) &MapOrbiter<ZBS>::orbit,
                    return_value_policy<manage_new_object>());
  orbiter_class.def(self_ns::str(self));

  {
    class_< MapOrbiter<EZBS> > orbiter_class("ParallelotopeOrbiter",no_init);
    orbiter_class.def("lower_evolve",&lower_evolve<EZBS,R>);
    orbiter_class.def("lower_reach",&MapOrbiter<EZBS>::lower_reach);
    orbiter_class.def("upper_evolve",&MapOrbiter<EZBS>::upper_evolve);
    orbiter_class.def("upper_reach",&MapOrbiter<EZBS>::upper_reach);
    orbiter_class.def("orbit",(Orbit<N,EZBS>*(MapOrbiter<EZBS>::*)(const F&,const Box<R>&,const N&)const) &MapOrbiter<EZBS>::orbit,
                      return_value_policy<manage_new_object>());
  }

}


template void export_orbiter<FloatPy>();
