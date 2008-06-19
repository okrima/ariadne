/***************************************************************************
 *            level_set.code.h
 *
 *  Copyright  2007  Alberto Casagrande, Pieter Collins
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
 
#include "level_set.h"
#include "linear_algebra/vector.h"
#include "geometry/euclidean_space.h"
#include "geometry/point.h"
#include "geometry/box.h"

namespace {

using namespace Ariadne;

template<class X1, class X2> inline
tribool less(const Vector<X1>& v1, const Vector<X2>& v2)
{
  tribool result=true;
  for(size_type i=0; i!=v1.size(); ++i) {
    result=result && (v1[0]<v2[0]);
    if(!result) { return false; }
  }
  return result;
}
  
}


namespace Ariadne {
    
template<class R>
LevelSet<R>::LevelSet(const FunctionInterface<R>& f)
  : _function_ptr(f.clone()) 
{ 
}


template<class R>
LevelSet<R>::~LevelSet() 
{ 
}


template<class R>
LevelSet<R>* 
LevelSet<R>::clone() const 
{
  return new LevelSet<R>(*this->_function_ptr); 
}

template<class R>
dimension_type 
LevelSet<R>::dimension() const 
{
  return this->_function_ptr->argument_size(); 
}

template<class R>
EuclideanSpace
LevelSet<R>::space() const 
{
  return EuclideanSpace(this->dimension());
}


template<class R>
tribool 
LevelSet<R>::separates(const Point<A>& pt1, const Point<A>& pt2) const 
{
  if(this->number_of_constraints()!=1) {
    return false;
  }
  const FunctionInterface<R>& f=*this->_function_ptr;
  A z=0;
  Vector<A> v1=pt1.position_vector();
  A fv1=f(pt1.position_vector())[0];
  A fv2=f(pt2.position_vector())[0];
  tribool less1=z<fv1;
  tribool less2=z<fv2;
  return less1 xor less2;
}



template<class R>
tribool 
LevelSet<R>::contains(const Point<R>& pt) const 
{
  const FunctionInterface<R>& f=*this->_function_ptr;
  Vector<A> v=pt.position_vector();
  Vector<R> z=Vector<R>::zero(this->number_of_constraints());
  Vector<A> fv=f(v);
  return !(::less(z,fv) || ::less(fv,z));
}


template<class R>
tribool 
LevelSet<R>::superset(const Box<R>& r) const 
{
  return false;
}


template<class R>
tribool 
LevelSet<R>::intersects(const Box<R>& r) const 
{
  return !this->disjoint(r);
}


template<class R>
tribool 
LevelSet<R>::disjoint(const Box<R>& r) const 
{
  const FunctionInterface<R>& f=*this->_function_ptr;
  Vector<A> v=r.position_vectors();
  Vector<R> z=Vector<R>::zero(this->number_of_constraints());
  Vector<A> fv=f(v);
  return ::less(fv,z) || ::less(z,fv);
}

template<class R>
tribool 
LevelSet<R>::subset(const Box<R>& r) const 
{
  return indeterminate;
}

template<class R>
tribool 
LevelSet<R>::bounded() const 
{
  return indeterminate;
}      

template<class R>
Box<R> 
LevelSet<R>::bounding_box() const 
{
  throw UnboundedSet("LevelSet::bounding_box(): cannot be computed in general case");
}      

template<class R>
std::ostream& 
LevelSet<R>::write(std::ostream& os) const 
{
  return os << "LevelSet( function=" << *this->_function_ptr << " )";
}



template<class R>
size_type
LevelSet<R>::number_of_constraints() const
{
  return this->_function_ptr->result_size();
}


template<class R>
const FunctionInterface<R>&
LevelSet<R>::function() const 
{
  return *this->_function_ptr;
}


template<class R>
Point<typename LevelSet<R>::A> 
LevelSet<R>::function(const Point<A>& pt) const 
{
  const FunctionInterface<R>& f=*this->_function_ptr;
  Vector<A> v=pt.position_vector();
  return Point<A>(f(v));
}

} // namespace Ariadne
