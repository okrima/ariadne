/***************************************************************************
 *            orthotope.inline.h
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
 
#include "geometry/rectangle_expression.h"
#include "geometry/rectangle.h"

namespace Ariadne {
  namespace Geometry {
  
    template<class RC,class RG>
    class OrthotopeVerticesIterator 
      : public boost::iterator_facade<OrthotopeVerticesIterator<RC,RG>,
                                      Point<typename Numeric::traits<RC,RG>::arithmetic_type>,
                                      boost::forward_traversal_tag,
                                      Point<typename Numeric::traits<RC,RG>::arithmetic_type> const&,
                                      Point<typename Numeric::traits<RC,RG>::arithmetic_type> const*
                                     >
    {
      friend class Orthotope<RC,RG>;
      typedef typename Numeric::traits<RC,RG>::arithmetic_type F;
      const Orthotope<RC,RG>* _z; long unsigned int _i; bool _parity; Point<F> _pt;
     public:
      OrthotopeVerticesIterator(const Orthotope<RC,RG>& z, bool end) 
        : _z(&z), _i(end ? (1u<<(z.number_of_generators()-1))*3 : 0), _parity(0), _pt(z.centre())
      {
        if(end) { return; }
        for(dimension_type i=0; i!=z.dimension(); ++i) {
          for(dimension_type j=0; j!=z.number_of_generators(); ++j) {
            this->_pt[i]-=z.generators()(i,j); } }
      }
      bool equal(const OrthotopeVerticesIterator<RC,RG>& other) const {
        //std::cerr << "OrthotopeVerticesIterator<RC,RG>::equal" << std::endl;
        return this->_i==other._i && this->_z==other._z; }
      const Point<F>& dereference() const { 
        //std::cerr << "OrthotopeVerticesIterator<RC,RG>::dereference" << std::endl;
        return this->_pt; }
      void increment() { 
        //std::cerr << "OrthotopeVerticesIterator<RC,RG>::increment" << std::endl;
        uint j=0; uint m=1; if(this->_parity) { while(!(m&(this->_i))) { ++j; m*=2u; } ++j; m*=2u; }
        this->_parity=!this->_parity;
        if(j==this->_z->number_of_generators()) { this->_i+=m; return; }
        if(m&(this->_i)) { this->_pt=this->_pt-RC(2)*this->_z->generator(j); this->_i-=m; }
        else { this->_pt=this->_pt+RC(2)*this->_z->generator(j); this->_i+=m; }
      }

      std::ostream& write(std::ostream& os) const;
    };

    template<class RC,class RG> inline 
    std::ostream& operator<<(std::ostream& os, const OrthotopeVerticesIterator<RC,RG>& iter) {
      return iter.write(os);
    }

  } // namespace Geometry
} // namespace Ariadne




namespace Ariadne {
    
    
template<class RC,class RG> inline
Geometry::Orthotope<RC,RG>::Orthotope()
  : _centre(0),_generators(0,0)
{ 
}


template<class RC,class RG> inline
Geometry::Orthotope<RC,RG>::Orthotope(dimension_type d)
  : _centre(d),_generators(d,0)
{ 
}


template<class RC,class RG> inline
Geometry::Orthotope<RC,RG>::Orthotope(dimension_type d, size_type m)
  : _centre(d),_generators(d,m)
{ 
}



template<class RC,class RG> template<class E> inline
Geometry::Orthotope<RC,RG>::Orthotope(const RectangleExpression<E>& re)
  : _centre(re().dimension()),
    _generators(re().dimension(),re().dimension())
{
  (*this)=re;
}




template<class RC,class RG> inline
Geometry::Orthotope<RC,RG>::Orthotope(const Orthotope<RC,RG>& z)
  : _centre(z._centre), 
    _generators(z._generators)
{ 
}


template<class RC,class RG> inline
Geometry::Orthotope<RC,RG>& 
Geometry::Orthotope<RC,RG>::operator=(const Orthotope<RC,RG>& z) {
  if(this != &z) {
    this->_centre = z._centre;
    this->_generators = z._generators;
  }
  return *this;
}


template<class RC,class RG> template<class RC1,class RG1> inline
Geometry::Orthotope<RC,RG>::Orthotope(const Orthotope<RC1,RG1>& z)
  : _centre(z.centre()), 
    _generators(z.generators())
{ 
}


template<class RC,class RG> template<class RC1,class RG1> inline
Geometry::Orthotope<RC,RG>& 
Geometry::Orthotope<RC,RG>::operator=(const Orthotope<RC1,RG1>& z) {
  this->_centre = z.centre();
  this->_generators = z.generators();
  return *this;
}


template<class RC,class RG> template<class RRC, class RRG> inline
Geometry::Orthotope<RC,RG>& 
Geometry::Orthotope<RC,RG>::operator=(const Parallelotope<RRC,RRG>& p) 
{
  this->_centre = p.centre();
  this->_generators = p.generators();
  return *this;
}



template<class RC,class RG> template<class E> inline      
Geometry::Orthotope<RC,RG>& 
Geometry::Orthotope<RC,RG>::operator=(const RectangleExpression<E>& re) 
{
  typedef typename Orthotope<RC,RG>::real_type R;
  
  this->resize(re().dimension(),re().dimension());
  R zero=static_cast<R>(0);
  R two=static_cast<R>(2);
  const E& r=re();
  
  Point<RC>& c=this->_centre;
  LinearAlgebra::Matrix<RG>& g=this->_generators;
  for(size_type i=0; i!=r.dimension(); ++i) {
    c[i]=med_approx(r.lower_bound(i),r.upper_bound(i));
    for(size_type j=0; j!=r.dimension(); ++j) {
      g(i,j)=zero;
    }
    g(i,i)=div_up(sub_up(r.upper_bound(i),r.lower_bound(i)),two);
  }
  return *this;
}



template<class RC,class RG> inline
const Geometry::Point<RC>&
Geometry::Orthotope<RC,RG>::centre() const 
{ 
  return this->_centre;
}


template<class RC,class RG> inline
const RC& 
Geometry::Orthotope<RC,RG>::centre(size_type i) const
{
  return this->_centre[i];
}


template<class RC,class RG> inline
const LinearAlgebra::Matrix<RG>& 
Geometry::Orthotope<RC,RG>::generators() const
{
  return this->_generators;
}

template<class RC,class RG> inline
LinearAlgebra::Vector<RG> 
Geometry::Orthotope<RC,RG>::generator(size_type n) const
{
  return this->_generators.column(n);
}
template<class RC,class RG> inline
const RG& 
Geometry::Orthotope<RC,RG>::generators(size_type i, size_type j) const
{
  return this->_generators(i,j);
  
}  



template<class RC,class RG> inline
void 
Geometry::Orthotope<RC,RG>::resize(dimension_type d, size_type m) 
{ 
  this->_centre.resize(d);
  this->_generators.resize(d,m);
}






template<class RC,class RG> inline
dimension_type 
Geometry::Orthotope<RC,RG>::dimension() const 
{
  return this->_centre.dimension();
}


template<class RC,class RG> inline
size_type 
Geometry::Orthotope<RC,RG>::number_of_generators() const 
{
  return this->_generators.number_of_columns();
}


template<class RC,class RG> inline
tribool 
Geometry::Orthotope<RC,RG>::empty() const 
{ 
  return false; 
}


template<class RC,class RG> inline
tribool 
Geometry::Orthotope<RC,RG>::bounded() const 
{ 
  return true; 
}


template<class RC,class RG> inline
typename Geometry::Orthotope<RC,RG>::R
Geometry::Orthotope<RC,RG>::radius() const 
{
  return this->bounding_box().radius();
}



template<class RC,class RG> inline
Geometry::Orthotope<RC,RG> 
Geometry::operator+(const Orthotope<RC,RG>& z, const LinearAlgebra::Matrix<RG>& A) 
{
  return Orthotope<RC,RG>(z.centre(),concatenate_columns(z.generators(),A));
}



template<class RC,class RG> inline 
std::ostream& 
Geometry::operator<<(std::ostream& os, const Orthotope<RC,RG>& z) 
{
  return z.write(os);
}

template<class RC,class RG> inline
std::istream& 
Geometry::operator>>(std::istream& is, Orthotope<RC,RG>& z) 
{
  return z.read(is);
}










} // namespace Ariadne