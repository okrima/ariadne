/***************************************************************************
 *            python/export_matrix.cc
 *
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

#include <utility>  //for std::pair

#include "numeric/rational.h"
#include "numeric/interval.h"
#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"
#include "linear_algebra/matrix_function.h"

#include "python/name.h"
#include "python/operators.h"
#include "python/float.h"
#include "python/read_matrix.h"

using namespace Ariadne;
using namespace Ariadne::Python;

#include <boost/python.hpp>
#include <boost/python/detail/api_placeholder.hpp>
using namespace boost::python;

template<class T1, class T2>
boost::python::tuple
make_tuple(const std::pair<T1,T2>& p)
{
  return make_tuple(p.first,p.second);
}

template<class X>
std::string
__str__(const Matrix<X>& A)
{
  std::stringstream ss;
  ss << A;
  return ss.str();
}

template<class X>
std::string
__repr__(const Matrix<X>& A)
{
  std::stringstream ss;
  ss << "Matrix(";
  for(size_type i=0; i!=A.number_of_rows(); ++i) {
    ss << (i==0 ? '[' : ',');
    for(size_type j=0; j!=A.number_of_columns(); ++j) {
      ss << (j==0 ? '[' : ',') << A(i,j);
    }
    ss << ']';
  }
  ss << "])";
  return ss.str();
}


template<class X>  
Matrix<X>*
make_matrix(const boost::python::object& obj) 
{
  Matrix<X>* A=new Matrix<X>;
  read_matrix(*A,obj);
  return A;
}

template<class R, class A> 
Matrix<R> 
matrix_exp(const Matrix<A>& mx) 
{
  return exp(mx); 
}

template<class Mx>  
typename Mx::value_type
__getitem__(const Mx& M, boost::python::object obj) 
{
  tuple index=extract<tuple>(obj);
  int i=extract<int>(index[0]);
  int j=extract<int>(index[1]);
  return M(i,j);
}

template<class Mx, class A>  
void __setitem__(Mx& M, boost::python::object obj, const A& x) {
  tuple index=extract<tuple>(obj);
  int i=extract<int>(index[0]);
  int j=extract<int>(index[1]);
  M(i,j)=x;
}


template<class R1,class R2> inline 
Vector<typename traits<R1,R2>::arithmetic_type> 
matrix_vector_solve(const Matrix<R1>& A, const Vector<R2>& v) {
  typedef typename traits<R1,R2>::arithmetic_type F;
  return solve(static_cast<const Matrix<F>&>(A),static_cast<const Vector<F>&>(v));
}

template<class R> inline 
Matrix<typename traits<R>::arithmetic_type> 
matrix_inverse(const Matrix<R>& A) {
  return inverse(A);
}

template<class R> inline 
Matrix<R> 
matrix_transpose(const Matrix<R>& A) {
  return  transpose(A);
}

template<class R> inline 
boost::python::tuple
matrix_qr_approx(const Matrix<R>& A) {
  return make_tuple(qr_approx(A));
}

template<class R>
void export_matrix()
{
  typedef Interval<R> I;
  typedef Vector<R> Vec;
  typedef Vector< Interval<R> > IVec;
  typedef Matrix<R> Mx;
  typedef Matrix< Interval<R> > IMx;

  class_< Matrix<R> >(python_name<R>("Matrix").c_str(),no_init)
    .def("__init__", make_constructor(&make_matrix<R>) )
    .def(init<int,int>())
    .def(init<Mx>())
    //.def(init<std::string>())
    .def("__getitem__",&__getitem__<Mx>)
    .def("__setitem__",&__setitem__<Mx,R>)
    .def("__setitem__",&__setitem__<Mx,double>)
    .def("__neg__",&__neg__<Mx,Mx>)
    .def("__add__",&__add__<IMx,Mx,Mx>)
    .def("__add__",&__add__<IMx,Mx,IMx>)
    .def("__sub__",&__sub__<IMx,Mx,Mx>)
    .def("__sub__",&__sub__<IMx,Mx,IMx>)
    .def("__rmul__",&__mul__<IMx,Mx,int,Mx,R>)
    .def("__rmul__",&__mul__<IMx,Mx,double,Mx,R>)
    .def("__rmul__",&__mul__<IMx,Mx,R>)
    .def("__rmul__",&__mul__<IMx,Mx,I>)
    .def("__mul__",&__mul__<IMx,Mx,int,Mx,R>)
    .def("__mul__",&__mul__<IMx,Mx,double,Mx,R>)
    .def("__mul__",&__mul__<IMx,Mx,R>)
    .def("__mul__",&__mul__<IMx,Mx,I>)
    .def("__mul__",&__mul__<IVec,Mx,Vec>)
    .def("__mul__",&__mul__<IVec,Mx,IVec>)
    .def("__mul__",&__mul__<IMx,Mx,Mx>)
    .def("__mul__",&__mul__<IMx,Mx,IMx>)
    .def("__div__",&__div__<IMx,Mx,int,Mx,R>)
    .def("__div__",&__div__<IMx,Mx,double,Mx,R>)
    .def("__div__",&__div__<IMx,Mx,R>)
    .def("__div__",&__div__<IMx,Mx,I>)
    .def("__str__",&__str__<R>)
    .def("__repr__",&__repr__<R>)
  ;
  def("solve",(Vector<I>(*)(const Matrix<R>&,const Vector<R>&))&solve);
  def("solve",(Vector<I>(*)(const Matrix<R>&,const Vector<I>&))&solve);
  def("transpose",(Matrix<R>(*)(const Matrix<R>&))&transpose);
  def("inverse",(Matrix<I>(*)(const Matrix<R>&))&inverse);

}

 
template<>
void export_matrix<Rational>() 
{
  typedef Rational R;
  typedef Vector<R> Vec;
  typedef Matrix<R> Mx;
  
  class_<Mx>(python_name<R>("Matrix").c_str(),no_init)
    .def("__init__", make_constructor(&make_matrix<R>) )
    .def(init<int,int>())
    .def(init<Mx>())
    //.def(init<std::string>())

    .def("__getitem__",&__getitem__<Mx>)
    .def("__setitem__",&__setitem__<Mx,R>)
    .def("__setitem__",&__setitem__<Mx,double>)

    .def("__neg__",&__neg__<Mx,Mx>)
    .def("__add__",&__add__<Mx,Mx,Mx>)
    .def("__sub__",&__sub__<Mx,Mx,Mx>)
    .def("__rmul__",&__mul__<Mx,Mx,int,Mx,R>)
    .def("__rmul__",&__mul__<Mx,Mx,double,Mx,R>)
    .def("__rmul__",&__mul__<Mx,Mx,R>)
    .def("__mul__",&__mul__<Mx,Mx,int,Mx,R>)
    .def("__mul__",&__mul__<Mx,Mx,double,Mx,R>)
    .def("__mul__",&__mul__<Mx,Mx,R>)
    .def("__mul__",&__mul__<Vec,Mx,Vec>)
    .def("__mul__",&__mul__<Mx,Mx,Mx>)
    .def("__div__",&__div__<Mx,Mx,int,Mx,R>)
    .def("__div__",&__div__<Mx,Mx,double,Mx,R>)
    .def("__div__",&__div__<Mx,Mx,R>)
    .def("__str__",&__str__<R>)
    .def("__repr__",&__repr__<R>)
  ;
  
  def("solve",(Vector<R>(*)(const Matrix<R>&,const Vector<R>&))&solve);
  def("transpose",(Matrix<R>(*)(const Matrix<R>&))&transpose);
  def("inverse",(Matrix<R>(*)(const Matrix<R>&))&inverse);
}

template<class R>
void export_interval_matrix() 
{
  typedef Interval<R> I;
  typedef Vector<R> Vec;
  typedef Matrix<R> Mx;
  typedef Vector< Interval<R> > IVec;
  typedef Matrix< Interval<R> > IMx;
  
  class_<IMx>(python_name<R>("IntervalMatrix").c_str(),no_init)
    .def("__init__", make_constructor(&make_matrix<I>) )
    .def(init<int,int>())
    //.def(init<std::string>())
    .def(init<Mx>())
    .def(init<IMx>())

    .def("__getitem__",&__getitem__<IMx>)
    .def("__setitem__",&__setitem__<IMx,I>)
    .def("__setitem__",&__setitem__<IMx,R>)
    .def("__setitem__",&__setitem__<IMx,double>)

    .def("__neg__",&__neg__<IMx,IMx>)
    .def("__add__",&__add__<IMx,IMx,Mx>)
    .def("__add__",&__add__<IMx,IMx,IMx>)
    .def("__sub__",&__sub__<IMx,IMx,Mx>)
    .def("__sub__",&__sub__<IMx,IMx,IMx>)
    .def("__rmul__",&__mul__<IMx,IMx,int,IMx,R>)
    .def("__rmul__",&__mul__<IMx,IMx,double,IMx,R>)
    .def("__rmul__",&__mul__<IMx,IMx,R>)
    .def("__rmul__",&__mul__<IMx,IMx,I>)
    .def("__mul__",&__mul__<IMx,IMx,int,IMx,R>)
    .def("__mul__",&__mul__<IMx,IMx,double,IMx,R>)
    .def("__mul__",&__mul__<IMx,IMx,R>)
    .def("__mul__",&__mul__<IMx,IMx,I>)
    .def("__mul__",&__mul__<IVec,IMx,Vec>)
    .def("__mul__",&__mul__<IVec,IMx,IVec>)
    .def("__mul__",&__mul__<IMx,IMx,Mx>)
    .def("__mul__",&__mul__<IMx,IMx,IMx>)
    .def("__div__",&__div__<IMx,IMx,int,IMx,R>)
    .def("__div__",&__div__<IMx,IMx,double,IMx,R>)
    .def("__div__",&__div__<IMx,IMx,R>)
    .def("__div__",&__div__<IMx,IMx,I>)
    .def("__str__",&__str__<I>)
    .def("__repr__",&__repr__<I>)
  ;
  
  def("solve",(Vector<I>(*)(const Matrix<R>&,const Vector<I>&))&solve);
  def("solve",(Vector<I>(*)(const Matrix<I>&,const Vector<R>&))&solve);
  def("solve",(Vector<I>(*)(const Matrix<I>&,const Vector<I>&))&solve);
  def("transpose",(Matrix<I>(*)(const Matrix<I>&))&transpose);
  def("inverse",(Matrix<I>(*)(const Matrix<I>&))&inverse);
  def("exp",(Matrix<I>(*)(const Matrix<R>&))&matrix_exp<I,R>);
  def("exp",(Matrix<I>(*)(const Matrix<I>&))&matrix_exp<I,I>);

  def("solve_approx",(Vector<R>(*)(const Matrix<R>&,const Vector<R>&))&solve_approx);
  def("inverse_approx",(Matrix<R>(*)(const Matrix<R>&))&inverse_approx);
  def("qr_approx",&matrix_qr_approx<R>);

  def("midpoint",(Mx(*)(const IMx&))&midpoint);
  def("encloses",(bool(*)(const IMx&,const Mx&))&encloses);
  def("refines",(bool(*)(const IMx&,const IMx&))&refines);

}

template void export_matrix<FloatPy>();
template void export_matrix<Rational>();

template void export_interval_matrix<FloatPy>();
