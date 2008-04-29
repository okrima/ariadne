/***************************************************************************
 *            affine_integrator.code.h
 *
 *  Copyright  2006-7  Alberto Casagrande, Pieter Collins
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
 
//#define DEBUG

#include <iosfwd>
#include <string>
#include <sstream>
#include <algorithm>
#include <typeinfo>
#include <cassert>

#include <list>
#include <set>
#include <vector>
#include <valarray>

#include "affine_integrator.h"

#include "base/array.h"

#include "numeric/interval.h"

#include "linear_algebra/vector.h"

#include "linear_algebra/matrix.h"
#include "linear_algebra/matrix_function.h"

#include "function/affine_model.h"

#include "geometry/box.h"
#include "geometry/zonotope.h"

#include "system/exceptions.h"
#include "system/vector_field.h"
#include "system/affine_vector_field.h"

#include "output/logging.h"

namespace Ariadne {

namespace Evaluation { static int& verbosity = integrator_verbosity; }

using Numeric::Interval;
using Geometry::Point;
using Geometry::Box;
using Geometry::Zonotope;
using Geometry::UniformErrorTag;
using Geometry::IntervalTag;

template<class R> 
R
Evaluation::gexp_up(const R& x, uint k)
{
  using namespace Numeric;
  
  R result=div_up(static_cast<R>(1),R(Numeric::fac<Integer>(k)));
  uint n=k;
  R term=result;
  while(add_down(result,term)==result) {
    n=n+1;
    term=div_up(mul_up(result,term),static_cast<R>(n));
    result=add_up(result,term);
  }
  return result;
}

template<class R> 
LinearAlgebra::Vector< Numeric::Interval<R> > 
Evaluation::gexp(const LinearAlgebra::Matrix< Numeric::Interval<R> >& A, 
                 const LinearAlgebra::Vector< Numeric::Interval<R> >& b, 
                 const Numeric::Interval<R>& t, 
                 const uint& k)
{
  using namespace Numeric;
  using namespace LinearAlgebra;
  
  //FIXME: Make number of steps depend on precision
  static const int MAX_STEPS=12;
  typedef Interval<R> I;
  
  int ns=MAX_STEPS;
  
  Matrix<R> id=LinearAlgebra::Matrix<R>::identity(A.number_of_rows());
  
  Vector<I> result=b/R(Numeric::fac<Integer>(k));
  Vector<I> term=result;
  
  for(uint n=k+1; n!=k+ns; ++n) {
    term=(A*term)*(t/static_cast<R>(n));
    result=result+term;
  }
  
  R err=mul_up(pow_up(I(norm(A)*t).upper(),ns),gexp_up(I(norm(A)*t).upper(),k+ns));
  I ierr=err*I(-1,1);
  result+=Vector<I>(result.size(),ierr);
  
  if(Evaluation::verbosity>7) { 
    std::clog << "A=" << A << ",  t=" << t << ",  k=" << k << "\n" 
              << "gexp(A,t,k)=" << result << ", err=" << err << std::endl; 
  }
  
  return result;
}

template<class R> 
LinearAlgebra::Matrix< Numeric::Interval<R> > 
Evaluation::gexp(
                 const LinearAlgebra::Matrix< Numeric::Interval<R> >& A, 
                 const Numeric::Interval<R>& t, 
                 const uint& k)
{
  using namespace Numeric;
  using namespace LinearAlgebra;
  
  //FIXME: Make number of steps depend on precision
  static const int MAX_STEPS=12;
  typedef Interval<R> I;
  
  int ns=MAX_STEPS;
  
  Matrix<R> id=LinearAlgebra::Matrix<R>::identity(A.number_of_rows());
  Matrix<I> result=id/static_cast<R>(Numeric::fac<Integer>(k));
  Matrix<I> term=result;
  
  for(uint n=k+1; n!=k+ns; ++n) {
    term=(A*term)*(t/static_cast<R>(n));
    result=result+term;
  }
  
  R err=mul_up(pow_up(I(norm(A)*t).upper(),ns),gexp_up(I(norm(A)*t).upper(),k+ns));
  I ierr=err*I(-1,1);
  result+=Matrix<I>(result.number_of_rows(),result.number_of_columns(),&ierr,0,0);
  
  if(Evaluation::verbosity>7) { 
    std::clog << "A=" << A << ",  t=" << t << ",  k=" << k << "\n" 
              << "gexp(A,t,k)=" << result << ", err=" << err << std::endl; 
  }
  
  return result;
}


template<class R>
Evaluation::AffineIntegrator< Geometry::Zonotope<R> >::
AffineIntegrator()
{
}



template<class R>
Evaluation::AffineIntegrator< Geometry::Zonotope<R> >*
Evaluation::AffineIntegrator< Geometry::Zonotope<R> >::
clone() const
{
  return new AffineIntegrator< Zonotope<R> >();
}


template<class R>
const System::AffineVectorField<R>*
Evaluation::AffineIntegrator< Geometry::Zonotope<R> >::
cast(const System::VectorField<R>* vf) const
{
 if(dynamic_cast<Function::AffineFunction<R>*>(&vf->function())) {
    return static_cast<const System::AffineVectorField<R>*>(vf);
  } else {
    return 0;
  }
}


template<class R>
std::pair< Numeric::Rational, Geometry::Box<R> >
Evaluation::AffineIntegrator< Geometry::Zonotope<R> >::
flow_bounds(const System::VectorField<R>& f, 
            const Geometry::Zonotope<R>& z,
            const Numeric::Rational& t) const
{
  // Since we don't need a bounding box for an affine integrator, we can just return the entire space.
  return std::make_pair(t,Geometry::Box<R>::entire_space(z.dimension())); 
}


template<class R> 
Geometry::Point< Numeric::Interval<R> > 
Evaluation::AffineIntegrator< Geometry::Zonotope<R> >::
flow_step(const System::VectorField<R>& vf,
          const Geometry::Point<I>& p,
          const Numeric::Rational& h,
          const Geometry::Box<R>& bb) const
{
  const System::AffineVectorField<R>* avf=this->cast(&vf);
  if(!avf) {
    ARIADNE_THROW(std::runtime_error,"AffineIntegrator::flow_step(...)","vector_field is not affine.");
  }
  return this->flow_step(*avf,p,h);
}
     
 
template<class R> 
LinearAlgebra::Matrix< Numeric::Interval<R> > 
Evaluation::AffineIntegrator< Geometry::Zonotope<R> >::
flow_step_jacobian(const System::VectorField<R>& vf,
                   const Point<I>& p,
                   const Numeric::Rational& h,
                   const Box<R>& bb) const
{
  const System::AffineVectorField<R>* avf=this->cast(&vf);
  if(!avf) {
    ARIADNE_THROW(std::runtime_error,"AffineIntegrator::flow_step_jacobian(...)","vector_field is not affine.");
  }
  return this->flow_step_jacobian(*avf,p,h);
}

     

template<class R>
Geometry::Zonotope<R>
Evaluation::AffineIntegrator< Geometry::Zonotope<R> >::
integration_step(const System::AffineVectorField<R>& affine_vector_field, 
                 const Geometry::Zonotope<R>& initial_set, 
                 const Numeric::Rational& step_size) const
{
  ARIADNE_LOG(6,"AffineIntegrator::integration_step(AffineVectorField,Zonotope<Interval>,Time) const\n");
  
  using namespace LinearAlgebra;
  
  const System::AffineVectorField<R>& vf=affine_vector_field;
  const Geometry::Zonotope<R>& z=initial_set;
  I h=step_size;
  
  ARIADNE_LOG(9,"vf="<<vf<<"\n");
  ARIADNE_LOG(7,"z="<<z<<"\n");
  ARIADNE_LOG(7,"h="<<h<<"\n");
  ARIADNE_LOG(7,"vf(c)="<<vf(z.centre())<<", h="<<h<<"\n");
  
  // Use the formula x(t) = x0 + h * P * ( A * x0 + b ) 
  //                      = D * x0 + h * P * b
  // where P = gexp(A*h,1) = sum (Ah)^n / (n+1)!
  const Matrix<I>& A=vf.A();
  const Vector<I>& b=vf.b();
  
  Matrix<I> iP=gexp(A,h,1);
  Matrix<I> iDf=iP*(h*A)+Matrix<R>::identity(vf.dimension());

  ARIADNE_LOG(9,"approximating derivative=" << iDf << "\n");
  ARIADNE_LOG(9,"approximating twist=" << iP << "\n");
  
  Vector<I> iv1=(iDf*z.centre().position_vector());
  ARIADNE_LOG(9,"iv1="<<iv1<<"\n");
  Vector<I> iv2=h*(iP*b);
  ARIADNE_LOG(9,"iv2="<<iv2<<"\n");
  Vector<I> icv=iv1+iv2;
  Point<I> ic(icv);
  ARIADNE_LOG(9,"ic="<<icv<<"\n");
  
  Function::AffineModel<R> model(initial_set.bounding_box(),initial_set.centre(),ic,iDf);
  return Geometry::apply(model,initial_set);
}




template<class R>
Geometry::Zonotope<R> 
Evaluation::AffineIntegrator< Geometry::Zonotope<R> >::
reachability_step(const System::AffineVectorField<R>& vector_field, 
                  const Geometry::Zonotope<R>& initial_set, 
                  const Numeric::Rational& step_size) const
{
  using namespace Numeric;
  using namespace LinearAlgebra;
  using namespace Geometry;
  using namespace System;
  
  ARIADNE_LOG(6,"AffineIntegrator::reachability_step(AffineVectorField,Zonotope<Interval>,Time)\n");
  
  ARIADNE_CHECK_EQUAL_DIMENSIONS(vector_field,initial_set,"AffineIntegrator::reachability_step(AffineVectorField,Zonotope<Interval>,Time)");
  
  const AffineVectorField<R>& avf(vector_field);
  Zonotope<R> ez=initial_set;
  const size_type n=vector_field.dimension();
  const Matrix<R> id=Matrix<R>::identity(n);
  const Rational hc=step_size/2;
  // FIXME: There should be no need to convert to time_type / Rational
  Rational hh(step_size/2);
  
  const Matrix<I>& A=avf.A();
  const Vector<I>& b=avf.b();
  
  // No change of step size for affine integrator
  ez=this->integration_step(avf,ez,hc);
  
  /* Use centre c, generators G and t*(Ac+b), and error (exp(At)-I)Ge+inv(A)(exp(At)-I-At)(Ac+b)
   * 
   */
  const I ihh=hh;
  const Point<I>& c=ez.centre();
  const Matrix<R>& G=ez.generators();
  Vector<I> Acpb=A*c.position_vector()+b;
  //Acpb=Acpb+b;
  R nrmA=norm(A).upper();
  R nrmAh=I(nrmA*ihh).upper();
  R err=I(gexp_up(nrmAh,1)*nrmA*ihh*norm(G)+gexp_up(nrmAh,2)*ihh*ihh*nrmA*norm(Acpb).upper()).upper();
  Vector<I> errv=Vector<I>(avf.dimension(),Interval<R>(-err,err));
  Vector<I> itg=ihh*Acpb;
  Vector<R> tg=midpoint(itg);
  Point<I> nc=c+errv+(itg-tg);
  return Zonotope<R>(nc,concatenate_columns(G,tg));
}



template<class R>
Geometry::Zonotope<R> 
Evaluation::AffineIntegrator< Geometry::Zonotope<R> >::
integration_step(const System::VectorField<R>& vector_field, 
                 const Geometry::Zonotope<R>& initial_set, 
                 const Numeric::Rational& step_size, 
                 const Geometry::Box<R>& bounding_set) const
{
  const System::AffineVectorField<R>* affine_vector_field=this->cast(&vector_field);
  if(!affine_vector_field) {
    ARIADNE_THROW(std::runtime_error,"AffineIntegrator::integration_step(...)","vector_field is not affine.");
  }
  return this->integration_step(*affine_vector_field,initial_set,step_size);
}


template<class R>
Geometry::Zonotope<R> 
Evaluation::AffineIntegrator< Geometry::Zonotope<R> >::
reachability_step(const System::VectorField<R>& vector_field, 
                  const Geometry::Zonotope<R>& initial_set, 
                  const Numeric::Rational& step_size, 
                  const Geometry::Box<R>& bounding_set) const
{
  const System::AffineVectorField<R>* affine_vector_field=this->cast(&vector_field);
  if(!affine_vector_field) {
    ARIADNE_THROW(std::runtime_error,"AffineIntegrator::reachability_step(...)","vector_field is not affine.");
  }
  return this->reachability_step(*affine_vector_field,initial_set,step_size);
}






template<class R> 
Geometry::Point< Numeric::Interval<R> > 
Evaluation::AffineIntegrator<  Geometry::Zonotope<R> >::
flow_step(const System::AffineVectorField<R>& avf,
          const Geometry::Point<I>& p,
          const Numeric::Rational& h) const
{
  ARIADNE_LOG(6,"AffineIntegrator::flow_step(AffineVectorField,Point<Interval>,Time) const\n");
  const LinearAlgebra::Matrix<I>& A=avf.A();
  const LinearAlgebra::Vector<I>& b=avf.b();
  return Point<I>(gexp(A,I(h),0)*p.position_vector() + gexp(A,I(h),1)*b);
}
     
 
template<class R> 
LinearAlgebra::Matrix< Numeric::Interval<R> > 
Evaluation::AffineIntegrator< Geometry::Zonotope<R> >::
flow_step_jacobian(const System::AffineVectorField<R>& avf,
                   const Geometry::Point<I>& p,
                   const Numeric::Rational& h) const
{
  ARIADNE_LOG(6,"AffineIntegrator::flow_step_jacobian(AffineVectorField,Point<Interval>,Time) const\n");
  const LinearAlgebra::Matrix<I>& A=avf.A();
  return gexp(A,I(h),0);
}




template<class R>
std::ostream&
Evaluation::AffineIntegrator< Geometry::Zonotope<R> >::
write(std::ostream& os) const
{
  return os << "AffineIntegrator( )";
}


}
