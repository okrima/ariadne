/***************************************************************************
 *            approximate_taylor_model.code.h
 *
 *  Copyright 2008  Pieter Collins
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

#include "base/exceptions.h"
#include "numeric/approximate_float.h"
#include "numeric/interval.h"
#include "linear_algebra/slice.h"
#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"
#include "differentiation/multi_index.h"
#include "differentiation/sparse_differential.h"
#include "function/approximate_taylor_model.h"
#include "function/function_interface.h"
#include "output/latexstream.h"


namespace Ariadne {

template<class R>
ApproximateTaylorModel<R>::ApproximateTaylorModel() 
  : _domain(), 
    _centre(),
    _expansion()
{
}


template<class R>
ApproximateTaylorModel<R>::ApproximateTaylorModel(uint rs, uint as, ushort o, ushort s) 
  : _domain(as,I(-1,1)),
    _centre(as),
    _expansion(rs,as,o)
{
}


template<class R>
ApproximateTaylorModel<R>::ApproximateTaylorModel(const Vector<I>& d, 
                                                  const Vector<A>& c, 
                                                  const SparseDifferentialVector<A>& e)
  : _domain(d),
    _centre(c),
    _expansion(e)
{
  assert(d.size()==e.argument_size());
  assert(c.size()==e.argument_size());
  for(uint i=0; i!=e.size(); ++i) {
    assert(e[i].argument_size()==e[0].argument_size());
  }
}




template<class R>
ApproximateTaylorModel<R>::ApproximateTaylorModel(const Vector<I>& d, const Vector<A>& c, 
                                                  ushort o, ushort s,
                                                  const FunctionInterface<R>& f)
  : _domain(d),
    _centre(c),
    _expansion(f.expansion(Vector<A>(c),o))
{
  assert(d.size()==f.argument_size());
  assert(c.size()==f.argument_size());
}




template<class R>
ApproximateTaylorModel<R>
ApproximateTaylorModel<R>::zero(uint rs, uint as)  
{
  return ApproximateTaylorModel<R>(rs,as,0,0);
}



template<class R>
ApproximateTaylorModel<R>
ApproximateTaylorModel<R>::one(uint as)  
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}



template<class R>
ApproximateTaylorModel<R>
ApproximateTaylorModel<R>::constant(uint as, const R& c)  
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}



template<class R>
Vector< Interval<R> >
ApproximateTaylorModel<R>::domain() const
{ 
  return this->_domain; 
}


template<class R>
Vector<R>
ApproximateTaylorModel<R>::centre() const
{ 
  Vector<R> r(this->_centre.size()); 
  for(uint i=0; i!=this->_centre.size(); ++i) {
    r[i]=this->_centre[i]._value; 
  }
  return r;
}


template<class R>
Vector< Interval<R> >
ApproximateTaylorModel<R>::range() const
{ 
  //  return Vector<I>(this->_expansion.value());
}



template<class R>
const SparseDifferentialVector<typename ApproximateTaylorModel<R>::A>&
ApproximateTaylorModel<R>::expansion() const
{ 
  return this->_expansion; 
}


template<class R>
uint 
ApproximateTaylorModel<R>::argument_size() const
{ 
  return this->_expansion.argument_size(); 
}


template<class R>
uint 
ApproximateTaylorModel<R>::result_size() const 
{ 
  return this->_expansion.result_size();
}


template<class R>
ushort 
ApproximateTaylorModel<R>::order() const 
{
  return this->_expansion.degree();
}
      

template<class R>
ushort 
ApproximateTaylorModel<R>::smoothness() const 
{ 
  return this->_expansion.degree();
}
      

template<class R> template<class X>
array< array<X> >
ApproximateTaylorModel<R>::_powers(const Vector<X>& v) const
{
  array< array<X> > powers(this->argument_size(), array<X>(this->order()+1));
  for(uint i=0; i!=this->argument_size(); ++i) {
    powers[i][0]=1;
    if(this->order()>=1) {
      powers[i][1]=v(i);
      if(this->order()>=2) {
        powers[i][2]=pow(v(i),2);
        for(uint j=3; j<=this->order(); ++j) {
          powers[i][j]=powers[i][2]*powers[i][j-2];
        }
      }
    }
  }
  return powers;
}


template<class R> 
Vector<typename ApproximateTaylorModel<R>::A> 
ApproximateTaylorModel<R>::evaluate(const Vector<A>& x) const
{
  assert(this->argument_size()==x.size());

  // TODO: Make this more efficient
  Vector<A> result(this->result_size());
  Vector<A> w=Vector<A>(x)-this->_centre;
  for(MultiIndex j(this->argument_size()); j.degree()<=this->order(); ++j) {
    A wa=1;
    for(uint k=0; k!=j.number_of_variables(); ++k) {
      wa*=pow(w[k],int(j[k]));
    }
    for(uint i=0; i!=this->result_size(); ++i) {
      const SparseDifferential<A>& cd=this->expansion()[i];
      result[i]+=cd[j]*wa;
    }
  }

  return result;
}


template<class R> 
Vector< Interval<R> > 
ApproximateTaylorModel<R>::evaluate(const Vector<I>& x) const
{
  assert(this->argument_size()==x.size());

  // TODO: Make this more efficient
  Vector<I> w=x-this->centre();
  Vector<I> result(this->result_size());
  for(MultiIndex j(this->argument_size()); j.degree()<=this->order(); ++j) {
    I wa=1;
    for(uint k=0; k!=j.number_of_variables(); ++k) {
      wa*=pow(w[k],int(j[k]));
    }
    for(uint i=0; i!=this->result_size(); ++i) {
      const SparseDifferential<A>& cd=this->expansion()[i];
      result[i]+=cd[j]*wa;
    }
  }

  return result;
}




template<class R> 
Matrix<typename ApproximateTaylorModel<R>::A> 
ApproximateTaylorModel<R>::jacobian(const Vector<A>& x) const
{
  assert(this->argument_size()==x.size());
  Matrix<A> J(this->result_size(),this->argument_size());
  Vector<A> w=x-this->_centre;
  array< array<A> > powers=this->_powers(w);
  for(uint j=0; j!=this->argument_size(); ++j) {
    for(MultiIndex m(this->argument_size()); m.degree()<=this->order(); ++m) {
      MultiIndex n=m;
      int c=n[j];
      if(c!=0) {
        n.decrement_index(j);
        A a=c;
        for(uint k=0; k!=this->argument_size(); ++k) {
          a*=powers[k][n[k]];
        }
        for(uint i=0; i!=this->result_size(); ++i) {
          A& x=J[i][j];
          x+=a*this->_expansion[i][m];
        }
      }
    }
  }  
  return J;
}


template<class R> 
Matrix< Interval<R> > 
ApproximateTaylorModel<R>::jacobian(const Vector<I>& x) const
{
  Matrix<I> J(this->result_size(),this->argument_size());
  Vector<I> w=x-this->centre();
  array< array<I> > powers=this->_powers(w);

  for(uint j=0; j!=this->argument_size(); ++j) {
    for(MultiIndex m(this->argument_size()); m.degree()<=this->order(); ++m) {
      MultiIndex n=m;
      int c=n[j];
      if(c!=0) {
        n.decrement_index(j);
        I a=c;
        for(uint k=0; k!=this->argument_size(); ++k) {
          a*=powers[k][n[k]];
        }
        for(uint i=0; i!=this->result_size(); ++i) {
          I amim=a*this->_expansion[i][m];
          J[i][j]+=amim;
        }
      }
    }
  }
  
  return J;
}


template<class R> 
ApproximateTaylorModel<R>
ApproximateTaylorModel<R>::identity(const Vector<I>& d, uint o)
{
  uint n=d.size();
  SparseDifferentialVector<A> ide(n,n,o);
  for(uint i=0; i!=n; ++i) {
    ide[i]=SparseDifferential<A>::variable(n,o,0.0,i);
  }
  return ApproximateTaylorModel<R>(d,midpoint(d),ide);
}                     


template<class R> 
ApproximateTaylorModel<R> 
project(const ApproximateTaylorModel<R>& f, Range rng)
{
  return ApproximateTaylorModel<R>(f.domain(),f.centre(),project(f.expansion(),rng));
}

template<class R> 
ApproximateTaylorModel<R> 
join(const ApproximateTaylorModel<R>& f, const ApproximateTaylorModel<R>& g)
{
  assert(f.domain()==g.domain());
  assert(f.centre()==g.centre());

  SparseDifferentialVector<R> he(f.result_size()+g.result_size(),f.argument_size(),f.order());
  project(he,range(0,f.result_size()))=f.expansion();
  project(he,range(f.result_size(),f.result_size()+g.result_size()))=g.expansion();

  return ApproximateTaylorModel<R> (f.domain(),f.centre(),he);
}


/*

ApproximateTaylorModel<R> 
compose(const ApproximateTaylorModel<R>& f,
        const AffineTransformation<R>& g)
{
  assert(f.centre()==g.b());
  SparseDifferentialVector<R> const& fe=f.expansion();
  Vector<R> const& b=g.b();
  Matrix<R> const& A=g.A();

  SparseDifferentialVector<R> he;

  // Check that the matrix has zero or one unit in each row
  for(uint i=0; i!=A.row_size(); ++i) {
    uint ones=0;
    for(uint j=0; i!=A.column_size(); ++i) {
      if(A[i][j]==1) { ++ones; }
      else { assert(A[i][j]==0); }
    }
    assert(ones<=1);
  }
  
  assert(false); // not called
}

ApproximateTaylorModel<R> 
compose(const AffineTransformation<R>& f,
        const ApproximateTaylorModel<R>& g)
{
  assert(f.centre()==g.expansion().value());
  Vector<R> const& c=f.centre();
  Vector<R> const& b=f.b();
  Matrix<R> const& A=f.A();
  SparseDifferentialVector<R> const& ge=g.expansion();

  SparseDifferentialVector<R> he(f.result_size(),g.argument_size(),g.order());
  for(uint i=0; i!=f.result_size(); ++i) {
    he[i]+=b[i];
    for(uint j=0; j!=f.argument_size(); ++j) {
      he[i]+=A[i][j]*ge[j];
    }
  }
  return ApproximateTaylorModel<R>(g.domain(),g.centre(),he);
}

*/



template<class R> 
ApproximateTaylorModel<R> 
compose(const ApproximateTaylorModel<R>& f,
        const ApproximateTaylorModel<R>& g)
{
  //std::cerr<<__PRETTY_FUNCTION__<<std::endl;
  assert(f.argument_size()==g.result_size());
  typedef Interval<R> I;
  SparseDifferentialVector<R> const& fe=f.expansion();
  //std::cerr<<"fe="<<fe<<std::endl;
  SparseDifferentialVector<R> const& ge=g.expansion();
  //std::cerr<<"ge="<<ge<<std::endl;
  Vector<R> tr=-f.centre()+ge.value();
  //std::cerr<<"tr="<<tr<<std::endl;
  SparseDifferentialVector<R> fet=translate(fe,tr);
  //std::cerr<<"fet="<<fet<<std::endl;
  Vector<I> hd=g.domain();
  Vector<R> hc=g.centre();
  //std::cerr<<"he="<<std::flush;
  SparseDifferentialVector<R> he=compose(fet,ge);
  // FIXME: Change domain
  //std::cerr<<he<<std::endl;
  return ApproximateTaylorModel<R>(hd,hc,he);
}


template<class R> 
ApproximateTaylorModel<R> 
recentre(const ApproximateTaylorModel<R>& f,
         const Vector< Interval<R> >& d,
         const Vector<R>& c)
{
}


template<class R> 
ApproximateTaylorModel<R> 
inverse(const ApproximateTaylorModel<R>& f,
        const Vector<R>& x)
{
  typedef Interval<R> I;
  Vector<R> tr=f.centre()-x;
  SparseDifferentialVector<R> fet=translate(f.expansion(),tr);
  Vector<I> gd=f.domain();
  Vector<R> gc=fet.value();
  SparseDifferentialVector<R> ge=inverse(fet);
  // FIXME: Change domain
  return ApproximateTaylorModel<R>(gd,gc,ge);
}


template<class R> 
ApproximateTaylorModel<R>
implicit(const ApproximateTaylorModel<R>& f)
{
  typedef Interval<R> I;
  assert(f.argument_size()>f.result_size());
  uint m=f.argument_size(); 
  uint n=f.result_size();

  array<uint> p(n);
  for(uint i=0; i!=n; ++i) { p[i]=i+m-n; }
  //std::cerr << "p=" << p << std::endl;

  // Construct the Taylor model for g(y)=f(c,y)
  Interval<R> gd=project(f.domain(),range(m-n,m));
  //std::cerr << "gd=" << gd << std::endl;
  Vector<R> gc=project(f.centre(),range(m-n,m));
  //std::cerr << "gc=" << gc << std::endl;
  SparseDifferentialVector<R> ge=restrict(f.expansion(),p);
  //std::cerr << "ge=" << ge << std::endl;
  ApproximateTaylorModel<R> g(gd,gc,ge);

  Vector<R> z(n);
  Vector<I> iv = solve(g,z);
  Vector<R> v = midpoint(iv);
  //std::cerr<<"iv="<<iv<<std::endl;
  //std::cerr<<"v="<<v<<std::endl;
  Vector<R> t(m);
  project(t,range(m-n,m))=v;
  //std::cerr<<"t="<<t<<std::endl;

  SparseDifferentialVector<R> fe=f.expansion();
  //std::cerr<<"fe="<<fe<<std::endl;
  SparseDifferentialVector<R> fet=translate(fe,t);
  //std::cerr<<"fet="<<fet<<std::endl;
  fet.set_value(Vector<R>(fe.result_size(),0.0));
  //std::cerr<<"fet="<<fet<<std::endl;

  SparseDifferentialVector<R> he=implicit(fet);
  //std::cerr<<"he="<<he<<std::endl;
  he+=v;
  //std::cerr<<"he="<<he<<std::endl;
  Vector<I> hd=project(f.domain(),range(0,m-n));
  Vector<R> hc=project(f.centre(),range(0,m-n));
  return ApproximateTaylorModel<R>(hd,hc,he);
}

// If \f$f:\R^m\rightarrow \R^n\f$, compute the function
// \f$ h:\R^{m-n}\rightarrow \R^n\f$ satisfying 
// \f$ f(x,h(x)) = f(c) \f$ for all \f$x\f$.
// The function \f$ h\f$ is centred at
// \f$ c=(c_1,\ldots,c_{m-n}) \f$
// and has value \f$h(x)=(c_{m-n+1},\ldots,c_n)\f$.
template<class R> 
ApproximateTaylorModel<R>
implicit1(const ApproximateTaylorModel<R>& f,
          const Vector<R>& c)
{
  typedef Interval<R> I;
  uint m=f.argument_size(); 
  uint n=f.result_size();
  Vector<R> v=f.evaluate(c);
  Vector<R> hc=project(c,range(0,m-n));
  Vector<I> hd=project(f.domain(),range(0,m-n));
  
}

// If \f$f:\R^m\rightarrow \R^n\f$, compute the function
// \f$ h:\R^{m-n}\rightarrow \R^n\f$ satisfying 
// \f$ f(x,h(x)) = z \f$ for all \f$c\f$.
// The function \f$ h\f$ is centred at
// \f$ c=(c_1,\ldots,c_{m-n}) \f$.
//
// Precondition: there is a unique solution 
// of f(c,y)=z$.
template<class R> 
ApproximateTaylorModel<R>
implicit2(const ApproximateTaylorModel<R>& f,
          const Vector<R>& c,
          const Vector<R>& z)
{
  typedef Interval<R> I;
  // Solve for the value of h(c)
  uint m=f.argument_size(); 
  uint n=f.result_size();
  Vector<R> v=f.evaluate(c);
  Matrix<R> A=f.expansion().jacobian();
  //A=A[slice(0,n)][slice(m-n,n)];
  Matrix<R> B=inverse(A);
  Vector<R> y(n);
  Vector<R> tr=y-f.centre();
  
  SparseDifferentialVector<R> fet=translate(f.expansion(),tr)-z;
  SparseDifferentialVector<R> he=implicit(fet)+y;
  
  Vector<I> hd(m-n,I(-inf(),+inf()));
  Vector<R> hc=c;
  return ApproximateTaylorModel<R>(hd,hc,he);
}

template<class R> 
ApproximateTaylorModel<R> 
flow(const ApproximateTaylorModel<R>& vf, uint ox)
{
  typedef Interval<R> I;

  Vector< SparseSeries< SparseDifferential<R> > > y=flow(vf.expansion(),vf.centre(),ox);
  
  uint n=vf.result_size();
  SparseDifferentialVector<R> x(n,n+1,ox);
  for(MultiIndex jx(n+1); jx.degree()<=ox; ++jx) {
    MultiIndex jy(n);
    for(uint k=0; k!=n; ++k) { 
      jy.set(k,jx[k]); 
    }
    uint jt=jx[n];
    for(uint i=0; i!=n; ++i) {
      x[i][jx]=y[i][jt][jy];
    }
  }

  I h(-0.1,0.1);
  Vector<I> d=join(vf.domain(),h);
  Vector<R> c=join(vf.centre(),0.0);
  return ApproximateTaylorModel<R>(d,c,x);
}


template<class R> 
ApproximateTaylorModel<R> 
integrate(const ApproximateTaylorModel<R>& vf, const R& h, uint ox)
{
  typedef Interval<R> I;
  assert(vf.result_size()==vf.argument_size());
  uint n=vf.result_size();
  SparseDifferentialVector<R> xe=vector_variable(n,ox,Vector<R>(n+1,0.0));
  SparseDifferential<R> te=scalar_constant(n,ox,h);
  ApproximateTaylorModel<R> f=flow(vf,ox);
  Vector<I> xtd=join(vf.domain(),I(-1,1));
  Vector<R> xtc=join(vf.centre(),0.0);
  ApproximateTaylorModel<R> xt(xtd,xtc,join(xe,te));
  return compose(f,xt);
}


// Find the set of values such that \f$f(x)=y\f$
// assuming a unique solution in 
template<class R> 
Vector< Interval<R> >
solve(const ApproximateTaylorModel<R>& f,
      const Vector<R>& y)
{
  typedef Interval<R> I;
  Vector<I> x=f.centre();
  //Vector<I> x=f.domain();
  Matrix<R> J,Jinv;
  Vector<I> nx,fx,fxmy;
  Vector<R> m;
  for(uint i=0; i!=6; ++i) {
    //std::cerr << "  x[" << i << "]="<<x <<" y="<<y<<"\n";
    m=midpoint(x);
    J=f.jacobian(m);
    Jinv=inverse(J);
    fx=f.evaluate(x);
    fxmy=fx-y;
    nx=m-Jinv*Vector<I>(fx-y);
    if(disjoint(x,nx)) {
      x=nx;
    } else {
      x=intersection(x,nx);
    }
    //std::cerr<<"    m="<<m<<" J="<<J<<" Jinv="<<Jinv<<" fx="<<fx<<" fx-y="<<fxmy<<" nx="<<nx<<"\n";
  }
  return x;
}



// Compute the hitting map of the flow of vf(x)
// with the guard condition g(x)=0
template<class R> 
ApproximateTaylorModel<R>
hitting(const ApproximateTaylorModel<R>& vf,
        const ApproximateTaylorModel<R>& g)
{
  typedef Interval<R> I;
  assert(vf.result_size()==vf.argument_size());
  assert(g.argument_size()==vf.result_size());
  assert(g.result_size()==1);
  assert(g.expansion().degree()==vf.expansion().degree());
  uint n=vf.result_size();
  uint d=vf.order();
  //std::cerr<<"vf="<<vf<<std::endl; // n,n
  //std::cerr<<"g="<<g<<std::endl; // 1,n+1
  //Differential<R> t=variable(n+1,d,0.0,n); // 1,n+1
  //std::cerr<<"t="<<t<<std::endl;
  ApproximateTaylorModel<R> f=flow(vf,d); // n,n+1
  //std::cerr<<"f="<<f<<std::endl;
  ApproximateTaylorModel<R> gf=compose(g,f); // 1,n+1
  //std::cerr<<"gf="<<gf<<std::endl;
  ApproximateTaylorModel<R> ht=implicit(gf); // 1,n
  //std::cerr<<"ht="<<ht<<std::endl;

  SparseDifferentialVector<R> xe=vector_variable(n,d,vf.centre());
  //std::cerr<<"  xe="<<xe<<std::endl;
  ApproximateTaylorModel<R> hxt(ht.domain(),ht.centre(),join(xe,ht.expansion())); // n+1,n
  //std::cerr<<"hxt="<<hxt<<std::endl;
  
  ApproximateTaylorModel<R> hy=compose(f,hxt); // n,n
  //std::cerr<<"hy="<<hy<<std::endl;
  return join(hy,ht);
}



template<class R> 
ApproximateTaylorModel<R> 
ApproximateTaylorModel<R>::truncate(const Vector<I>& domain, const Vector<R>& centre,
                                    ushort order, ushort smoothness) const
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}



template<class R> 
std::ostream&
operator<<(std::ostream& os, const ApproximateTaylorModel<R>& tm)
{
  os << "ApproximateTaylorModel(\n";
  os << "  domain=" << tm.domain() << ",\n" << std::flush;
  os << "  centre=" << tm.centre() << ",\n" << std::flush;
  os << "  expansion=" << tm.expansion() << ",\n" << std::flush;
  os << ")\n";
  return os;
}


template<class R> 
latexstream&
operator<<(latexstream& texs, const ApproximateTaylorModel<R>& p)
{
  texs << "%ApproximateTaylorModel<R>\n";
  texs << "\\ensuremath{\n";
  texs << "\\left( \\begin{array}{c}\n";
  char var='x';
  for(uint i=0; i!=p.result_size(); ++i) {
    bool first = true;
    if(i!=0) { texs << "\\\\"; }
    for(MultiIndex j(p.argument_size()); j.degree()<=p.order(); ++j) {
      const R& a=p.expansion()[i][j];
      if(a!=0) {
        if(first) { first=false; }
        else { if(a>0) { texs << '+'; } }
        if(a==1) { if(j.degree()==0) { texs << a; } }
        else if(a==-1) { if(j.degree()==0) { texs << a; } else { texs << '-'; } }
        else { texs << a << ' '; }
        for(uint k=0; k!=p.argument_size(); ++k) {
          if(j[k]!=0) {
            texs << var << "_{ " << k+1 << "}";
            if(j[k]!=1) {
              texs << "^{" << j[k] << "}";
            }
          }
          texs << " ";
        }
      } 
    }
    texs << "\n"; 
  }
  texs << "\\end{array}\\right)\n}\n";
  return texs;
}


} // namespace Ariadne