/***************************************************************************
 *            taylor_model.cc
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
 *  GNU Library General Public License for more detai1ls.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#include "numeric/numeric.h"
#include "config.h"

#include <iomanip>
#include <limits>

#include "numeric/rounding.h"
#include "numeric/numeric.h"
#include "algebra/vector.h"
#include "algebra/matrix.h"
#include "algebra/expansion.h"
#include "algebra/series.h"
#include "algebra/differential.h"
#include "function/taylor_model.h"
#include "function/taylor_series.h"
#include "function/function.h"
#include "utility/exceptions.h"

#include "algebra/algebra_mixin.tcc"

#define VOLATILE ;
#include "algebra/multi_index-noaliasing.h"
#include "function/function_mixin.h"
#include "algebra/vector.h"

namespace Ariadne {

namespace {

Bool operator<(const MultiIndex& a1, const MultiIndex& a2) {
    return reverse_lexicographic_less(a1,a2); }

} // namespace


const double em=2.2204460492503131e-16;
const double ec=em/2;

/*

typedef Float ErrorFloat;
typedef Float ApproxFloat;

Void axpy(Float& te, ApproxFloat& r, const Float& a, const Float& x, const Float& y) {
    VOLATILE Float& xv=const_cast<VOLATILE Float&>(x);
    VOLATILE Float mxv=-x;
    set_rounding_upward();
    VOLATILE Float u=a*xv+y;
    VOLATILE Float ml=a*mxv-y;
    te+=(u+ml);
    set_rounding_to_nearest();
    r=a*xv+y;
}

Void add_op(Float& te, ApproxFloat& r, const Float& x, const Float& y) {
    VOLATILE Float& xv=const_cast<VOLATILE Float&>(x);
    VOLATILE Float mxv=-x;
    set_rounding_upward();
    VOLATILE Float u=xv+y;
    VOLATILE Float ml=mxv-y;
    te+=(u+ml);
    set_rounding_to_nearest();
    r=xv+y;
}

Void sub_op(Float& te, ApproxFloat& r, const Float& x, const Float& y) {
    VOLATILE Float& xv=const_cast<VOLATILE Float&>(x);
    VOLATILE Float mxv=-x;
    set_rounding_upward();
    VOLATILE Float u=xv-y;
    VOLATILE Float ml=mxv+y;
    te+=(u+ml);
    set_rounding_to_nearest();
    r=xv-y;
}

Void mul_op(Float& te, ApproxFloat& r, const Float& s, const Float& x) {
    VOLATILE Float& xv=const_cast<VOLATILE Float&>(x);
    VOLATILE Float mxv=-x;
    set_rounding_upward();
    VOLATILE Float u=xv*s;
    VOLATILE Float ml=mxv*s;
    te+=(u+ml);
    set_rounding_to_nearest();
    r=xv*s;
}

Void mul_op(Float& te, ApproxFloat& r, const Float& sl, const Float& sm, const Float& su, const Float& x) {
    VOLATILE Float& xv=const_cast<VOLATILE Float&>(x);
    VOLATILE Float mxv=-x;
    set_rounding_upward();
    VOLATILE Float u,ml;
    if(x>=0) {
        u=xv*su;
        ml=mxv*sl;
    } else {
        u=xv*sl;
        ml=mxv*su;
    }
    te+=(u+ml);
    set_rounding_to_nearest();
    r=xv*sm;
}

*/



Vector<ValidatedNumber> unscale(const Vector<ValidatedNumber>& x, const Vector<ExactInterval>& d) {
    Vector<ValidatedNumber> r(x);
    for(Nat i=0; i!=r.size(); ++i) {
        if(d[i].lower().raw()==d[i].upper().raw()) {
            if(x[i]==d[i].midpoint()) {
                r[i]=ValidatedNumber(0.0,0.0);
            } else {
                r[i]=ValidatedNumber(-inf,+inf);
            }
        } else {
            //r[i]=(2*r[i]-add_ivl(d[i].lower().raw(),d[i].upper().raw()))/sub_ivl(d[i].upper().raw(),d[i].lower().raw());
            r[i]=(2*r[i]-(d[i].lower()+d[i].upper()))/(d[i].upper()-d[i].lower());
        }
    }
    return r;
}



/* */

ValidatedTaylorModel::TaylorModel()
    : _expansion(0), _error(0), _sweeper()
{ }


ValidatedTaylorModel::TaylorModel(Nat as, Sweeper swp)
    : _expansion(as), _error(0), _sweeper(swp)
{
}

ValidatedTaylorModel::TaylorModel(const Expansion<CoefficientType>& f, const ErrorType& e, Sweeper swp)
    : _expansion(f), _error(e), _sweeper(swp)
{
    this->unique_sort();
    this->sweep();
}

ValidatedTaylorModel::TaylorModel(const Expansion<Float>& f, const Float& e, Sweeper swp)
    : TaylorModel(reinterpret_cast<Expansion<CoefficientType>const&>(f),reinterpret_cast<ErrorType const&>(e),swp)
{
}


ValidatedTaylorModel
ValidatedTaylorModel::create() const
{
    return ValidatedTaylorModel(this->argument_size(),this->_sweeper);
}

ValidatedTaylorModel
ValidatedTaylorModel::create_zero() const
{
    return ValidatedTaylorModel(this->argument_size(),this->_sweeper);
}

ValidatedTaylorModel
ValidatedTaylorModel::create_ball(ErrorType e) const
{
    ARIADNE_PRECONDITION(e>=0);
    ValidatedTaylorModel r(this->argument_size(),this->_sweeper);
    r.set_error(e);
    return r;
}

Void
ValidatedTaylorModel::swap(ValidatedTaylorModel& tm)
{
    this->_expansion.swap(tm._expansion);
    std::swap(this->_error,tm._error);
    std::swap(this->_sweeper,tm._sweeper);
}

Void
ValidatedTaylorModel::clear()
{
    this->_expansion.clear();
    this->_error=0u;
}

Nat
ValidatedTaylorModel::degree() const
{
    uchar deg=0u;
    for(ValidatedTaylorModel::ConstIterator iter=this->begin(); iter!=this->end(); ++iter) {
        deg=std::max(deg,iter->key().degree());
    }
    return deg;
}


ValidatedTaylorModel&
ValidatedTaylorModel::operator=(const ValidatedNumber& c)
{
    this->_expansion.clear();
    ExactFloat m=c.midpoint();
    if(m!=0) {
        this->_expansion.append(MultiIndex::zero(this->argument_size()),m);
    }
    this->_error=c.radius();
    return *this;
}


namespace { // Internal code for arithmetic

// Inplace negation
Void _neg(ValidatedTaylorModel& r)
{
    for(ValidatedTaylorModel::Iterator iter=r.begin(); iter!=r.end(); ++iter) {
        iter->data().raw()=-iter->data().raw();
    }
}

inline Void _scal_exact(ValidatedTaylorModel& r, const Float& c)
{
    // Operation can be performed exactly
    RawFloat pc=static_cast<RawFloat>(c);
    for(ValidatedTaylorModel::Iterator iter=r.begin(); iter!=r.end(); ++iter) {
        iter->data().raw()*=pc;
    }
    r.error().raw()*=abs(pc);
    return;
}




inline Void _scal_approx(ValidatedTaylorModel& r, const Float& c)
{
    // General case with error analysis
    Float& re=r.error().raw();
    set_rounding_upward();
    Float u,ml;
    Float te=0; // Twice the maximum accumulated error
    Float pc=c;
    Float mc=-c;
    for(ValidatedTaylorModel::ConstIterator riter=r.begin(); riter!=r.end(); ++riter) {
        const Float& rv=riter->data().raw();
        u=rv*pc;
        ml=rv*mc;
        te+=(u+ml);
    }
    re*=abs(c);
    re+=te/2;

    set_rounding_to_nearest();
    Float m=c;
    for(ValidatedTaylorModel::Iterator riter=r.begin(); riter!=r.end(); ++riter) {
        riter->data().raw()*=m;
    }
    return;
}



inline Void _scal_approx4(ValidatedTaylorModel& r, const Float& c)
{
    // General case with error analysis
    Float& re=r.error().raw();
    set_rounding_upward();
    Float u,ml;
    //Float u,ml;
    Float te=0; // Twice the maximum accumulated error
    Float pc=c;
    Float mc=-pc;
    for(ValidatedTaylorModel::ConstIterator riter=r.begin(); riter!=r.end(); ++riter) {
        const Float& rv=riter->data().raw();
        //volatile Float& rv=const_cast<volatile Float&>(riter->data().raw());
        u=rv*pc;
        ml=rv*mc;
        te+=u+ml;
    }
    re*=abs(pc);
    re+=te/2;

    set_rounding_to_nearest();
    Float m=c;
    for(ValidatedTaylorModel::Iterator riter=r.begin(); riter!=r.end(); ++riter) {
        riter->data().raw()*=m;
    }
    return;
}

inline Void _scal_approx3(ValidatedTaylorModel& r, const Float& c)
{
    // General case with error analysis
    double& re=internal_cast<double&>(r.error().raw());
    set_rounding_upward();
    volatile double u,ml;
    double te=0; // Twice the maximum accumulated error
    double pc=internal_cast<const double&>(c);
    double mc=-pc;
    for(ValidatedTaylorModel::ConstIterator riter=r.begin(); riter!=r.end(); ++riter) {
        const double& rv=internal_cast<const double&>(riter->data().raw());
        u=rv*pc;
        ml=rv*mc;
        te+=(u+ml);
    }
    re*=std::abs(pc);
    re+=te/2;

    set_rounding_to_nearest();
    Float m=c;
    for(ValidatedTaylorModel::Iterator riter=r.begin(); riter!=r.end(); ++riter) {
        riter->data().raw()*=m;
    }
    return;
}


// Version using only Float class directly
inline Void _scal_approx2(ValidatedTaylorModel& r, const Float& c)
{
    // General case with error analysis
    Float& re=r.error().raw();
    set_rounding_upward();
    Float u,ml;
    Float te=0; // Twice the maximum accumulated error
    Float pc=c;
    Float mc=-c;
    for(ValidatedTaylorModel::ConstIterator riter=r.begin(); riter!=r.end(); ++riter) {
        const Float& rv=riter->data().raw();
        u=mul_rnd(rv,pc);
        ml=mul_rnd(rv,mc);
        te+=(u+ml);
    }
    re*=abs(pc);
    re+=te/2;

    set_rounding_to_nearest();
    Float m=c;
    for(ValidatedTaylorModel::Iterator riter=r.begin(); riter!=r.end(); ++riter) {
        riter->data().raw()*=m;
    }
    return;
}

// Modified version to try to obtain optimal performance
inline Void _scal_approx1(ValidatedTaylorModel& rr, const Float& cc)
{
    // General case with error analysis
    Expansion<double>& r=reinterpret_cast<Expansion<double>&>(rr.expansion());
    double& re=reinterpret_cast<double&>(rr.error().raw());
    const double& c=reinterpret_cast<const double&>(cc);

    set_rounding_upward();
    volatile double u,ml;
    double te=0; // Twice the maximum accumulated error
    double pc=c;
    double mc=-c;
    for(Expansion<double>::ConstIterator riter=r.begin(); riter!=r.end(); ++riter) {
        const double& rv=riter->data();
        u=rv*pc;
        ml=rv*mc;
        te+=(u+ml);
    }
    re*=std::abs(pc);
    re+=te/2;

    set_rounding_to_nearest();
    double m=c;
    for(Expansion<double>::Iterator riter=r.begin(); riter!=r.end(); ++riter) {
        riter->data()*=m;
    }
    return;
}


// Original version except with double
inline Void _scal_approx0(ValidatedTaylorModel& rr, const Float& cc)
{
    Expansion<double>& r=reinterpret_cast<Expansion<double>&>(rr.expansion());
    volatile double& re=reinterpret_cast<volatile double&>(rr.error());
    volatile double c=reinterpret_cast<const double&>(cc);
    //const double& c=cc.v;

    set_rounding_upward();
    volatile double u,ml;
    double te=0; // Twice the maximum accumulated error
    double pc=c;
    double mc=-c;
    for(Expansion<double>::ConstIterator riter=r.begin(); riter!=r.end(); ++riter) {
        const double& rv=riter->data();
        u=rv*pc;
        ml=rv*mc;
        te+=(u+ml);
    }
    re*=std::abs(c);
    re+=te/2;

    set_rounding_to_nearest();
    double m=c;
    for(Expansion<double>::Iterator riter=r.begin(); riter!=r.end(); ++riter) {
        riter->data()*=m;
    }
    return;
}

Void _scal(ValidatedTaylorModel& r, const Float& c)
{
    // No measurable speedup in general case by avoiding checks
    _scal_approx4(r,c); return;
    if(c==0.0) { r.expansion().clear(); r.error().raw()=0; }
    else if(c==1.0) { }
    else if(c==2.0 || c==0.5) { _scal_exact(r,c); }
    else { _scal_approx(r,c); }
}


Void _scal(ValidatedTaylorModel& r, const ValidatedNumber& c)
{
    ARIADNE_ASSERT_MSG(c.lower().raw()<=c.upper().raw(),c);
    ARIADNE_ASSERT_MSG(r.error().raw()>=0,"r="<<r);

    if(r.error().raw()==inf) { r.expansion().clear(); return; }

    //std::cerr<<"ValidatedTaylorModel::scal(ValidatedNumber c) c="<<c<<std::endl;
    Float& re=r.error().raw();
    set_rounding_upward();
    VOLATILE Float u,ml;
    Float te=0; // Twice the maximum accumulated error
    Float cu=c.upper().raw();
    Float mcl=-c.lower().raw();
    for(ValidatedTaylorModel::ConstIterator riter=r.begin(); riter!=r.end(); ++riter) {
        const Float& rv=riter->data().raw();
        if(rv>=0) {
            u=rv*cu;
            ml=rv*mcl;
        } else {
            Float mrv=-rv;
            u=mrv*mcl;
            ml=mrv*cu;
        }
        te+=(u+ml);
    }
    re*=mag(c).raw();
    re+=te/2;

    set_rounding_to_nearest();
    Float m=(c.upper().raw()+c.lower().raw())/2;
    for(ValidatedTaylorModel::Iterator riter=r.begin(); riter!=r.end(); ++riter) {
        riter->data().raw()*=m;
    }

    ARIADNE_ASSERT(r.error().raw()>=0);
    return;
}

Void _scal2(ValidatedTaylorModel& r, const ValidatedNumber& c)
{
    //std::cerr<<"ValidatedTaylorModel::scal(ValidatedNumber c) c="<<c<<std::endl;
    Float& re=r.error().raw();
    VOLATILE Float u,ml;
    Float te=0; // Twice the maximum accumulated error
    VOLATILE Float cu=c.upper().raw();
    VOLATILE Float mcl=-c.lower().raw();
    set_rounding_to_nearest();
    Float cm=(c.lower().raw()+c.upper().raw())/2;
    for(ValidatedTaylorModel::Iterator riter=r.begin(); riter!=r.end(); ++riter) {
        Float& rv=riter->data().raw();
        set_rounding_upward();
        if(rv>=0) {
            u=rv*cu;
            ml=rv*mcl;
        } else {
            Float mrv=-rv;
            u=mrv*mcl;
            ml=mrv*cu;
        }
        te+=(u+ml);
        set_rounding_to_nearest();
        rv*=cm;
    }
    set_rounding_upward();
    re*=mag(c).raw();
    re+=te/2;
    set_rounding_to_nearest();

    ARIADNE_ASSERT(r.error().raw()>=0);
    return;
}


struct UnitMultiIndex { Nat argument_size; Nat unit_index; };


inline Void _incr(ValidatedTaylorModel& r, const MultiIndex& a)
{
    for(ValidatedTaylorModel::Iterator iter=r.begin(); iter!=r.end(); ++iter) {
        static_cast<MultiIndex&>(iter->key())+=a;
    }
}

inline Void _incr(ValidatedTaylorModel& r, Nat j)
{
    for(ValidatedTaylorModel::Iterator iter=r.begin(); iter!=r.end(); ++iter) {
        ++static_cast<MultiIndex&>(iter->key())[j];
    }
}


inline
Void _acc(ValidatedTaylorModel& r, const ExactFloat& c)
{
    // Compute self+=c
    ARIADNE_DEBUG_ASSERT(r.error().raw()>=0);
    if(c==0) { return; }
    if(r.expansion().empty()) {
        r.expansion().raw().append(MultiIndex(r.argument_size()),c.raw());
    } else if((r.end()-1)->key().degree()>0) {
        r.expansion().raw().append(MultiIndex(r.argument_size()),c.raw());
    } else {
        Float& rv=(r.end()-1)->data().raw();
        Float& re=r.error().raw();
        set_rounding_upward();
        VOLATILE Float rvu=rv+c.raw();
        VOLATILE Float mrvl=(-rv)-c.raw();
        //std::cerr<<"re="<<re.u<<" ";
        re+=(rvu+mrvl)/2;
        //std::cerr<<"nre="<<re.u<<"\n";
        set_rounding_to_nearest();
        rv+=c.raw();
    }
    ARIADNE_ASSERT_MSG(r.error().raw()>=0,"c="<<c<<" r="<<r);
    return;
}



inline
Void _acc(ValidatedTaylorModel& r, const ValidatedNumber& c)
{
    // Compute self+=c
    ARIADNE_ASSERT_MSG(r.error().raw()>=0,r);

    if(c.lower().raw()==-inf || c.upper().raw()==+inf) {
        r.clear();
        r.set_error(+infty);
        return;
    }

    if(c.lower().raw()==-c.upper().raw()) { // The midpoint of the interval is zero, so no need to change constant term
        set_rounding_upward();
        r.error().raw()+=c.upper().raw();
        set_rounding_to_nearest();
        return;
    }
    if(r.expansion().empty()) { // Append a constant term zero
        r.expansion().raw().append(MultiIndex(r.argument_size()),0.0);
    } else if((r.end()-1)->key().degree()>0) { // Append a constant term zero
        r.expansion().raw().append(MultiIndex(r.argument_size()),0.0);
    }

    Float& rv=(r.end()-1)->data().raw();
    Float& re=r.error().raw();
    set_rounding_upward();
    VOLATILE Float rvu=rv+c.upper().raw();
    VOLATILE Float mrvl=(-rv)-c.lower().raw();
    re+=(rvu+mrvl)/2;
    set_rounding_to_nearest();
    VOLATILE Float m=(c.upper().raw()+c.lower().raw())/2;
    rv+=m;
    ARIADNE_DEBUG_ASSERT(r.error().raw()>=0);
}



// Compute r=x+y, assuming r is empty
// First compute the errors and then at the end compute the expansion
// This saves rounding mode changes, but requires an extra loop
inline Void _add1(ValidatedTaylorModel& r, const ValidatedTaylorModel& x, const ValidatedTaylorModel& y)
{
    // Compute r=x+y, assuming r is empty
    set_rounding_upward();
    Float te=0.0;
    ValidatedTaylorModel::ConstIterator xiter=x.begin();
    ValidatedTaylorModel::ConstIterator yiter=y.begin();
    while(xiter!=x.end() && yiter!=y.end()) {
        if(xiter->key()<yiter->key()) {
            ++xiter;
        } else if(yiter->key()<xiter->key()) {
            ++yiter;
        } else {
            const Float& xv=xiter->data().raw();
            const Float& yv=yiter->data().raw();
            VOLATILE Float u=xv+yv;
            VOLATILE Float t=-xv;
            VOLATILE Float ml=t-yv;
            te+=(u+ml);
            ++xiter; ++yiter;
        }
    }
    r.error().raw()=x.error().raw();
    r.error().raw()+=y.error().raw();
    r.error().raw()+=(te/2);

    set_rounding_to_nearest();
    xiter=x.begin();
    yiter=y.begin();
    while(xiter!=x.end() && yiter!=y.end()) {
        if(xiter->key()<yiter->key()) {
            r.expansion().raw().append(xiter->key(),xiter->data().raw());
            ++xiter;
        } else if(yiter->key()<xiter->key()) {
            r.expansion().raw().append(yiter->key(),yiter->data().raw());
            ++yiter;
        } else {
            RawFloat c=xiter->data().raw()+yiter->data().raw();
            if(c!=0) { r.expansion().raw().append(xiter->key(),c); }
            ++xiter; ++yiter;
        }
    }
    while(xiter!=x.end()) {
        r.expansion().raw().append(xiter->key(),xiter->data().raw());
        ++xiter;
    }
    while(yiter!=y.end()) {
        r.expansion().raw().append(yiter->key(),yiter->data().raw());
        ++yiter;
    }
    ARIADNE_ASSERT(r.error().raw()>=0);

}

// Compute r=x+y, assuming r is empty.
// Use a rounding mode change every iteration, as this appears to be faster
//   than using two loops
// Use opposite rounding to compute difference of upward and downward roundings,
//   as this seems to be marginally faster than changing the rounding mode
inline Void _add2(ValidatedTaylorModel& r, const ValidatedTaylorModel& x, const ValidatedTaylorModel& y)
{
    set_rounding_upward();
    Float te=0.0;
    ValidatedTaylorModel::ConstIterator xiter=x.begin();
    ValidatedTaylorModel::ConstIterator yiter=y.begin();
    while(xiter!=x.end() && yiter!=y.end()) {
        if(xiter->key()<yiter->key()) {
            r.expansion().raw().append(xiter->key(),xiter->data().raw());
            ++xiter;
        } else if(yiter->key()<xiter->key()) {
            r.expansion().raw().append(yiter->key(),yiter->data().raw());
            ++yiter;
        } else {
            assert(xiter->key()==yiter->key());
            VOLATILE Float& xv=const_cast<Float&>(static_cast<const Float&>(xiter->data().raw()));
            VOLATILE Float& yv=const_cast<Float&>(static_cast<const Float&>(yiter->data().raw()));
            set_rounding_upward();
            VOLATILE Float u=xv+yv;
            VOLATILE Float ml=-xv; ml-=yv;
            te+=(u+ml);
            set_rounding_to_nearest();
            Float c=xiter->data().raw()+yiter->data().raw();
            if(c!=0) { r.expansion().raw().append(xiter->key(),xiter->data().raw()+yiter->data().raw()); }
            ++xiter; ++yiter;
        }
    }

    while(xiter!=x.end()) {
        r.expansion().raw().append(xiter->key(),xiter->data().raw());
        ++xiter;
    }
    while(yiter!=y.end()) {
        r.expansion().raw().append(yiter->key(),yiter->data().raw());
        ++yiter;
    }

    set_rounding_upward();
    r.error().raw()=x.error().raw();
    r.error().raw()+=y.error().raw();
    r.error().raw()+=(te/2);
    set_rounding_to_nearest();
    ARIADNE_ASSERT(r.error().raw()>=0);

}

inline Void _add(ValidatedTaylorModel& r, const ValidatedTaylorModel& x, const ValidatedTaylorModel& y)
{
    _add2(r,x,y);
    ARIADNE_ASSERT(r.error().raw()>=0);
}

inline Void _acc(ValidatedTaylorModel& r, const ValidatedTaylorModel& x)
{
    ValidatedTaylorModel s(r.argument_size(),r.sweeper()); _add(s,r,x); s.swap(r);
    ARIADNE_ASSERT(r.error().raw()>=0);
}

inline Void _sub(ValidatedTaylorModel& r, const ValidatedTaylorModel& x, const ValidatedTaylorModel& y)
{
    // Compute r=x+y, assuming r is empty
    set_rounding_upward();
    Float te=0.0;
    ValidatedTaylorModel::ConstIterator xiter=x.begin();
    ValidatedTaylorModel::ConstIterator yiter=y.begin();
    while(xiter!=x.end() && yiter!=y.end()) {
        if(xiter->key()<yiter->key()) {
            r.expansion().raw().append(xiter->key(),xiter->data().raw());
            ++xiter;
        } else if(yiter->key()<xiter->key()) {
            r.expansion().raw().append(yiter->key(),-yiter->data().raw());
            ++yiter;
        } else {
            set_rounding_upward();
            const Float& xv=xiter->data().raw();
            const Float& yv=yiter->data().raw();
            VOLATILE Float u=xv-yv;
            VOLATILE Float t=-xv;
            VOLATILE Float ml=t+yv;
            te+=(u+ml);
            set_rounding_to_nearest();
            r.expansion().raw().append(xiter->key(),xiter->data().raw()-yiter->data().raw());
            ++xiter; ++yiter;
        }
    }

    while(xiter!=x.end()) {
        r.expansion().raw().append(xiter->key(),xiter->data().raw());
        ++xiter;
    }
    while(yiter!=y.end()) {
        r.expansion().raw().append(yiter->key(),-yiter->data().raw());
        ++yiter;
    }

    set_rounding_upward();
    r.error().raw()=x.error().raw();
    r.error().raw()+=y.error().raw();
    r.error().raw()+=(te/2);
    set_rounding_to_nearest();

    ARIADNE_ASSERT(r.error().raw()>=0);
}

inline Void _sma(ValidatedTaylorModel& r, const ValidatedTaylorModel& x, const Float& c, const ValidatedTaylorModel& y)
{
    if(c==1.0) { _add(r,x,y); return; }
    if(c==-1.0) { _sub(r,x,y); return; }

    // Compute r=x+y, assuming r is empty
    set_rounding_upward();
    Float te=0.0;
    VOLATILE Float u,ml;
    Float mc=-c;

    ValidatedTaylorModel::ConstIterator xiter=x.begin();
    ValidatedTaylorModel::ConstIterator yiter=y.begin();
    while(xiter!=x.end() && yiter!=y.end()) {
        if(xiter->key()<yiter->key()) {
            r.expansion().raw().append(xiter->key(),xiter->data().raw());
            ++xiter;
        } else if(yiter->key()<xiter->key()) {
            set_rounding_upward();
            const Float& yv=yiter->data().raw();
            u=c*yv;
            ml=mc*yv;
            te+=(u+ml);
            set_rounding_to_nearest();
            r.expansion().raw().append(yiter->key(),c*yiter->data().raw());
            ++yiter;
        } else {
            set_rounding_upward();
            const Float& xv=xiter->data().raw();
            const Float& yv=yiter->data().raw();
            u=xv+c*yv;
            ml=mc*yv-xv;
            te+=(u+ml);
            set_rounding_to_nearest();
            r.expansion().raw().append(xiter->key(),xiter->data().raw()+c*yiter->data().raw());
            ++xiter; ++yiter;
        }
    }

    while(xiter!=x.end()) {
        r.expansion().raw().append(xiter->key(),xiter->data().raw());
        ++xiter;
    }
    while(yiter!=y.end()) {
        set_rounding_upward();
        const Float& yv=yiter->data().raw();
        VOLATILE Float u=c*yv;
        VOLATILE Float ml=mc*yv;
        te+=(u+ml);
        set_rounding_to_nearest();
        r.expansion().raw().append(yiter->key(),c*yiter->data().raw());
        ++yiter;
    }

    set_rounding_upward();
    r.error().raw()=x.error().raw();
    r.error().raw()+=y.error().raw();
    r.error().raw()+=(te/2);
    set_rounding_to_nearest();

    ARIADNE_ASSERT(r.error().raw()>=0);
}

inline Void _sma(ValidatedTaylorModel& r, const ValidatedTaylorModel& x, const ValidatedNumber& c, const ValidatedTaylorModel& y)
{
    if(c.lower().raw()==c.upper().raw()) {
        _sma(r,x,c.lower().raw(),y); return;
    }

    ARIADNE_ASSERT_MSG(c.lower().raw()<=c.upper().raw(),c);
    ARIADNE_ASSERT_MSG(x.error().raw()>=0,"x="<<x);
    ARIADNE_ASSERT_MSG(y.error().raw()>=0,"y="<<y);

    //std::cerr<<"ValidatedTaylorModel::scal(ValidatedNumber c) c="<<c<<std::endl;
    set_rounding_upward();
    VOLATILE Float u,ml,myv;
    Float te=0; // Twice the maximum accumulated error
    Float cu=c.upper().raw();
    Float mcl=-c.lower().raw();
    Float cm=c.midpoint().raw();

    // Compute r=x+y, assuming r is empty
    set_rounding_upward();
    ValidatedTaylorModel::ConstIterator xiter=x.begin();
    ValidatedTaylorModel::ConstIterator yiter=y.begin();
    while(xiter!=x.end() && yiter!=y.end()) {
        if(xiter->key()<yiter->key()) {
            r.expansion().raw().append(xiter->key(),xiter->data().raw());
            ++xiter;
        } else if(yiter->key()<xiter->key()) {
            set_rounding_upward();
            const Float& yv=yiter->data().raw();
            if(yv>=0.0) {
                u=cu*yv;
                ml=mcl*yv;
            } else {
                myv = -yv;
                u = mcl*myv;
                ml = cu*myv;
            }
            te+=(u+ml);
            set_rounding_to_nearest();
            r.expansion().raw().append(yiter->key(),cm*yiter->data().raw());
            ++yiter;
        } else {
            set_rounding_upward();
            const Float& xv=xiter->data().raw();
            const Float& yv=yiter->data().raw();
            if(yv>0.0) {
                u=cu*yv+xv;
                ml=mcl*yv-xv;
            } else {
                myv=-yv;
                u=mcl*myv+xv;
                ml=cu*myv-xv;
            }
            te+=(u+ml);
            set_rounding_to_nearest();
            r.expansion().raw().append(xiter->key(),xiter->data().raw()+cm*yiter->data().raw());
            ++xiter; ++yiter;
        }
    }

    while(xiter!=x.end()) {
        r.expansion().raw().append(xiter->key(),xiter->data().raw());
        ++xiter;
    }
    while(yiter!=y.end()) {
        set_rounding_upward();
        const Float& yv=yiter->data().raw();
        if(yv>0.0) {
            u=cu*yv;
            ml=mcl*yv;
        } else {
            myv = -yv;
            u=mcl*myv;
            ml=cu*myv;
        }
        te+=(u+ml);
        set_rounding_to_nearest();
        r.expansion().raw().append(yiter->key(),cm*yiter->data().raw());
        ++yiter;
    }

    set_rounding_upward();
    r.error().raw()=x.error().raw();
    r.error().raw()+=y.error().raw();
    r.error().raw()+=(te/2);
    set_rounding_to_nearest();

    ARIADNE_ASSERT(r.error().raw()>=0);
}

struct Ivl { Float u; Float ml; };

/*
inline Void _mul1(ValidatedTaylorModel& r, const ValidatedTaylorModel& x, const ValidatedTaylorModel& y)
{
    // Compute r+=x*y
    typedef ValidatedTaylorModel::ConstIterator ConstIterator;
    typedef std::map<MultiIndex,Ivl>::ConstIterator ivl_const_iterator;
    Float& re=r.error().raw();
    std::map<MultiIndex,Ivl> z;

    set_rounding_upward();
    for(ConstIterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
        for(ConstIterator yiter=y.begin(); yiter!=y.end(); ++yiter) {
            const Float& xv=xiter->data().raw();;
            const Float& yv=yiter->data().raw();;
            Ivl& zv=z[xiter->key()+yiter->key()];
            zv.u+=xv*yv;
            VOLATILE Float t=-xv;
            zv.ml+=t*yv;
        }
    }

    for(ConstIterator riter=r.begin(); riter!=r.end(); ++riter) {
        Ivl& zv=z[riter->key()];
        const Float& rv=riter->data().raw();
        zv.u+=rv; zv.ml-=rv;
    }

    VOLATILE Float te=0;
    for(ivl_const_iterator ziter=z.begin(); ziter!=z.end(); ++ziter) {
        const Ivl& zv=ziter->second;
        te+=(zv.u+zv.ml);
    }
    te/=2;

    Float xs=0;
    for(ConstIterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
        xs+=abs(xiter->data().raw());
    }

    Float ys=0;
    for(ConstIterator yiter=y.begin(); yiter!=y.end(); ++yiter) {
        ys+=abs(yiter->data().raw());
    }

    const Float& xe=x.error().raw();
    const Float& ye=y.error().raw();

    re+=xs*ye+ys*xe+te+xe*ye;

    set_rounding_to_nearest();
    for(ivl_const_iterator ziter=z.begin(); ziter!=z.end(); ++ziter) {
        const Ivl& zv=ziter->second;
        r[ziter->first]=(zv.u-zv.ml)/2;
    }

    return;
}
*/


// Compute r+=x*y
// Compute monomial-by-monomial in y
// Avoid changing rounding mode
inline
Void _mul2(ValidatedTaylorModel& r, const ValidatedTaylorModel& x, const ValidatedTaylorModel& y)
{
    ValidatedTaylorModel t(x.argument_size(),r.sweeper());
    ValidatedTaylorModel s(x.argument_size(),r.sweeper());
    for(ValidatedTaylorModel::ConstIterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
        set_rounding_upward();
        Float te=0.0;
        for(ValidatedTaylorModel::ConstIterator yiter=y.begin(); yiter!=y.end(); ++yiter) {
            const Float& xv=xiter->data().raw();
            const Float& yv=yiter->data().raw();
            VOLATILE Float u=xv*yv;
            VOLATILE Float ml=-xv; ml=ml*yv;
            te+=(u+ml);
        }
        t.error().raw()=te/2;
        set_rounding_to_nearest();
        for(ValidatedTaylorModel::ConstIterator yiter=y.begin(); yiter!=y.end(); ++yiter) {
            t.expansion().raw().append(xiter->key(),yiter->key(),xiter->data().raw()*yiter->data().raw());
        }
        _add2(s,r,t);
        r.expansion().swap(s.expansion());
        r.error().raw()=s.error().raw();
        s.expansion().clear();
        s.error().raw()=0.0;
        t.expansion().clear();
        t.error().raw()=0.0;
    }

    set_rounding_upward();
    Float xs=0;
    for(ValidatedTaylorModel::ConstIterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
        xs+=abs(xiter->data().raw());
    }

    Float ys=0;
    for(ValidatedTaylorModel::ConstIterator yiter=y.begin(); yiter!=y.end(); ++yiter) {
        ys+=abs(yiter->data().raw());
    }

    Float& re=r.error().raw();
    const Float& xe=x.error().raw();
    const Float& ye=y.error().raw();
    re+=xs*ye+ys*xe+xe*ye;

    r.sweep();

    set_rounding_to_nearest();
    return;
}


// Compute r+=x*y
// Compute monomial-by-monomial in y
// Change the rounding mode to avoid iterating using opposite rounding
inline Void _mul3(ValidatedTaylorModel& r, const ValidatedTaylorModel& x, const ValidatedTaylorModel& y)
{
    ValidatedTaylorModel t(x.argument_size(),r.sweeper());
    for(ValidatedTaylorModel::ConstIterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
        VOLATILE Float pxv=xiter->data().raw();
        VOLATILE Float nxv=-pxv;
        VOLATILE Float te=0.0;
        for(ValidatedTaylorModel::ConstIterator yiter=y.begin(); yiter!=y.end(); ++yiter) {
            set_rounding_upward();
            const Float& yv=yiter->data().raw();
            te+=(pxv*yv)+(nxv*yv);
            set_rounding_to_nearest();
            t.expansion().raw().append(xiter->key(),yiter->key(),xiter->data().raw()*yiter->data().raw());
        }
        t.error().raw()=te/2;
        r+=t;
        //std::cerr<<"  t="<<t<<"\n  r="<<r<<std::endl;
        t.clear();
    }

    set_rounding_upward();
    Float xs=0;
    for(ValidatedTaylorModel::ConstIterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
        xs+=abs(xiter->data().raw());
    }

    Float ys=0;
    for(ValidatedTaylorModel::ConstIterator yiter=y.begin(); yiter!=y.end(); ++yiter) {
        ys+=abs(yiter->data().raw());
    }

    Float& re=r.error().raw();
    const Float& xe=x.error().raw();
    const Float& ye=y.error().raw();
    re+=xs*ye+ys*xe+xe*ye;

    r.sweep();

    set_rounding_to_nearest();
    return;
}

inline Void _add(MultiIndex& r, const MultiIndex& a1, const MultiIndex& a2) {
    for(MultiIndex::SizeType j=0; j!=r.word_size(); ++j) {
        r.word_at(j)=a1.word_at(j)+a2.word_at(j);
    }
}


// Compute r+=x*y
// Compute monomial-by-monomial in y
// Avoid changing rounding mode
inline
Void _mul4(ValidatedTaylorModel& r, const ValidatedTaylorModel& x, const ValidatedTaylorModel& y, const Sweeper& accuracy)
{
    const Nat as=r.argument_size();
    ValidatedTaylorModel t(as,r.sweeper());
    ValidatedTaylorModel s(as,r.sweeper());
    MultiIndex ta(as);
    Float tv;
    VOLATILE Float u;
    VOLATILE Float ml;
    for(ValidatedTaylorModel::ConstIterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
        Float tte=0.0; // trucation error
        Float tre=0.0; // roundoff error
        const MultiIndex& xa=xiter->key();
        VOLATILE Float& xv=const_cast<Float&>(static_cast<const Float&>(xiter->data().raw()));
        VOLATILE Float mxv=-xv;
        VOLATILE Float axv=abs(xv);
        for(ValidatedTaylorModel::ConstIterator yiter=y.begin(); yiter!=y.end(); ++yiter) {
            const MultiIndex& ya=yiter->key();
            VOLATILE Float& yv=const_cast<Float&>(static_cast<const Float&>(yiter->data().raw()));
            set_rounding_to_nearest();
            ta=xa+ya;
            tv=xv*yv;
            set_rounding_upward();
            if(accuracy.discard(ta,tv)) {
                tte+=axv*abs(yv);
            } else {
                t.expansion().raw().append(ta,tv);
                u=xv*yv;
                ml=mxv*yv;
                tre+=(u+ml);
            }
        }
        t.error().raw()=tte+(tte/2);

        _add2(s,r,t);
        r.expansion().swap(s.expansion());
        r.error().raw()=s.error().raw();
        s.expansion().clear();
        s.error().raw()=0.0;
        t.expansion().clear();
        t.error().raw()=0.0;
    }

    set_rounding_upward();
    Float xs=0;
    for(ValidatedTaylorModel::ConstIterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
        xs+=abs(xiter->data().raw());
    }

    Float ys=0;
    for(ValidatedTaylorModel::ConstIterator yiter=y.begin(); yiter!=y.end(); ++yiter) {
        ys+=abs(yiter->data().raw());
    }

    Float& re=r.error().raw();
    const Float& xe=x.error().raw();
    const Float& ye=y.error().raw();
    re+=xs*ye+ys*xe+xe*ye;

    set_rounding_to_nearest();
    return;
}

inline Void _mul(ValidatedTaylorModel& r, const ValidatedTaylorModel& x, const ValidatedTaylorModel& y) {
    //_mul2(r,x,y);
    _mul4(r,x,y,r.sweeper());
    ARIADNE_ASSERT(r.error().raw()>=0);
}



} // namespace

Void _mul_full(ValidatedTaylorModel& r, const ValidatedTaylorModel& x, const ValidatedTaylorModel& y)
{
    _mul2(r,x,y);
}

Void _mul_clear(ValidatedTaylorModel& r, const ValidatedTaylorModel& x, const ValidatedTaylorModel& y)
{
    _mul4(r,x,y,r.sweeper());
}


///////////////////////////////////////////////////////////////////////////////

// Inplace arithmetical operations for Algebra concept

Void ValidatedTaylorModel::iadd(const ValidatedNumber& c)
{
    _acc(*this,c);
    this->sweep();
    ARIADNE_ASSERT_MSG(this->error().raw()>=0,*this);
}

Void ValidatedTaylorModel::imul(const ValidatedNumber& c)
{
    _scal(*this,c);
    this->sweep();
    ARIADNE_ASSERT_MSG(this->error().raw()>=0,*this);
}

Void ValidatedTaylorModel::isma(const ValidatedNumber& c, const ValidatedTaylorModel& y)
{
    ValidatedTaylorModel& x=*this;
    ValidatedTaylorModel r=this->create();
    _sma(r,x,c,y);
    this->swap(r);
    this->sweep();
    ARIADNE_ASSERT_MSG(this->error().raw()>=0,*this);
}

Void ValidatedTaylorModel::ifma(const ValidatedTaylorModel& x1, const ValidatedTaylorModel& x2)
{
    _mul(*this,x1,x2);
    this->sweep();
    ARIADNE_ASSERT_MSG(this->error().raw()>=0,*this);
}



///////////////////////////////////////////////////////////////////////////////

// Truncation and error control


ValidatedTaylorModel&
ValidatedTaylorModel::unique_sort()
{
    this->_expansion.reverse_lexicographic_sort();

    ValidatedTaylorModel::ConstIterator advanced =this->begin();
    ValidatedTaylorModel::ConstIterator end =this->end();
    ValidatedTaylorModel::Iterator current=this->begin();
    Float te=0.0;
    while(advanced!=end) {
        current->key()=advanced->key();
        VOLATILE Float u=advanced->data().raw();
        VOLATILE Float ml=-advanced->data().raw();
        VOLATILE Float v=advanced->data().raw();
        ++advanced;
        while(advanced!=end && advanced->key()==current->key()) {
            const Float& xv=advanced->data().raw();
            set_rounding_upward();
            u+=xv;
            ml-=xv;
            set_rounding_to_nearest();
            v+=xv;
            ++advanced;
        }
        current->data().raw()=v;
        te+=(u+ml);
        ++current;
    }
    set_rounding_upward();
    this->error().raw()+=te;
    set_rounding_to_nearest();
    this->_expansion.resize(current-this->begin());

    return *this;
}


ValidatedTaylorModel&
ValidatedTaylorModel::sweep()
{
    this->_sweeper.sweep(this->_expansion.raw(),this->_error.raw());
    return *this;
}

ValidatedTaylorModel&
ValidatedTaylorModel::sweep(const SweeperInterface& sweeper)
{
    sweeper.sweep(this->_expansion.raw(),this->_error.raw());
    return *this;
}



ValidatedTaylorModel&
ValidatedTaylorModel::clobber()
{
    this->_error.raw()=0;
    return *this;
}



///////////////////////////////////////////////////////////////////////////////

// Accuracy control

Float
ValidatedTaylorModel::tolerance() const
{
    const ThresholdSweeper* ptr=dynamic_cast<const ThresholdSweeper*>(&static_cast<const SweeperInterface&>(this->_sweeper));
    if(ptr) {
        return ptr->sweep_threshold();
    } else {
        return std::numeric_limits<double>::epsilon();
    }
}



//////////////////////////////////////////////////////////////////////////////

// Exact functions (max, min, abs, neg) and arithmetical functions (sqr, pow)


ValidatedTaylorModel max(const ValidatedTaylorModel& x, const ValidatedTaylorModel& y) {
    UpperInterval xr=x.range();
    UpperInterval yr=y.range();
    if(xr.lower().raw()>=yr.upper().raw()) {
        return x;
    } else if(yr.lower().raw()>=xr.upper().raw()) {
        return y;
    } else {
        return ((x+y)+abs(x-y))/Float(2);
    }
}


ValidatedTaylorModel min(const ValidatedTaylorModel& x, const ValidatedTaylorModel& y) {
    return -max(-x,-y);
}

ValidatedTaylorModel abs(const ValidatedTaylorModel& x) {
    UpperInterval xr=x.range();
    if(xr.lower().raw()>=0.0) {
        return x;
    } else if(xr.upper().raw()<=0.0) {
        return -x;
    } else {
        // Use power series expansion $abs(x)=\sum_{i=0}^{7} p_i x^{2i} \pm e$ for $x\in[-1,+1]$ with
        // p=[0.0112167620474, 5.6963263292747541, -31.744583789655049, 100.43002481377681, -162.01366698662306, 127.45243493284417, -38.829743345344667] and e=0.035
        // TODO: Find more accurate and stable formula
        static const Nat n=7u;
        static const double p[n]={0.0112167620474, 5.6963263292747541, -31.744583789655049, 100.43002481377681, -162.01366698662306, 127.45243493284417, -38.829743345344667};
        static const double err=0.035;
        ValidatedTaylorModel r(x.argument_size(),x.sweeper());
        Float xmag=mag(xr).raw();
        ValidatedTaylorModel s=x/xmag;
        s=sqr(s);
        r=static_cast<ExactNumber>(p[n-1]);
        for(Nat i=0; i!=(n-1); ++i) {
            Nat j=(n-2)-i;
            r=s*r+p[j];
        }
        r+=ValidatedFloat(-err,+err);
        return r*xmag;
    }
}

ValidatedTaylorModel neg(const ValidatedTaylorModel& x) {
    ValidatedTaylorModel r(x.argument_size(),x.sweeper());
    for(ValidatedTaylorModel::ConstIterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
        r.expansion().raw().append(xiter->key(),-xiter->data().raw());
    }
    r.error().raw()=x.error().raw();
    return r;
}

//////////////////////////////////////////////////////////////////////////////

// Arithmetical functions (sqr, pow)

ValidatedTaylorModel sqr(const ValidatedTaylorModel& x) {
    ValidatedTaylorModel r=x*x;
    return r;
}

ValidatedTaylorModel pow(const ValidatedTaylorModel& x, Int n) {
    ValidatedTaylorModel r(x.argument_size(),x.sweeper()); r+=1;
    ValidatedTaylorModel p(x);
    while(n) {
        if(n%2) { r=r*p; }
        p=sqr(p);
        n/=2;
    }
    return r;
}



//////////////////////////////////////////////////////////////////////////////

// Basic function operators (domain, range, evaluate)

/* FIXME: The following code does not work with optimisation turned on using gcc */
UpperInterval _range1(const ValidatedTaylorModel& tm) {
    set_rounding_mode(upward);
    VOLATILE Float t=tm.error().raw();
    VOLATILE Float v=0.0;
    for(ValidatedTaylorModel::ConstIterator iter=tm.begin(); iter!=tm.end(); ++iter) {
        if(iter->key().degree()==0) {
            v=iter->data().raw();
        } else {
            t+=abs(iter->data().raw());
        }
    }
    set_rounding_mode(to_nearest);
    return v+UpperInterval(-t,t);
}


UpperInterval _range2(const ValidatedTaylorModel& tm) {
    UpperInterval r(-tm.error().raw(),+tm.error().raw());
    for(ValidatedTaylorModel::ConstIterator iter=tm.begin(); iter!=tm.end(); ++iter) {
        if(iter->key().degree()==0) {
            r+=iter->data().raw();
        } else {
            r+=abs(iter->data().raw())*UpperInterval(-1,1);
        }
    }
    return r;
}

// Compute the range by grouping all quadratic terms x[i]^2 with linear terms x[i]
// The range of ax^2+bx+c is a([-1,1]+b/2a)^2+(c-b^2/4a)
UpperInterval _range3(const ValidatedTaylorModel& tm) {
    const Nat as=tm.argument_size();
    Array<Float> linear_terms(as,0.0);
    Array<Float> quadratic_terms(as,0.0);
    UpperInterval r(-tm.error().raw(),+tm.error().raw());
    for(ValidatedTaylorModel::ConstIterator iter=tm.begin(); iter!=tm.end(); ++iter) {
        if(iter->key().degree()==0) {
            r+=iter->data().raw();
        } else if(iter->key().degree()==1) {
            for(Nat j=0; j!=tm.argument_size(); ++j) {
                if(iter->key()[j]==1) { linear_terms[j]=iter->data().raw(); break; }
            }
        } else if(iter->key().degree()==2) {
            for(Nat j=0; j!=tm.argument_size(); ++j) {
                if(iter->key()[j]==2) { quadratic_terms[j]=iter->data().raw(); break; }
                if(iter->key()[j]==1) { r+=abs(iter->data().raw())*UpperInterval(-1,1); break; }
            }
        } else {
            r+=abs(iter->data().raw())*UpperInterval(-1,1);
        }
    }
    // If the ratio b/a is very large, then roundoff error can cause a significant
    // additional error. We compute both |a|+|b| and a([-1,+1]+b/2a)-b^2/4a and take best bound
    for(Nat j=0; j!=as; ++j) {
        const Float& a=quadratic_terms[j];
        const Float& b=linear_terms[j];
        UpperInterval ql=abs(a)*UpperInterval(-1,1) + abs(b)*UpperInterval(-1,+1);
        UpperInterval qf=a*(sqr(UpperInterval(-1,+1)+UpperInterval(b)/(2*a)))-sqr(UpperInterval(b))/(4*a);
        r += intersection(ql,qf); // NOTE: ql must be the first term in case of NaN in qf
    }
    return r;
}



UpperInterval
ValidatedTaylorModel::range() const {
    return Ariadne::_range3(*this);
}


Vector<UnitInterval>
ValidatedTaylorModel::domain() const
{
    return Vector<UnitInterval>(this->argument_size(),UnitInterval());
}

ExactInterval
ValidatedTaylorModel::codomain() const
{
    return make_exact_interval(this->range());
}



//////////////////////////////////////////////////////////////////////////////

// Composition with power series

template<class X> class Series;
class TaylorSeries;


ValidatedTaylorModel
_compose(const TaylorSeries& ts, const ValidatedTaylorModel& tv, double eps)
{
    Sweeper threshold_sweeper(new ThresholdSweeper(eps));
    //std::cerr<<"_compose(TaylorSeries,ValidatedTaylorModel,Error)\n";
    //std::cerr<<"\n  ts="<<ts<<"\n  tv="<<tv<<"\n";
    Float& vref=const_cast<Float&>(tv.value().raw());
    Float vtmp=vref;
    vref=0.0;
    ValidatedTaylorModel r(tv.argument_size(),tv.sweeper());
    r+=ts.expansion[ts.expansion.size()-1];
    for(Nat i=1; i!=ts.expansion.size(); ++i) {
        //std::cerr<<"    r="<<r<<std::endl;
        r=r*tv;
        r+=ts.expansion[ts.expansion.size()-i-1];
        r.sweep(threshold_sweeper);
    }
    //std::cerr<<"    r="<<r<<std::endl;
    r+=ts.error;
    //std::cerr<<"    r="<<r<<std::endl;
    vref=vtmp;
    return r;
}

ValidatedTaylorModel
compose(const TaylorSeries& ts, const ValidatedTaylorModel& tm)
{
    return _compose(ts,tm,ec);
}


// Compose using the Taylor formula directly. The final term is the Taylor series computed
// over the range of the series. This method tends to suffer from blow-up of the
// truncation error
ValidatedTaylorModel
_compose1(const ValidatedSeriesFunctionPointer& fn, const ValidatedTaylorModel& tm, double eps)
{
    Sweeper threshold_sweeper(new ThresholdSweeper(eps));
    static const Nat DEGREE=18;
    static const double TRUNCATION_ERROR=1e-8;
    Nat d=DEGREE;
    ExactFloat c=tm.value();
    UpperInterval r=tm.range();
    Series<ValidatedFloat> centre_series=fn(d,c);
    Series<ValidatedFloat> range_series=fn(d,static_cast<ValidatedFloat>(r));

    Float truncation_error_estimate=(mag(range_series[d])*pow(mag(r-c),d)).raw();
    if(truncation_error_estimate>TRUNCATION_ERROR) {
        ARIADNE_WARN("Truncation error estimate "<<truncation_error_estimate
                     <<" is greater than maximum allowable truncation error "<<TRUNCATION_ERROR<<"\n");
    }

    ValidatedTaylorModel x=tm-c;
    ValidatedTaylorModel res(tm.argument_size(),tm.sweeper());
    res+=range_series[d];
    for(Nat i=0; i!=d; ++i) {
        //std::cerr<<"i="<<i<<" r="<<res<<"\n";
        res=centre_series[d-i-1]+x*res;
        res.sweep(threshold_sweeper);
    }
    //std::cerr<<"i="<<d<<" r="<<res<<"\n";
    return res;
}

// Compose using the Taylor formula with a constant truncation error. This method
// is usually better than _compose1 since there is no blow-up of the trunction
// error. The radius of convergence of this method is still quite low,
// typically only half of the radius of convergence of the power series itself
ValidatedTaylorModel
_compose2(const ValidatedSeriesFunctionPointer& fn, const ValidatedTaylorModel& tm, double eps)
{
    Sweeper threshold_sweeper(new ThresholdSweeper(eps));
    static const Nat DEGREE=20;
    static const Float TRUNCATION_ERROR=1e-8;
    Nat d=DEGREE;
    ExactFloat c=tm.value();
    UpperInterval r=tm.range();
    Series<ValidatedFloat> centre_series=fn(d,c);
    Series<ValidatedFloat> range_series=fn(d,static_cast<ValidatedFloat>(r));

    //std::cerr<<"c="<<c<<" r="<<r<<" r-c="<<r-c<<" e="<<mag(r-c)<<"\n";
    //std::cerr<<"cs[d]="<<centre_series[d]<<" rs[d]="<<range_series[d]<<"\n";
    //std::cerr<<"cs="<<centre_series<<"\nrs="<<range_series<<"\n";
    ErrorFloat truncation_error=mag(range_series[d]-centre_series[d])*pow(mag(r-c),d);
    //std::cerr<<"te="<<truncation_error<<"\n";
    if(truncation_error.raw()>TRUNCATION_ERROR) {
        ARIADNE_WARN("Truncation error estimate "<<truncation_error
                 <<" is greater than maximum allowable truncation error "<<TRUNCATION_ERROR<<"\n");
    }

    ValidatedTaylorModel x=tm-c;
    ValidatedTaylorModel res(tm.argument_size(),tm.sweeper());
    res+=centre_series[d];
    for(Nat i=0; i!=d; ++i) {
        res=centre_series[d-i-1]+x*res;
        res.sweep(threshold_sweeper);
    }
    res+=ValidatedFloat(-truncation_error,+truncation_error);
    return res;
}


// Compose using the Taylor formula with a constant truncation error. This method
// is usually better than _compose1 since there is no blow-up of the trunction
// error. This method is better than _compose2 since the truncation error is
// assumed at the ends of the intervals
ValidatedTaylorModel
_compose3(const ValidatedSeriesFunctionPointer& fn, const ValidatedTaylorModel& tm, Float eps)
{
    static const Nat DEGREE=20;
    static const Float TRUNCATION_ERROR=1e-8;
    Nat d=DEGREE;
    ExactFloat c=tm.value();
    UpperInterval r=tm.range();
    Series<ValidatedFloat> centre_series=fn(d,c);
    Series<ValidatedFloat> range_series=fn(d,static_cast<ValidatedFloat>(r));

    //std::cerr<<"c="<<c<<" r="<<r<<" r-c="<<r-c<<" e="<<mag(r-c)<<"\n";
    //std::cerr<<"cs[d]="<<centre_series[d]<<" rs[d]="<<range_series[d]<<"\n";
    //std::cerr<<"cs="<<centre_series<<"\nrs="<<range_series<<"\n";
    ValidatedFloat se=range_series[d]-centre_series[d];
    ValidatedFloat e=static_cast<ValidatedFloat>(r)-c;
    ValidatedFloat p=pow(e,d-1);
    p=ValidatedFloat(-p.lower().raw()*e.lower().raw(),p.upper().raw()*e.upper().raw());
    //std::cerr<<"se="<<se<<" e="<<e<<" p="<<p<<std::endl;
    // FIXME: Here we assume the dth derivative of f is monotone increasing
    Float truncation_error=max(se.lower().raw()*p.lower().raw(),se.upper().raw()*p.upper().raw());
    //std::cerr<<"te="<<truncation_error<<"\n";
    if(truncation_error>TRUNCATION_ERROR) {
        ARIADNE_WARN("Truncation error estimate "<<truncation_error
                 <<" is greater than maximum allowable truncation error "<<TRUNCATION_ERROR<<"\n");
    }

    ValidatedTaylorModel x=tm;
    ValidatedTaylorModel res(tm.argument_size(),tm.sweeper());
    res+=centre_series[d];
    for(Nat i=0; i!=d; ++i) {
        res=centre_series[d-i-1]+x*res;
        //res.sweep(eps);
    }
    res+=ValidatedFloat(-truncation_error,+truncation_error);
    return res;
}


ValidatedTaylorModel
_compose(const ValidatedSeriesFunctionPointer& fn, const ValidatedTaylorModel& tm, Float eps) {
    //std::cerr<<"_compose(SeriesFunction,ValidatedTaylorModel,Error)\n";
    return _compose3(fn,tm,eps);
}




///////////////////////////////////////////////////////////////////////////////

// Algebraic and trancendental functions
//   bounded domain (rec,sqrt,log,tan)
//   unbounded domain (exp,sin,cos)


ValidatedTaylorModel sqrt(const ValidatedTaylorModel& x) {
    //std::cerr<<"rec(ValidatedTaylorModel)\n";
    // Use a special routine to minimise errors
    // Given range [rl,ru], rescale by constant a such that rl/a=1-d; ru/a=1+d
    UpperInterval r=x.range();

    if(r.lower().raw()<=0) {
        ARIADNE_THROW(DomainException,"sqrt",x.range());
    }

    assert(r.lower().raw()>0);
    Float a=(r.lower().raw()+r.upper().raw())/2;
    set_rounding_upward();
    Float eps=(r.upper().raw()-r.lower().raw())/(r.upper().raw()+r.lower().raw());
    set_rounding_to_nearest();
    assert(eps<1);
    Nat d=integer_cast<Int>((log((1-eps)*x.tolerance())/log(eps)+1));
    //std::cerr<<"x="<<x<<std::endl;
    //std::cerr<<"x/a="<<x/a<<" a="<<a<<std::endl;
    ValidatedTaylorModel y=(x/a)-1.0;
    //std::cerr<<"y="<<y<<std::endl;
    ValidatedTaylorModel z(x.argument_size(),x.sweeper());
    Series<ValidatedFloat> sqrt_series=Series<ValidatedFloat>::sqrt(d,ValidatedFloat(1));
    //std::cerr<<"sqrt_series="<<sqrt_series<<std::endl;
    //std::cerr<<"y="<<y<<std::endl;
    z+=sqrt_series[d-1];
    for(Nat i=0; i!=d; ++i) {
        z=sqrt_series[d-i-1] + z * y;
        z.sweep();
        //std::cerr<<"z="<<z<<std::endl;
    }
    Float trunc_err=(pow(eps,d)/(1-eps))*mag(sqrt_series[d]).raw();
    //std::cerr<<"te="<<trunc_err<<" te*[-1,+1]="<<trunc_err*ValidatedFloat(-1,1)<<std::endl;
    z.error().raw()+=trunc_err;
    //std::cerr<<"z="<<z<<std::endl;
    ValidatedFloat sqrta=sqrt(ValidatedFloat(a));
    //std::cerr<<"sqrt(a)="<<sqrta<<std::endl;
    z*=sqrt(ValidatedFloat(a));
    //std::cerr<<"z="<<z<<std::endl;
    return z;
}

ValidatedTaylorModel rec(const ValidatedTaylorModel& x) {
    //std::cerr<<"rec(ValidatedTaylorModel)\n";
    // Use a special routine to minimise errors
    // Given range [rl,ru], rescale by constant a such that rl/a=1-d; ru/a=1+d
    ValidatedFloat r=static_cast<ValidatedFloat>(x.range());
    if(r.upper().raw()>=0 && r.lower().raw()<=0) {
        ARIADNE_THROW(DivideByZeroException,"rec(ValidatedTaylorModel x)","x="<<x<<", x.range()="<<x.range());
    }
    Float a=(r.lower().raw()+r.upper().raw())/2;
    set_rounding_upward();
    Float eps=abs((r.upper().raw()-r.lower().raw())/(r.upper().raw()+r.lower().raw()));
    set_rounding_to_nearest();
    assert(eps<1);

    Nat d=integer_cast<Nat>((log((1-eps)*x.tolerance())/log(eps))+1);

    ValidatedTaylorModel y=1-(x/a);
    ValidatedTaylorModel z(x.argument_size(),x.sweeper());
    z+=Float(d%2?-1:+1);
    //std::cerr<<"  y="<<y<<"\n";
    //std::cerr<<"  z="<<z<<"\n";
    for(Nat i=0; i!=d; ++i) {
        z=1.0 + z * y;
        //std::cerr<<"  z="<<z<<"\n";
    }

    // Compute the truncation error
    Float te=pow(eps,d)/(1-eps);
    //std::cerr<<"  te="<<te<<"\n";
    set_rounding_upward();
    Float nze=te+z.error().raw();
    //std::cerr<<"  nze="<<nze<<"\n";
    set_rounding_to_nearest();
    z.set_error(nze);
    //std::cerr<<"  z="<<z<<"\n";
    z/=a;
    //std::cerr<<"  z="<<z<<"\n";
    return z;
}

ValidatedTaylorModel log(const ValidatedTaylorModel& x) {
    // Use a special routine to minimise errors
    // Given range [rl,ru], rescale by constant a such that rl/a=1-d; ru/a=1+d
    UpperInterval r=x.range();
    if(r.lower().raw()<=0) {
        ARIADNE_THROW(DomainException,"sqrt",x.range());
    }
    Float a=(r.lower().raw()+r.upper().raw())/2;
    set_rounding_upward();
    Float eps=(r.upper().raw()-r.lower().raw())/(r.upper().raw()+r.lower().raw());
    set_rounding_to_nearest();
    assert(eps<1);
    Nat d=integer_cast<Nat>((log((1-eps)*x.tolerance())/log(eps)+1));
    ValidatedTaylorModel y=x/a-1;
    ValidatedTaylorModel z(x.argument_size(),x.sweeper());
    z+=Float(d%2?-1:+1)/d;
    for(Nat i=1; i!=d; ++i) {
        z=Float((d-i)%2?+1:-1)/(d-i) + z * y;
        z.sweep();
    }
    z=z*y;
    z.sweep();
    Float trunc_err=pow(eps,d)/(1-eps)/d;
    return z+log(ValidatedFloat(a))+ValidatedFloat(-trunc_err,+trunc_err);
}

static const Nat MAXIMUM_DEGREE = 18;

// Use special code to utilise exp(ax+b)=exp(x)^a*exp(b)
ValidatedTaylorModel exp(const ValidatedTaylorModel& x) {
    // FIXME: Truncation error may be incorrect

    // Scale to unit interval
    ValidatedTaylorModel y=x;
    Float xval=x.value().raw();
    y-=xval;
    Float xrad=mag(y.range()).raw();
    Nat sfp=0; // A number such that 2^sfp>rad(x.range())
    while(Float(1<<sfp)<xrad) { ++sfp; }
    Float sf=1.0/(1<<sfp);
    _scal_exact(y,sf);
    Float yrad=xrad*sf;

    Nat degree=MAXIMUM_DEGREE;

    // Since x is in unit domain, truncation error is no worse than maximum omitted term, i.e. xr/fac(d+1)
    ValidatedTaylorModel res(x.argument_size(),x.sweeper());
    res.set_error(pow_up(yrad,degree+1));
    for(Nat i=0; i!=degree; ++i) {
        res/=(degree-i);
        res=y*res+1.0;
    }

    // Square r a total of sfp times
    ValidatedTaylorModel square(x.argument_size(),x.sweeper());
    for(Nat i=0; i!=sfp; ++i) {
        _mul(square,res,res);
        res.swap(square);
        square.clear();
    }

    // Multiply by exp(xv)
    res*=Ariadne::exp(ValidatedFloat(xval));

    return res;
    //return _compose(&Series<ValidatedFloat>::exp,x,x.tolerance());
}

// Use special code to utilise sin(x+2pi)=sin(x)
// and that the power series is of the form x*f(x^2)
ValidatedTaylorModel sin(const ValidatedTaylorModel& x) {
    // FIXME: Truncation error may be incorrect
    ValidatedTaylorModel y(x.argument_size(),x.sweeper());
    ValidatedTaylorModel s(x.argument_size(),x.sweeper());
    ValidatedTaylorModel r(x.argument_size(),x.sweeper());
    ValidatedTaylorModel t(x.argument_size(),x.sweeper());

    ExactNumber one(1.0);

    Float two_pi_approx=2*pi_approx;
    Int n=integer_cast<Int>(floor(x.value().raw()/two_pi_approx + 0.5));
    y=x-(n*2*pi_val);

    if(y.error().raw()>two_pi_approx/2) {
        r.error().raw()=one.raw();
    } else {
        _mul(s,y,y);

        Int d=(MAXIMUM_DEGREE+3)/2;
        Float srad=mag(s.range()).raw();
        Float truncation_error=pow_up(srad,d+1)*rec_fac_up((d+1)*2);

        // Compute x(1-y/6+y^2/120-y^3/5040+... = x(1-y/6*(1-y/20*(1-y/42*...)
        r=one;
        for(Int i=0; i!=d; ++i) {
            r/=Float(-2*(d-i)*(2*(d-i)+1));
            _mul(t,s,r); r.swap(t); t.clear();
            r+=one;
        }
        _mul(t,y,r); r.swap(t);

        r.error().raw()+=truncation_error;
    }

    return r;
}

// Use special code to utilise sin(x+2pi)=sin(x)
// and that the power series is of the form f(x^2)
ValidatedTaylorModel cos(const ValidatedTaylorModel& x) {
    // FIXME: Truncation error may be incorrect
    ValidatedTaylorModel y(x.argument_size(),x.sweeper());
    ValidatedTaylorModel s(x.argument_size(),x.sweeper());
    ValidatedTaylorModel r(x.argument_size(),x.sweeper());
    ValidatedTaylorModel t(x.argument_size(),x.sweeper());

    Float two_pi=2*pi_approx;
    Int n=integer_cast<Int>(floor(x.value().raw()/two_pi + 0.5));

    y=x-2*n*pi_val;
    ExactNumber one(1.0);

    if(y.error().raw()>two_pi/2) {
        r.error().raw()=one.raw();
    } else {
        _mul(s,y,y);

        Int d=(MAXIMUM_DEGREE+3)/2;
        Float srad=mag(s.range()).raw();
        Float truncation_error=pow_up(srad,d+1)*rec_fac_up((d+1)*2);

        // Compute 1-y/2+y^2/24-y^3/720+... = (1-y/2*(1-y/12*(1-y/30*...)
        r=one;
        for(Int i=0; i!=d; ++i) {
            r/=double(-2*(d-i)*(2*(d-i)-1));
            _mul(t,s,r); r.swap(t); t.clear();
            r+=one;
        }

        r.error().raw()+=truncation_error;
    }

    return r;
}

ValidatedTaylorModel tan(const ValidatedTaylorModel& x) {
    return sin(x)*rec(cos(x));
}

ValidatedTaylorModel asin(const ValidatedTaylorModel& x) {
    static const Nat DEG=18;
    return compose(TaylorSeries(DEG,&Series<ValidatedNumber>::asin,
                                x.value(),static_cast<ValidatedNumber>(x.range())),x);
}

ValidatedTaylorModel acos(const ValidatedTaylorModel& x) {
    static const Nat DEG=18;
    return compose(TaylorSeries(DEG,&Series<ValidatedNumber>::acos,
                                x.value(),static_cast<ValidatedNumber>(x.range())),x);
}

ValidatedTaylorModel atan(const ValidatedTaylorModel& x) {
    static const Nat DEG=18;
    return compose(TaylorSeries(DEG,&Series<ValidatedNumber>::atan,
                                x.value(),static_cast<ValidatedNumber>(x.range())),x);
}



///////////////////////////////////////////////////////////////////////////////

// Inplace function operators (rescale, restrict, antidifferentiate)


ValidatedTaylorModel& ValidatedTaylorModel::rescale(const ExactInterval& ocd, const ExactInterval& ncd)
{
    ValidatedTaylorModel& x=*this;

    const ExactFloat a=ocd.lower(); const ExactFloat b=ocd.upper();
    const ExactFloat c=ncd.lower(); const ExactFloat d=ncd.upper();

    // Scale the interval [a,b] onto [c,d]
    // The function is given by x:-> alpha*x+beta where
    // alpha=(d-c)/(b-a) and beta=(cb-ad)/(b-a)
    // If a==b, return a ValidatedTaylorModel with unbounded error

    ARIADNE_ASSERT_MSG(ocd.radius().lower().raw()>=0,"Illegal scaling from interval "<<ocd<<" with zero radius to interval "<<ncd);
    if(ocd.lower()==ocd.upper()) {
        x.clear();
        x.set_error(+inf);
    } else {
        ValidatedNumber tmp=1/(b-a);
        ValidatedNumber alpha=(d-c)*tmp;
        ValidatedNumber beta=(b*c-a*d)*tmp;
        x*=alpha;
        x+=beta;
    }

    return x;
}



// Replace the kth variable x[k] by a*x[k]+b.
ValidatedTaylorModel preaffine(const ValidatedTaylorModel& tm, Nat k, const ValidatedNumber& a, const ValidatedNumber& b) {
    Nat d=tm.degree();
    Nat as=tm.argument_size();
    Sweeper swp=tm.sweeper();

    ValidatedTaylorModel r(as,swp);
    r.set_error(tm.error().raw());

    // Create a temporary TaylorModels containing just terms x[k]^i
    Array<ValidatedTaylorModel> atm(d+1,ValidatedTaylorModel(as,swp));
    for(ValidatedTaylorModel::ConstIterator iter=tm.begin(); iter!=tm.end(); ++iter) {
        MultiIndex a=iter->key();
        const Float& c=iter->data().raw();
        Nat ak=a[k];
        a[k]=0;
        atm[ak].expansion().raw().append(a,c);
    }

    ValidatedTaylorModel xk=ValidatedTaylorModel::variable(as,k,swp);

    for(Nat i=0; i<=d; ++i) {
        for(Nat j=i; j<=d; ++j) {
            ValidatedNumber c=(Ariadne::bin(j,i)*Ariadne::pow(a,i)*Ariadne::pow(b,j-i));
            r+=c*atm[j];
            atm[j]*=xk;
        }
    }
    return r;
}


ValidatedTaylorModel discard(const ValidatedTaylorModel& tm, Array<Nat>& discarded_variables) {
    return recondition(tm,discarded_variables,0u,0u);
}

ValidatedTaylorModel recondition(const ValidatedTaylorModel& tm, Array<Nat>& discarded_variables, Nat number_of_error_variables) {
    return recondition(tm,discarded_variables,number_of_error_variables,number_of_error_variables);
}

ValidatedTaylorModel recondition(const ValidatedTaylorModel& tm, Array<Nat>& discarded_variables, Nat number_of_error_variables, Nat index_of_error)
{
    for(Nat i=0; i!=discarded_variables.size()-1; ++i) {
        ARIADNE_PRECONDITION(discarded_variables[i]<discarded_variables[i+1]);
    }
    ARIADNE_PRECONDITION(discarded_variables[discarded_variables.size()-1]<tm.argument_size());
    ARIADNE_PRECONDITION(index_of_error<=number_of_error_variables);

    const Nat number_of_discarded_variables = discarded_variables.size();
    const Nat number_of_kept_variables = tm.argument_size() - discarded_variables.size();

    // Make an Array of the variables to be kept
    Array<Nat> kept_variables(number_of_kept_variables);
    Nat kd=0; Nat kk=0;
    for(Nat j=0; j!=tm.argument_size(); ++j) {
        if(kd==number_of_discarded_variables || j!=discarded_variables[kd]) {
            kept_variables[kk]=j; ++kk;
        } else {
            ++kd;
        }
    }

    // Construct result and reserve memory
    ValidatedTaylorModel r(number_of_kept_variables+number_of_error_variables,tm.sweeper());
    r.expansion().reserve(tm.number_of_nonzeros()+1u);
    MultiIndex a(number_of_kept_variables+number_of_error_variables);

    // Set the uniform error of the original model
    // If index_of_error == number_of_error_variables, then the error is kept as a uniform error bound
    Float* error_ptr;
    if(number_of_error_variables==index_of_error) {
        error_ptr = &r.error().raw();
    } else {
        a[number_of_kept_variables+index_of_error]=1;
        r.expansion().raw().append(a,tm.error().raw());
        a[number_of_kept_variables+index_of_error]=0;
        error_ptr = &r.begin()->data().raw();
    }

    set_rounding_upward();
    for(ValidatedTaylorModel::ConstIterator iter=tm.begin(); iter!=tm.end(); ++iter) {
        Bool keep=true;
        for(Nat k=0; k!=number_of_discarded_variables; ++k) {
            if(iter->key()[discarded_variables[k]]!=0) {
                *error_ptr = add_rnd(*error_ptr,abs(iter->data().raw()));
                keep=false;
                break;
            }
        }
        if(keep) {
            for(Nat k=0; k!=number_of_kept_variables; ++k) {
                a[k]=iter->key()[kept_variables[k]];
            }
            r.expansion().raw().append(a,iter->data().raw());
        }
    }
    set_rounding_to_nearest();

    return r;
}


ValidatedTaylorModel restrict(const ValidatedTaylorModel& tm, Nat k, const ExactInterval& nd) {
    ARIADNE_ASSERT(k<tm.argument_size());
    ARIADNE_ASSERT(nd.lower().raw()>=-1 && nd.upper().raw()<=+1);
    if(nd.lower().raw()==-1 && nd.upper().raw()==1) {
        return tm;
    } else if(nd.lower().raw()==-1 && nd.upper().raw()==0) {
        return split(tm,k,false);
    } else if(nd.lower().raw()==0 && nd.upper().raw()==1) {
        return split(tm,k,true);
    } else if(nd.lower().raw()==-0.5 && nd.upper().raw()==0.5) {
        return split(tm,k,indeterminate);
    } else {
        ValidatedNumber a=(nd.upper()-nd.lower())/2;
        ValidatedNumber b=(nd.lower()+nd.upper())/2;
        return preaffine(tm,k,a,b);
    }
}


ValidatedTaylorModel& ValidatedTaylorModel::restrict(const Vector<ExactInterval>& nd)
{
    ValidatedTaylorModel& x=*this;
    ARIADNE_ASSERT(x.argument_size()==nd.size());
    const Nat as=x.argument_size();
    const Nat d=x.expansion().degree();
    if(d==0) { return x; }

    Array< Array<ValidatedNumber> > sf(as,Array<ValidatedNumber>(d+1u));
    for(Nat j=0; j!=as; ++j) {
        sf[j][0]=1.0;
        sf[j][1]=(nd[j].upper()-nd[j].lower())/2;
        for(Nat k=2; k<=d; ++k) {
            sf[j][k]=sf[j][k-1]*sf[j][1];
        }
    }

    // TODO: separate into roundoff computation and value computation
    for(Iterator iter=x.begin(); iter!=x.end(); ++iter) {
        for(Nat j=0; j!=as; ++j) {
            Float& c=iter->data().raw();
            UpperInterval ci=UpperInterval(c);
            for(Nat j=0; j!=as; ++j) {
                ci*=sf[j][iter->key()[j]];
            }
            iter->data().raw()=ci.centre().raw();
            x._error.raw()=add_up(x._error.raw(),mag(ci-iter->data()).raw());
        }
    }

    return x;
}

inline double div_rnd(double x, unsigned int m) { return (volatile double&)x/(volatile unsigned int&)m; }
inline double div_rnd(double x, int n) { return (volatile double&)x/(volatile int&)n; }

// Compute antiderivative by first computing the error and then the individual terms
// Seems to have problems with setting the rounding mode; with certain compiler flags
// sometimes the truncation error is not set
Void _antidifferentiate1(ValidatedTaylorModel& x, Nat k)
{
    ARIADNE_ASSERT(k<x.argument_size());

    //std::cerr<<"xe="<<xe<<"\n";
    set_rounding_mode(upward);
    VOLATILE Float tre=0; // Twice the maximum accumulated roundoff error
    for(ValidatedTaylorModel::ConstIterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
        const Nat c=xiter->key()[k]+1;
        Float xv=xiter->data().raw();
        Float mxv=-xv;
        assert(c>0);
        //VOLATILE Float ml,u;
        if(xv>=0) {
            //VOLATILE Float u=xv/c;
            //VOLATILE Float ml=mxv/c;
            tre+=add_rnd(div_rnd(xv,c),div_rnd(mxv,c));
        } else {
            tre+=add_rnd(div_rnd(mxv,c),div_rnd(xv,c));
            //VOLATILE Float u=mxv/c;
            //VOLATILE Float ml=xv/c;
            //te+=(u+ml);
        }
        //std::cerr<<"  te="<<te;;
    }
    //std::cerr<<"  te="<<tre;;
    x.error().raw()+=(tre/2);
    //xe+=te/2;
    //std::cerr<<"xe="<<xe<<"\n";

    set_rounding_to_nearest();
    for(ValidatedTaylorModel::Iterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
        MultiIndex& a=xiter->key();
        Float& v=xiter->data().raw();

        // Use the following code as opposed to ++a[k]; c=a[k];
        // as this seems to be safer against compiler optimisations through MultiIndexElementReference
        a.increment(k);
        const Int& c=a.at(k);
        assert(c>0);
        v/=c;
    }
}

// Compute antiderivative by computing term-by-term, switching the rounding mode
// May be slow if no assembly rounding mode switching is available.
Void _antidifferentiate2(ValidatedTaylorModel& x, Nat k)
{
    ARIADNE_ASSERT(k<x.argument_size());

    //std::cerr<<"xe="<<xe<<"\n";
    set_rounding_mode(upward);
    VOLATILE Float tre=0; // Twice the maximum accumulated roundoff error
    VOLATILE Float u,ml;
    Nat c;
    for(ValidatedTaylorModel::Iterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
        MultiIndex& xa=xiter->key();
        Float& xv=xiter->data().raw();
        Float mxv=-xv;
        c=xa.at(k)+1;

        set_rounding_upward();
        if(xv>=0) {
            u=xv/c;
            ml=mxv/c;
        } else {
            u=mxv/c;
            ml=xv/c;
        }
        tre+=(u+ml);

        set_rounding_to_nearest();
        xa.increment(k);
        xv/=c;
    }

    set_rounding_upward();
    x.error().raw()+=(tre/2);
    set_rounding_to_nearest();
}


ValidatedTaylorModel& ValidatedTaylorModel::antidifferentiate(Nat k)
{
    _antidifferentiate2(*this,k); return *this;
}

// Compute derivative by computing term-by-term, switching the rounding mode
// May be slow if no assembly rounding mode switching is available.


ValidatedTaylorModel derivative(const ValidatedTaylorModel& x, Nat k)
{
    ARIADNE_ASSERT(k<x.argument_size());
    ValidatedTaylorModel r(x.argument_size(),x.sweeper());

    //std::cerr<<"xe="<<xe<<"\n";
    VOLATILE Float tre=0; // Twice the maximum accumulated roundoff error
    Float xv;
    VOLATILE Float u,ml,mxv;
    MultiIndex a(x.argument_size());
    Nat c;
    for(ValidatedTaylorModel::Iterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
        c=xiter->key().at(k);
        if(c!=0) {
            a=xiter->key();
            xv=xiter->data().raw();
            mxv=-xv;

            if(c!=1 && c!=2 && c!=4) {
                set_rounding_upward();
                u=xv*c;
                ml=mxv*c;
                tre+=(u+ml);
                set_rounding_to_nearest();
            }
            --a[k];
            xv*=c;
            r.expansion().raw().append(a,xv);
        }
    }
    r.error().raw()=(tre/2);

    return r;
}



///////////////////////////////////////////////////////////////////////////////

// Scalar function operators (evaluate, split, unscale, embed, restrict)
// and predicates (refines)

ValidatedNumber
evaluate(const ValidatedTaylorModel& tm, const Vector<ValidatedNumber>& x)
{
    return horner_evaluate(tm.expansion(),x)+ValidatedNumber(-tm.error().raw(),+tm.error().raw());
}

Vector<ValidatedNumber>
evaluate(const Vector<ValidatedTaylorModel>& tv, const Vector<ValidatedNumber>& x)
{
    Vector<ValidatedNumber> r(tv.size());
    for(Nat i=0; i!=r.size(); ++i) {
        r[i]=evaluate(tv[i],x);
    }
    return r;
}


template<class T> class Powers {
  public:
    explicit Powers(const T& t) { _values.push_back(t*0+1); _values.push_back(t); }
    explicit Powers(const T& z, const T& t) { _values.push_back(z); _values.push_back(t); }
    const T& operator[](Nat i) const { while(_values.size()<=i) { _values.push_back(_values[1]*_values.back()); } return _values[i]; }
  private:
    mutable std::vector<T> _values;
};

template<> class Powers<ValidatedNumber> {
  public:
    explicit Powers(const ValidatedNumber& t) { _values.push_back(ValidatedNumber(1)); _values.push_back(t); }
    explicit Powers(const ValidatedNumber& z, const ValidatedNumber& t) { _values.push_back(z); _values.push_back(t); }
    const ValidatedNumber& operator[](Nat i) const {
        while(_values.size()<=i) {
            if(_values.size()%2==0) { _values.push_back(sqr(_values[_values.size()/2])); }
            else { _values.push_back(_values[1]*_values.back()); } }
        return _values[i]; }
  private:
    mutable std::vector<ValidatedNumber> _values;
};


ValidatedTaylorModel substitute(const ValidatedTaylorModel& x, Nat k, const ValidatedTaylorModel& s) {
    ARIADNE_ASSERT(x.argument_size()==s.argument_size()+1u);
    const Nat n=s.argument_size();
    Sweeper swp=s.sweeper();
    Vector<ValidatedTaylorModel> y(n+1,ValidatedTaylorModel::zero(n,swp));
    for(Nat i=0; i!=n; ++i) {
        y[i]=ValidatedTaylorModel::variable(n,i,swp);
    }
    y[n]=s;
    return compose(x,y);
}

Vector<ValidatedTaylorModel> substitute(const Vector<ValidatedTaylorModel>& x, Nat k, const ValidatedTaylorModel& s) {
    const Nat n=s.argument_size();
    Sweeper swp=s.sweeper();
    Vector<ValidatedTaylorModel> y(n+1,ValidatedTaylorModel::zero(n,swp));
    for(Nat i=0; i!=n; ++i) {
        y[i]=ValidatedTaylorModel::variable(n,i,swp);
    }
    y[n]=s;
    return compose(x,y);
}


ValidatedTaylorModel
partial_evaluate(const ValidatedTaylorModel& x, Nat k, ExactNumber c)
{
    return partial_evaluate(x,k,ValidatedNumber(c));
}


ValidatedTaylorModel
partial_evaluate(const ValidatedTaylorModel& x, Nat k, ValidatedNumber c)
{
    ValidatedTaylorModel r(x.argument_size()-1,x.sweeper());
    MultiIndex ra(r.argument_size());
    if(c==ValidatedNumber(0)) {
        for(ValidatedTaylorModel::ConstIterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
            const MultiIndex& xa=xiter->key();
            MultiIndex::IndexType xak=xa[k];
            if(xak==0) {
                const Float& xv=xiter->data().raw();
                for(Nat i=0; i!=k; ++i) { ra[i]=xa[i]; }
                for(Nat i=k; i!=ra.size(); ++i) { ra[i]=xa[i+1]; }
                r.expansion().raw().append(ra,xv);
            }
        }
        r.set_error(x.error().raw());
    } else if(c==ValidatedNumber(1)) {
        ValidatedTaylorModel s(x.argument_size()-1,x.sweeper());
        Array<ValidatedTaylorModel> p(x.degree()+1,ValidatedTaylorModel(x.argument_size()-1,x.sweeper()));

        for(ValidatedTaylorModel::ConstIterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
            const MultiIndex& xa=xiter->key();
            const Float& xv=xiter->data().raw();
            MultiIndex::IndexType xak=xa[k];
            for(Nat i=0; i!=k; ++i) { ra[i]=xa[i]; }
            for(Nat i=k; i!=ra.size(); ++i) { ra[i]=xa[i+1]; }
            assert(ra.degree()+xak==xa.degree());
            p[xak].expansion().raw().append(ra,xv);
        }

        r=p[0];
        r.set_error(x.error().raw());
        for(Nat i=1; i!=p.size(); ++i) {
            _add(s,r,p[i]);
            r.swap(s);
            s.clear();
        }
    } else {
        ValidatedTaylorModel s(x.argument_size()-1,x.sweeper());
        Array<ValidatedTaylorModel> p(x.degree()+1,ValidatedTaylorModel(x.argument_size()-1,x.sweeper()));

        Powers<ValidatedNumber> cpowers(c);

        for(ValidatedTaylorModel::ConstIterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
            const MultiIndex& xa=xiter->key();
            const Float& xv=xiter->data().raw();
            MultiIndex::IndexType xak=xa[k];
            for(Nat i=0; i!=k; ++i) { ra[i]=xa[i]; }
            for(Nat i=k; i!=ra.size(); ++i) { ra[i]=xa[i+1]; }
            assert(ra.degree()+xak==xa.degree());
            p[xak].expansion().raw().append(ra,xv);
        }
        for(Nat i=1; i!=p.size(); ++i) {
            p[i]*=cpowers[i];
        }

        r=p[0];
        r.set_error(x.error().raw());
        for(Nat i=1; i!=p.size(); ++i) {
            _add(s,r,p[i]);
            r.swap(s);
            s.clear();
        }
    }

    return r;
}

Vector<ValidatedTaylorModel>
partial_evaluate(const Vector<ValidatedTaylorModel>& tv, Nat k, ExactNumber c)
{
    // FIXME: Fails if tv.size()==0
    Vector<ValidatedTaylorModel> r(tv.size(),ValidatedTaylorModel::zero(tv.zero_element().argument_size(),tv.zero_element().sweeper()));
    for(Nat i=0; i!=r.size(); ++i) {
        r[i]=partial_evaluate(tv[i],k,c);
    }
    return r;
}

Vector<ValidatedTaylorModel>
partial_evaluate(const Vector<ValidatedTaylorModel>& tv, Nat k, ValidatedNumber c)
{
    Vector<ValidatedTaylorModel> r(tv.size(),ValidatedTaylorModel::zero(tv.zero_element().argument_size(),tv.zero_element().sweeper()));
    for(Nat i=0; i!=r.size(); ++i) {
        r[i]=partial_evaluate(tv[i],k,c);
    }
    return r;
}


Pair<ValidatedTaylorModel,ValidatedTaylorModel>
split(const ValidatedTaylorModel& tv, Nat j)
{
    // TODO: improve efficiency of implementation
    Nat as=tv.argument_size();
    Sweeper swp=tv.sweeper();
    Vector<ValidatedTaylorModel> s=ValidatedTaylorModel::variables(as,swp);
    s[j]=ValidatedTaylorModel::scaling(as,j,ExactInterval(-1,0),swp);
    ValidatedTaylorModel r1=compose(tv,s);
    r1.set_sweeper(tv.sweeper());
    s[j]=ValidatedTaylorModel::scaling(as,j,ExactInterval(0,+1),swp);
    ValidatedTaylorModel r2=compose(tv,s);
    r2.set_sweeper(tv.sweeper());
    return make_pair(r1,r2);
}




ValidatedTaylorModel
_split1(const ValidatedTaylorModel& tm, Nat k, Tribool b)
{
    const Nat deg=tm.degree();
    const Nat as=tm.argument_size();
    Sweeper swp=tm.sweeper();

    ValidatedTaylorModel r(tm);

    // Divide all coefficients by 2^a[k]
    // This can be done exactly
    for(ValidatedTaylorModel::Iterator iter=r.begin(); iter!=r.end(); ++iter) {
        const uchar ak=iter->key()[k];
        Float& c=iter->data().raw();
        c/=(1<<ak);
    }

    Int tr=( definitely(b) ? +1 : possibly(b) ? 0 : -1 );

    if(tr==0) { return r; }

    // Replace x[k] with x[k]+tr

    // Split variables by degree in x[k]
    Array<ValidatedTaylorModel> ary(deg+1,ValidatedTaylorModel(as,swp));
    for(ValidatedTaylorModel::ConstIterator iter=r.begin(); iter!=r.end(); ++iter) {
        MultiIndex a=iter->key();
        const Float& c=iter->data().raw();
        uchar ak=a[k];
        a[k]=0u;
        ary[ak].expansion().raw().append(a,c);
    }

    Float re=r.error().raw();
    r.clear();
    r.set_error(re);

    for(Nat i=0; i<=deg; ++i) {
        for(Nat j=i; j<=deg; ++j) {
            Int sf=bin(j,i);
            if(tr==-1 && (j-i)%2==1) { sf=-sf; }
            r+=ary[j]*sf;
            for(ValidatedTaylorModel::Iterator iter=ary[j].begin(); iter!=ary[j].end(); ++iter) {
                ++iter->key()[k];
            }
         }
    }

    return r;
}


ValidatedTaylorModel
split(const ValidatedTaylorModel& tm, Nat j, Tribool b)
{
    return _split1(tm,j,b);
}

ValidatedTaylorModel
unscale(const ValidatedTaylorModel& tv, const ExactInterval& ivl)
{
    // Scale tv so that the interval ivl maps into [-1,1]
    // The result is given by  (tv-c)*s where c is the centre
    // and s the reciprocal of the radius of ivl

    // If the radius of ivl is less than a constant, then
    // there are three possible policies. We can either map to zero or to
    // everything, or map a constant to zero and other models to everything.
    // The motivation for mapping to zero is that the domain is
    // restricted to the point and this is
    // The motivation for mapping to everything is that any function on the
    // resulting interval should be independent of the unneeded component

    auto& l=ivl.lower().raw();
    auto& u=ivl.upper().raw();
    ARIADNE_ASSERT_MSG(l<=u,"Cannot unscale ValidatedTaylorModel "<<tv<<" from empty interval "<<ivl);

    if(l==u) {
        ValidatedTaylorModel r=tv.create();
        r+=static_cast<ValidatedFloat>(ivl);
        // Uncomment out line below to make unscaling to a singleton interval undefined
        //r.set_error(+inf);
        return r;
    } else {
        ValidatedTaylorModel r=tv;
        ValidatedNumber c=ValidatedNumber(l/2)+ValidatedNumber(u/2);
        ValidatedNumber s=2/(ValidatedNumber(u)-ValidatedNumber(l));
        r-=c;
        r*=s;

        return r;
    }
}


ValidatedTaylorModel
rescale(const ValidatedTaylorModel& tv, const ExactInterval& ivl)
{
    // Scale tv so that the interval [-1,1] maps into ivl
    // The result is given by  (tv*s)+c where c is the centre
    // and r the radius of ivl
    const ExactFloat& l=ivl.lower();
    const ExactFloat& u=ivl.upper();

    ValidatedTaylorModel r=tv;
    ValidatedNumber c=(l+u)/2;
    ValidatedNumber s=(u+l)/2;
    r*=s;
    r+=c;

    return r;
}



ErrorFloat ValidatedTaylorModel::radius() const {
    set_rounding_mode(upward);
    Float r=this->error().raw();
    for(ValidatedTaylorModel::ConstIterator iter=this->begin(); iter!=this->end(); ++iter) {
        if(iter->key().degree()!=0) {
            r+=abs(iter->data().raw());
        }
    }
    set_rounding_mode(to_nearest);
    return ErrorFloat(r);
}

NormType ValidatedTaylorModel::norm() const {
    set_rounding_mode(upward);
    Float r=this->error().raw();
    for(ValidatedTaylorModel::ConstIterator iter=this->begin(); iter!=this->end(); ++iter) {
        r+=abs(iter->data().raw());
    }
    set_rounding_mode(to_nearest);
    return ErrorFloat(r);
}

NormType norm(const ValidatedTaylorModel& tv) {
    return tv.norm();
}

Bool
refines(const ValidatedTaylorModel& tv1, const ValidatedTaylorModel& tv2)
{
    ARIADNE_ASSERT(tv1.argument_size()==tv2.argument_size());
    ValidatedTaylorModel d=tv2;
    d.error().raw()=0.0;
    d-=tv1;
    return norm(d).raw() <= tv2.error().raw();
}


Bool
disjoint(const ValidatedTaylorModel& tv1, const ValidatedTaylorModel& tv2)
{
    ARIADNE_ASSERT(tv1.argument_size()==tv2.argument_size());
    Float tv1e=tv1.error().raw();
    Float tv2e=tv2.error().raw();
    const_cast<ValidatedTaylorModel&>(tv1).error().raw()=0.0;
    const_cast<ValidatedTaylorModel&>(tv2).error().raw()=0.0;
    ValidatedTaylorModel d=tv2-tv1;
    const_cast<ValidatedTaylorModel&>(tv1).error().raw()=tv1e;
    const_cast<ValidatedTaylorModel&>(tv2).error().raw()=tv2e;
    //std::cerr<<"\ntv1="<<tv1<<"\ntv2="<<tv2<<"\nd="<<d
    //         <<"\n|d|="<<norm(d)<<" |e|="<<add_up(tv1.error().raw(),tv2.error().raw())<<"\n\n";
    return norm(d).raw()>add_up(tv1.error().raw(),tv2.error().raw());
}


ValidatedTaylorModel embed(const ValidatedTaylorModel& x, Nat as)
{
    return ValidatedTaylorModel(embed(0u,x.expansion(),as),x.error(),x.sweeper());
}

ValidatedTaylorModel embed(Nat as, const ValidatedTaylorModel& x)
{
    return ValidatedTaylorModel(embed(as,x.expansion(),0u),x.error(),x.sweeper());
}


///////////////////////////////////////////////////////////////////////////////

// Input/output operators

OutputStream&
operator<<(OutputStream& os, const ValidatedTaylorModel& tm) {
    // Set the variable names to be 'parameter' s0,s1,..
    Array<StringType> variable_names(tm.argument_size());
    for(Nat j=0; j!=tm.argument_size(); ++j) {
        StringStream sstr;
        sstr << 's' << j;
        variable_names[j]=sstr.str();
    }

    //os << "ValidatedTaylorModel";
    os << "TM["<<tm.argument_size()<<"](";
    Expansion<RawFloat> e=tm.expansion().raw();
    e.graded_sort();
    e.write(os,variable_names);
    return os << "+/-" << tm.error().raw() << ")";
}


///////////////////////////////////////////////////////////////////////////////

// Vector-valued named constructors


Vector<ValidatedTaylorModel> ValidatedTaylorModel::zeros(Nat rs, Nat as, Sweeper swp)
{
    Vector<ValidatedTaylorModel> result(rs,ValidatedTaylorModel::zero(as,swp));
    return result;
}

Vector<ValidatedTaylorModel> ValidatedTaylorModel::constants(Nat as, const Vector<ExactNumber>& c, Sweeper swp)
{
    Vector<ValidatedTaylorModel> result(c.size(),ValidatedTaylorModel::zero(as,swp));
    for(Nat i=0; i!=c.size(); ++i) {
        result[i]=ValidatedTaylorModel::constant(as,c[i],swp);
    }
    return result;
}

Vector<ValidatedTaylorModel> ValidatedTaylorModel::constants(Nat as, const Vector<ValidatedNumber>& c, Sweeper swp)
{
    Vector<ValidatedTaylorModel> result(c.size(),ValidatedTaylorModel::zero(as,swp));
    for(Nat i=0; i!=c.size(); ++i) {
        result[i]=ValidatedTaylorModel::constant(as,c[i],swp);
    }
    return result;
}

Vector<ValidatedTaylorModel> ValidatedTaylorModel::variables(Nat as, Sweeper swp)
{
    Vector<ValidatedTaylorModel> result(as,ValidatedTaylorModel::zero(as,swp));
    for(Nat i=0; i!=as; ++i) { result[i]=ValidatedTaylorModel::variable(as,i,swp); }
    return result;
}

Vector<ValidatedTaylorModel> ValidatedTaylorModel::scalings(const Vector<ExactInterval>& d, Sweeper swp)
{
    Vector<ValidatedTaylorModel> result(d.size(),ValidatedTaylorModel::zero(d.size(),swp));
    for(Nat i=0; i!=d.size(); ++i) {
        result[i]=ValidatedTaylorModel::scaling(d.size(),i,d[i],swp);
    }
    return result;
}

Vector<ValidatedTaylorModel> ValidatedTaylorModel::unscalings(const Vector<ExactInterval>& cd, Sweeper swp)
{
    Vector<ValidatedTaylorModel> result(cd.size(),ValidatedTaylorModel::zero(cd.size(),swp));
    for(Nat i=0; i!=cd.size(); ++i) {
        result[i]=ValidatedTaylorModel::unscaling(cd.size(),i,cd[i],swp);
    }
    return result;
}



///////////////////////////////////////////////////////////////////////////////

// Vector-valued versions of scalar operators




Vector<ValidatedTaylorModel> embed(const Vector<ValidatedTaylorModel>& x, Nat as)
{
    Vector<ValidatedTaylorModel> r(x.size());
    for(Nat i=0; i!=x.size(); ++i) {
        r[i]=embed(x[i],as);
    }
    return r;
}

Vector<ValidatedTaylorModel> embed(Nat as, const Vector<ValidatedTaylorModel>& x)
{
    Vector<ValidatedTaylorModel> r(x.size());
    for(Nat i=0; i!=x.size(); ++i) {
        r[i]=embed(as,x[i]);
    }
    return r;
}

Vector<ValidatedTaylorModel> combine(const Vector<ValidatedTaylorModel>& x1, const Vector<ValidatedTaylorModel>& x2) {
    return join(embed(x1,x2.zero_element().argument_size()),embed(x1.zero_element().argument_size(),x2));
}

Vector<ValidatedTaylorModel> combine(const Vector<ValidatedTaylorModel>& x1, const ValidatedTaylorModel& x2) {
    return join(embed(x1,x2.argument_size()),embed(x1.zero_element().argument_size(),x2));
}

Vector<ValidatedTaylorModel> combine(const ValidatedTaylorModel& x1, const Vector<ValidatedTaylorModel>& x2) {
    return join(embed(x1,x2.zero_element().argument_size()),embed(x1.argument_size(),x2));
}

Vector<ValidatedTaylorModel> combine(const ValidatedTaylorModel& x1, const ValidatedTaylorModel& x2) {
    return join(embed(x1,x2.argument_size()),embed(x1.argument_size(),x2));
}

Bool
refines(const Vector<ValidatedTaylorModel>& tv1, const Vector<ValidatedTaylorModel>& tv2)
{
    ARIADNE_ASSERT(tv1.size()==tv2.size());
    for(Nat i=0; i!=tv1.size(); ++i) {
        if(!refines(tv1[i],tv2[i])) { return false; }
    }
    return true;
}

Bool
disjoint(const Vector<ValidatedTaylorModel>& tv1, const Vector<ValidatedTaylorModel>& tv2)
{
    ARIADNE_ASSERT(tv1.size()==tv2.size());
    for(Nat i=0; i!=tv1.size(); ++i) {
        if(disjoint(tv1[i],tv2[i])) { return true; }
    }
    return false;
}

Pair< Vector<ValidatedTaylorModel>, Vector<ValidatedTaylorModel> >
split(const Vector<ValidatedTaylorModel>& tv, Nat j)
{
    Vector<ValidatedTaylorModel> r1(tv.size());
    Vector<ValidatedTaylorModel> r2(tv.size());
    for(Nat i=0; i!=tv.size(); ++i) {
        make_lpair(r1[i],r2[i])=split(tv[i],j);
    }
    return make_pair(r1,r2);
}

Vector<ValidatedTaylorModel>
split(const Vector<ValidatedTaylorModel>& tv, Nat j, Bool h)
{
    Vector<ValidatedTaylorModel> r(tv.size());
    for(Nat i=0; i!=tv.size(); ++i) {
        r[i]=split(tv[i],j,h);
    }
    return r;
}

Vector<ValidatedTaylorModel>
unscale(const Vector<ValidatedTaylorModel>& tvs, const Vector<ExactInterval>& ivls)
{
    Vector<ValidatedTaylorModel> r(tvs.size());
    for(Nat i=0; i!=r.size(); ++i) {
        r[i]=unscale(tvs[i],ivls[i]);
    }
    return r;
}

ValidatedTaylorModel
rescale(const ValidatedTaylorModel& tv, const ExactInterval& old_codom, const ExactInterval& new_codom)
{
    ValidatedTaylorModel res(tv); res.rescale(old_codom,new_codom); return res;
}

ValidatedTaylorModel&
rescale(ValidatedTaylorModel& tv, const ExactInterval& old_codom, const ExactInterval& new_codom)
{
    tv.rescale(old_codom,new_codom); return tv;
}

Vector<ValidatedTaylorModel>
rescale(const Vector<ValidatedTaylorModel>& tvs, const Vector<ExactInterval>& old_codom, const Vector<ExactInterval>& new_codom)
{
    Vector<ValidatedTaylorModel> r(tvs.size());
    for(Nat i=0; i!=r.size(); ++i) {
        r[i]=rescale(tvs[i],old_codom[i],new_codom[i]);
    }
    return r;
}

Vector<ValidatedTaylorModel>&
rescale(Vector<ValidatedTaylorModel>& tvs, const Vector<ExactInterval>& old_codom, const Vector<ExactInterval>& new_codom)
{
    for(Nat i=0; i!=tvs.size(); ++i) {
        rescale(tvs[i],old_codom[i],new_codom[i]);
    }
    return tvs;
}

Vector<ValidatedTaylorModel>
scale(const Vector<ValidatedTaylorModel>& tvs, const Vector<ExactInterval>& new_codom)
{
    Vector<ValidatedTaylorModel> r(tvs);
    for(Nat i=0; i!=tvs.size(); ++i) {
        r[i].rescale(ExactInterval(-1,+1),new_codom[i]);
    }
    return r;
}

ValidatedTaylorModel antiderivative(const ValidatedTaylorModel& x, Nat k) {
    ValidatedTaylorModel r(x);
    r.antidifferentiate(k);
    return r;
}

Vector<ValidatedTaylorModel> antiderivative(const Vector<ValidatedTaylorModel>& x, Nat k) {
    Vector<ValidatedTaylorModel> r(x);
    for(Nat i=0; i!=x.size(); ++i) {
        r[i].antidifferentiate(k);
    }
    return r;
}

ValidatedTaylorModel antiderivative(const ValidatedTaylorModel& x, const ExactInterval& dk, Nat k) {
    ValidatedNumber dkr=static_cast<ValidatedNumber>(rad_ivl(dk.lower().raw(),dk.upper().raw()));
    ValidatedTaylorModel r(x);
    r.antidifferentiate(k);
    r*=dkr;
    return r;
}

Vector<ValidatedTaylorModel> antiderivative(const Vector<ValidatedTaylorModel>& x, const ExactInterval& dk, Nat k) {
    Vector<ValidatedTaylorModel> r(x.size());
    for(Nat i=0; i!=x.size(); ++i) {
        r[i]=antiderivative(x[i],dk,k);
    }
    return r;
}

ValidatedTaylorModel intersection(const ValidatedTaylorModel& x, const ValidatedTaylorModel& y) {
    ValidatedTaylorModel r(x.argument_size(),x.sweeper());
    Float twice_max_error=0.0;

    const Float& xe=x.error().raw();
    const Float& ye=y.error().raw();
    VOLATILE Float rv,xv,yv,xu,yu,mxl,myl,u,ml;
    //const MultiIndex* aptr;
    MultiIndex a;

    ValidatedTaylorModel::ConstIterator xiter=x.begin();
    ValidatedTaylorModel::ConstIterator yiter=y.begin();
    while(xiter!=x.end() || yiter!=y.end()) {
        // Can't use const MultiIndex& here since references change as the iterators change
        // We would need to use a smart reference
        //const MultiIndex& xa=xiter->key();
        //const MultiIndex& ya=yiter->key();
        if(xiter==x.end()) {
            a=yiter->key();
            yv=yiter->data().raw();
            xv=0.0;
            ++yiter;
        } else if(yiter==y.end()) {
            a=xiter->key();
            xv=xiter->data().raw();
            yv=0.0;
            ++xiter;
        } else if(xiter->key()==yiter->key()) {
            a=xiter->key();
            xv=xiter->data().raw();
            yv=yiter->data().raw();
            ++xiter;
            ++yiter;
        } else if(xiter->key()<yiter->key()) {
            a=xiter->key();
            xv=xiter->data().raw();
            yv=0.0;
            ++xiter;
        } else { // xa>ya
            a=yiter->key();
            yv=yiter->data().raw();
            xv=0.0;
            ++yiter;
        }

        set_rounding_upward();
        xu=xv+xe;
        yu=yv+ye;
        mxl=xe-xv;
        myl=ye-yv;
        u=min(xu,yu);
        ml=min(mxl,myl);
        if(u+ml<0) {
            ARIADNE_THROW(IntersectionException,"intersection(ValidatedTaylorModel,ValidatedTaylorModel)",x<<" and "<<y<<" are disjoint.");
        }
        twice_max_error=max(twice_max_error,(u+ml));

        set_rounding_to_nearest();
        xu=xv+xe;
        yu=yv+ye;
        mxl=xe-xv;
        myl=ye-yv;
        u=min(xu,yu);
        ml=min(mxl,myl);
        rv=(u-ml)/2;
        if(rv!=0.0) { r.expansion().raw().append(a,(u-ml)/2); }
    }

    r.error().raw()=twice_max_error/2;

    return r;
}

Vector<UpperInterval> ranges(const Vector<ValidatedTaylorModel>& f) {
    Vector<UpperInterval> r(f.size()); for(Nat i=0; i!=f.size(); ++i) { r[i]=f[i].range(); } return r;
}

inline Vector<ValidatedTaylorModel>& clobber(Vector<ValidatedTaylorModel>& h) {
    for(Nat i=0; i!=h.size(); ++i) { h[i].set_error(0u); } return h; }

inline Vector<ErrorFloat> errors(const Vector<ValidatedTaylorModel>& h) {
    Vector<ErrorFloat> e(h.size()); for(Nat i=0; i!=h.size(); ++i) { e[i]=h[i].error(); } return e; }

inline Vector<ErrorFloat> norms(const Vector<ValidatedTaylorModel>& h) {
    Vector<ErrorFloat> r(h.size()); for(Nat i=0; i!=h.size(); ++i) { r[i]=norm(h[i]); } return r; }

ErrorFloat norm(const Vector<ValidatedTaylorModel>& h) {
    return norm(norms(h));
}

///////////////////////////////////////////////////////////////////////////////

// Jacobian matrices

// Compute the Jacobian over an arbitrary domain
Matrix<ValidatedNumber>
jacobian(const Vector<ValidatedTaylorModel>& f, const Vector<ValidatedNumber>& x)
{
    Vector< Differential<ValidatedNumber> > dx=Differential<ValidatedNumber>::variables(1u,x);
    Vector< Differential<ValidatedNumber> > df(f.size(),x.size(),1u);
    for(Nat i=0; i!=f.size(); ++i) {
        df[i]=evaluate(f[i].expansion(),dx);
    }
    return jacobian(df);
}

// Compute the Jacobian over an arbitrary domain
Matrix<ValidatedNumber>
jacobian2(const Vector<ValidatedTaylorModel>& f, const Vector<ValidatedNumber>& x)
{
    Vector< Differential<ValidatedNumber> > dx(x.size());
    for(Nat i=0; i!=x.size()-f.size(); ++i) {
        dx[i]=Differential<ValidatedNumber>::constant(f.size(),1u,x[i]); }
    for(Nat i=0; i!=f.size(); ++i) {
        Nat j=i+(x.size()-f.size());
        dx[j]=Differential<ValidatedNumber>::variable(f.size(),1u,x[j],i); }
    Vector< Differential<ValidatedNumber> > df(f.size());
    for(Nat i=0; i!=f.size(); ++i) {
        df[i]=evaluate(f[i].expansion(),dx);
    }
    Matrix<ValidatedNumber> J=jacobian(df);
    return J;
}

/*
// Compute the Jacobian over an arbitrary domain
Matrix<ValidatedNumber>
jacobian(const Vector<ValidatedTaylorModel>& f, const Vector<ValidatedNumber>& d)
{
    Nat rs=f.size();
    Nat as=f.zero_element().argument_size();
    Matrix<ValidatedNumber> J(rs,as);
    for(Nat i=0; i!=rs; ++i) {
        for(ValidatedTaylorModel::ConstIterator iter=f[i].begin(); iter!=f[i].end(); ++iter) {
            const MultiIndex& a=iter->key();
            const Float& x=iter->data().raw();
            for(Nat k=0; k!=as; ++k) {
                const Nat c=a[k];
                if(c>0) {
                    if(iter->key().degree()==1) {
                        J[i][k]+=x;
                    }
                    else {
                        ValidatedNumber p(-x,+x);
                        p*=Float(c);
                        for(Nat l=0; l!=as; ++l) {
                            if(l==k) { if(a[l]>1) { p*=pow(d[l],a[l]-1); } }
                            else { if(a[k]>0) { p*=pow(d[k],a[k]); } }
                        }
                        J[i][k]+=p;
                    }
                }
            }
        }
    }
    return J;
}
*/

// Compute the Jacobian over the unit domain
Matrix<Float>
jacobian_value(const Vector<ValidatedTaylorModel>& f)
{
    Nat rs=f.size();
    Nat as=f.zero_element().argument_size();
    Matrix<Float> J(rs,as);
    MultiIndex a(as);
    for(Nat i=0; i!=rs; ++i) {
        for(Nat j=0; j!=as; ++j) {
            a[j]=1; const Float x=f[i][a].raw(); J[i][j]=x; a[j]=0;
        }
    }
    return J;
}

// Compute the Jacobian over the unit domain
Matrix<Float>
jacobian2_value(const Vector<ValidatedTaylorModel>& f)
{
    const Nat rs=f.size();
    const Nat fas=f.zero_element().argument_size();
    const Nat has=fas-rs;
    Matrix<Float> J(rs,rs);
    MultiIndex a(fas);
    for(Nat i=0; i!=rs; ++i) {
        for(Nat j=0; j!=rs; ++j) {
            a[has+j]=1; const Float x=f[i][a].raw(); J[i][j]=x; a[has+j]=0;
        }
    }
    return J;
}



// Compute the Jacobian over the unit domain
Matrix<UpperInterval>
jacobian_range(const Vector<ValidatedTaylorModel>& f)
{
    Nat rs=f.size();
    Nat as=f.zero_element().argument_size();
    Matrix<UpperInterval> J(rs,as);
    for(Nat i=0; i!=rs; ++i) {
        for(ValidatedTaylorModel::ConstIterator iter=f[i].begin(); iter!=f[i].end(); ++iter) {
            for(Nat k=0; k!=as; ++k) {
                const Nat c=iter->key()[k];
                if(c>0) {
                    const Float& x=iter->data().raw();
                    if(iter->key().degree()==1) { J[i][k]+=x; }
                    else { J[i][k]+=UpperInterval(-1,1)*x*c; }
                    //std::cerr<<"  J="<<J<<" i="<<i<<" a="<<iter->key()<<" k="<<k<<" c="<<c<<" x="<<x<<std::endl;
                }
            }
        }
    }
    return J;
}

// Compute the Jacobian over the unit domain
Matrix<UpperInterval>
jacobian2_range(const Vector<ValidatedTaylorModel>& f)
{
    Nat rs=f.size();
    Nat fas=f.zero_element().argument_size();
    Nat has=fas-rs;
    Matrix<UpperInterval> J(rs,rs);
    for(Nat i=0; i!=rs; ++i) {
        for(ValidatedTaylorModel::ConstIterator iter=f[i].begin(); iter!=f[i].end(); ++iter) {
            for(Nat k=0; k!=rs; ++k) {
                const Nat c=iter->key()[has+k];
                if(c>0) {
                    const Float& x=iter->data().raw();
                    if(iter->key().degree()==1) { J[i][k]+=x; }
                    else { J[i][k]+=UpperInterval(-1,1)*x*c; }
                    //std::cerr<<"  J="<<J<<" i="<<i<<" a="<<iter->key()<<" k="<<k<<" c="<<c<<" x="<<x<<std::endl;
                }
            }
        }
    }
    return J;
}


ValidatedNumber
jacobian2(const ValidatedTaylorModel& f, const Vector<ValidatedNumber>& x)
{
    return jacobian2(Vector<ValidatedTaylorModel>(1u,f),x)[0][0];
}

Float
jacobian2_value(const ValidatedTaylorModel& f)
{
    return jacobian2_value(Vector<ValidatedTaylorModel>(1u,f))[0][0];
}

UpperInterval jacobian2_range(const ValidatedTaylorModel& f) {
    return jacobian2_range(Vector<ValidatedTaylorModel>(1u,f))[0][0];
}

///////////////////////////////////////////////////////////////////////////////

// Vector operators (evaluate, compose




// Compose by computing each and every term individually without caching
// Easy to implement, but far too slow
inline
Vector<ValidatedTaylorModel>
_compose1(const Vector<ValidatedTaylorModel>& x,
          const Vector<ValidatedTaylorModel>& ys)
{
    //std::cerr<<"compose1"<<std::endl;
    ARIADNE_ASSERT(x.size()>0);
    ARIADNE_ASSERT(ys.size()==x.zero_element().argument_size());
    for(Nat i=0; i!=x.size(); ++i) { ARIADNE_ASSERT(x[i].argument_size()==x.zero_element().argument_size()); }
    for(Nat i=0; i!=ys.size(); ++i) { ARIADNE_ASSERT_MSG(ys[i].argument_size()==ys.zero_element().argument_size(),"ys="<<ys); }

    Nat as=ys.zero_element().argument_size();
    Sweeper swp=ys.zero_element().sweeper();

    Vector<ValidatedTaylorModel> r(x.size(),ValidatedTaylorModel(as,swp));
    ValidatedTaylorModel t(as,swp);
    for(Nat i=0; i!=x.size(); ++i) {
        r[i].set_error(x[i].error().raw());
        for(ValidatedTaylorModel::ConstIterator iter=x[i].begin(); iter!=x[i].end(); ++iter) {
            t=static_cast<ExactNumber>(iter->data().raw());
            for(Nat j=0; j!=iter->key().size(); ++j) {
                ValidatedTaylorModel p=pow(ys[j],iter->key()[j]);
                t=t*p;
            }
            r[i]+=t;
        }
    }

    return r;
}


inline
Vector<ValidatedTaylorModel>
_compose2(const Vector<ValidatedTaylorModel>& x,
          const Vector<ValidatedTaylorModel>& ys)
{
    //std::cerr<<"compose2"<<std::endl;
    Nat yrs=ys.size();
    Nat xas=ys.size();
    Nat as=ys.zero_element().argument_size();
    Sweeper sweeper=ys.zero_element().sweeper();

    Array<uchar> max_power(ys.size());
    for(Nat j=0; j!=ys.size(); ++j) { max_power[j]=1; }

    for(Nat i=0; i!=x.size(); ++i) {
        for(ValidatedTaylorModel::ConstIterator iter=x[i].begin(); iter!=x[i].end(); ++iter) {
            assert(xas==iter->key().size());
            for(Nat j=0; j!=iter->key().size(); ++j) {
                max_power[j]=max(max_power[j],iter->key()[j]);
            }
        }
    }

    uchar max_max_power = 1;
    for(Nat j=0; j!=ys.size(); ++j) {
        max_max_power = max(max_max_power,max_power[j]);
    }

    Array< Array< ValidatedTaylorModel > > powers(yrs, Array<ValidatedTaylorModel>(max_max_power+1,ValidatedTaylorModel::zero(as,sweeper)));
    for(Nat j=0; j!=yrs; ++j) {
        powers[j][0]=ys[j]*0;
        powers[j][1]=ys[j];
        for(Nat k=2; k!=powers[j].size(); ++k) {
            powers[j][k]=powers[j][k/2]*powers[j][(k+1)/2];
        }
    }

    Vector<ValidatedTaylorModel> r(x.size(),ValidatedTaylorModel(as,sweeper));
    ValidatedTaylorModel t(as,sweeper);
    MultiIndex a;
    Float c;
    for(Nat i=0; i!=x.size(); ++i) {
        r[i].set_error(x[i].error().raw());
        for(ValidatedTaylorModel::ConstIterator iter=x[i].begin(); iter!=x[i].end(); ++iter) {
            a=iter->key();
            c=iter->data().raw();
            t=static_cast<ExactNumber>(c);
            for(Nat j=0; j!=a.size(); ++j) {
                if(a[j]>0) {
                    t=t*powers[j][a[j]];
                }
            }
            r[i]+=t;
        }
    }

    return r;
}

Vector<ValidatedTaylorModel>
_compose(const Vector<ValidatedTaylorModel>& x,
         const Vector<ValidatedTaylorModel>& ys)
{
    ARIADNE_ASSERT_MSG(x.size()>0,"x="<<x<<", ys="<<ys);
    ARIADNE_ASSERT_MSG(ys.size()==x[0].argument_size(),"x="<<x<<", ys="<<ys);
    for(Nat i=1; i!=x.size(); ++i) { ARIADNE_ASSERT(x[i].argument_size()==x[0].argument_size()); }
    for(Nat i=1; i!=ys.size(); ++i) { ARIADNE_ASSERT_MSG(ys[i].argument_size()==ys[0].argument_size(),"ys="<<ys); }

    return _compose2(x,ys);
}


Vector<ValidatedTaylorModel>
unchecked_compose(const Vector<ValidatedTaylorModel>& x,
                  const Vector<ValidatedTaylorModel>& ys)
{
    return _compose(x,ys);
}

ValidatedTaylorModel
unchecked_compose(const ValidatedTaylorModel& x,
                  const Vector<ValidatedTaylorModel>& ys)
{
    return _compose(Vector<ValidatedTaylorModel>(1u,x),ys)[0];
}

Vector<ValidatedTaylorModel>
compose(const Vector<ValidatedTaylorModel>& x,
        const Vector<ValidatedTaylorModel>& ys)
{
    return _compose(x,ys);
}

ValidatedTaylorModel
compose(const ValidatedTaylorModel& x,
        const Vector<ValidatedTaylorModel>& ys)
{
    return Ariadne::power_evaluate(x.expansion(),ys)+ValidatedNumber(-x.error().raw(),+x.error().raw());
    return _compose(Vector<ValidatedTaylorModel>(1u,x),ys)[0];
}

ValidatedTaylorModel
compose(const ValidatedTaylorModel& x,
        const Vector<ExactInterval>& d,
        const Vector<ValidatedTaylorModel>& y)
{
    return _compose(Vector<ValidatedTaylorModel>(1u,x),unscale(y,d))[0];
}

Vector<ValidatedTaylorModel>
compose(const Vector<ValidatedTaylorModel>& x,
        const Vector<ExactInterval>& d,
        const Vector<ValidatedTaylorModel>& y)
{
    return _compose(x,unscale(y,d));
}

Vector<ValidatedTaylorModel>
unchecked_compose(const Vector<ValidatedTaylorModel>& x,
        const Vector<ExactInterval>& d,
        const Vector<ValidatedTaylorModel>& y)
{
    return _compose(x,unscale(y,d));
}



template TaylorModel<ValidatedNumber> rec(const TaylorModel<ValidatedNumber>&);
template TaylorModel<ValidatedNumber> sqrt(const TaylorModel<ValidatedNumber>&);
template TaylorModel<ValidatedNumber> exp(const TaylorModel<ValidatedNumber>&);
template TaylorModel<ValidatedNumber> log(const TaylorModel<ValidatedNumber>&);
template TaylorModel<ValidatedNumber> sin(const TaylorModel<ValidatedNumber>&);
template TaylorModel<ValidatedNumber> cos(const TaylorModel<ValidatedNumber>&);
template TaylorModel<ValidatedNumber> tan(const TaylorModel<ValidatedNumber>&);







ApproximateTaylorModel::TaylorModel(Nat as)
    : _expansion(as), _sweeper()
{
}

ApproximateTaylorModel::TaylorModel(Nat as, Sweeper swp)
    : _expansion(as), _sweeper(swp)
{
}

Void ApproximateTaylorModel::iadd(const ApproximateFloat& c)
{
    // Compute self+=c
    ApproximateTaylorModel& x=*this;
    if(c==0) { return; }
    if(x._expansion.empty()) {
        x._expansion.append(MultiIndex(x.argument_size()),c);
    } else if((x._expansion.end()-1)->key().degree()>0) {
        x._expansion.append(MultiIndex(x.argument_size()),c);
    } else {
        ApproximateFloat& rv=(x._expansion.end()-1)->data();
        rv+=c;
    }
}

Void ApproximateTaylorModel::imul(const ApproximateFloat& c)
{
    // Compute self*=c
    if(c==0) { this->clear(); return; }
    if(c==1) { return; }
    for(ExpansionType::Iterator iter=this->_expansion.begin(); iter!=this->_expansion.end(); ++iter) {
        iter->data() *= c;
    }
}



Void ApproximateTaylorModel::isma(const ApproximateFloat& c, const ApproximateTaylorModel& y)
{
    const ApproximateTaylorModel& x=*this;
    ARIADNE_ASSERT_MSG(x.argument_size()==y.argument_size(),"x="<<x<<", y="<<y);
    ApproximateTaylorModel r(x.argument_size());

    ConstIterator xiter=x._expansion.begin();
    ConstIterator yiter=y._expansion.begin();
    while(xiter!=x._expansion.end() && yiter!=y._expansion.end()) {
        if(xiter->key()<yiter->key()) {
            r._expansion.append(xiter->key(),xiter->data());
            ++xiter;
        } else if(yiter->key()<xiter->key()) {
            r._expansion.append(yiter->key(),c*yiter->data());
            ++yiter;
        } else {
            r._expansion.append(xiter->key(),xiter->data()+c*yiter->data());
            ++xiter; ++yiter;
        }
    }
    while(xiter!=x._expansion.end()) {
        r._expansion.append(xiter->key(),xiter->data());
        ++xiter;
    }
    while(yiter!=y._expansion.end()) {
        r._expansion.append(yiter->key(),c*yiter->data());
        ++yiter;
    }

    this->swap(r);
}


Void ApproximateTaylorModel::ifma(const ApproximateTaylorModel& x, const ApproximateTaylorModel& y)
{
    ApproximateTaylorModel& r(*this);
    ARIADNE_ASSERT_MSG(x.argument_size()==y.argument_size(),"x="<<x<<",y="<<y);

    ApproximateTaylorModel t(x.argument_size());
    MultiIndex sa;
    for(ExpansionType::ConstIterator xiter=x._expansion.begin(); xiter!=x._expansion.end(); ++xiter) {
        ExpansionType::ConstIterator yiter=y._expansion.begin();
        ExpansionType::ConstIterator riter=r._expansion.begin();
        while(riter!=r._expansion.end() && yiter!=y._expansion.end()) {
            const MultiIndex& ra=riter->key();
            sa=xiter->key()+yiter->key();
            if(sa==ra) {
                t._expansion.append(sa,riter->data()+xiter->data()*yiter->data());
                ++riter; ++yiter;
            } else if(sa<ra) {
                t._expansion.append(sa,xiter->data()*yiter->data());
                ++yiter;
            } else {
                t._expansion.append(ra,riter->data());
                ++riter;
            }
        }
        while(riter!=r._expansion.end()) {
            t._expansion.append(riter->key(),riter->data());
            ++riter;
        }
        while(yiter!=y._expansion.end()) {
            t._expansion.append(xiter->key(),yiter->key(),xiter->data()*yiter->data());
            ++yiter;

        }

        r._expansion.swap(t._expansion);
        t._expansion.clear();
    }
}




ErrorType ApproximateTaylorModel::norm() const {
    return ExactFloat(((*this)[MultiIndex(this->argument_size())]).raw());
}

CoefficientType ApproximateTaylorModel::average() const {
    return ExactFloat(((*this)[MultiIndex(this->argument_size())]).raw());
}

ErrorType ApproximateTaylorModel::radius() const {
    ApproximateFloat r=0.0;
    for(Expansion<ApproximateFloat>::ConstIterator iter=this->_expansion.begin(); iter!=this->_expansion.end(); ++iter) {
        if(iter->key().degree()!=0) {
            r+=abs(iter->data());
        }
    }
    return ErrorType(r.raw());
}

UpperInterval ApproximateTaylorModel::range() const {
    ErrorType rad=this->radius(); return this->average()+UpperInterval(-rad.raw(),+rad.raw());
}

Float ApproximateTaylorModel::tolerance() const {
    return dynamic_cast<const ThresholdSweeper&>(static_cast<const SweeperInterface&>(this->_sweeper)).sweep_threshold();
}


OutputStream& ApproximateTaylorModel::write(OutputStream& os) const {
    return os << "TM["<<this->argument_size()<<"](" << this->_expansion << ")";
}


template TaylorModel<ApproximateNumber> neg(const TaylorModel<ApproximateNumber>&);
template TaylorModel<ApproximateNumber> rec(const TaylorModel<ApproximateNumber>&);
template TaylorModel<ApproximateNumber> sqr(const TaylorModel<ApproximateNumber>&);
template TaylorModel<ApproximateNumber> pow(const TaylorModel<ApproximateNumber>&, Int);
template TaylorModel<ApproximateNumber> sqrt(const TaylorModel<ApproximateNumber>&);
template TaylorModel<ApproximateNumber> exp(const TaylorModel<ApproximateNumber>&);
template TaylorModel<ApproximateNumber> log(const TaylorModel<ApproximateNumber>&);
template TaylorModel<ApproximateNumber> sin(const TaylorModel<ApproximateNumber>&);
template TaylorModel<ApproximateNumber> cos(const TaylorModel<ApproximateNumber>&);
template TaylorModel<ApproximateNumber> tan(const TaylorModel<ApproximateNumber>&);


OutputStream& ApproximateTaylorModel::str(OutputStream& os) const {
    return os << this->_expansion;
}
} //namespace Ariadne

