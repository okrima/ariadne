/***************************************************************************
 *            interval.cc
 *
 *  Copyright 2008-10  Pieter Collins
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

#include <iostream>
#include <iomanip>
#include <cassert>

#include "config.h"
#include "macros.h"
#include "integer.h"
#include "float.h"
#include "rational.h"
#include "interval.h"


namespace Ariadne {

namespace {
static const double pi_up=3.1415926535897936;
static const double pi_approx=3.1415926535897931;
static const double pi_down=3.1415926535897931;
inline double _add_down(volatile double x, volatile double y) { set_rounding_downward(); return x+y; }
inline double _add_up(volatile double x, volatile double y) { set_rounding_upward(); return x+y; }
inline double _sub_down(volatile double x, volatile double y) { set_rounding_downward(); return x-y; }
inline double _sub_up(volatile double x, volatile double y) { set_rounding_upward(); return x-y; }
inline double _mul_down(volatile double x, volatile double y) { set_rounding_downward(); return x*y; }
inline double _mul_up(volatile double x, volatile double y) { set_rounding_upward(); return x*y; }
inline double _div_down(volatile double x, volatile double y) { set_rounding_downward(); return x/y; }
inline double _div_up(volatile double x, volatile double y) { set_rounding_upward(); return x/y; }
inline Float cos_down(Float x) { set_rounding_downward(); Float y=cos_rnd(x); return y; }
inline Float cos_up(Float x) { set_rounding_upward(); Float y=cos_rnd(x); return y; }
} // namespace




Interval widen(Interval x)
{
    rounding_mode_t rm=get_rounding_mode();
    const double& xl=internal_cast<const double&>(x.lower());
    const double& xu=internal_cast<const double&>(x.upper());
    const double m=std::numeric_limits<float>::min();
    set_rounding_upward();
    volatile double wu=xu+m;
    volatile double mwl=-xl+m;
    volatile double wl=-mwl;
    set_rounding_mode(rm);
    assert(wl<xl); assert(wu>xu);
    return Interval(wl,wu);
}

Interval trunc(Interval x)
{

    rounding_mode_t rm=get_rounding_mode();
    const double& xl=internal_cast<const double&>(x.lower());
    const double& xu=internal_cast<const double&>(x.upper());
    // Use machine epsilon instead of minimum to move away from zero
    const float fm=std::numeric_limits<float>::epsilon();
    volatile float tu=xu;
    if(tu<xu) { set_rounding_upward(); tu+=fm; }
    volatile float tl=xl;
    if(tl>xl) { set_rounding_downward(); tl-=fm; }
    set_rounding_mode(rm);
    assert(tl<=xl); assert(tu>=xu);
    return Interval(double(tl),double(tu));
}

Interval trunc(Interval x, uint n)
{
    Interval e=Interval(std::pow(2.0,52-n));
    Interval y=x+e;
    return y-e;
}

Interval rec(Interval i)
{
    volatile double& il=internal_cast<volatile double&>(i.lower());
    volatile double& iu=internal_cast<volatile double&>(i.upper());
    volatile double rl,ru;
    if(il>0 || iu<0) {
        rounding_mode_t rnd=get_rounding_mode();
        rl=_div_down(1.0,iu);
        ru=_div_up(1.0,il);
        set_rounding_mode(rnd);
    } else {
        rl=-std::numeric_limits<double>::infinity();
        ru=+std::numeric_limits<double>::infinity();
        ARIADNE_THROW(DivideByZeroException,"Interval rec(Interval ivl)","ivl="<<i);
    }
    return Interval(rl,ru);
}


Interval mul(Interval i1, Interval i2)
{
    volatile double& i1l=internal_cast<volatile double&>(i1.lower());
    volatile double& i1u=internal_cast<volatile double&>(i1.upper());
    volatile double& i2l=internal_cast<volatile double&>(i2.lower());
    volatile double& i2u=internal_cast<volatile double&>(i2.upper());
    volatile double rl,ru;
    rounding_mode_t rnd=get_rounding_mode();
    if(i1l>=0) {
        if(i2l>=0) {
            rl=_mul_down(i1l,i2l); ru=_mul_up(i1u,i2u);
        } else if(i2u<=0) {
            rl=_mul_down(i1u,i2l); ru=_mul_up(i1l,i2u);
        } else {
            rl=_mul_down(i1u,i2l); ru=_mul_up(i1u,i2u);
        }
    }
    else if(i1u<=0) {
        if(i2l>=0) {
            rl=_mul_down(i1l,i2u); ru=_mul_up(i1u,i2l);
        } else if(i2u<=0) {
            rl=_mul_down(i1u,i2u); ru=_mul_up(i1l,i2l);
        } else {
            rl=_mul_down(i1l,i2u); ru=_mul_up(i1l,i2l);
        }
    } else {
        if(i2l>=0) {
            rl=_mul_down(i1l,i2u); ru=_mul_up(i1u,i2u);
        } else if(i2u<=0) {
            rl=_mul_down(i1u,i2l); ru=_mul_up(i1l,i2l);
        } else {
            set_rounding_mode(downward);
            rl=std::min(i1u*i2l,i1l*i2u);
            set_rounding_mode(upward);
            ru=std::max(i1l*i2l,i1u*i2u);
        }
    }
    set_rounding_mode(rnd);
    return Interval(rl,ru);
}


Interval mul(Interval i1, Float x2)
{
    rounding_mode_t rnd=get_rounding_mode();
    volatile double& i1l=internal_cast<volatile double&>(i1.lower());
    volatile double& i1u=internal_cast<volatile double&>(i1.upper());
    volatile double& x2v=internal_cast<volatile double&>(x2);
    volatile double rl,ru;
    if(x2>=0) {
        rl=_mul_down(i1l,x2v); ru=_mul_up(i1u,x2v);
    } else {
        rl=_mul_down(i1u,x2v); ru=_mul_up(i1l,x2v);
    }
    set_rounding_mode(rnd);
    return Interval(rl,ru);
}


Interval mul(Float x1, Interval i2)
{
    return mul(i2,x1);
}


Interval div(Interval i1, Interval i2)
{
    rounding_mode_t rnd=get_rounding_mode();
    volatile double& i1l=internal_cast<volatile double&>(i1.lower());
    volatile double& i1u=internal_cast<volatile double&>(i1.upper());
    volatile double& i2l=internal_cast<volatile double&>(i2.lower());
    volatile double& i2u=internal_cast<volatile double&>(i2.upper());
    volatile double rl,ru;
    if(i2l>0) {
        if(i1l>=0) {
            rl=_div_down(i1l,i2u); ru=_div_up(i1u,i2l);
        } else if(i1u<=0) {
            rl=_div_down(i1l,i2l); ru=_div_up(i1u,i2u);
        } else {
            rl=_div_down(i1l,i2l); ru=_div_up(i1u,i2l);
        }
    }
    else if(i2u<0) {
        if(i1l>=0) {
            rl=_div_down(i1u,i2u); ru=_div_up(i1l,i2l);
        } else if(i1u<=0) {
            rl=_div_down(i1u,i2l); ru=_div_up(i1l,i2u);
        } else {
            rl=_div_down(i1u,i2u); ru=_div_up(i1l,i2u);
        }
    }
    else {
        // ARIADNE_THROW(DivideByZeroException,"Interval div(Interval ivl1, Interval ivl2)","ivl1="<<i1<<", ivl2="<<i2);
        rl=-std::numeric_limits<double>::infinity();
        ru=+std::numeric_limits<double>::infinity();
    }
    set_rounding_mode(rnd);
    return Interval(rl,ru);
}



Interval div(Interval i1, Float x2)
{
    rounding_mode_t rnd=get_rounding_mode();
    volatile double& i1l=internal_cast<volatile double&>(i1.lower());
    volatile double& i1u=internal_cast<volatile double&>(i1.upper());
    volatile double& x2v=internal_cast<volatile double&>(x2);
    volatile double rl,ru;
    if(x2v>0) {
        rl=_div_down(i1l,x2v); ru=_div_up(i1u,x2v);
    } else if(x2v<0) {
        rl=_div_down(i1u,x2v); ru=_div_up(i1l,x2v);
    } else {
        rl=-std::numeric_limits<double>::infinity();
        ru=+std::numeric_limits<double>::infinity();
    }
    set_rounding_mode(rnd);
    return Interval(rl,ru);
}


Interval div(Float x1, Interval i2)
{
    rounding_mode_t rnd=get_rounding_mode();
    volatile double& x1v=internal_cast<volatile double&>(x1);
    volatile double& i2l=internal_cast<volatile double&>(i2.lower());
    volatile double& i2u=internal_cast<volatile double&>(i2.upper());
    volatile double rl,ru;
    if(i2l<=0 && i2u>=0) {
        ARIADNE_THROW(DivideByZeroException,"Interval div(Float x1, Interval ivl2)","x1="<<x1<<", ivl2="<<i2);
        rl=-std::numeric_limits<double>::infinity();
        ru=+std::numeric_limits<double>::infinity();
    } else if(x1v>=0) {
        rl=_div_down(x1v,i2u); ru=_div_up(x1v,i2l);
    } else {
        rl=_div_down(x1v,i2l); ru=_div_up(x1v,i2u);
    }
    set_rounding_mode(rnd);
    return Interval(rl,ru);
}

Interval sqr(Interval i)
{
    rounding_mode_t rnd=get_rounding_mode();
    volatile double& il=internal_cast<volatile double&>(i.lower());
    volatile double& iu=internal_cast<volatile double&>(i.upper());
    volatile double rl,ru;
    if(il>0.0) {
        set_rounding_mode(downward);
        rl=il*il;
        set_rounding_mode(upward);
        ru=iu*iu;
    } else if(iu<0.0) {
        set_rounding_mode(downward);
        rl=iu*iu;
        set_rounding_mode(upward);
        ru=il*il;
    } else {
        rl=0.0;
        set_rounding_mode(upward);
        volatile double ru1=il*il;
        volatile double ru2=iu*iu;
        ru=max(ru1,ru2);
    }
    set_rounding_mode(rnd);
    return Interval(rl,ru);
}




Interval pow(Interval i, int n)
{
    if(n<0) { return pow(rec(i),uint(-n)); }
    else return pow(i,uint(n));
}

Interval pow(Interval i, uint m)
{
    rounding_mode_t rnd = get_rounding_mode();
    const Interval& nvi=i;
    if(m%2==0) { i=abs(nvi); }
    set_rounding_mode(downward);
    Float rl=pow_rnd(i.lower(),m);
    set_rounding_mode(upward);
    Float ru=pow_rnd(i.upper(),m);
    set_rounding_mode(rnd);
    return Interval(rl,ru);
}



Interval sqrt(Interval i)
{
    rounding_mode_t rnd = get_rounding_mode();
    set_rounding_downward();
    Float rl=sqrt_rnd(i.lower());
    set_rounding_upward();
    Float ru=sqrt_rnd(i.upper());
    set_rounding_mode(rnd);
    return Interval(rl,ru);
}

Interval exp(Interval i)
{
    rounding_mode_t rnd = get_rounding_mode();
    set_rounding_downward();
    Float rl=exp_rnd(i.lower());
    set_rounding_upward();
    Float ru=exp_rnd(i.upper());
    set_rounding_mode(rnd);
    return Interval(rl,ru);
}

Interval log(Interval i)
{
    rounding_mode_t rnd = get_rounding_mode();
    set_rounding_downward();
    Float rl=log_rnd(i.lower());
    set_rounding_upward();
    Float ru=log_rnd(i.upper());
    set_rounding_mode(rnd);
    return Interval(rl,ru);
}


template<> Interval pi<Interval>()
{
    return Interval(pi_down,pi_up);
}


Interval sin(Interval i)
{
    return cos(i-pi<Interval>()/2);
}

Interval cos(Interval i)
{
    ARIADNE_ASSERT(i.lower()<=i.upper());

    rounding_mode_t rnd = get_rounding_mode();

    static const Interval pi(pi_down,pi_up);
    if(i.radius()>2*pi_down) { return Interval(-1.0,+1.0); }

    Float n=floor(i.lower()/(2*pi_approx)+0.5);
    i=i-2*n*pi;

    ARIADNE_ASSERT(i.lower()>=-pi_up);
    ARIADNE_ASSERT(i.lower()<=pi_up);

    Float rl,ru;
    if(i.lower()<=0.0) {
        if(i.upper()<=0.0) { rl=cos_down(i.lower()); ru=cos_up(i.upper()); }
        else if(i.upper()<=pi_down) { rl=cos_down(max(-i.lower(),i.upper())); ru=+1.0; }
        else { return Interval(-1.0,+1.0); }
    } else if(i.lower()<=pi_up) {
        if(i.upper()<=pi_down) { rl=cos_down(i.upper()); ru=cos_up(i.lower()); }
        else if(i.upper()<=2*pi_down) { rl=-1.0; ru=cos_up(min(i.lower(),sub_down(2*pi_down,i.upper()))); }
        else { return Interval(-1.0,+1.0); }
    } else {
        assert(false);
    }

    set_rounding_mode(rnd);

    return Interval(rl,ru);
}

Interval tan(Interval i)
{
    ARIADNE_NOT_IMPLEMENTED;
}

Interval asin(Interval i)
{
    ARIADNE_NOT_IMPLEMENTED;
}

Interval acos(Interval i)
{
    ARIADNE_NOT_IMPLEMENTED;
}

Interval atan(Interval i)
{
    ARIADNE_NOT_IMPLEMENTED;
}



#ifdef HAVE_GMPXX_H


Interval::Interval(const Rational& q) : l(q.get_d()), u(l) {
    rounding_mode_t rounding_mode=get_rounding_mode();
    set_rounding_mode(downward);
    while(l>q) { l-=std::numeric_limits<double>::min(); }
    set_rounding_mode(upward);
    while(u<q) { u+=std::numeric_limits<double>::min(); }
    set_rounding_mode(rounding_mode);
}

Interval::Interval(const Rational& ql, const Rational& qu) : l(ql.get_d()), u(qu.get_d())  {
    rounding_mode_t rounding_mode=get_rounding_mode();
    set_rounding_mode(downward);
    while(l>ql) { l-=std::numeric_limits<double>::min(); }
    set_rounding_mode(upward);
    while(u<qu) { u+=std::numeric_limits<double>::min(); }
    set_rounding_mode(rounding_mode);
}

Interval& Interval::operator=(const Rational& q) {
    return *this = Interval(q);
}


#endif // HAVE_GMPXX_H


std::ostream&
operator<<(std::ostream& os, const Interval& ivl)
{
    if(ivl.lower()==ivl.upper()) { return os << ivl.lower(); }
    return os << '{' << ivl.lower() << ':' << ivl.upper() << '}';
}

/*
std::ostream&
operator<<(std::ostream& os, const Interval& ivl)
{
    return os << '[' << ivl.l << ':' << ivl.u << ']';
}
*/

/*
std::ostream&
operator<<(std::ostream& os, const Interval& ivl)
{
    if(ivl.lower()==ivl.upper()) {
        return os << std::setprecision(18) << ivl.lower();
    }

    std::stringstream iss,uss;
    iss << std::setprecision(18) << ivl.lower();
    uss << std::setprecision(18) << ivl.upper();

    std::string lstr,ustr;
    iss >> lstr; uss >> ustr;

    // Test if one endpoint is an integer and the other is not
    // If this is the case, append ".0" to to integer value
    if( (lstr.find('.')==lstr.size()) xor (ustr.find('.')==lstr.size()) ) {
        if(lstr.find('.')==lstr.size()) {
            lstr+=".0";
        } else {
            ustr+=".0";
        }
    }

    // Write common head
    uint i;
    for(i=0; (i<std::min(lstr.size(),ustr.size()) && lstr[i]==ustr[i]); ++i) {
        os << lstr[i];
    }

    os << "[";
    if(i==lstr.size()) {
        os << "0";
    }
    for(uint li=i; li != lstr.size(); ++li) {
        os << lstr[li];
    }
    os << ":";
    if(i==ustr.size()) {
        os << "0";
    }
    for(uint ui=i; ui != ustr.size(); ++ui) {
        os << ustr[ui];
    }
    os << "]";
    return os;

}
*/
std::istream&
operator>>(std::istream& is, Interval& ivl)
{
    double l,u;
    char cl,cm,cr;
    is >> cl >> l >> cm >> u >> cr;
    ARIADNE_ASSERT(is);
    ARIADNE_ASSERT(cl=='[' || cl=='(');
    ARIADNE_ASSERT(cm==':' || cm==',' || cm==';');
    ARIADNE_ASSERT(cr==']' || cr==')');
    ivl.set(l,u);
    return is;
}



} // namespace Ariadne
