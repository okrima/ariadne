/***************************************************************************
 *            affine_model.cc
 *
 *  Copyright 2009  Pieter Collins
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

#include "algebra/vector.h"
#include "function/function.h"
#include "function/taylor_model.h"
#include "function/affine_model.h"

#include "function/affine.h"
#include "function/taylor_function.h"
#include "algebra/vector.h"

namespace Ariadne {

ValidatedAffineModel operator+(const ValidatedAffineModel& a1, const ValidatedAffineModel& a2) {
    ARIADNE_ASSERT_MSG(a1.argument_size()==a2.argument_size(),"a1="<<a1<<" a2="<<a2);
    Nat n=a1.argument_size();
    ValidatedAffineModel r(n);
    r=CoefficientType( a1.value().raw()+a2.value().raw() );
    for(Nat i=0; i!=n; ++i) {
        r[i]=CoefficientType(a1[i].raw()+a2[i].raw());
    }

    set_rounding_upward();

    RawFloat te=0.0;
    for(Nat j=0; j!=n; ++j) {
        RawFloat mrjl = (-a1.gradient(j).raw())-a2.gradient(j).raw();
        RawFloat  rju = ( a1.gradient(j).raw())+a2.gradient(j).raw();
        te+=(rju+mrjl);
    }
    RawFloat mrl = (-a1.value().raw())-a2.value().raw();
    RawFloat  ru = ( a1.value().raw())+a2.value().raw();
    te += (ru+mrl);

    RawFloat re=0.0;
    r.set_error( ErrorType(te/2 + (a1.error().raw()+a2.error().raw()) ) );

    set_rounding_to_nearest();

    return r;
}

ValidatedAffineModel operator+(const ValidatedNumber& c, const ValidatedAffineModel& a) {
    ValidatedAffineModel r=a;
    RawFloat cm=c.midpoint().raw();
    r.set_value( static_cast<CoefficientType>( cm + a.value().raw() ) );

    set_rounding_upward();

    RawFloat mrl = (-a.value().raw())-cm;
    RawFloat  ru = ( a.value().raw())+cm;
    RawFloat te = (ru+mrl)/2;

    r.set_error( ErrorType( a.error().raw() + max(c.upper().raw()-cm,cm-c.lower().raw()) + te) );

    set_rounding_to_nearest();

    return r;
}

ValidatedAffineModel operator+(const ValidatedAffineModel& a, const ValidatedNumber& c) {
    return c+a;
}

ValidatedAffineModel operator*(const ValidatedNumber& c, const ValidatedAffineModel& a) {
    Nat n=a.argument_size();
    RawFloat cm=c.midpoint().raw();
    ValidatedAffineModel r(n);
    r=CoefficientType(a.value().raw()*cm);
    for(Nat i=0; i!=n; ++i) {
        r[i].raw()=a[i].raw()*cm;
    }

    set_rounding_upward();

    RawFloat te=0.0;
    for(Nat j=0; j!=n; ++j) {
        RawFloat mca=(-cm)*a.gradient(j).raw();
        RawFloat ca= cm*a.gradient(j).raw();
        te+=(ca+mca);
    }
    RawFloat mca=(-cm)*a.value().raw();
    RawFloat ca= cm*a.value().raw();

    RawFloat re=0.0;
    if(c.lower()!=c.upper()) {
        RawFloat ce=max(c.upper().raw()-cm,cm-c.lower().raw());
        for(Nat j=0; j!=n; ++j) {
            re+=abs(a.gradient(j).raw()*ce);
        }
    }

    r.set_error(ErrorType(abs(cm)*a.error().raw() + ((ca+mca) + te)/2 + re));

    set_rounding_to_nearest();

    return r;
}

ValidatedAffineModel operator*(const ValidatedAffineModel& a, const ValidatedNumber& c) {
    return c*a;
}

ValidatedAffineModel affine_model(const ValidatedAffine& a) {
    ValidatedAffineModel am(a.argument_size());
    am = CoefficientType( midpoint(a.value()).raw() );
    for(Nat j=0; j!=a.argument_size(); ++j) {
        am[j] = midpoint(a[j]);
    }
    set_rounding_upward();
    RawFloat e = 0.0;
    for(Nat j=0; j!=a.argument_size(); ++j) {
        e += max(a.gradient(j).upper().raw()-am.gradient(j).raw(),am.gradient(j).raw()-a.gradient(j).lower().raw());
    }
    e += max(a.value().upper().raw()-am.value().raw(),am.value().raw()-a.value().lower().raw());
    am.set_error(ErrorType(e));
    set_rounding_to_nearest();
    return am;
}

ValidatedAffineModel affine_model(const ValidatedTaylorModel& taylor_model) {
    ValidatedAffineModel affine_model(taylor_model.argument_size());

    rounding_mode_t rnd=get_rounding_mode();
    set_rounding_upward();
    for(ValidatedTaylorModel::ConstIterator iter=taylor_model.begin(); iter!=taylor_model.end(); ++iter) {
        if(iter->key().degree()>=2) {
            affine_model.set_error(abs(iter->data()+affine_model.error()));
        } else if(iter->key().degree()==1) {
            for(Nat i=0; i!=taylor_model.argument_size(); ++i) {
                if(iter->key()[i]==1) {
                    affine_model.set_gradient(i,iter->data());
                    break;
                }
            }
        } else {
            affine_model.set_value(iter->data());
        }
    }
    affine_model.set_error(taylor_model.error()+affine_model.error());
    set_rounding_mode(rnd);

    return affine_model;
}

Vector< ValidatedAffineModel > affine_models(const Vector< ValidatedTaylorModel >& models)
{
    Vector< ValidatedAffineModel > result(models.size(),ValidatedAffineModel(models.size()==0?0u:models[0].argument_size()));
    for(Nat i=0; i!=result.size(); ++i) { result[i]=affine_model(models[i]); }
    return result;
}

ValidatedAffineModel affine_model(const ExactBox& domain, const ValidatedScalarFunction& function)
{
    return affine_model(ScalarTaylorFunction(domain,function,AffineSweeper()).model());
}

Vector< ValidatedAffineModel > affine_models(const ExactBox& domain, const ValidatedVectorFunction& function)
{
    return affine_models(VectorTaylorFunction(domain,function,AffineSweeper()).models());
}

OutputStream& operator<<(OutputStream& os, const ValidatedAffineModel& f)
{
    os << f.value().raw();
    for(Nat j=0; j!=f.argument_size(); ++j) {
        if(f.gradient(j)>0.0) { os << "+"; }
        if(f.gradient(j)!=0.0) { os << f.gradient(j) << "*x" << j; }
    }
    return os << "+/-" << f.error();

}

} //namespace Ariadne

