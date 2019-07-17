/***************************************************************************
 *            float_bounds.hpp
 *
 *  Copyright 2008-17  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

/*! \file float_bounds.hpp
 *  \brief Floating-point bounds for real numbers.
 */

#ifndef ARIADNE_FLOAT_BOUNDS_HPP
#define ARIADNE_FLOAT_BOUNDS_HPP

#include "../utility/macros.hpp"

#include "number.decl.hpp"
#include "float.decl.hpp"

#include "float_operations.hpp"
#include "float_factory.hpp"

namespace Ariadne {

template<class F> struct NumericTraits<Bounds<F>> {
    typedef ValidatedNumber GenericType;
    typedef Bounds<F> OppositeType;
    typedef PositiveBounds<F> PositiveType;
    typedef ValidatedKleenean LessType;
    typedef ValidatedKleenean EqualsType;
};

//! \ingroup NumericModule
//! \brief Validated bounds on a number with floating-point endpoints supporting outwardly-rounded arithmetic.
//! \details
//! Note that direct construction from a floating-point number is prohibited, since <c>%FloatDPBounds(3.3)</c> would the singleton interval \f$[3.2999999999999998224,3.2999999999999998224]\f$ (the constant is first interpreted by the C++ compiler to give a C++ \c double, whereas <c>%FloatDPBounds(3.3_decimal)</c> yields the interval \f$[3.2999999999999998224,3.3000000000000002665]\f$ enclosing \f$3.3\f$.
//!
//! Comparison tests on \c FloatBounds use the idea that an interval represents a single number with an unknown value.
//! Hence the result is of type \c ValidatedKleenean, which can take values { \c True, \c False, \c Indeterminate }.
//! Hence a test \f$[\underline{x},\overline{x}]\leq [\underline{y},\overline{y}]\f$ returns \c True if \f$\overline{x}\leq \underline{y}\f$, since in this case \f$x\leq x\f$ whenever \f$x_1\in[\underline{x},\overline{x}]\f$ and \f$y\in[\underline{y},\overline{y}]\f$, \c False if \f$\underline{x}>\overline{y}\f$, since in this case we know \f$x>y\f$, and \c Indeterminate otherwise, since in this case we can find \f$x,y\f$ making the result either true or false.
//! In the case of equality, the comparison \f$[\underline{x},\overline{x}]\f$==\f$[\underline{y},\overline{y}]\f$ only returns \c True if both intervals are singletons, since otherwise we can find values making the result either true of false.
//! To test equality of representation, use \c same(x,y)
//!
//! To obtain the lower and upper bounds of the possible values, use \c x.lower() and \c x.upper().
//! To obtain a best estimate of the value, use \c x.value(), which has an error at most \a x.error().
//! If \f$v\f$ and \f$e\f$ are the returned value and error for the bounds \f$[l,u]\f$, then it is guaranteed that \f$v-e\leq l\f$ and \f$v+e\geq u\f$ in exact arithmetic.
//!
//! To test if the bounds contain a number , use \c models(FloatBounds,FloatValue), and to test if bounds are inconsistent use \c inconsistent(x,y), and to test if \c x provides a better approximation, use \c refines(x,y).
//! \sa Real, FloatDP, FloatMP, FloatValue, FloatBall, FloatUpperBound, FloatLowerBound, FloatApproximation.
//!
//! \par Python interface
//!
//! In the Python interface, %Ariadne validated bounds can be constructed from Python literals of the form \c {a:b} or (deprecated) \c [a,b] .
//! The former is preferred, as it cannot be confused with literals for other classes such as Vector and Array types.
//! Automatic conversion is used to convert FloatBounds literals of the form \c {a,b} to an FloatBounds in functions.
//!
//! Care must be taken when defining intervals using floating-point coefficients, since values are first converted to the nearest
//! representable value by the Python interpreter. <br><br>
//! \code
//!   FloatDPBounds({1.1:2.3}) # Create the interval [1.1000000000000001, 2.2999999999999998]
//!   FloatDPBounds({2.5:4.25}) # Create the interval [2.5, 4.25], which can be represented exactly
//!   FloatDPBounds([2.5,4.25]) # Alternative syntax for creating the interval [2.5, 4.25]
//! \endcode
template<class F> class Bounds
    : public DefineConcreteGenericOperators<Bounds<F>>
    , public DefineFieldOperators<Bounds<F>>
    , public DefineComparisonOperators<Bounds<F>,LessTrait<Bounds<F>>,EqualsTrait<Bounds<F>>>
{
  protected:
    typedef ValidatedTag P; typedef typename F::RoundingModeType RND; typedef typename F::PrecisionType PR;
  public:
    typedef P Paradigm;
    typedef Bounds<F> NumericType;
    typedef Number<P> GenericType;
    typedef F RawType;
    typedef PR PrecisionType;
    typedef PR PropertiesType;
  public:
    Bounds<F>() : _l(0.0), _u(0.0) { }
    explicit Bounds<F>(PrecisionType pr) : _l(0.0,pr), _u(0.0,pr) { }
    explicit Bounds<F>(RawType const& v) : _l(v), _u(v) { }
    explicit Bounds<F>(RawType const& l, RawType const& u) : _l(l), _u(u) { }
    Bounds<F>(LowerBound<F> const& lower, UpperBound<F> const& upper);
    Bounds<F>(LowerBound<F> const& lower, ValidatedUpperNumber const& upper);
    Bounds<F>(ValidatedLowerNumber const& lower, UpperBound<F> const& upper);
    Bounds<F>(ValidatedLowerNumber const& lower, ValidatedUpperNumber const& upper, PR pr);
    template<class N1, class N2, EnableIf<And<IsBuiltinIntegral<N1>,IsBuiltinIntegral<N2>>> = dummy> Bounds<F>(N1 n1, N2 n2, PR pr) : _l(n1,pr), _u(n2,pr) { }
    Bounds<F>(ExactDouble const& dl, ExactDouble const& du, PrecisionType pr) : _l(dl,pr), _u(du,pr) { }
    Bounds<F>(Dyadic const& wl, Dyadic const& wu, PrecisionType pr) : _l(wl,down,pr), _u(wu,up,pr) { }
    Bounds<F>(Rational const& ql, Rational const& qu, PrecisionType pr) : _l(ql,down,pr), _u(qu,up,pr) { }

    template<class FF, EnableIf<IsConstructible<F,FF,RND,PR>> =dummy>
        Bounds<F>(Bounds<FF> const& x, PR pr) : Bounds<F>(F(x.lower_raw(),downward,pr),F(x.upper_raw(),upward,pr)) { }
    Bounds<F>(const Bounds<F>& x, PR pr) : _l(x._l,down,pr), _u(x._u,up,pr) { }

    template<class N, EnableIf<IsBuiltinIntegral<N>> = dummy> Bounds<F>(N n, PR pr) : Bounds<F>(ExactDouble(n),pr) { }
    Bounds<F>(ExactDouble const& d, PR pr) : _l(d,pr), _u(d,pr) { }
        Bounds<F>(TwoExp const& t, PR pr) : _l(t,pr), _u(t,pr) { }
        Bounds<F>(const Integer& z, PR pr) : _l(z,down,pr), _u(z,up,pr) { }
        Bounds<F>(const Dyadic& w, PR pr) : _l(w,down,pr), _u(w,up,pr) { }
        Bounds<F>(const Decimal& d, PR pr) : _l(d,down,pr), _u(d,up,pr) { }
        Bounds<F>(const Rational& q, PR pr) : _l(q,down,pr), _u(q,up,pr) { }
        Bounds<F>(const Real& x, PR pr);
    Bounds<F>(const ValidatedNumber& y, PR pr);


    template<class FE> Bounds<F>(Ball<F,FE> const& x);
    Bounds<F>(Value<F> const& x);

        Bounds<F>& operator=(const Value<F>& x) { return *this=Bounds<F>(x); }
    Bounds<F>& operator=(const ValidatedNumber& y) { return *this=Bounds<F>(y,this->precision()); }
    template<class N, EnableIf<IsBuiltinIntegral<N>> = dummy> Bounds<F>& operator=(N n) { return *this=ValidatedNumber(n); }

    operator ValidatedNumber () const;

    Bounds<F> create(const ValidatedNumber& y) const { return Bounds<F>(y,this->precision()); }

    LowerBound<F> const lower() const;
    UpperBound<F> const upper() const;
    Value<F> const value() const;
    Error<F> const error() const;

    friend Value<F> value(Bounds<F> const& x) { return x.value(); }
    friend Error<F> error(Bounds<F> const& x) { return x.error(); }


    RawType const& lower_raw() const { return _l; }
    RawType const& upper_raw() const { return _u; }
    RawType const value_raw() const { return hlf(add(near,_l,_u)); }
    RawType const error_raw() const { RawType v=value_raw(); return max(sub(up,_u,v),sub(up,v,_l)); }
    double get_d() const { return value_raw().get_d(); }

    PrecisionType precision() const { ARIADNE_DEBUG_ASSERT(_l.precision()==_u.precision()); return _u.precision(); }
    PropertiesType properties() const { ARIADNE_DEBUG_ASSERT(_l.precision()==_u.precision()); return _u.precision(); }
    GenericType generic() const { return this->operator GenericType(); }

    Bounds<F> pm(Error<F> const& e) const;

    // DEPRECATED
    explicit operator RawType () const { return value_raw(); }
    friend Approximation<F> round(Approximation<F> const& x);
    friend Value<F> midpoint(Bounds<F> const& x);
  public:
    friend Bounds<F> max(Bounds<F> const& x1, Bounds<F> const& x2) {
        return Bounds<F>(max(x1.lower_raw(),x2.lower_raw()),max(x1.upper_raw(),x2.upper_raw())); }
    friend Bounds<F> min(Bounds<F> const& x1, Bounds<F> const& x2) {
        return Bounds<F>(min(x1.lower_raw(),x2.lower_raw()),min(x1.upper_raw(),x2.upper_raw())); }
    friend Bounds<F> abs(Bounds<F> const& x) {
        if(x.lower_raw()>=0) { return Bounds<F>(x.lower_raw(),x.upper_raw());}
        else if(x.upper_raw()<=0) { return Bounds<F>(neg(x.upper_raw()),neg(x.lower_raw())); }
        else { return Bounds<F>(F(0.0,x.precision()),max(neg(x.lower_raw()),x.upper_raw())); } }
    friend PositiveLowerBound<F> mig(Bounds<F> const& x) {
        return PositiveLowerBound<F>(max(0,max(x._l,neg(x._u)))); }
    friend PositiveUpperBound<F> mag(Bounds<F> const& x) {
        return PositiveUpperBound<F>(max(neg(x._l),x._u)); }
    friend ValidatedKleenean sgn(Bounds<F> const& x) {
        if (x._l>0) { return true; } else if (x._u<0) { return false; } else { return indeterminate; } }

    friend Bounds<F> nul(Bounds<F> const& x) {
        return Bounds<F>(nul(x._l),nul(x._u)); }
    friend Bounds<F> pos(Bounds<F> const& x) {
        return Bounds<F>(pos(x._l),pos(x._u)); }
    friend Bounds<F> neg(Bounds<F> const& x) {
        return Bounds<F>(neg(x._u),neg(x._l)); }
    friend Bounds<F> hlf(Bounds<F> const& x) {
        return Bounds<F>(hlf(x._l),hlf(x._u)); }
    friend Bounds<F> sqr(Bounds<F> const& x) {
        if(x._l>0.0) { return Bounds<F>(mul(down,x._l,x._l),mul(up,x._u,x._u)); }
        else if(x._u<0.0) { return Bounds<F>(mul(down,x._u,x._u),mul(up,x._l,x._l)); }
        else { return Bounds<F>(nul(x._l),max(mul(up,x._l,x._l),mul(up,x._u,x._u))); } }
    friend Bounds<F> rec(Bounds<F> const& x) {
        if(x._l>0 || x._u<0) {  return Bounds<F>(rec(down,x._u),rec(up,x._l)); }
    //ARIADNE_THROW(DivideByZeroException,"FloatBounds rec(FloatBounds x)","x="<<x);
        else { F inf_=F::inf(x.precision()); return Bounds<F>(-inf_,+inf_); } }

    friend Bounds<F> add(Bounds<F> const& x1, Bounds<F> const& x2) {
        return Bounds<F>(add(down,x1._l,x2._l),add(up,x1._u,x2._u)); }
    friend Bounds<F> sub(Bounds<F> const& x1, Bounds<F> const& x2) {
        return Bounds<F>(sub(down,x1._l,x2._u),sub(up,x1._u,x2._l)); }
    friend Bounds<F> mul(Bounds<F> const& x1, Bounds<F> const& x2) {
        return Operations<Bounds<F>>::_mul(x1,x2); }
    friend Bounds<F> div(Bounds<F> const& x1, Bounds<F> const& x2) {
        return Operations<Bounds<F>>::_div(x1,x2); }
    friend Bounds<F> fma(Bounds<F> const& x1, Bounds<F> const& x2, Bounds<F> const& x3) {
        return Operations<Bounds<F>>::_fma(x1,x2,x3); }
    friend Bounds<F> pow(Bounds<F> const& x, Int n) {
        if(n<0) { return pow(rec(x),Nat(-n)); } else return pow(x,Nat(n));}
    friend Bounds<F> pow(Bounds<F> const& x, Nat m) {
        Bounds<F> y = (m%2==0) ? abs(x) : x;  Int n=static_cast<Int>(m); return Bounds<F>(pow(down,y._l,n),pow(up,y._u,n)); }

    friend Bounds<F> sqrt(Bounds<F> const& x) {
        return Bounds<F>(sqrt(down,x.lower_raw()),sqrt(up,x.upper_raw())); }
    friend Bounds<F> exp(Bounds<F> const& x) {
        return Bounds<F>(exp(down,x.lower_raw()),exp(up,x.upper_raw())); }
    friend Bounds<F> log(Bounds<F> const& x) {
        return Bounds<F>(log(down,x.lower_raw()),log(up,x.upper_raw())); }
    friend Bounds<F> sin(Bounds<F> const& x) {
        return Operations<Bounds<F>>::_sin(x); }
    friend Bounds<F> cos(Bounds<F> const& x) {
        return Operations<Bounds<F>>::_cos(x); }
    friend Bounds<F> tan(Bounds<F> const& x) {
        return mul(sin(x),rec(cos(x))); }
    friend Bounds<F> asin(Bounds<F> const& x) {
        ARIADNE_NOT_IMPLEMENTED; }
    friend Bounds<F> acos(Bounds<F> const& x) {
        ARIADNE_NOT_IMPLEMENTED; }
    friend Bounds<F> atan(Bounds<F> const& x) {
        return Bounds<F>(atan(down,x._l),atan(up,x._u)); }

    //! \related Bounds<F> \brief Strict greater-than comparison operator. Tests equality of represented real-point value.
    friend LogicalType<ValidatedTag> eq(Bounds<F> const& x1, Bounds<F> const& x2) {
        if(x1.upper_raw()<x2.lower_raw() || x1.lower_raw()>x2.upper_raw()) { return false; }
        else if(x1.lower_raw()==x2.upper_raw() && x1.upper_raw() == x2.lower_raw()) { return true; }
        else { return indeterminate; } }
    //! \related Bounds<F> \brief Strict greater-than comparison operator. Tests equality of represented real-point value.
    friend LogicalType<ValidatedTag> lt(Bounds<F> const& x1, Bounds<F> const& x2) {
        if(x1.upper_raw()< x2.lower_raw()) { return true; }
        else if(x1.lower_raw()>=x2.upper_raw()) { return false; }
        else { return indeterminate; } }


    friend Bounds<F> round(Bounds<F> const& x) {
        return Bounds<F>(round(x.lower_raw()),round(x.upper_raw())); }
    friend Bounds<F> widen(Bounds<F> const& x) {
        const F m=std::numeric_limits<float>::min(); return Bounds<F>(sub(down,x._l,m),add(up,x._u,m)); }
    friend Bounds<F> narrow(Bounds<F> const& x) {
        const F m=std::numeric_limits<float>::min(); return Bounds<F>(add(up,x._l,m),add(down,x._u,m)); }
    friend Bounds<F> trunc(Bounds<F> const& x) {
        return Operations<Bounds<F>>::_trunc(x); }
    friend Bounds<F> trunc(Bounds<F> const& x, Nat n) {
        return Operations<Bounds<F>>::_trunc(x,n); }

    friend Integer integer_cast(Bounds<F> const& x) {
        return Integer(static_cast<int>(x.value_raw().get_d())); }

    friend auto is_zero(Bounds<F> const& x) -> LogicalType<ValidatedTag> {
        if(x.lower_raw()>0.0 || x.upper_raw()<0.0) { return false; }
        else if(x.lower_raw()==0.0 && x.upper_raw()==0.0) { return true; }
        else { return indeterminate; } }
    friend auto is_positive(Bounds<F> const& x) -> LogicalType<ValidatedTag> {
        if(x.lower_raw()>=0.0) { return true; } else if(x.upper_raw()<0.0) { return false; } else { return indeterminate; } }

    friend Bool same(Bounds<F> const& x1, Bounds<F> const& x2) {
        return x1._l==x2._l && x1._u==x2._u; }
    friend Bool models(Bounds<F> const& x1, Value<F> const& x2) {
        return x1._l<=x2._v && x1._u >= x2._v; }
    friend Bool consistent(Bounds<F> const& x1, Bounds<F> const& x2) {
        return x1._l<=x2._u && x1._u >= x2._l; }
    friend Bool inconsistent(Bounds<F> const& x1, Bounds<F> const& x2) {
        return x1._l>x2._u || x1._u < x2._l; }
    friend Bool refines(Bounds<F> const& x1, Bounds<F> const& x2) {
        return x1._l>=x2._l && x1._u <= x2._u; }
    friend Bounds<F> refinement(Bounds<F> const& x1, Bounds<F> const& x2) {
        return Bounds<F>(max(x1._l,x2._l),min(x1._u,x2._u)); }

        // FIXME: Added functionality
        friend Bool same(Bounds<F> const& x, Bounds<Dyadic> const& w);
        friend Bool same(Bounds<F> const& x, Dyadic const& w) { return x._l==w && x._u==w; }
        friend Bool models(Bounds<F> const& x, Dyadic const& w) { return x._l<=w && w<=x._u; }
        friend Bounds<F> round(Bounds<F> const&);

    friend OutputStream& operator<<(OutputStream& os, Bounds<F> const& x) { return Operations<Bounds<F>>::_write(os,x); }
    friend InputStream& operator>>(InputStream& is, Bounds<F>& x) { return Operations<Bounds<F>>::_read(is,x); }
  public:
    static Nat output_places;
    static Void set_output_places(Nat p) { output_places=p; }
  private: public:
    RawType _l, _u;
};

template<class PR> Bounds(ValidatedNumber, PR) -> Bounds<RawFloatType<PR>>;
template<class PR> Bounds(ValidatedLowerNumber, ValidatedUpperNumber, PR) -> Bounds<RawFloatType<PR>>;
template<class F> Bounds(F,F) -> Bounds<F>;


template<class F> inline FloatFactory<PrecisionType<F>> factory(Bounds<F> const& flt) { return FloatFactory<PrecisionType<F>>(flt.precision()); }
template<class PR> inline FloatBounds<PR> FloatFactory<PR>::create(Number<ValidatedTag> const& y) { return FloatBounds<PR>(y,_pr); }
template<class PR> inline FloatBounds<PR> FloatFactory<PR>::create(Number<EffectiveTag> const& y) { return FloatBounds<PR>(y,_pr); }
template<class PR> inline FloatBounds<PR> FloatFactory<PR>::create(Number<ExactTag> const& y) { return FloatBounds<PR>(y,_pr); }
template<class PR> inline FloatBounds<PR> FloatFactory<PR>::create(Real const& y) { return FloatBounds<PR>(y,_pr); }
template<class PR> inline FloatBounds<PR> FloatFactory<PR>::create(Rational const& y) { return FloatBounds<PR>(y,_pr); }
template<class PR> inline FloatBounds<PR> FloatFactory<PR>::create(Dyadic const& y) { return FloatBounds<PR>(y,_pr); }
template<class PR> inline FloatBounds<PR> FloatFactory<PR>::create(Integer const& y) { return FloatBounds<PR>(y,_pr); }


template<class F> class Positive<Bounds<F>> : public Bounds<F>
    , public DeclarePositiveFloatOperations<PositiveBounds<F>>
{
    using typename Bounds<F>::PR;
  public:
    Positive<Bounds<F>>() : Bounds<F>() { }
    explicit Positive<Bounds<F>>(PR const& pr) : Bounds<F>(pr) { }
    template<class M, EnableIf<IsBuiltinUnsignedIntegral<M>> =dummy> Positive<Bounds<F>>(M m, PR pr) : Bounds<F>(m,pr) { }
    explicit Positive<Bounds<F>>(F const& x) : Bounds<F>(x) { }
    explicit Positive<Bounds<F>>(F const& l, F const& u) : Bounds<F>(l,u) { }
    explicit Positive<Bounds<F>>(Bounds<F> const& x) : Bounds<F>(x) { }
    Positive<Bounds<F>>(Positive<LowerBound<F>> const& xl, Positive<UpperBound<F>> const& xu) : Bounds<F>(xl,xu) { }
  public:
    Positive<Value<F>> value() const { return cast_positive(this->Bounds<F>::value()); }
    Positive<LowerBound<F>> lower() const { return cast_positive(this->Bounds<F>::lower()); }
    Positive<UpperBound<F>> upper() const { return cast_positive(this->Bounds<F>::upper()); }
  public:
    friend PositiveBounds<F> nul(PositiveBounds<F> const& x) {
        return PositiveBounds<F>(nul(x._l),nul(x._u)); }
    friend PositiveBounds<F> hlf(PositiveBounds<F> const& x) {
        return PositiveBounds<F>(hlf(x._l),hlf(x._u)); }
    friend PositiveBounds<F> sqr(PositiveBounds<F> const& x) {
        return PositiveBounds<F>(sqr(down,x._l),sqr(up,x._u)); }
    friend PositiveBounds<F> rec(PositiveBounds<F> const& x) {
        return PositiveBounds<F>(rec(down,x._u),rec(up,x._l)); }
    friend PositiveBounds<F> add(PositiveBounds<F> const& x1, PositiveBounds<F> const& x2) {
        return PositiveBounds<F>(add(down,x1._l,x2._l),add(up,x1._u,x2._u)); }
    friend PositiveBounds<F> mul(PositiveBounds<F> const& x1, PositiveBounds<F> const& x2) {
        return PositiveBounds<F>(mul(down,x1._l,x2._l),mul(up,x1._u,x2._u)); }
    friend PositiveBounds<F> div(PositiveBounds<F> const& x1, PositiveBounds<F> const& x2) {
        return PositiveBounds<F>(div(down,x1._l,x2._u),div(up,x1._u,x2._l)); }
    friend PositiveBounds<F> pow(PositiveBounds<F> const& x, Nat m) {
        return PositiveBounds<F>(pow(down,x._l,static_cast<Int>(m)),pow(up,x._u,static_cast<Int>(m))); }
    friend PositiveBounds<F> sqrt(PositiveBounds<F> const& x) {
        return PositiveBounds<F>(sqrt(down,x._l),sqrt(up,x._u)); }
    friend PositiveBounds<F> atan(PositiveBounds<F> const& x) {
        return PositiveBounds<F>(atan(down,x._l),atan(up,x._u)); }
    friend PositiveBounds<F> max(PositiveBounds<F> const& x1, PositiveBounds<F> const& x2) {
        return PositiveBounds<F>(max(down,x1._l,x2._l),max(up,x1._u,x2._u)); }
    friend PositiveBounds<F> max(PositiveBounds<F> const& x1, Bounds<F> const& x2) {
        return PositiveBounds<F>(max(down,x1._l,x2._l),max(up,x1._u,x2._u)); }
    friend PositiveBounds<F> max(Bounds<F> const& x1, PositiveBounds<F> const& x2) {
        return PositiveBounds<F>(max(down,x1._l,x2._l),max(up,x1._u,x2._u)); }
    friend PositiveBounds<F> min(PositiveBounds<F> const& x1, PositiveBounds<F> const& x2) {
        return PositiveBounds<F>(min(down,x1._l,x2._l),min(up,x2._u,x2._u)); }
    friend PositiveBounds<F> abs(PositiveBounds<F> const& x) { return PositiveBounds<F>(abs(x._l),abs(x._u)); }

    friend Bounds<F> const& cast_unsigned(PositiveBounds<F> const& x) { return x; }
};

template<class F> inline PositiveBounds<F> cast_positive(Bounds<F> const& x) {
    return PositiveBounds<F>(x); }

extern template Ariadne::Nat Ariadne::Bounds<Ariadne::FloatDP>::output_places;
extern template Ariadne::Nat Ariadne::Bounds<Ariadne::FloatMP>::output_places;


template<class F> class Operations<Bounds<F>> {
    typedef PrecisionType<F> PR;
  public:
    static Bounds<F> _mul(Bounds<F> const& x1, Bounds<F> const& x2) {
        const F& x1l=x1._l; const F& x1u=x1._u;
        const F& x2l=x2._l; const F& x2u=x2._u;
        F rl,ru;
        if(x1l>=0) {
            if(x2l>=0) {
                rl=mul(down,x1l,x2l); ru=mul(up,x1u,x2u);
            } else if(x2u<=0) {
                rl=mul(down,x1u,x2l); ru=mul(up,x1l,x2u);
            } else {
                rl=mul(down,x1u,x2l); ru=mul(up,x1u,x2u);
            }
        }
        else if(x1u<=0) {
            if(x2l>=0) {
                rl=mul(down,x1l,x2u); ru=mul(up,x1u,x2l);
            } else if(x2u<=0) {
                rl=mul(down,x1u,x2u); ru=mul(up,x1l,x2l);
            } else {
                rl=mul(down,x1l,x2u); ru=mul(up,x1l,x2l);
            }
        } else {
            if(x2l>=0) {
                rl=mul(down,x1l,x2u); ru=mul(up,x1u,x2u);
            } else if(x2u<=0) {
                rl=mul(down,x1u,x2l); ru=mul(up,x1l,x2l);
            } else {
                rl=min(mul(down,x1u,x2l),mul(down,x1l,x2u));
                ru=max(mul(up,x1l,x2l),mul(up,x1u,x2u));
            }
        }
        return Bounds<F>(rl,ru);
    }

    static Bounds<F> _div(Bounds<F> const& x1, Bounds<F> const& x2) {
        const F& x1l=x1.lower_raw(); const F& x1u=x1.upper_raw();
        const F& x2l=x2.lower_raw(); const F& x2u=x2.upper_raw();
        F rl,ru;

        // IMPORTANT: Need to be careful when one of the bounds is 0, since if x2l=-0.0 and x1u>0, then x2l>=0 but x1u/x2l=-inf
        if(x2l>0) {
            if(x1l>=0) {
                rl=div(down,x1l,x2u); ru=div(up,x1u,x2l);
            } else if(x1u<=0) {
                rl=div(down,x1l,x2l); ru=div(up,x1u,x2u);
            } else {
                rl=div(down,x1l,x2l); ru=div(up,x1u,x2l);
            }
        }
        else if(x2u<0) {
            if(x1l>=0) {
                rl=div(down,x1u,x2u); ru=div(up,x1l,x2l);
            } else if(x1u<=0) {
                rl=div(down,x1u,x2l); ru=div(up,x1l,x2u);
            } else {
                rl=div(down,x1u,x2u); ru=div(up,x1l,x2u);
            }
        }
        else {
            //ARIADNE_THROW(DivideByZeroException,"FloatBounds div(FloatBounds x1, FloatBounds x2)","x1="<<x1<<", x2="<<x2);
            PR pr=max(x1.precision(),x2.precision());
            rl=-F::inf(pr);
            ru=+F::inf(pr);
        }
        return Bounds<F>(rl,ru);
    }

    static Bounds<F> _pi(PR pr) {
        return Bounds<F>(F::pi(down,pr),F::pi(up,pr));
    }

    static Bounds<F> _sin(Bounds<F> const& x) {
        return cos(x-hlf(_pi(x.precision())));
    }

    static Bounds<F> _cos(Bounds<F> const& x) {
        ARIADNE_ASSERT(x.lower_raw()<=x.upper_raw());
        typename F::RoundingModeType rnd = F::get_rounding_mode();
        PR prec=x.precision();

        const F one(1,prec);
        const Value<F> two(2,prec);
        const Bounds<F> pi_val=_pi(prec);
        const Bounds<F> two_pi_val=2*pi_val;
        if(x.error().raw()>two_pi_val.lower().raw()) { return Bounds<F>(-one,+one); }

        Value<F> n(round(div(near,x.value_raw(),(two_pi_val.value_raw()))));
        Bounds<F> y=x-two*(n*pi_val);

        ARIADNE_ASSERT(y.lower_raw()<=pi_val.upper_raw());
        ARIADNE_ASSERT(y.upper_raw()>=-pi_val.upper_raw());

        F rl,ru;
        if(y.lower_raw()<=-pi_val.lower_raw()) {
            if(y.upper_raw()<=0.0) { rl=-one; ru=cos(up,y.upper_raw()); }
            else { rl=-one; ru=+one; }
        } else if(y.lower_raw()<=0.0) {
            if(y.upper_raw()<=0.0) { rl=cos(down,y.lower_raw()); ru=cos(up,y.upper_raw()); }
            else if(y.upper_raw()<=pi_val.lower_raw()) { rl=cos(down,max(-y.lower_raw(),y.upper_raw())); ru=+one; }
            else { rl=-one; ru=+one; }
        } else if(y.lower_raw()<=pi_val.upper_raw()) {
            if(y.upper_raw()<=pi_val.lower_raw()) { rl=cos(down,y.upper_raw()); ru=cos(up,y.lower_raw()); }
            else if(y.upper_raw()<=two_pi_val.lower_raw()) { rl=-one; ru=cos(up,min(y.lower_raw(),sub(down,two_pi_val.lower_raw(),y.upper_raw()))); }
            else { rl=-one; ru=+one; }
        } else {
            assert(false);
        }

        F::set_rounding_mode(rnd);
        return Bounds<F>(rl,ru);
    }

    static Bounds<F> _trunc(Bounds<F> const& x) {
        typename F::RoundingModeType rm=F::get_rounding_mode();
        const double& xl=x.lower_raw().get_d();
        const double& xu=x.upper_raw().get_d();
        // Use machine epsilon instead of minimum to move away from zero
        const float fm=std::numeric_limits<float>::epsilon();
        volatile float tu=xu;
        if(tu<xu) { F::set_rounding_upward(); tu+=fm; }
        volatile float tl=xl;
        if(tl>xl) { F::set_rounding_downward(); tl-=fm; }
        F::set_rounding_mode(rm);
        assert(tl<=xl); assert(tu>=xu);
        return Bounds<F>(double(tl),double(tu));
    }

    static Bounds<F> _trunc(Bounds<F> const& x, Nat n) {
        Bounds<F> _e=Bounds<F>(std::pow(2.0,52-(Int)n));
        Bounds<F> y=x+_e; return y-_e;
    }

    // Mixed BoundedTag - ExactTag operations
    static Bounds<F> _add(Bounds<F> const& x1, Value<F> const& x2) {
        return Bounds<F>(add(down,x1._l,x2._v),add(up,x1._u,x2._v));
    }

    static Bounds<F> _add(Value<F> const& x1, Bounds<F> const& x2) {
        return Bounds<F>(add(down,x1._v,x2._l),add(down,x1._v,x2._u));
    }

    static Bounds<F> _sub(Bounds<F> const& x1, Value<F> const& x2) {
        return Bounds<F>(sub(down,x1._l,x2._v),sub(up,x1._u,x2._v));
    }

    static Bounds<F> _sub(Value<F> const& x1, Bounds<F> const& x2) {
        return Bounds<F>(sub(down,x1._v,x2._u),sub(up,x1._v,x2._l));
    }

    static Bounds<F> _mul(Bounds<F> const& x1, Value<F> const& x2) {
        const F& x1l=x1.lower_raw(); const F& x1u=x1.upper_raw();
        const F& x2v=x2.raw();
        F rl,ru;
        if(x2v>=0.0) {
            rl=mul(down,x1l,x2v); ru=mul(up,x1u,x2v);
        } else {
            rl=mul(down,x1u,x2v); ru=mul(up,x1l,x2v);
        }
        return Bounds<F>(rl,ru);
    }

    static Bounds<F> _mul(Value<F> const& x1, Bounds<F> const& x2) {
        const F& x1v=x1.raw();
        const F& x2l=x2.lower_raw(); const F& x2u=x2.upper_raw();
        F rl,ru;
        if(x1v>=0.0) {
            rl=mul(down,x1v,x2l); ru=mul(up,x1v,x2u);
        } else {
            rl=mul(down,x1v,x2u); ru=mul(up,x1v,x2l);
        }
        return Bounds<F>(rl,ru);
    }

    static Bounds<F> _div(Bounds<F> const& x1, Value<F> const& x2) {
        const F& x1l=x1.lower_raw();
        const F& x1u=x1.upper_raw();
        const F& x2v=x2.raw();
        F rl,ru;
        if(x2v>0) {
            rl=div(down,x1l,x2v); ru=div(up,x1u,x2v);
        } else if(x2v<0) {
            rl=div(down,x1u,x2v); ru=div(up,x1l,x2v);
        } else {
            //ARIADNE_THROW(DivideByZeroException,"FloatBounds div(FloatBounds const& x1, FloatValue x2)","x1="<<x1<<", x2="<<x2);
            auto pr=min(x1.precision(),x2.precision());
            rl=-F::inf(pr);
            ru=+F::inf(pr);
        }
        return Bounds<F>(rl,ru);
    }

    static Bounds<F> _div(Value<F> const& x1, Bounds<F> const& x2) {
        const F& x1v=x1.raw();
        const F& i2l=x2.lower_raw();
        const F& i2u=x2.upper_raw();
        F rl,ru;
        if(i2l<=0 && i2u>=0) {
            //ARIADNE_THROW(DivideByZeroException,"FloatBounds div(FloatValue const& x1, FloatBounds x2)","x1="<<x1<<", x2="<<x2);
            auto pr=min(x1.precision(),x2.precision());
            rl=-F::inf(pr);
            ru=+F::inf(pr);
        } else if(x1v>=0) {
            rl=div(down,x1v,i2u); ru=div(up,x1v,i2l);
        } else {
            rl=div(down,x1v,i2l); ru=div(up,x1v,i2u);
        }
        return Bounds<F>(rl,ru);
    }

    static OutputStream& _write(OutputStream& os, const Bounds<F>& x) {
        os << '{';
        write(os,x.lower().raw(),Bounds<F>::output_places,downward);
        os << ':';
        write(os,x.upper().raw(),Bounds<F>::output_places,upward);
        os << '}';
        return os;

    }

    static InputStream& _read(InputStream& is, Bounds<F>& x) {
        char cl,cm,cr;
        F _l,_u;
        auto rnd=F::get_rounding_mode();
        is >> cl;
        F::set_rounding_downward();
        is >> _l;
        is >> cm;
        F::set_rounding_upward();
        is >> _u;
        is >> cr;
        F::set_rounding_mode(rnd);
        ARIADNE_ASSERT(not is.fail());
        ARIADNE_ASSERT(cl=='[' || cl=='(');
        ARIADNE_ASSERT(cm==':' || cm==',' || cm==';');
        ARIADNE_ASSERT(cr==']' || cr==')');
        x._l=_l; x._u=_u;
        return is;
    }
};


}


#endif
