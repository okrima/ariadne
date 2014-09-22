/***************************************************************************
 *            float-validated.h
 *
 *  Copyright 2008-14  Pieter Collins
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

/*! \file float-validated.h
 *  \brief Floating-point lower and upper bounds for a number.
 */
#ifndef ARIADNE_FLOAT_VALIDATED_H
#define ARIADNE_FLOAT_VALIDATED_H

#ifdef HAVE_GMPXX_H
#include <gmpxx.h>
#endif // HAVE_GMPXX_H

#include <iostream>
#include <cassert>

#include "tribool.h"
#include "rounding.h"
#include "float.h"
#include "float-exact.h"

// Simplifying typedefs for unsigned types
typedef unsigned int uint;

namespace Ariadne {

// Forward declarations
class Float;
class ValidatedFloat;
class ExactFloat;

class Real;

class Integer;
class Rational;
class Dyadic;
class Decimal;

class LowerFloat;
class UpperFloat;

//! \ingroup NumericModule
//! \brief Floating-point lower bounds for real numbers.
class LowerFloat {
  public:
    typedef LowerFloat NumericType;
  public:
    //! \brief Default constructor yields a lower bound of \a 0.
    LowerFloat() : l(0.0) { }
    //! \brief Convert from a builtin double-precision floating-point value.
    LowerFloat(double x) : l(x) { }
    //! \brief Explicitly construct from a raw floating-point value.
    explicit LowerFloat(Float x) : l(x) { }
    //! \brief Convert from floating-point bounds on a number.
    inline LowerFloat(const ValidatedFloat& x);
    //! \brief Convert from a floating-point number with an exact representation.
    LowerFloat(const ExactFloat& x) : l(x.value()) { }
    //! \brief Explicitly convert to the raw floating-point value.
    explicit operator Float const& () const { return l; }
    //! \brief Get the raw value.
    Float const& value() const { return l; }
    //! \brief Get the value to double-precision.
    double get_d() const { return l.get_d(); }
    friend UpperFloat operator-(LowerFloat);
    friend LowerFloat operator+(LowerFloat, LowerFloat);
    friend LowerFloat operator-(LowerFloat, UpperFloat);
    friend UpperFloat operator-(UpperFloat, LowerFloat);
    friend std::ostream& operator<<(std::ostream& os, LowerFloat);
  private:
    Float l;
};

//! \ingroup NumericModule
//! \brief Floating-point upper bounds for real numbers.
class UpperFloat {
  public:
    typedef UpperFloat NumericType;
  public:
    //! \brief Default constructor yields an upper bound of \a 0.
    UpperFloat() : u(0.0) { }
    //! \brief Convert from a builtin double-precision floating-point value.
    UpperFloat(double x) : u(x) { }
    //! \brief Explicitly construct from a raw floating-point value.
    UpperFloat(Float x) : u(x) { }
    //! \brief Convert from floating-point bounds on a number.
    inline UpperFloat(const ValidatedFloat& x);
    //! \brief Convert from a floating-point number with an exact representation.
    UpperFloat(const ExactFloat& x) : u(x.value()) { }
    //! \brief Explicitly convert to the raw floating-point value.
    explicit operator Float const& () const { return u; }
    //! \brief Get the raw value.
    Float const& value() const { return u; }
    //! \brief Get the value to double-precision.
    double get_d() const { return u.get_d(); }
    friend LowerFloat operator-(UpperFloat);
    friend UpperFloat operator+(UpperFloat, UpperFloat);
    friend UpperFloat operator-(UpperFloat, LowerFloat);
    friend LowerFloat operator-(LowerFloat, UpperFloat);
    friend std::ostream& operator<<(std::ostream& os, UpperFloat);
  private:
    Float u;
};

//! \ingroup NumericModule
//! \brief Intervals with floating-point endpoints supporting outwardly-rounded arithmetic.
//! \details
//! Note that <c>%ValidatedFloat(3.3)</c> yields the singleton interval \f$[3.2999999999999998224,3.2999999999999998224]\f$ (the constant is first interpreted by the C++ compiler to give a C++ \c double, whereas <c>%ValidatedFloat("3.3")</c> yields the interval \f$[3.2999999999999998224,3.3000000000000002665]\f$ enclosing \f$3.3\f$.
//!
//! Comparison tests on \c ValidatedFloat use the idea that an interval represents a single number with an unknown value.
//! Hence the result is of type \c tribool, which can take values { \c True, \c False, \c Indeterminate }.
//! Hence a test \f$[l_1,u_1]\leq [l_2,u_2]\f$ returns \c True if \f$u_1\leq u_2\f$, since in this case \f$x_1\leq x_2\f$ whenever \f$x_1\in[l_1,u_2]\f$ and \f$x_2\in[l_2,u_2]\f$, \c False if \f$l_1>u_2\f$, since in this case we know \f$x_1>x_2\f$, and \c Indeterminate otherwise, since in this case we can find \f$x_1,x_2\f$ making the result either true or false.
//! In the case of equality, the comparison \f$[l_1,u_1]\f$==\f$[l_2,u_2]\f$ only returns \c True if both intervals are singletons, since otherwise we can find values making the result either true of false.
//!
//! To obtain the lower and upper bounds of an interval, use \c ivl.lower() and \c ivl.upper().
//! To obtain the midpoint and radius, use \c ivl.midpoint() and \c ivl.radius().
//! Alternatives \c midpoint(ivl) and \c radius(ivl) are also provided.
//! Note that \c midpoint and \c radius return approximations to the true midpoint and radius of the interval. If \f$m\f$ and \f$r\f$ are the returned midpoint and radius of the interval \f$[l,u]\f$, the using exact arithmetic, we guarentee \f$m-r\leq l\f$ and \f$m+r\geq u\f$
//!
//! To test if an interval contains a point or another interval, use \c encloses(ValidatedFloat,Float) or \c encloses(ValidatedFloat,ValidatedFloat).
//! The test \c refines(ValidatedFloat,ValidatedFloat) can also be used.
//! \sa Float
//!
//! \par Python interface
//!
//! In the Python interface, %Ariadne intervals can be constructed from Python literals of the form \c {a:b} or (deprecated) \c [a,b] .
//! The former is preferred, as it cannot be confused with literals for other classes such as Vector and Array types.
//! Automatic conversion is used to convert ValidatedFloat literals of the form \c {a,b} to an ValidatedFloat in functions.
//!
//! Care must be taken when defining intervals using floating-point coefficients, since values are first converted to the nearest
//! representable value by the Python interpreter. <br><br>
//! \code
//!   ValidatedFloat({1.1:2.3}) # Create the interval [1.1000000000000001, 2.2999999999999998]
//!   ValidatedFloat({2.5:4.25}) # Create the interval [2.5, 4.25], which can be represented exactly
//!   ValidatedFloat([2.5,4.25]) # Alternative syntax for creating the interval [2.5, 4.25]
//! \endcode
class ValidatedFloat {
  public:
    typedef ValidatedFloat NumericType;
  public:
    //! \brief Default constructor yields the singleton zero interval \a [0,0].
    ValidatedFloat() : l(0.0), u(0.0) { }
    ValidatedFloat(uint m) : l(m), u(m) { }
    ValidatedFloat(int n) : l(n), u(n) { }
    //! \brief Convert from a builtin double-precision floating-point value. Yields the singleton interval \a [x,x].
    ValidatedFloat(double x) : l(x), u(x) { }
    //! \brief Create from a floating-point value. Yields the singleton interval \a [x,x].
    //! Cannot be used in conversions since the \c %ValidatedFloat class provides stronger accuracy guarantees than the \c %Float class.
    explicit ValidatedFloat(const Float& x) : l(x), u(x) { }
    //! \brief Copy constructor.
    ValidatedFloat(const ValidatedFloat& i) : l(i.l), u(i.u) { }
    //! \brief Convert from a general real number. Yields an interval containing the exact value.
    ValidatedFloat(const Real& x);
    //! \brief Convert from a floating-point number with an exact representation.
    ValidatedFloat(const ExactFloat& x) : l(x.value()), u(x.value()) { }
    //! \brief Convert from a dyadic number.
    ValidatedFloat(const Dyadic& x);
    //! \brief Convert from a decimal number.
    ValidatedFloat(const Decimal& x);

    //! \brief Convert to a floating-point approximation.
    operator Float () const { return this->midpoint(); }

    //! \brief Create from explicitly given lower and upper bounds. Yields the interval \a [lower,upper].
    ValidatedFloat(double lower, double upper) : l(lower), u(upper) { }
    //! \brief Create from explicitly given lower and upper bounds. Yields the interval \a [lower,upper].
    ValidatedFloat(const Float& lower, const Float& upper) : l(lower), u(upper) { }
        // ARIADNE_ASSERT_MSG(lower<=upper, "lower = "<<lower<<", upper ="<<upper);
    //! \brief Create from explicitly given lower and upper bounds. Yields the interval \a [lower,upper].
    ValidatedFloat(const LowerFloat& lower, const UpperFloat& upper) : l(lower.value()), u(upper.value()) { }
        // ARIADNE_ASSERT_MSG(lower<=upper, "lower = "<<lower<<", upper ="<<upper);
#ifdef HAVE_GMPXX_H
    ValidatedFloat(const Integer& z);
    ValidatedFloat(const Rational& q);
    ValidatedFloat& operator=(const Rational& q);
    ValidatedFloat(const Rational& lower, const Rational& upper);
#endif // HAVE_GMPXX_H

    ValidatedFloat& operator=(uint m) { l=m; u=m; return *this; }
    ValidatedFloat& operator=(int n) { l=n; u=n; return *this; }
    ValidatedFloat& operator=(double c) { l=c; u=c; return *this; }
    ValidatedFloat& operator=(const Float& x) { l=x; u=x; return *this; }
    ValidatedFloat& operator=(const Real& x);
    ValidatedFloat& operator=(const ExactFloat& x) { l=x.value(); u=x.value(); return *this; };

    //! \brief The lower bound of the interval.
    LowerFloat lower() const { return LowerFloat(u); }
    //! \brief The upper bound of the interval.
    UpperFloat upper() const { return UpperFloat(u); }

    //! \brief The lower bound of the interval.
    const Float& lower_value() const { return l; }
    //! \brief The upper bound of the interval.
    const Float& upper_value() const { return u; }
    //! \brief An approximation to the midpoint of the interval.
    const Float midpoint() const { return half_exact(add_approx(l,u)); }
    //! \brief An over-approximation to the radius of the interval.
    const Float radius() const { return half_exact(sub_up(u,l)); }
    //! \brief An over-approximation to the width of the interval.
    const Float width() const { return sub_up(u,l); }

    //! \brief Tests if the interval is empty.
    bool empty() const { return l>u; }
    //! \brief Tests if the interval is a singleton.
    bool singleton() const { return l==u; }

    //! \brief Sets the interval to a "canonical" empty interval \a [1,0].
    void set_empty() { l=+std::numeric_limits< double >::infinity(); u=-std::numeric_limits< double >::infinity(); }
    void set_lower(const Float& lower) { l=lower; }
        // ARIADNE_ASSERT(lower<=this->u);
    void set_upper(const Float& upper) { u=upper; }
        // ARIADNE_ASSERT(this->l<=upper);
    void set(const Float& lower, const Float& upper) { l=lower; u=upper; }
        // ARIADNE_ASSERT(lower<=upper);
  public:
    //! \brief Extract a double-precision point approximation to the value represented by the interval.
    double get_d() const { return (this->l.get_d()+this->u.get_d())/2; }
    static uint output_precision;
    static void set_output_precision(uint p) { output_precision=p; }
  private:
    Float l, u;
};

std::ostream& operator<<(std::ostream& os, const ValidatedFloat& ivl);

inline LowerFloat::LowerFloat(ValidatedFloat const& i) : l(i.lower_value()) {
}

inline UpperFloat::UpperFloat(ValidatedFloat const& i) : u(i.upper_value()) {
}

inline Float midpoint(ValidatedFloat i) {
    return half_exact(add_approx(i.lower_value(),i.upper_value()));
}

inline Float radius(ValidatedFloat i) {
    return half_exact(sub_up(i.upper_value(),i.lower_value()));
}

inline Float width(ValidatedFloat i) {
    return sub_up(i.upper_value(),i.lower_value());
}

//! \related ValidatedFloat \brief Test if the intervals are equal (as sets).
inline bool equal(ValidatedFloat i1, ValidatedFloat i2) {
    //std::cerr<<"equal(i1,i2) with i1="<<i1<<"; i2="<<i2<<std::endl;
    return i1.lower_value()==i2.lower_value() && i1.upper_value()==i2.upper_value();
}

// An interval one ulp wider
//! \related ValidatedFloat \brief An interval containing the given interval in its interior.
ValidatedFloat widen(ValidatedFloat i);
//! \related ValidatedFloat \brief An interval contained in the interior of the given interval.
ValidatedFloat narrow(ValidatedFloat i);

// Over-approximate by an interval with float coefficients
//! \related ValidatedFloat \brief Over-approximate the interval by one using builtin single-precision floating-point values as endpoints.
ValidatedFloat trunc(ValidatedFloat);
ValidatedFloat trunc(ValidatedFloat, uint eps);

//! \related ValidatedFloat \brief The interval of possible maximum values. Yields the interval between \c i1.upper() and \c i2.upper().
inline ValidatedFloat max(ValidatedFloat i1,ValidatedFloat i2);
//! \related ValidatedFloat \brief The interval of possible minimum values. Yields the interval between \c i1.lower() and \c i2.lower().
inline ValidatedFloat min(ValidatedFloat,ValidatedFloat);
//! \related ValidatedFloat \brief The interval of possible absolute values. Yields \f$\{ |x| \mid x\in I\}\f$.
inline ValidatedFloat abs(ValidatedFloat);

//! \related ValidatedFloat \brief Unary plus function. Yields the identity \f$I=\{+x | x\in I\}\f$.
inline ValidatedFloat pos(ValidatedFloat i);
//! \related ValidatedFloat \brief Unary negation function. Yields the exact interval \f$\{-x | x\in I\}\f$.
inline ValidatedFloat neg(ValidatedFloat i);
//! \related ValidatedFloat \brief Unary square function. Yields an over-approximation to \f$\{ x^2 \mid x\in I\}\f$.
//! Note that if \a I contains positive and negative values, \c sqr(I) is tighter than \c I*I .
ValidatedFloat sqr(ValidatedFloat i);
//! \related ValidatedFloat \brief Unary reciprocal function. Yields an over-approximation to \f$\{ 1/x \mid x\in I\}\f$.
//! Yields \f$[-\infty,+\infty]\f$ if \a I contains \a 0 in its interior, and an interval containing \f$[1/u,+\infty]\f$ if \a I=[0,u] .
ValidatedFloat rec(ValidatedFloat i);

//! \related ValidatedFloat \brief Binary addition function. Yields an over-approximation to \f$\{ x_1+x_2 \mid x_1\in I_1 \wedge x_2\in I_2\}\f$.
inline ValidatedFloat add(ValidatedFloat, ValidatedFloat);
//! \related ValidatedFloat \brief Subtraction function. Yields an over-approximation to \f$\{ x_1-x_2 \mid x_1\in I_1 \wedge x_2\in I_2\}\f$.
inline ValidatedFloat sub(ValidatedFloat, ValidatedFloat);
//! \related ValidatedFloat \brief Binary multiplication function. Yields an over-approximation to \f$\{ x_1\times x_2 \mid x_1\in I_1 \wedge x_2\in I_2\}\f$.
ValidatedFloat mul(ValidatedFloat, ValidatedFloat);
//! \related ValidatedFloat \brief Division function. Yields an over-approximation to \f$\{ x_1 \div x_2 \mid x_1\in I_1 \wedge x_2\in I_2\}\f$.
ValidatedFloat div(ValidatedFloat, ValidatedFloat);

inline ValidatedFloat add(ValidatedFloat, ExactFloat);
inline ValidatedFloat add(ExactFloat, ValidatedFloat);
inline ValidatedFloat sub(ValidatedFloat, ExactFloat);
inline ValidatedFloat sub(ExactFloat, ValidatedFloat);
ValidatedFloat mul(ValidatedFloat, ExactFloat);
ValidatedFloat mul(ExactFloat,ValidatedFloat);
ValidatedFloat div(ValidatedFloat, ExactFloat);
ValidatedFloat div(ExactFloat, ValidatedFloat);

inline ValidatedFloat rec(ExactFloat);
inline ValidatedFloat add(ExactFloat, ExactFloat);
inline ValidatedFloat sub(ExactFloat, ExactFloat);
inline ValidatedFloat mul(ExactFloat, ExactFloat);
inline ValidatedFloat div(ExactFloat, ExactFloat);

//! \related ValidatedFloat \brief Positive integer power function. Yields an over-approximation to \f$\{ x^m \mid x\in I\}\f$.
ValidatedFloat pow(ValidatedFloat i, uint m);
//! \related ValidatedFloat \brief %Integer power function. Yields an over-approximation to \f$\{ x^n \mid x\in I\}\f$.
ValidatedFloat pow(ValidatedFloat i, int n);

//! \related ValidatedFloat \brief Square-root function. Yields an over-approximation to \f$\{ \sqrt{x} \mid x\in I\}\f$.
//! Requires \c I.lower()>=0 .
ValidatedFloat sqrt(ValidatedFloat);
//! \related ValidatedFloat \brief Exponential function. Yields an over-approximation to \f$\{ \exp{x} \mid x\in I\}\f$.
ValidatedFloat exp(ValidatedFloat);
//! \related ValidatedFloat \brief Natural logarithm function. Yields an over-approximation to \f$\{ \log{x} \mid x\in I\}\f$.
//! Requires \c I.lower()>0 .
ValidatedFloat log(ValidatedFloat);

//! \related ValidatedFloat \brief Sine function. Yields an over-approximation to \f$\{ \sin{x} \mid x\in I\}\f$.
ValidatedFloat sin(ValidatedFloat);
//! \related ValidatedFloat \brief Cosine function. Yields an over-approximation to \f$\{ \sin{x} \mid x\in I\}\f$.
ValidatedFloat cos(ValidatedFloat);
//! \related ValidatedFloat \brief Tangent function. Yields an over-approximation to \f$\{ \sin{x} \mid x\in I\}\f$.
ValidatedFloat tan(ValidatedFloat);
ValidatedFloat asin(ValidatedFloat);
ValidatedFloat acos(ValidatedFloat);
ValidatedFloat atan(ValidatedFloat);


//! \related ValidatedFloat \brief The magnitude of the interval \a I. Yields \f$ \max\{ |x|\,\mid\,x\in I \}\f$.
inline Float mag(ValidatedFloat i) { return max(abs(i.lower_value()),abs(i.upper_value())); }
//! \related ValidatedFloat \brief The mignitude of the interval \a I. Yields \f$ \min\{ |x|\,\mid\,x\in I \}\f$.
inline Float mig(ValidatedFloat i) { return min(abs(i.lower_value()),abs(i.upper_value())); }

inline ValidatedFloat max(ValidatedFloat i1, ValidatedFloat i2)
{
    return ValidatedFloat(max(i1.lower_value(),i2.lower_value()),max(i1.upper_value(),i2.upper_value()));
}

inline ValidatedFloat min(ValidatedFloat i1, ValidatedFloat i2)
{
    return ValidatedFloat(min(i1.lower_value(),i2.lower_value()),min(i1.upper_value(),i2.upper_value()));
}


inline ValidatedFloat abs(ValidatedFloat i)
{
    if(i.lower_value()>=0) {
        return ValidatedFloat(i.lower_value(),i.upper_value());
    } else if(i.upper_value()<=0) {
        return ValidatedFloat(neg(i.upper_value()),neg(i.lower_value()));
    } else {
        return ValidatedFloat(static_cast<Float>(0.0),max(neg(i.lower_value()),i.upper_value()));
    }
}

inline ValidatedFloat pos(ValidatedFloat i)
{
    return ValidatedFloat(+i.lower_value(),+i.upper_value());
}

inline ExactFloat pos(ExactFloat x);

inline ValidatedFloat neg(ValidatedFloat i)
{
    return ValidatedFloat(-i.upper_value(),-i.lower_value());
}

inline ExactFloat neg(ExactFloat x);

inline ValidatedFloat half(ValidatedFloat i) {
    return ValidatedFloat(half(i.lower_value()),half(i.upper_value()));
}

inline ValidatedFloat sqr(ExactFloat x)
{
    rounding_mode_t rnd=get_rounding_mode();
    volatile double xv=internal_cast<volatile double&>(x.value());
    set_rounding_mode(downward);
    volatile double rl=xv*xv;
    set_rounding_mode(upward);
    volatile double ru=xv*xv;
    set_rounding_mode(rnd);
    return ValidatedFloat(rl,ru);
}

inline ValidatedFloat rec(ExactFloat x)
{
    rounding_mode_t rnd=get_rounding_mode();
    volatile double xv=internal_cast<volatile double&>(x.value());
    set_rounding_mode(downward);
    volatile double rl=1.0/xv;
    set_rounding_mode(upward);
    volatile double ru=1.0/xv;
    set_rounding_mode(rnd);
    return ValidatedFloat(rl,ru);
}



inline ValidatedFloat add(ValidatedFloat i1, ValidatedFloat i2)
{
    rounding_mode_t rnd=get_rounding_mode();
    volatile double i1l=internal_cast<volatile double&>(i1.lower_value());
    volatile double i1u=internal_cast<volatile double&>(i1.upper_value());
    volatile double i2l=internal_cast<volatile double&>(i2.lower_value());
    volatile double i2u=internal_cast<volatile double&>(i2.upper_value());
    set_rounding_mode(downward);
    volatile double rl=i1l+i2l;
    set_rounding_mode(upward);
    volatile double ru=i1u+i2u;
    set_rounding_mode(rnd);
    return ValidatedFloat(rl,ru);
}

inline ValidatedFloat add(ValidatedFloat i1, ExactFloat x2)
{
    rounding_mode_t rnd=get_rounding_mode();
    volatile double i1l=internal_cast<volatile double&>(i1.lower_value());
    volatile double i1u=internal_cast<volatile double&>(i1.upper_value());
    volatile double x2v=internal_cast<volatile double&>(x2.value());
    set_rounding_mode(downward);
    volatile double rl=i1l+x2v;
    set_rounding_mode(upward);
    volatile double ru=i1u+x2v;
    set_rounding_mode(rnd);
    return ValidatedFloat(rl,ru);
}

inline ValidatedFloat add(ExactFloat x1, ValidatedFloat i2)
{
    return add(i2,x1);
}

inline ValidatedFloat add(ExactFloat x1, ExactFloat x2)
{
    rounding_mode_t rnd=get_rounding_mode();
    volatile double x1v=internal_cast<volatile double&>(x1.value());
    volatile double x2v=internal_cast<volatile double&>(x2.value());
    set_rounding_mode(downward);
    volatile double rl=x1v+x2v;
    set_rounding_mode(upward);
    volatile double ru=x1v+x2v;
    set_rounding_mode(rnd);
    return ValidatedFloat(rl,ru);
}

inline ValidatedFloat sub(ValidatedFloat i1, ValidatedFloat i2)
{
    rounding_mode_t rnd=get_rounding_mode();
    volatile double i1l=internal_cast<volatile double&>(i1.lower_value());
    volatile double i1u=internal_cast<volatile double&>(i1.upper_value());
    volatile double i2l=internal_cast<volatile double&>(i2.lower_value());
    volatile double i2u=internal_cast<volatile double&>(i2.upper_value());
    set_rounding_mode(downward);
    volatile double rl=i1l-i2u;
    set_rounding_mode(upward);
    volatile double ru=i1u-i2l;
    set_rounding_mode(rnd);
    return ValidatedFloat(rl,ru);
}

inline ValidatedFloat sub(ValidatedFloat i1, ExactFloat x2)
{
    rounding_mode_t rnd=get_rounding_mode();
    volatile double i1l=internal_cast<volatile double&>(i1.lower_value());
    volatile double i1u=internal_cast<volatile double&>(i1.upper_value());
    volatile double x2v=internal_cast<volatile double&>(x2.value());
    set_rounding_mode(downward);
    volatile double rl=i1l-x2v;
    set_rounding_mode(upward);
    volatile double ru=i1u-x2v;
    set_rounding_mode(rnd);
    return ValidatedFloat(rl,ru);
}

inline ValidatedFloat sub(ExactFloat x1, ValidatedFloat i2)
{
    rounding_mode_t rnd=get_rounding_mode();
    volatile double x1v=internal_cast<volatile double&>(x1.value());
    volatile double i2l=internal_cast<volatile double&>(i2.lower_value());
    volatile double i2u=internal_cast<volatile double&>(i2.upper_value());
    set_rounding_mode(downward);
    volatile double rl=x1v-i2u;
    set_rounding_mode(upward);
    volatile double ru=x1v-i2l;
    set_rounding_mode(rnd);
    return ValidatedFloat(rl,ru);
}

inline ValidatedFloat sub(ExactFloat x1, ExactFloat x2)
{
    rounding_mode_t rnd=get_rounding_mode();
    volatile double x1v=internal_cast<volatile double&>(x1.value());
    volatile double x2v=internal_cast<volatile double&>(x2.value());
    set_rounding_mode(downward);
    volatile double rl=x1v-x2v;
    set_rounding_mode(upward);
    volatile double ru=x1v-x2v;
    set_rounding_mode(rnd);
    return ValidatedFloat(rl,ru);
}

inline ValidatedFloat mul(ExactFloat x1, ExactFloat x2)
{
    rounding_mode_t rnd=get_rounding_mode();
    volatile double x1v=internal_cast<volatile double&>(x1.value());
    volatile double x2v=internal_cast<volatile double&>(x2.value());
    set_rounding_mode(downward);
    volatile double rl=x1v*x2v;
    set_rounding_mode(upward);
    volatile double ru=x1v*x2v;
    set_rounding_mode(rnd);
    return ValidatedFloat(rl,ru);
}

inline ValidatedFloat div(ExactFloat x1, ExactFloat x2)
{
    rounding_mode_t rnd=get_rounding_mode();
    volatile double x1v=internal_cast<volatile double&>(x1.value());
    volatile double x2v=internal_cast<volatile double&>(x2.value());
    set_rounding_mode(downward);
    volatile double rl=x1v/x2v;
    set_rounding_mode(upward);
    volatile double ru=x1v/x2v;
    set_rounding_mode(rnd);
    return ValidatedFloat(rl,ru);
}

inline ValidatedFloat pow(ExactFloat x1, int n2)
{
    return pow(ValidatedFloat(x1),n2);
}

inline ValidatedFloat med(ExactFloat x1, ExactFloat x2)
{
    return add(half(x1),half(x2));
}

inline ValidatedFloat rad(ExactFloat x1, ExactFloat x2)
{
    return sub(half(x2),half(x1));
}

inline ValidatedFloat med(ValidatedFloat x);

inline ValidatedFloat rad(ValidatedFloat i);

//! \related ValidatedFloat \brief Unary plus operator. Should be implemented exactly and yield \f$\{ +x \mid x\in I\}\f$.
inline ValidatedFloat operator+(const ValidatedFloat& i) { return pos(i); }
//! \related ValidatedFloat \brief Unary negation operator. Should be implemented exactly and yield \f$\{ -x \mid x\in I\}\f$.
inline ValidatedFloat operator-(const ValidatedFloat& i) { return neg(i); }
//! \related ValidatedFloat \brief Binary addition operator. Guaranteed to yield an over approximation to \f$\{ x_1+x_2 \mid x_1\in I_1 \wedge x_2\in I_2\}\f$.
inline ValidatedFloat operator+(const ValidatedFloat& i1, const ValidatedFloat& i2) { return add(i1,i2); }
//! \related ValidatedFloat \brief Binary addition operator. Guaranteed to yield an over approximation to \f$\{ x_1-x_2 \mid x_1\in I_1 \wedge x_2\in I_2\}\f$.
inline ValidatedFloat operator-(const ValidatedFloat& i1, const ValidatedFloat& i2) { return sub(i1,i2); }
//! \related ValidatedFloat \brief Binary addition operator. Guaranteed to yield an over approximation to \f$\{ x_1*x_2 \mid x_1\in I_1 \wedge x_2\in I_2\}\f$.
inline ValidatedFloat operator*(const ValidatedFloat& i1, const ValidatedFloat& i2) { return mul(i1,i2); }
//! \related ValidatedFloat \brief Binary addition operator. Guaranteed to yield an over approximation to \f$\{ x_1/x_2 \mid x_1\in I_1 \wedge x_2\in I_2\}\f$. Yields \f$[-\infty,+\infty]\f$ if \f$0\in I_2\f$.
inline ValidatedFloat operator/(const ValidatedFloat& i1, const ValidatedFloat& i2) { return div(i1,i2); };

//! \related ValidatedFloat \brief Inplace addition operator.
inline ValidatedFloat& operator+=(ValidatedFloat& i1, const ValidatedFloat& i2) { i1=add(i1,i2); return i1; }
//! \related ValidatedFloat \brief Inplace subtraction operator.
inline ValidatedFloat& operator-=(ValidatedFloat& i1, const ValidatedFloat& i2) { i1=sub(i1,i2); return i1; }
//! \related ValidatedFloat \brief Inplace multiplication operator.
inline ValidatedFloat& operator*=(ValidatedFloat& i1, const ValidatedFloat& i2) { i1=mul(i1,i2); return i1; }
//! \related ValidatedFloat \brief Inplace division operator.
inline ValidatedFloat& operator/=(ValidatedFloat& i1, const ValidatedFloat& i2) { i1=div(i1,i2); return i1; }

inline ValidatedFloat operator+(const ValidatedFloat& i1, const ExactFloat& x2) { return add(i1,x2); }
inline ValidatedFloat operator-(const ValidatedFloat& i1, const ExactFloat& x2) { return sub(i1,x2); }
inline ValidatedFloat operator*(const ValidatedFloat& i1, const ExactFloat& x2) { return mul(i1,x2); }
inline ValidatedFloat operator/(const ValidatedFloat& i1, const ExactFloat& x2) { return div(i1,x2); }
inline ValidatedFloat operator+(const ExactFloat& x1, const ValidatedFloat& i2) { return add(i2,x1); }
inline ValidatedFloat operator-(const ExactFloat& x1, const ValidatedFloat& i2) { return sub(x1,i2); }
inline ValidatedFloat operator*(const ExactFloat& x1, const ValidatedFloat& i2) { return mul(i2,x1); }
inline ValidatedFloat operator/(const ExactFloat& x1, const ValidatedFloat& i2) { return div(x1,i2); }

inline ValidatedFloat operator+(const ExactFloat& x1, const ExactFloat& x2) { return add(x1,x2); }
inline ValidatedFloat operator-(const ExactFloat& x1, const ExactFloat& x2) { return sub(x1,x2); }
inline ValidatedFloat operator*(const ExactFloat& x1, const ExactFloat& x2) { return mul(x1,x2); }
inline ValidatedFloat operator/(const ExactFloat& x1, const ExactFloat& x2) { return div(x1,x2); }
inline ValidatedFloat pow(const ExactFloat& x, int n) { return pow(ValidatedFloat(x),n); }

inline ValidatedFloat& operator+=(ValidatedFloat& i1, const ExactFloat& x2) { i1=add(i1,x2); return i1; }
inline ValidatedFloat& operator-=(ValidatedFloat& i1, const ExactFloat& x2) { i1=sub(i1,x2); return i1; }
inline ValidatedFloat& operator*=(ValidatedFloat& i1, const ExactFloat& x2) { i1=mul(i1,x2); return i1; }
inline ValidatedFloat& operator/=(ValidatedFloat& i1, const ExactFloat& x2) { i1=div(i1,x2); return i1; }

//inline ValidatedFloat operator/(const ValidatedFloat& i1, int n2) { return div(i1,Float(n2)); }
//inline ValidatedFloat operator/(const ValidatedFloat& i1, double x2) { return div(i1,Float(x2)); }

// Standard equality operators
//! \related ValidatedFloat \brief Equality operator. Tests equality of intervals as geometric objects, so \c [0,1]==[0,1] returns \c true.
inline bool operator==(const ValidatedFloat& i1, const ValidatedFloat& i2) { return i1.lower_value()==i2.lower_value() && i1.upper_value()==i2.upper_value(); }
//! \related ValidatedFloat \brief Inequality operator. Tests equality of intervals as geometric objects, so \c [0,2]!=[1,3] returns \c true,
//! even though the intervals possibly represent the same exact real value.
inline bool operator!=(const ValidatedFloat& i1, const ValidatedFloat& i2) { return i1.lower_value()!=i2.lower_value() || i1.upper_value()!=i2.upper_value(); }

// Boost-style tribool (in)equality operators
//inline tribool operator==(const ValidatedFloat& i1, const ValidatedFloat& i2) {
//  if(i1.lower_value()>i2.upper_value() || i1.upper_value()<i2.lower_value()) { return false; } else if(i1.lower_value()==i2.upper_value() && i1.upper_value()==i2.lower_value()) { return true; } else { return indeterminate; } }
//inline tribool operator!=(const ValidatedFloat& i1, const ValidatedFloat& i2) { return !(i1==i2); }

//! \related ValidatedFloat \brief Equality operator. Tests equality of represented real-point value.
//! Hence \c [0.0,2.0]==1.0 yields \c indeterminate since the interval may represent a real number other than \c 1.0 .
inline tribool operator==(const ValidatedFloat& i1, const ExactFloat& x2) {
    if(i1.upper_value()<x2.value() || i1.lower_value()>x2.value() ) { return false; }
    else if(i1.lower_value()==x2.value()  && i1.upper_value()==x2.value() ) { return true; }
    else { return indeterminate; }
}

//! \related ValidatedFloat \brief Equality operator. Tests equality of represented real-point value.
//! Hence \c [0.0,2.0]!=1.0 yields \c indeterminate since the interval may represent a real number equal to \c 1.0 .
inline tribool operator!=(const ValidatedFloat& i1, const ExactFloat& x2) {
    if(i1.upper_value()<x2.value()  || i1.lower_value()>x2.value() ) { return true; }
    else if(i1.lower_value()==x2.value()  && i1.upper_value()==x2.value() ) { return false; }
    else { return indeterminate; }
}

inline tribool operator> (const ValidatedFloat& i1, const ExactFloat& x2) {
    if(i1.lower_value()> x2.value() ) { return true; }
    else if(i1.upper_value()<=x2.value() ) { return false; }
    else { return indeterminate; }
}

inline tribool operator< (const ValidatedFloat& i1, const ExactFloat& x2) {
    if(i1.upper_value()< x2.value() ) { return true; }
    else if(i1.lower_value()>=x2.value() ) { return false; }
    else { return indeterminate; }
}

inline tribool operator>=(const ValidatedFloat& i1, const ExactFloat& x2) {
    if(i1.lower_value()>=x2.value() ) { return true; }
    else if(i1.upper_value()< x2.value() ) { return false; }
    else { return indeterminate; }
}

inline tribool operator<=(const ValidatedFloat& i1, const ExactFloat& x2) {
    if(i1.upper_value()<=x2.value() ) { return true; }
    else if(i1.lower_value()> x2.value() ) { return false; }
    else { return indeterminate; }
}

/*
inline tribool operator==(const ValidatedFloat& i1, double x2) { return i1==static_cast<ExactFloat>(x2); }
inline tribool operator!=(const ValidatedFloat& i1, double x2) { return i1!=static_cast<ExactFloat>(x2); }
inline tribool operator<=(const ValidatedFloat& i1, double x2) { return i1<=static_cast<ExactFloat>(x2); }
inline tribool operator>=(const ValidatedFloat& i1, double x2) { return i1>=static_cast<ExactFloat>(x2); }
inline tribool operator< (const ValidatedFloat& i1, double x2) { return i1< static_cast<ExactFloat>(x2); }
inline tribool operator> (const ValidatedFloat& i1, double x2) { return i1> static_cast<ExactFloat>(x2); }
*/


//! \related ValidatedFloat \brief Strict greater-than comparison operator. Tests equality of represented real-point value.
//! Hence \c [1.0,3.0]>[0.0,2.0] yields \c indeterminate since the first interval could represent the number 1.25 and the second 1.75.
inline tribool operator> (ValidatedFloat i1, ValidatedFloat i2) {
    if(i1.lower_value()> i2.upper_value()) { return true; }
    else if(i1.upper_value()<=i2.lower_value()) { return false; }
    else { return indeterminate; }
}

//! \related ValidatedFloat \brief Strict greater-than comparison operator. Tests equality of represented real-point value.
inline tribool operator< (ValidatedFloat i1, ValidatedFloat i2) {
    if(i1.upper_value()< i2.lower_value()) { return true; }
    else if(i1.lower_value()>=i2.upper_value()) { return false; }
    else { return indeterminate; }
}

//! \related ValidatedFloat \brief Strict greater-than comparison operator. Tests equality of represented real-point value.
inline tribool operator>=(ValidatedFloat i1, ValidatedFloat i2) {
    if(i1.lower_value()>=i2.upper_value()) { return true; }
    else if(i1.upper_value()< i2.lower_value()) { return false; }
    else { return indeterminate; }
}

//! \related ValidatedFloat \brief Strict greater-than comparison operator. Tests equality of represented real-point value.
inline tribool operator<=(ValidatedFloat i1, ValidatedFloat i2) {
    if(i1.upper_value()<=i2.lower_value()) { return true; }
    else if(i1.lower_value()> i2.upper_value()) { return false; }
    else { return indeterminate; }
}

#ifdef ARIADNE_ENABLE_SERIALIZATION
  template<class A> void serialize(A& a, ValidatedFloat& ivl, const uint version) {
    a & ivl.lower_value() & ivl.upper_value(); }
#endif

std::ostream& operator<<(std::ostream&, const ValidatedFloat&);
std::istream& operator>>(std::istream&, ValidatedFloat&);


} // namespace Ariadne

#endif