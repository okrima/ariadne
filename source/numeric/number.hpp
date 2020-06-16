/***************************************************************************
 *            numeric/number.hpp
 *
 *  Copyright  2013-20  Pieter Collins
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

/*! \file numeric/number.hpp
 *  \brief Generic numbers
 */



#ifndef ARIADNE_NUMBER_HPP
#define ARIADNE_NUMBER_HPP

#include "../utility/handle.hpp"
#include "../numeric/paradigm.hpp"
#include "../utility/prototype.hpp"

#include "logical.decl.hpp"
#include "number.decl.hpp"
#include "float.decl.hpp"

#include "number_interface.hpp"

#include "arithmetic.hpp"
#include "integer.hpp"
#include "rational.hpp"
#include "real.hpp"

#include "number_interface.hpp"

namespace Ariadne {

/************ Number *********************************************************/

template<class X> struct IsNumber;

class NumberInterface;

template<class P> class Number;
template<class P> struct IsNumber<Number<P>> : True { };

struct DispatchException : public std::runtime_error {
    using std::runtime_error::runtime_error;
};

class DeclareNumberOperators {
    friend ApproximateNumber operator+(ApproximateNumber const&, ApproximateNumber const&);
    friend ApproximateNumber operator-(ApproximateNumber const&, ApproximateNumber const&);
    friend ApproximateNumber operator*(ApproximateNumber const&, ApproximateNumber const&);
    friend ApproximateNumber operator/(ApproximateNumber const&, ApproximateNumber const&);

    friend ValidatedLowerNumber operator+(ValidatedLowerNumber const& y1, ValidatedLowerNumber const& y2);
    friend ValidatedLowerNumber operator-(ValidatedLowerNumber const& y1, ValidatedUpperNumber const& y2);
    friend ValidatedUpperNumber operator+(ValidatedUpperNumber const& y1, ValidatedUpperNumber const& y2);
    friend ValidatedUpperNumber operator-(ValidatedUpperNumber const& y1, ValidatedLowerNumber const& y2);

    template<GenericNumber N, BuiltinFloatingPoint D> friend decltype(auto)
    operator+(N const& y1, D const& d2) { return y1+Number<ApproximateTag>(d2); }
    template<GenericNumber N, BuiltinFloatingPoint D> friend decltype(auto)
    operator-(N const& y1, D const& d2) { return y1-Number<ApproximateTag>(d2); }
    template<GenericNumber N, BuiltinFloatingPoint D> friend decltype(auto)
    operator*(N const& y1, D const& d2) { return y1*Number<ApproximateTag>(d2); }
    template<GenericNumber N, BuiltinFloatingPoint D> friend decltype(auto)
    operator/(N const& y1, D const& d2) { return y1/Number<ApproximateTag>(d2); }
    template<GenericNumber N, BuiltinFloatingPoint D> friend decltype(auto)
    operator+(D const& d1, N const& y2) { return Number<ApproximateTag>(d1)+y2; }
    template<GenericNumber N, BuiltinFloatingPoint D> friend decltype(auto)
    operator-(D const& d1, N const& y2) { return Number<ApproximateTag>(d1)-y2; }
    template<GenericNumber N, BuiltinFloatingPoint D> friend decltype(auto)
    operator*(D const& d1, N const& y2) { return Number<ApproximateTag>(d1)*y2; }
    template<GenericNumber N, BuiltinFloatingPoint D> friend decltype(auto)
    operator/(D const& d1, N const& y2) { return Number<ApproximateTag>(d1)/y2; }

    template<ConcreteNumber R, class P> friend decltype(auto)
    operator+(R const& r1, Number<P> const& y2) { return Number<Paradigm<R>>(r1)+y2; }
    template<ConcreteNumber R, class P> friend decltype(auto)
    operator+(Number<P> const& y1, R const& r2) { return y1+Number<Paradigm<R>>(r2); }
};


template<class X, class P> concept ConvertibleBuiltinFloatingPointToNumber
    = Same<P,ApproximateTag> and BuiltinFloatingPoint<X>;
template<class X, class P> concept ConvertibleViaRealToNumber
    = WeakerThan<P,ParadigmTag<X>> and Convertible<X,Real>;
template<class X, class P> concept ConvertibleViaNumberToNumber
    = WeakerThan<P,ParadigmTag<X>> and (not Convertible<X,Real>) and Convertible<X,Number<ParadigmTag<X>>>;

//! \ingroup NumericModule
//! \brief Generic numbers with computational paradigm \a P, which may be %EffectiveTag, %ValidatedTag, %UpperTag, %LowerTag or %ApproximateTag.
template<class P> class Number
    : public Handle<NumberInterface>
    , public DeclareNumberOperators
{
    static_assert(IsParadigm<P>,"P must be a paradigm");
    template<class PP> friend class Number;
    template<class X> static const bool IsGettableAs = ANumber<X> and WeakerThan<typename X::Paradigm,P> and (not Same<typename X::Paradigm,ExactTag>);
  private:
    typedef Opposite<P> NP;
    typedef Weaker<P,NP> SP;
    typedef Widen<P> WP;
    //friend class DispatchGenericField<Number<P>>;
  public:
    typedef NumberInterface Interface;
    typedef P Paradigm;
    typedef Number<P> NumericType;
  public:
    using Handle<Interface>::Handle;
    explicit Number(NumberInterface* p) : Handle<NumberInterface>(p) { }
  private:
    explicit Number(Handle<NumberInterface> h) : Handle<NumberInterface>(h) { }
    Handle<NumberInterface> handle() const { return *this; }
  public:

    Number() : Number(Integer(0)) { }

    // Construct from a Number of a weaker paradigm
    template<class SP> requires WeakerThan<P,SP> Number(const Number<SP>& y) : Number<P>(y.handle()) { }

    // Construct from a builtin integer
    template<BuiltinIntegral N> Number(const N& n) : Number<P>(Integer(n)) { }
    // Construct from a builtin floating-point number
    template<ConvertibleBuiltinFloatingPointToNumber<P> X> Number(const X& x) : Number<P>(Dyadic(x)) { }

    // Construct from a type which is convertible to Real.
    template<ConvertibleViaRealToNumber<P> X> Number<P>(X const & x) : Number<P>(x.operator Number<ParadigmTag<X>>()) { }

    // Construct from a type which is convertible to another Number type.
    // TODO: Decide conversion properties from concrete type to Number<P>
    template<ConvertibleViaNumberToNumber<P> X>
        explicit Number<P>(X const & x) : Number<P>(x.operator Number<ParadigmTag<X>>()) { }

    //! \brief Get the value of the number as a double-precision floating-point type
    template<class WP> requires WeakerThan<WP,P>
    FloatType<WP,DoublePrecision> get(WP par, DoublePrecision const& prec) const { return pointer()->_get(WP(),prec); }
    //! \brief Get the value of the number as a multiple-precision floating-point type
    template<class WP> requires WeakerThan<WP,P>
    FloatType<WP,MultiplePrecision> get(WP par, MultiplePrecision const& prec) const { return pointer()->_get(WP(),prec); }
    //! \brief Get the value of the number as a double-precision floating-point type
    template<class PR> requires SameAs<PR,DoublePrecision>
    FloatType<P,PR> get(PR pr) const { return pointer()->_get(P(),pr); }
    //! \brief Get the value of the number as a multiple-precision floating-point type
    template<class PR> requires SameAs<PR,MultiplePrecision>
    FloatType<P,PR> get(PR pr) const { return pointer()->_get(P(),pr); }

    //! \brief Get the value of the number as a multiple-precision floating-point ball.
    template<class WP> requires WeakerThan<WP,P> and SameAs<WP,MetricTag>
    FloatType<WP,DoublePrecision> get(WP par, DoublePrecision const& prec, DoublePrecision const& errprec) const { return pointer()->_get(WP(),prec,errprec); }
    //! \brief Get the value of the number as a multiple-precision floating-point ball, with double-precision error.
    template<class WP> requires WeakerThan<WP,P> and SameAs<WP,MetricTag>
    FloatType<WP,MultiplePrecision,DoublePrecision> get(WP par, MultiplePrecision const& prec, DoublePrecision const& errprec) const { return pointer()->_get(WP(),prec,errprec); }
    //! \brief Get the value of the number as a multiple-precision floating-point ball.
    template<class WP> requires WeakerThan<WP,P> and SameAs<WP,MetricTag>
    FloatType<WP,MultiplePrecision,MultiplePrecision> get(WP par, MultiplePrecision const& prec, MultiplePrecision const& errprec) const { return pointer()->_get(WP(),prec,errprec); }

    //! \brief Get the value of the number as a double-precision floating-point type
    FloatType<P,DoublePrecision> get() const { return pointer()->_get(WP()); }
    //! \brief Get the value of the number as a double-precision floating-point type
    FloatType<P,DoublePrecision> get(DoublePrecision const& prec) const { return pointer()->_get(WP(),prec); }
    //! \brief Get the value of the number as a multiple-precision floating-point type
    FloatType<P,MultiplePrecision> get(MultiplePrecision const& prec) const { return pointer()->_get(WP(),prec); }

    template<class X> X extract() const;

    friend Number<P> operator+(Number<P> const& y) { return pos(y); }
    friend Number<NP> operator-(Number<P> const& y) { return neg(y); }
    friend Number<P> operator+(Number<P> const& y1, Number<P> const& y2) { return add(y1,y2); }
    friend Number<P> operator-(Number<P> const& y1, Number<NP> const& y2) { return sub(y1,y2); }
    friend Number<P> operator*(Number<P> const& y1, Number<P> const& y2) { return mul(y1,y2); }
    friend Number<P> operator/(Number<P> const& y1, Number<NP> const& y2) { return div(y1,y2); }
    friend Number<P>& operator+=(Number<P>& y1, Number<P> const& y2) { return y1=y1+y2; }
    friend Number<P>& operator-=(Number<P>& y1, Number<NP> const& y2) { return y1=y1-y2; }
    friend Number<P>& operator*=(Number<P>& y1, Number<P> const& y2) { return y1=y1*y2; }
    friend Number<P>& operator/=(Number<P>& y1, Number<NP> const& y2) { return y1=y1/y2; }

    friend Number<P> nul(Number<P> const& y) { return _apply<P>(Nul(),y); }
    friend Number<P> pos(Number<P> const& y) { return _apply<P>(Pos(),y); }
    friend Number<NP> neg(Number<P> const& y) { return _apply<NP>(Neg(),y); }
    friend Number<P> sqr(Number<P> const& y) { return _apply<P>(Sqr(),y); }
    friend Number<P> hlf(Number<P> const& y) { return _apply<P>(Hlf(),y); }
    friend Number<NP> rec(Number<P> const& y) { return _apply<NP>(Rec(),y); }
    friend Number<P> add(Number<P> const& y1, Number<P> const& y2) { return _apply<P>(Add(),y1,y2); }
    friend Number<P> sub(Number<P> const& y1, Number<NP> const& y2) { return _apply<P>(Sub(),y1,y2); }
    friend Number<P> mul(Number<P> const& y1, Number<P> const& y2) { return _apply<P>(Mul(),y1,y2); }
    friend Number<P> div(Number<P> const& y1, Number<NP> const& y2) { return _apply<P>(Div(),y1,y2); }
    friend Number<P> sqrt(Number<P> const& y) { return _apply<P>(Sqrt(),y); }
    friend Number<P> exp(Number<P> const& y) { return _apply<P>(Exp(),y); }
    friend Number<P> log(Number<P> const& y) { return _apply<P>(Log(),y); }
    friend Number<SP> sin(Number<P> const& y) { return _apply<SP>(Sin(),y); }
    friend Number<SP> cos(Number<P> const& y) { return _apply<SP>(Cos(),y); }
    friend Number<P> tan(Number<P> const& y) { return _apply<P>(Tan(),y); }
    friend Number<SP> asin(Number<P> const& y) { return _apply<SP>(Asin(),y); }
    friend Number<SP> acos(Number<P> const& y) { return _apply<SP>(Acos(),y); }
    friend Number<P> atan(Number<P> const& y) { return _apply<P>(Atan(),y); }

    friend Number<P> pow(Number<P> const& y, Nat m) { return _apply<P>(Pow(),y,m); }
    friend Number<SP> pow(Number<P> const& y, Int n) { return _apply<P>(Pow(),y,n); }

    friend Positive<Number<SP>> abs(Number<P> const& y) { return Positive<Number<SP>>(_apply<P>(Abs(),y)); }
    friend Number<P> max(Number<P> const& y1, Number<P> const& y2) { return _apply<P>(Max(),y1,y2); }
    friend Number<P> min(Number<P> const& y1, Number<P> const& y2) { return _apply<P>(Min(),y1,y2); }


    friend LogicalType<Equality<P>> operator==(Number<P> const& y1, Number<NP> const& y2) {
        return LogicalType<Equality<P>>(y1.ref()._equals(y2.ref())); }
    friend LogicalType<LessThan<P>> operator< (Number<P> const& y1, Number<Negated<P>> const& y2) {
        return LogicalType<LessThan<P>>(y1.ref()._less(y2.ref())); }
    friend LogicalType<LessThan<Negated<P>>> operator> (Number<P> const& y1, Number<Negated<P>> const& y2) { return (y2<y1); }
    friend LogicalType<Negated<Equality<P>>> operator!=(Number<P> const& y1, Number<Negated<P>> const& y2) { return !(y1==y2); }
    friend LogicalType<LessThan<P>> operator<=(Number<P> const& y1, Number<Negated<P>> const& y2) { return !(y1>y2); }
    friend LogicalType<LessThan<Negated<P>>> operator>=(Number<P> const& y1, Number<Negated<P>> const& y2) { return !(y1<y2); }

    String class_name() const { return this->ref()._class_name(); }

    friend OutputStream& operator<<(OutputStream& os, Number<P> const& y) { return y.ref()._write(os); }
  private:
    template<class RP, class OP> static Number<RP> _apply(OP op, Number<P> const& y) { return Number<RP>(y.ref()._apply(op)); }
    template<class RP, class OP> static Number<RP> _apply(OP op, Number<P> const& y, Int n) { return Number<RP>(y.ref()._apply(op,n)); }
    template<class RP, class OP, class P1, class P2> static Number<RP> _apply(OP op, Number<P1> const& y1, Number<P2> const& y2) {
        return Number<RP>(y1.ref()._apply(op,&y2.ref())); }
};



template<class P> Positive<Number<P>> cast_positive(Number<P> y);

class PositiveNumberOperations {
    friend Positive<ValidatedUpperNumber> mag(Positive<ValidatedUpperNumber> const& y);
};

template<class P> class Positive<Number<P>> : public Number<P>, public PositiveNumberOperations {
    friend Number<P> const& unsign(Positive<Number<P>> const& y) { return y; }
  public:
    Positive<Number<P>>() : Number<P>() { }
    explicit Positive<Number<P>>(Number<P> const& y) : Number<P>(y) { }
    template<class N> requires Constructible<ExactNumber,N> Positive<Number<P>>(const Positive<N>& n) : Number<P>(ExactNumber(static_cast<N const&>(n))) { }
    template<class N> requires Constructible<ExactNumber,N> Positive<Number<P>>(const N& n) : Number<P>(ExactNumber(n)) { }
    explicit operator Number<P> () const { return *this; }

    friend Positive<Number<P>> operator+(Positive<Number<P>> const& y1, Positive<Number<P>> const& y2) { return cast_positive(add(unsign(y1),unsign(y2))); }
    friend Positive<Number<P>> operator*(Positive<Number<P>> const& y1, Positive<Number<P>> const& y2) { return cast_positive(mul(unsign(y1),unsign(y2))); }
    friend Positive<Number<P>> add(Positive<Number<P>> const& y1, Positive<Number<P>> const& y2) { return cast_positive(add(unsign(y1),unsign(y2))); }
    friend Positive<Number<P>> mul(Positive<Number<P>> const& y1, Positive<Number<P>> const& y2) { return cast_positive(mul(unsign(y1),unsign(y2))); }
    friend Positive<Number<P>> max(Positive<Number<P>> const& y1, Positive<Number<P>> const& y2) { return cast_positive(max(unsign(y1),unsign(y2))); }
    friend Positive<Number<P>> abs(Positive<Number<P>> const& y) { return y; }

};
template<class P> Positive<Number<P>> cast_positive(Number<P> y) { return Positive<Number<P>>(y); }
template<class P> Positive<ExactNumber> cast_exact(Positive<Number<P>> const& y) { return cast_positive(cast_exact(unsign(y))); }

template<class X> decltype(auto) cast_generic(X const& x) { return x.generic(); }

} // namespace Ariadne

#endif
