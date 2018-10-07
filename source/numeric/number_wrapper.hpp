/***************************************************************************
 *            numeric/number_wrapper.hpp
 *
 *  Copyright 2013-17  Pieter Collins
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

/*! \file numeric/number_wrapper.hpp
 *  \brief
 */



#ifndef ARIADNE_NUMBER_WRAPPER_HPP
#define ARIADNE_NUMBER_WRAPPER_HPP

#include "../utility/module.hpp"
#include "../utility/dispatching.hpp"
#include "../numeric/paradigm.hpp"

#include "number_interface.hpp"

#include "number.hpp"
#include "logical.hpp"
#include "floatdp.hpp"
#include "floatmp.hpp"
#include "float-user.hpp"

#include "../symbolic/templates.hpp"
#include "../numeric/operators.hpp"

namespace Ariadne {


template<class X> using NumberMixin = Mixin<X,NumberInterface>;
template<class X> using NumberWrapper = Wrapper<X,NumberInterface>;

template<class X> using UpperNumberMixin = Mixin<X,UpperNumberInterface>;
template<class X> using UpperNumberWrapper = Wrapper<X,UpperNumberInterface>;

inline OutputStream& operator<<(OutputStream& os, Comparison c) {
    return os << ( (c==Comparison::EQUAL) ? "EQUAL" : (c==Comparison::LESS) ? "LESS" : "GREATER" );
}
inline OutputStream& operator<<(OutputStream& os, Sign s) {
    return os << ( (s==Sign::ZERO) ? "ZERO" : (s==Sign::NEGATIVE) ? "NEGATIVE" : "POSITIVE" );
}

// FIXME: Should test for other potential infinities
inline Comparison cmp(NumberInterface const& y1, NumberInterface const& y2) {
    Comparison res;
    FloatDPValue const* x1=extract<FloatDPValue>(&y1);
    FloatDPValue const* x2=extract<FloatDPValue>(&y2);
    if(x1) {
        if(x2) { res= cmp(ExactDouble(x1->raw().get_d()),ExactDouble(x2->raw().get_d())); }
        else { res= cmp(ExactDouble(x1->raw().get_d()),y2._get_q()); }
    } else {
        if(x2) { res= cmp(y1._get_q(),ExactDouble(x2->raw().get_d())); }
        else { res= cmp(y1._get_q(),y2._get_q()); }
    }
    return res;
}



template<class OP> inline NumberInterface* default_apply(OP op, NumberInterface const* yp1, NumberInterface const* yp2) {
    Handle<NumberInterface> y1(const_cast<NumberInterface*>(yp1)->shared_from_this());
    Handle<NumberInterface> y2(const_cast<NumberInterface*>(yp2)->shared_from_this());
    String yc1=yp1->_class_name(); String yc2=yp2->_class_name();    
    ARIADNE_THROW(DispatchException,op<<"(Number y1, Number y2) with y1="<<*yp1<<", y2="<<*yp2,"No dispatch for "<<op<<"("<<yc1<<", "<<yc2<<")");
}

template<class OP> inline UpperNumberInterface* default_apply(OP op, UpperNumberInterface const* yp1, UpperNumberInterface const* yp2) {
    Handle<UpperNumberInterface> y1(const_cast<UpperNumberInterface*>(yp1)->shared_from_this());
    Handle<UpperNumberInterface> y2(const_cast<UpperNumberInterface*>(yp2)->shared_from_this());
    String yc1=yp1->_class_name(); String yc2=yp2->_class_name();    
    ARIADNE_THROW(DispatchException,op<<"(UpperNumber y1, UpperNumber y2) with y1="<<*yp1<<", y2="<<*yp2,"No dispatch for "<<op<<"("<<yc1<<", "<<yc2<<")");
}



template<class X> struct HasUpper {
    template<class XX, class=decltype(declval<XX>().upper())> static std::true_type test(int);
    template<class XX> static std::false_type test(...);
    static const bool value = decltype(test<X>(1))::value;
};

template<class X> class NumberGetterMixin : public virtual NumberInterface {
  public:
  //  operator X const& () const { return static_cast<Mixin<X,NumberInterface>const&>(*this); }
    static X const& _cast(NumberGetterMixin<X> const& self) { return static_cast<Mixin<X,NumberInterface>const&>(self); }
    static X& _cast(NumberGetterMixin<X>& self) { return static_cast<Mixin<X,NumberInterface>&>(self); }

    typedef Paradigm<X> P;
    friend class Number<P>;

    virtual NumberInterface* _copy() const override { return new NumberWrapper<X>(_cast(*this)); }
    virtual NumberInterface* _move() override { return new NumberWrapper<X>(std::move(_cast(*this))); }

    virtual UpperNumberInterface* _upper() const override { 
        if constexpr (HasUpper<X>::value) { typedef RemoveConst<decltype(declval<X>().upper())> UX; return new UpperNumberWrapper<UX>(_cast(*this).upper()); }
        else { return new UpperNumberWrapper<X>(_cast(*this)); } }
        
    // FIXME: Proper comparisons for ExactNumber.
    virtual LogicalValue _apply(Equal, NumberInterface const* y) const override {
        if (this->_paradigm() == ParadigmCode::EXACT && y->_paradigm() == ParadigmCode::EXACT) {
            return LogicalValue( cmp(*this,*y)==Comparison::EQUAL ? LogicalValue::TRUE : LogicalValue::FALSE ); }
        if (this->_paradigm() == ParadigmCode::VALIDATED && y->_paradigm() == ParadigmCode::VALIDATED) {
            return LogicalValue(this->_get(OrderTag(),dp) == y->_get(OrderTag(),dp)); }
        else {
            return LogicalValue(this->_get(ApproximateTag(),dp) == y->_get(ApproximateTag(),dp)); }
    }
    virtual LogicalValue _apply(Less, NumberInterface const* y) const override {
        if (this->_paradigm() == ParadigmCode::EXACT && y->_paradigm() == ParadigmCode::EXACT) {
            return LogicalValue( cmp(*this,*y)==Comparison::LESS ? LogicalValue::TRUE : LogicalValue::FALSE ); }
        else if (this->_paradigm() == ParadigmCode::VALIDATED && y->_paradigm() == ParadigmCode::VALIDATED) {
            return LogicalValue(this->_get(OrderTag(),dp) < y->_get(OrderTag(),dp)); }
        else {
            return LogicalValue(this->_get(ApproximateTag(),dp) < y->_get(ApproximateTag(),dp));
        }
    }

    virtual Rational _get_q() const override {
        return this->_get_as<Rational>(); }

    virtual FloatDPBall _get(MetricTag,DoublePrecision pr,DoublePrecision pre) const override {
        return this->_get_as<FloatDPBall>(pr,pre); }
    virtual FloatDPBounds _get(OrderTag,DoublePrecision pr) const override {
        return this->_get_as<FloatDPBounds>(pr); }
    virtual FloatDPUpperBound _get(UpperTag,DoublePrecision pr) const override {
        return this->_get_as<FloatDPUpperBound>(pr); }
    virtual FloatDPLowerBound _get(LowerTag,DoublePrecision pr) const override {
        return this->_get_as<FloatDPLowerBound>(pr); }
    virtual FloatDPApproximation _get(ApproximateTag,DoublePrecision pr) const override {
        return this->_get_as<FloatDPApproximation>(pr); }
    virtual FloatMPDPBall _get(MetricTag,MultiplePrecision pr, DoublePrecision pre) const override {
        return this->_get_as<FloatMPDPBall>(pr,pre); }
    virtual FloatMPBall _get(MetricTag, MultiplePrecision pr, MultiplePrecision pre) const override {
        return this->_get_as<FloatMPBall>(pr,pre); }
    virtual FloatMPBounds _get(OrderTag, MultiplePrecision pr) const override {
        return this->_get_as<FloatMPBounds>(pr); }
    virtual FloatMPUpperBound _get(UpperTag, MultiplePrecision pr) const override {
        return this->_get_as<FloatMPUpperBound>(pr); }
    virtual FloatMPLowerBound _get(LowerTag, MultiplePrecision pr) const override {
        return this->_get_as<FloatMPLowerBound>(pr); }
    virtual FloatMPApproximation _get(ApproximateTag, MultiplePrecision pr) const override {
        return this->_get_as<FloatMPApproximation>(pr); }

    virtual ParadigmCode _paradigm() const override { return P::code(); }
    virtual String _class_name() const override { return class_name<X>(); }
    virtual OutputStream& _write(OutputStream& os) const override { return os << _cast(*this); }

  private:
    template<class R, class... PRS, EnableIf<IsConstructible<R,X,PRS...>> = dummy>
        inline R _get_as(PRS... prs) const { return R(_cast(*this),prs...); }
    template<class R, class... PRS, DisableIf<IsConstructible<R,X,PRS...>> = dummy>
        inline R _get_as(PRS... prs) const { 
            std::cerr<<"Warning: Cannot convert " << _cast(*this) << " of type " << this->_class_name() << " to " << class_name<R>();
//            if constexpr(std::tuple_size<Tuple<PRS...>>::value==0) { std::cerr << " with precision " << std::get<0>(prs...); }
            //std::cerr<< " with precision " << pr << " and error precision " << pre << "\n"; 
            std::cerr<< "\n"; throw ParadigmError(); }
};

template<class X> struct DispatchingTraits { typedef Aware<X> AwareOfTypes; };
template<class X> using Awares = typename DispatchingTraits<X>::AwareOfTypes;

template<class X> class Mixin<X,NumberInterface>
    : public UnaryOperationMixin<X,NumberInterface,UnaryOperator>
    , public BinaryOperationMixin<X,NumberInterface,BinaryOperator>
//    , public BinaryOperationMixin<X,LogicalValue,LogicalOperator, NumberInterface>
    , public BinaryOperableMixin<X,NumberInterface,BinaryOperator,Awares<X>>
    , public NumberGetterMixin<X>
{
  public:
    operator X const& () const { return static_cast<Wrapper<X,NumberInterface>const&>(*this); }
    operator X& () { return static_cast<Wrapper<X,NumberInterface>&>(*this); }
};


template<class X> class Wrapper<X,NumberInterface>
    : public X, public Mixin<X,NumberInterface>
{
    static_assert(Not<IsSame<X,Handle<NumberInterface>>>::value,"X must be a concrete number, not a handle");
    static_assert(Not<IsSame<X,Number<Paradigm<X>>>>::value,"X must be a concrete number, not a generic number");
  public:
    Wrapper<X,NumberInterface>(const X& a) : X(a) { }
    Wrapper<X,NumberInterface>(X&& a) : X(std::forward<X>(a)) { }
};


template<class X> class UpperNumberGetterMixin : public virtual UpperNumberInterface {
    typedef Paradigm<X> P;
    static X const& _cast(UpperNumberGetterMixin<X> const& self) { return static_cast<UpperNumberMixin<X>const&>(self); }

  public:
    virtual UpperNumberInterface* _copy() const override { return new UpperNumberWrapper<X>(_cast(*this)); }
    virtual UpperNumberInterface* _move() override { return new UpperNumberWrapper<X>(std::move(_cast(*this))); }

    virtual FloatDPUpperBound _get(DoublePrecision pr) const override { return this->_get_as<FloatDPUpperBound>(pr); }
    virtual FloatMPUpperBound _get(MultiplePrecision pr) const override { return this->_get_as<FloatMPUpperBound>(pr); }

    virtual ParadigmCode _paradigm() const override { assert(false); }
    virtual String _class_name() const override { return class_name<X>(); }
    virtual OutputStream& _write(OutputStream& os) const override { return os << _cast(*this); }

  private:
    template<class R, class... PRS, EnableIf<IsConstructible<R,X,PRS...>> = dummy>
        inline R _get_as(PRS... prs) const { return R(_cast(*this),prs...); }
    template<class R, DisableIf<IsConstructible<R,X>> = dummy>
        inline R _get_as() const { std::cerr<<"Warning: Cannot convert " << _cast(*this) << " of type " << this->_class_name() << " to " << class_name<R>() << "\n"; throw ParadigmError(); }
    template<class R, class PR, DisableIf<IsConstructible<R,X,PR>> = dummy>
        inline R _get_as(PR pr) const { std::cerr<<"Warning: Cannot convert " << _cast(*this) << " of type " << this->_class_name() << " to " << class_name<R>() << " with precision " << pr << "\n"; throw ParadigmError(); }
};

template<class X> class Mixin<X,UpperNumberInterface>
    : public UnaryOperationMixin<X,UpperNumberInterface,MonotoneUnaryOperator>
    , public BinaryOperationMixin<X,UpperNumberInterface,MonotoneBinaryOperator>
//    , public BinaryOperationMixin<X,LogicalValue,LogicalOperator, UpperNumberInterface>
    , public BinaryOperableMixin<X,UpperNumberInterface,MonotoneBinaryOperator,Awares<X>>
    , public UpperNumberGetterMixin<X>
{
  public:
    operator X const& () const { return static_cast<UpperNumberWrapper<X>const&>(*this); }
    operator X& () { return static_cast<UpperNumberWrapper<X>&>(*this); }
};

template<class X> class Wrapper<X,UpperNumberInterface>
    : public X, public Mixin<X,UpperNumberInterface>
{
    static_assert(Not<IsSame<X,Handle<UpperNumberInterface>>>::value,"X must be a concrete number, not a handle");
    static_assert(Not<IsSame<X,UpperNumber<Paradigm<X>>>>::value,"X must be a concrete number, not a generic number");
  public:
    Wrapper(const X& a) : X(a) { }
    Wrapper(X&& a) : X(std::forward<X>(a)) { }
};


} // namespace Ariadne

#endif /* ARIADNE_NUMBER_WRAPPER_HPP */
