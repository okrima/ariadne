/***************************************************************************
 *            utility/dispatching.hpp
 *
 *  Copyright 2017-18  Pieter Collins
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

/*! \file utility/dispatching.hpp
 *  \brief Classes for facilitating double-dispatching
 */



#ifndef ARIADNE_DISPATCHING_HPP
#define ARIADNE_DISPATCHING_HPP

namespace Ariadne {

template<class... AWS> struct Aware;

template<class X, class I> class Mixin;
template<class X, class I> class Wrapper;

template<class X, class I> inline X const* extract(I const* y) {
     return dynamic_cast<Wrapper<X,I>const*>(y);
}


template<class X, class I, class OP, class J=I> struct UnaryOperationMixin;
template<class X, class I, class OP, class J=I> struct BinaryOperationMixin;
template<class X, class I, class OP, class Y> struct BinaryOperable;

template<class X, class I, class OP> inline X const& _upcast(UnaryOperationMixin<X,I,OP> const& y) { 
    return static_cast<Mixin<X,I> const&>(y); }
template<class X, class I, class OP> inline X const& _upcast(BinaryOperationMixin<X,I,OP> const& y) { 
    return static_cast<Mixin<X,I> const&>(y); }
template<class X, class I, class OP, class Y> inline X const& _upcast(BinaryOperable<X,I,OP,Y> const& y) { 
    return static_cast<Mixin<X,I> const&>(y); }
    
template<class I, class R> inline I* _make_wrapper(R&& r, I*) { return new Wrapper<R,I>(r); }

template<class I, class R> inline I* _make_wrapper(R&& r) { I* p=nullptr; return _make_wrapper<I,R>(std::forward<R>(r),p); }


//------------ Single-dispatching code -----------------------

template<class X, class I, class OP, class J> struct UnaryOperationMixin : public virtual J {
    virtual I* _apply(OP op) const final { return _make_wrapper<I>(op(_upcast<X>(*this))); }
};

//------------ Double-dispatching code -----------------------

template<class I, class OP, class Y> struct BinaryOperableInterface {
    virtual ~BinaryOperableInterface() = default;
    virtual I* _apply_left(OP op, Y const& y) const = 0;
    virtual I* _apply_right(OP op, Y const& y) const = 0;
};
template<class X, class I, class OP, class Y> struct BinaryOperable : virtual BinaryOperableInterface<I,OP,Y> {
    virtual I* _apply_left(OP op, Y const& other) const { return _make_wrapper<I>(op(_upcast<X>(*this),other)); }
    virtual I* _apply_right(OP op, Y const& other) const { return _make_wrapper<I>(op(other,_upcast<X>(*this))); }
};
template<class X, class I, class OP, class AW> struct BinaryOperableMixin;
template<class X, class I, class OP, class Y, class... YS> struct BinaryOperableMixin<X,I,OP,Aware<Y,YS...>>
    : BinaryOperable<X,I,OP,Y>, BinaryOperableMixin<X,I,OP,Aware<YS...>> { };
template<class X, class I, class OP> struct BinaryOperableMixin<X,I,OP,Aware<>> { };

/*
template<class I, class OP> struct BinaryOperationInterface {
    virtual ~BinaryOperationInterface() = default;
    virtual I* _apply(OP op, I const* other) const = 0;
    virtual I* _rapply(OP op, I const* other) const = 0;
};
*/

template<class X, class I, class OP, class J> struct BinaryOperationMixin : public virtual J {
    virtual I* _apply(OP op, I const* other) const final;
    virtual I* _rapply(OP op, I const* other) const final;
};

template<class OP, class I> I* default_apply(OP op, I const* arg1, I const* arg2) {
    String ic = class_name<I>(); String arg1c=arg1->_class_name(); String arg2c=arg2->_class_name();    
    ARIADNE_THROW(DispatchException,op<<"(" << ic << " arg1, " << ic << " arg2) with arg1="<<*arg1<<", arg2="<<*arg2,"No dispatch for "<<op<<"("<<arg1c<<", "<<arg2c<<")");
}

template<class X, class I, class OP, class J> inline I* BinaryOperationMixin<X,I,OP,J>::_apply(OP op, I const* other) const {
    auto aware_other=dynamic_cast<BinaryOperableInterface<I,OP,X>const*>(other);
    if(aware_other) { X const& self = _upcast<X>(*this); return aware_other->_apply_right(op,self); }
    else { return other->_rapply(op,this); }
}
template<class X, class I, class OP, class J> inline I* BinaryOperationMixin<X,I,OP,J>::_rapply(OP op, I const* other) const {
    auto aware_other=dynamic_cast<BinaryOperableInterface<I,OP,X>const*>(other);
    if(aware_other) { X const& self = _upcast<X>(*this); return aware_other->_apply_left(op,self); }
    else { return default_apply(op,other,this); }
}


} // namespace Ariadne

#endif /* ARIADNE_NUMBER_WRAPPER_HPP */
