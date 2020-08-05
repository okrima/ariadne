/***************************************************************************
 *            function/function_patch_mixin.hpp
 *
 *  Copyright  2011-20  Pieter Collins
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

/*! \file function/function_patch_mixin.hpp
 *  \brief Mixin for functions on bounded sets.
 */

#ifndef ARIADNE_FUNCTION_PATCH_MIXIN_HPP
#define ARIADNE_FUNCTION_PATCH_MIXIN_HPP

#include "../numeric/number.decl.hpp"
#include "../function/function.decl.hpp"

#include "../numeric/operators.hpp"
#include "../numeric/numeric.hpp"

#include "../algebra/vector.hpp"
#include "../algebra/matrix.hpp"
#include "../algebra/range.hpp"
#include "../algebra/operations.hpp"
#include "../algebra/algebra_interface.hpp"
#include "../algebra/algebra_mixin.hpp"

#include "../function/domain.hpp"
#include "../function/function_interface.hpp"
#include "../function/function_patch_interface.hpp"
#include "../function/function_mixin.hpp"


namespace Ariadne {

template<class FM, class P, class D, class C> class FunctionPatchMixin;
template<class FM, class P, class D> using ScalarFunctionPatchMixin = FunctionPatchMixin<FM,P,D,IntervalDomainType>;
template<class FM, class P, class D> using VectorFunctionPatchMixin = FunctionPatchMixin<FM,P,D,BoxDomainType>;

template<class FM, class P, class D> class FunctionPatchMixin<FM,P,D,IntervalDomainType>
    : public virtual ScalarFunctionPatchInterface<P,D>
    , public ScalarFunctionMixin<FM,P,D>
    , public ElementaryAlgebraMixin<FM,Number<P>>
{
    using C=IntervalDomainType;
    typedef D DomainType;
    typedef C CodomainType;
    typedef PositiveValidatedUpperNumber NormType;
    typedef ElementTraits<C>::BoundedRangeType RangeType;
    typedef Number<P> X;
  public:
    DomainType const domain() const override { return static_cast<FM const&>(*this).domain(); }
    CodomainType const codomain() const override { return static_cast<FM const&>(*this).codomain(); }

    ScalarFunctionPatchInterface<P,D>* _clone() const override {
        return new FM(static_cast<const FM&>(*this)); }

    ScalarFunctionPatchInterface<P,D>* _create_copy() const override {
        return new FM(static_cast<const FM&>(*this)); }
    ScalarFunctionPatchInterface<P,D>* _create_zero() const override {
        return new FM(factory(static_cast<FM const&>(*this)).create_zero()); }
    ScalarFunctionPatchInterface<P,D>* _create_constant(Number<P> const& c) const override {
        return new FM(factory(static_cast<FM const&>(*this)).create_constant(c)); }

    RangeType const _generic_range() const override {
        return RangeType(static_cast<const FM&>(*this).range()); }
    NormType const _generic_norm() const override {
        return static_cast<const FM&>(*this).norm(); }
    ScalarFunctionPatchInterface<P,D>* _antiderivative(SizeType j) const override {
        return new FM(antiderivative(static_cast<const FM&>(*this),j)); }
    ScalarFunctionPatchInterface<P,D>* _antiderivative(SizeType j, Number<P> c) const override {
        return new FM(antiderivative(static_cast<const FM&>(*this),j,c)); }
     ScalarFunctionPatchInterface<P,D>* _restriction(const BoxDomainType& d) const override {
        return new FM(restriction(static_cast<const FM&>(*this),d)); }
    Number<P> _unchecked_evaluate(const Vector<Number<P>>& x) const override {
        return unchecked_evaluate(static_cast<const FM&>(*this),x); }
    ScalarFunctionPatchInterface<P,D>* _partial_evaluate(SizeType j, const Number<P>& c) const override {
        return heap_copy(partial_evaluate(static_cast<const FM&>(*this),j,c)); }
    ScalarFunctionPatchInterface<P,D>* _embed(const BoxDomainType& d1, const BoxDomainType& d2) const override {
        return new FM(embed(d1,static_cast<const FM&>(*this),d2)); }

    OutputStream& _write(OutputStream& os) const override {
        return os << static_cast<FM const&>(*this); }
};


template<class FM, class P, class D> class FunctionPatchMixin<FM,P,D,BoxDomainType>
    : public virtual VectorFunctionPatchInterface<P,D>
    , public VectorFunctionMixin<FM,P,D>
{
    using C=BoxDomainType;
    typedef D DomainType;
    typedef C CodomainType;
    typedef PositiveValidatedUpperNumber NormType;
    typedef ElementTraits<C>::BoundedRangeType RangeType;
    typedef typename Element<FM>::Type ScalarMultivariateFunctionType;
  public:
    DomainType const domain() const override {
        return static_cast<FM const&>(*this).domain(); }
    CodomainType const codomain() const override {
        return static_cast<FM const&>(*this).codomain(); }

    virtual VectorFunctionPatchInterface<P,D>* _clone() const override { return new FM(static_cast<const FM&>(*this)); }
    virtual Void _set(SizeType i, const ScalarFunctionPatchInterface<P,D>& sf) override {
        if(!dynamic_cast<const typename FM::ScalarMultivariateFunctionType*>(&sf)) {
            ARIADNE_FAIL_MSG("Cannot set element of VectorMultivariateFunctionPatch "<<*this<<" to "<<sf<<"\n"); }
        static_cast<FM&>(*this).FM::set(i,dynamic_cast<const ScalarMultivariateFunctionType&>(sf)); }
    virtual VectorFunctionPatchInterface<P,D>* _derivative(SizeType j) const override {
        ARIADNE_NOT_IMPLEMENTED; }
    RangeType const _generic_range() const override {
        return RangeType(static_cast<const FM&>(*this).range()); }
    NormType const _generic_norm() const override {
         return static_cast<const FM&>(*this).norm(); }
    VectorFunctionPatchInterface<P,D>* _embed(const BoxDomainType& d1, const BoxDomainType& d2) const override {
        return heap_copy(embed(d1,static_cast<const FM&>(*this),d2)); }
    VectorFunctionPatchInterface<P,D>* _restriction(const BoxDomainType& d) const override {
        return new FM(restriction(static_cast<const FM&>(*this),d)); }
    Void _adjoin(const ScalarFunctionPatchInterface<P,D>& f) override {
        static_cast<FM&>(*this).FM::adjoin(dynamic_cast<const ScalarMultivariateFunctionType&>(f)); }
    VectorFunctionPatchInterface<P,D>* _join(const VectorFunctionPatchInterface<P,D>& f) const override {
        return heap_copy(join(static_cast<const FM&>(*this),dynamic_cast<const FM&>(f))); }
    VectorFunctionPatchInterface<P,D>* _combine(const VectorFunctionPatchInterface<P,D>& f) const override {
        return heap_copy(combine(static_cast<const FM&>(*this),dynamic_cast<const FM&>(f))); }
    Vector<Number<P>> _unchecked_evaluate(const Vector<Number<P>>& x) const override {
        return unchecked_evaluate(static_cast<const FM&>(*this),x); }
/*
    ScalarFunctionPatchInterface<P,D>* _compose(const ScalarMultivariateFunctionInterface<P>& f) const override {
        return heap_copy(compose(f,static_cast<const FM&>(*this))); }
    VectorFunctionPatchInterface<P,D>* _compose(const VectorMultivariateFunctionInterface<P>& f) const override {
        return heap_copy(compose(f,static_cast<const FM&>(*this))); }
    ScalarFunctionPatchInterface<P,D>* _unchecked_compose(const ScalarMultivariateFunctionInterface<P>& f) const override {
        return heap_copy(unchecked_compose(dynamic_cast<const ScalarMultivariateFunctionType&>(f),static_cast<const FM&>(*this))); }
    VectorFunctionPatchInterface<P,D>* _unchecked_compose(const VectorMultivariateFunctionInterface<P>& f) const override {
        return heap_copy(unchecked_compose(dynamic_cast<const FM&>(f),static_cast<const FM&>(*this))); }
*/
    VectorFunctionPatchInterface<P,D>* _partial_evaluate(SizeType j, const Number<P>& c) const override {
        return heap_copy(partial_evaluate(static_cast<const FM&>(*this),j,c)); }
};


template<class FCTRY, class P> class FunctionPatchFactoryMixin
    : public FunctionPatchFactoryInterface<P>
{
    typedef BoxDomainType VD;
    typedef IntervalDomainType SD;
  public:
    typedef VD VectorDomainType;
    typedef SD ScalarDomainType;

    virtual FunctionPatchFactoryInterface<P>* clone() const override { return new FCTRY(this->upcast()); }
    virtual OutputStream& _write(OutputStream& os) const override { return os << this->upcast(); }
  private:
    template<class T> static inline T* heap_move(T&& t) { return new T(std::forward<T>(t)); }
    inline FCTRY const& upcast() const { return static_cast<FCTRY const&>(*this); }
    virtual Number<P> _create(const Number<P>& number) const override {
        return this->upcast().create(number); }
    virtual ScalarFunctionPatchInterface<P,VD>* _create(const VectorDomainType& domain, const ScalarMultivariateFunctionInterface<P>& function) const override {
        return heap_move(this->upcast().create(domain,function)); };
    virtual VectorFunctionPatchInterface<P,VD>* _create(const VectorDomainType& domain, const VectorMultivariateFunctionInterface<P>& function) const override {
        return heap_move(this->upcast().create(domain,function)); };

    virtual ScalarFunctionPatchInterface<P,VD>* _create_zero(const VectorDomainType& domain) const override {
        return heap_move(this->upcast().create_zero(domain)); };
    virtual ScalarFunctionPatchInterface<P,VD>* _create_constant(const VectorDomainType& domain, const Number<P>& value) const override {
        return heap_move(this->upcast().create_constant(domain,value)); };
    virtual ScalarFunctionPatchInterface<P,VD>* _create_coordinate(const VectorDomainType& domain, SizeType index) const override {
        return heap_move(this->upcast().create_coordinate(domain,index)); };
    virtual VectorFunctionPatchInterface<P,VD>* _create_zeros(SizeType rsize, const VectorDomainType& domain) const override {
        return heap_move(this->upcast().create_zeros(rsize,domain)); };
    virtual VectorFunctionPatchInterface<P,VD>* _create_constants(const VectorDomainType& domain, const Vector<Number<P>>& values) const override {
        return heap_move(this->upcast().create_constants(domain,values)); };
    virtual VectorFunctionPatchInterface<P,VD>* _create_projection(const VectorDomainType& domain, Range indices) const override {
        return heap_move(this->upcast().create_projection(domain,indices)); };
    virtual VectorFunctionPatchInterface<P,VD>* _create_identity(const VectorDomainType& domain) const override {
        return heap_move(this->upcast().create_identity(domain)); };
};


} // namespace Ariadne

#endif
