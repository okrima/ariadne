/***************************************************************************
 *            check_number.h
 *
 *  Copyright 2009-14  Pieter Collins
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

#include "test.h"
#include "test/utility.h"
#include "expression/operators.h"

using namespace Ariadne;

namespace Ariadne {

template<typename T1, typename T2, typename = Fallback> struct HasPow : False { };
template<typename T1, typename T2> struct HasPow<T1,T2,EnableIf<True, decltype(pow(declval<T1&>(),declval<T2&>()),Fallback())>> : True { };

//template<typename Op, typename T, typename = Fallback> struct HasUnaryOperator : False { };
//template<typename Op, typename T> struct HasUnaryOperator<Op,T,EnableIf<True,decltype(Op()(declval<T&>()),Fallback())>> : True { };
//template<typename Op, typename T1, typename T2, typename = Fallback>  struct HasBinaryOperator : False { };
//template<typename Op, typename T1, typename T2>  struct HasBinaryOperator<Op,T1,T2,EnableIf<True,decltype(Op()(declval<T1>(),declval<T2>()),Fallback())>> : True { };

template<class N> class CheckNumberConcepts {
  public:
    void check_signed_concept() {
        ARIADNE_TEST_STATIC_ASSERT(HasOperatorReturning<N,OperatorPositive,N>);
        ARIADNE_TEST_STATIC_ASSERT(HasOperatorReturning<N,OperatorNegative,N>);
        ARIADNE_TEST_STATIC_ASSERT(HasOperator<Pos,N,Return<N>>);
        ARIADNE_TEST_STATIC_ASSERT(HasOperator<Neg,N,Return<N>>);
    }

    void check_semiring_concept() {
        ARIADNE_TEST_STATIC_ASSERT(IsDefaultConstructible<N>);
        ARIADNE_TEST_STATIC_ASSERT(IsConstructible<N,uint>);
        ARIADNE_TEST_STATIC_ASSERT(HasOperator<OperatorPlus,N,N>);
        ARIADNE_TEST_STATIC_ASSERT(HasOperator<OperatorTimes,N,N>);
        ARIADNE_TEST_STATIC_ASSERT(HasOperator<Pow,N,uint>);
        ARIADNE_TEST_STATIC_ASSERT(HasOperator<Sqr,N>);
        ARIADNE_TEST_STATIC_ASSERT(HasOperator<Add,N,N>);
        ARIADNE_TEST_STATIC_ASSERT(HasOperator<Mul,N,N>);
    }

    void check_ring_concept() {
        check_semiring_concept();
        check_signed_concept();
        ARIADNE_TEST_STATIC_ASSERT(IsConstructible<N,int>);
        ARIADNE_TEST_STATIC_ASSERT(HasOperator<OperatorNegate,N>);
        ARIADNE_TEST_STATIC_ASSERT(HasOperator<OperatorMinus,N,N>);
        ARIADNE_TEST_STATIC_ASSERT(HasOperator<Neg,N>);
        ARIADNE_TEST_STATIC_ASSERT(HasOperator<Sub,N,N>);
    }

    void check_field_concept() {
        check_ring_concept();
        if (IsSame<Paradigm<N>,Approximate>::value) {
            ARIADNE_TEST_STATIC_ASSERT(IsConstructible<N,double>);
        }
        ARIADNE_TEST_STATIC_ASSERT(HasOperator<OperatorDivides,N,N>);
        ARIADNE_TEST_STATIC_ASSERT(HasOperator<Rec,N>);
        ARIADNE_TEST_STATIC_ASSERT(HasOperator<Div,N,N>);
        ARIADNE_TEST_STATIC_ASSERT(HasOperator<Pow,N,int>);
    };

    void check_transcendental_concept() {
        ARIADNE_TEST_STATIC_ASSERT(HasOperator<Sqrt,N>);
        ARIADNE_TEST_STATIC_ASSERT(HasOperator<Exp,N>);
        ARIADNE_TEST_STATIC_ASSERT(HasOperator<Log,N>);
        ARIADNE_TEST_STATIC_ASSERT(HasOperator<Sin,N>);
        ARIADNE_TEST_STATIC_ASSERT(HasOperator<Cos,N>);
        ARIADNE_TEST_STATIC_ASSERT(HasOperator<Tan,N>);
        ARIADNE_TEST_STATIC_ASSERT(HasOperator<Atan,N>);
    }

    void check_absolute_concept() {
        ARIADNE_TEST_STATIC_ASSERT(HasOperator<Abs,N>);
        ARIADNE_TEST_STATIC_ASSERT(HasOperator<Max,N,N>);
        ARIADNE_TEST_STATIC_ASSERT(HasOperator<Min,N,N>);
    }

    void check_equality_comparible_concept() {
        ARIADNE_TEST_STATIC_ASSERT(IsEqualityComparible<N,N>);
        ARIADNE_TEST_STATIC_ASSERT(HasOperator<Equal,N,N>);
        ARIADNE_TEST_STATIC_ASSERT(HasOperator<NotEqual,N,N>);
    }

    void check_order_concept() {
        ARIADNE_TEST_STATIC_ASSERT(IsLessThanCompartible<N,N>);
        ARIADNE_TEST_STATIC_ASSERT(HasOperator<Less,N,N>);
        ARIADNE_TEST_STATIC_ASSERT(HasOperator<LessEqual,N,N>);
        ARIADNE_TEST_STATIC_ASSERT(HasOperator<Greater,N,N>);
        ARIADNE_TEST_STATIC_ASSERT(HasOperator<GreaterEqual,N,N>);
    }

    void check_comparible_concept() {
        check_equality_comparible_concept();
        check_order_concept();
    }
};

} // namespace Ariadne
