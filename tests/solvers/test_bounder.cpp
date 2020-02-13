/***************************************************************************
 *            test_bounder.cpp
 *
 *  Copyright  2018-20  Luca Geretti
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

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "config.hpp"

#include "solvers/bounder.hpp"
#include "function/formula.hpp"
#include "function/taylor_model.hpp"

#include "../test.hpp"

using namespace Ariadne;
using namespace std;

class TestBounder
{
  private:
    std::unique_ptr<BounderInterface> bounder_ptr;
    EffectiveScalarMultivariateFunction v;
    EffectiveScalarMultivariateFunction x, y;
  public:
    TestBounder(const BounderInterface& b)
        : bounder_ptr(b.clone())
    {
        v=EffectiveScalarMultivariateFunction::coordinate(1,0);
        x=EffectiveScalarMultivariateFunction::coordinate(2,0);
        y=EffectiveScalarMultivariateFunction::coordinate(2,1);
    }

    Int test() {
        ARIADNE_TEST_CALL(test_print_name());
        ARIADNE_TEST_CALL(test_D_mismatch());
        ARIADNE_TEST_CALL(test_DA_mismatch());
        ARIADNE_TEST_CALL(test_suggested_step_acceptable());
        ARIADNE_TEST_CALL(test_suggested_step_not_acceptable());
        ARIADNE_TEST_CALL(test_unbounded());
        ARIADNE_TEST_CALL(test_logistic());
        ARIADNE_TEST_CALL(test_vanderpol());
        ARIADNE_TEST_CALL(test_time_invariant_without_parameter());
        ARIADNE_TEST_CALL(test_time_invariant_with_parameter());
        ARIADNE_TEST_CALL(test_time_variant_without_parameter());
        ARIADNE_TEST_CALL(test_time_variant_with_parameter());
        return 0;
    }

    Void test_print_name() {
        BounderInterface const& bounder=*bounder_ptr;
        ARIADNE_TEST_PRINT(bounder);
    }

    Void test_D_mismatch() {
        EffectiveVectorMultivariateFunction f={x,-y};
        ExactBoxType D={ExactIntervalType(-0.25,0.25)};
        StepSizeType hsug(0.25);

        ARIADNE_TEST_FAIL(bounder_ptr->compute(f,D,hsug));
    }

    Void test_DA_mismatch() {
        EffectiveVectorMultivariateFunction f={x,-y};
        ExactBoxType D={ExactIntervalType(-0.25,0.25),ExactIntervalType(-0.25,0.25)};
        ExactBoxType A={ExactIntervalType(-0.25,0.25)};
        StepSizeType hsug(0.25);

        ARIADNE_TEST_FAIL(bounder_ptr->compute(f,D,A,hsug));
    }

    Void test_suggested_step_acceptable() {
        EffectiveVectorMultivariateFunction f={x,-y};
        ExactBoxType dom={ExactIntervalType(-0.25,0.25),ExactIntervalType(-0.25,0.25)};
        StepSizeType hsug(0.25);

        StepSizeType h;
        UpperBoxType B;
        std::tie(h,B) = bounder_ptr->compute(f,dom,hsug);

        ARIADNE_TEST_PRINT(h);
        ARIADNE_TEST_PRINT(B);
        ARIADNE_TEST_EQUAL(h,hsug);
        ARIADNE_TEST_ASSERT(definitely(is_bounded(B)));
    }

    Void test_suggested_step_not_acceptable() {
        EffectiveVectorMultivariateFunction f={x,-y};
        ExactBoxType dom={ExactIntervalType(-0.25,0.25),ExactIntervalType(-0.25,0.25)};
        StepSizeType hsug(1.0);

        StepSizeType h;
        UpperBoxType B;
        std::tie(h,B) = bounder_ptr->compute(f,dom,hsug);

        ARIADNE_TEST_PRINT(h);
        ARIADNE_TEST_PRINT(B);
        ARIADNE_TEST_COMPARE(h,<,hsug);
        ARIADNE_TEST_ASSERT(definitely(is_bounded(B)));
    }

    Void test_unbounded() {
        // Attempt to find a bounding box for the flow
        ExpandContractBounderConfiguration& config=dynamic_cast<ExpandContractBounder&>(*bounder_ptr).configuration();
        auto minimum_step_size=config.minimum_step_size();

        EffectiveVectorMultivariateFunction f={sqr(v)};
        ExactBoxType dom={ExactIntervalType(1.0,1.0)};
        StepSizeType hsug(12.0);

        ARIADNE_TEST_PRINT(bounder_ptr->compute(f,dom,hsug));
        config.set_minimum_step_size(1.5_x);
        ARIADNE_TEST_THROWS(bounder_ptr->compute(f,dom,hsug),BoundingNotFoundException);
        config.set_minimum_step_size(minimum_step_size);
    }

    Void test_logistic() {
        EffectiveVectorMultivariateFunction f={v*(1-v)};

        // Clear that bounds are in [0.25:0.75]
        // For x in [a:b] with 0<=a<0.5<b, maximum f(x) is 0.25 at x=0.5.
        List<BoxDomainType> domains = { {{0.25,0.5}}, {{0.0,0.25}}, {{0.5,0.75}} };
        for (auto D : domains) {
            StepSizeType hsug(1.0);
            StepSizeType h;
            UpperBoxType B;
            for (int i=0; i!=5; ++i) {
                std::tie(h,B) = bounder_ptr->compute(f,D,hsug);
                auto DhB=make_tuple(D,h,B);
                ARIADNE_TEST_PRINT(DhB);
                hsug=hlf(hsug);
            }
        }
    }

    Void test_vanderpol() {
        RealConstant mu("mu",1.0_dec);
        EffectiveVectorMultivariateFunction f={y,mu*y*(1-sqr(x))-x};

        BoxDomainType D({{1.25,1.55},{2.35,2.45}});
        StepSizeType hsug(1.0);

        // Clear that bounds are in [1:3]
        // For x in [1:a] with a>2, maximum f(x) is 4 at x=2.

        StepSizeType h;
        UpperBoxType B;
        for (int i=0; i!=5; ++i) {
            std::tie(h,B) = bounder_ptr->compute(f,D,hsug);
            auto DhB=make_tuple(D,h,B);
            ARIADNE_TEST_PRINT(DhB);
            hsug=hlf(hsug);
        }
    }

    Void test_time_invariant_without_parameter() {
        EffectiveScalarMultivariateFunction x0=EffectiveScalarMultivariateFunction::coordinate(2,0);
        EffectiveScalarMultivariateFunction x1=EffectiveScalarMultivariateFunction::coordinate(2,1);

        EffectiveVectorMultivariateFunction f={x0,-x1};
        ExactBoxType D={ExactIntervalType(-0.25,0.25),ExactIntervalType(-0.25,0.25)};
        StepSizeType td(1.0);
        StepSizeType hsug(1.0);

        StepSizeType h;
        UpperBoxType B;

        std::tie(h,B) = bounder_ptr->compute(f,D,hsug);

        ARIADNE_TEST_COMPARE(h,<=,hsug);
        ARIADNE_TEST_ASSERT(definitely(is_bounded(B)));
    }

    Void test_time_invariant_with_parameter() {
        EffectiveScalarMultivariateFunction x0=EffectiveScalarMultivariateFunction::coordinate(3,0);
        EffectiveScalarMultivariateFunction x1=EffectiveScalarMultivariateFunction::coordinate(3,1);
        EffectiveScalarMultivariateFunction a0=EffectiveScalarMultivariateFunction::coordinate(3,2);

        EffectiveVectorMultivariateFunction f={x0+a0,-x1};
        ExactBoxType D={ExactIntervalType(-0.25,0.25),ExactIntervalType(-0.25,0.25)};
        ExactBoxType A={ExactIntervalType(-0.25,0.25)};
        StepSizeType td(1.0);
        StepSizeType hsug(1.0);

        StepSizeType h;
        UpperBoxType B;

        std::tie(h,B) = bounder_ptr->compute(f,D,A,hsug);

        ARIADNE_TEST_COMPARE(h,<=,hsug);
        ARIADNE_TEST_ASSERT(definitely(is_bounded(B)));
    }

    Void test_time_variant_without_parameter() {
        EffectiveScalarMultivariateFunction x0=EffectiveScalarMultivariateFunction::coordinate(3,0);
        EffectiveScalarMultivariateFunction x1=EffectiveScalarMultivariateFunction::coordinate(3,1);
        EffectiveScalarMultivariateFunction t=EffectiveScalarMultivariateFunction::coordinate(3,2);

        EffectiveVectorMultivariateFunction f={x0,-x1+2*t};
        ExactBoxType D={ExactIntervalType(-0.25,0.25),ExactIntervalType(-0.25,0.25)};
        StepSizeType td(1.0);
        StepSizeType hsug(1.0);

        StepSizeType h;
        UpperBoxType B;

        std::tie(h,B) = bounder_ptr->compute(f,D,td,hsug);

        ARIADNE_TEST_COMPARE(h,<=,hsug);
        ARIADNE_TEST_ASSERT(definitely(is_bounded(B)));
    }

    Void test_time_variant_with_parameter() {
        EffectiveScalarMultivariateFunction x0=EffectiveScalarMultivariateFunction::coordinate(4,0);
        EffectiveScalarMultivariateFunction x1=EffectiveScalarMultivariateFunction::coordinate(4,1);
        EffectiveScalarMultivariateFunction t=EffectiveScalarMultivariateFunction::coordinate(4,2);
        EffectiveScalarMultivariateFunction a0=EffectiveScalarMultivariateFunction::coordinate(4,3);

        EffectiveVectorMultivariateFunction f={x0+a0,-x1+2*t};
        ExactBoxType D={ExactIntervalType(-0.25,0.25),ExactIntervalType(-0.25,0.25)};
        ExactBoxType A={ExactIntervalType(-0.25,0.25)};
        StepSizeType td(1.0);
        StepSizeType hsug(1.0);

        StepSizeType h;
        UpperBoxType B;

        std::tie(h,B) = bounder_ptr->compute(f,D,td,A,hsug);

        ARIADNE_TEST_COMPARE(h,<=,hsug);
        ARIADNE_TEST_ASSERT(definitely(is_bounded(B)));
    }
};

Int main(Int argc, const char* argv[]) {

    EulerBounder euler_bounder;

    List<BounderHandle> bounders = { euler_bounder };

    for (BounderHandle bounder : bounders) {
        TestBounder(bounder).test(); }

    return ARIADNE_TEST_FAILURES;
}
