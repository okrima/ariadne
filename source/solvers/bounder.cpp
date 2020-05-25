/***************************************************************************
 *            solvers/bounder.cpp
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

#include "bounder.hpp"

#include "../function/formula.hpp"
#include "../function/taylor_model.hpp"

namespace Ariadne {



BounderConfiguration::BounderConfiguration() {
    _minimum_step_size=0;
    _maximum_step_size=Dyadic::inf();
    _lipschitz_tolerance=Dyadic::inf();
}

OutputStream& BounderConfiguration::_write(OutputStream& os) const {
    return os << "BounderConfiguration("
        "minimum_step_size="<<minimum_step_size()<<", "
        "maximum_step_size="<<maximum_step_size()<<", "
        "lipschitz_tolerance="<<lipschitz_tolerance()<<")";
}


BounderBase::BounderBase(BounderBaseConfiguration* config) : _configuration(config) { }

Pair<StepSizeType,UpperBoxType> BounderBase::compute(ValidatedVectorMultivariateFunction const& f, BoxDomainType const& D, StepSizeType const& hsug) const {
    return this->compute(f,D,BoxDomainType(0u),hsug);
}

Pair<StepSizeType,UpperBoxType> BounderBase::compute(ValidatedVectorMultivariateFunction const& f, BoxDomainType const& D, StepSizeType const& t, StepSizeType const& hsug) const {
    return this->compute(f,D,t,BoxDomainType(0u),hsug);
}


Pair<StepSizeType,UpperBoxType> BounderBase::compute(ValidatedVectorMultivariateFunction const& f, BoxDomainType const& D, BoxDomainType const& A, StepSizeType const& hsug) const {
    ARIADNE_PRECONDITION(f.result_size()==D.dimension());
    ARIADNE_PRECONDITION(f.argument_size()==D.dimension()+A.dimension());

    return this->compute(f,D,0,A,hsug);

}


ExpandContractBounderConfiguration::ExpandContractBounderConfiguration()
    : BounderBaseConfiguration()
{
    _expansion_steps=4;
    _refinement_steps=4;
    _domain_widening=1.25_x;
    _starting_widening=2;
    _expanding_widening=1.125_x;
}

OutputStream& ExpandContractBounderConfiguration::_write(OutputStream& os) const {
    return os << "ExpandContractBounderConfiguration("
        "minimum_step_size="<<minimum_step_size()<<", "
        "maximum_step_size="<<maximum_step_size()<<", "
        "lipschitz_tolerance="<<lipschitz_tolerance()<<", "
        "expansion_steps="<<expansion_steps()<<", "
        "refinement_steps="<<refinement_steps()<<", "
        "domain_widening="<<domain_widening()<<", "
        "starting_widening="<<starting_widening()<<", "
        "expanding_widening="<<expanding_widening()<<")";
}


ExpandContractBounder::ExpandContractBounder() : ExpandContractBounder(new ExpandContractBounderConfiguration()) { }

ExpandContractBounder::ExpandContractBounder(ExpandContractBounderConfiguration* config) : BounderBase(config) { }


Pair<StepSizeType,UpperBoxType> ExpandContractBounder::compute(ValidatedVectorMultivariateFunction const& f, BoxDomainType const& D, StepSizeType const& t, BoxDomainType const& A, StepSizeType const& hsug) const {
    ARIADNE_PRECONDITION(f.result_size()==D.dimension());
    ARIADNE_PRECONDITION(f.argument_size()==D.dimension()+1u+A.dimension() or (t==0 and f.argument_size()==D.dimension()+A.dimension()));

    auto config=this->configuration();

    StepSizeType h=hsug;

    h=min(hsug,config.maximum_step_size());

    // Apply condition that h*L should not be higher than lipschitz_tolerance
    if(is_finite(config.lipschitz_tolerance())) {
        FloatDPUpperBound lipschitz = norm(f.jacobian(Vector<FloatDPBounds>(cast_singleton(product(D,to_time_bounds(t,t+h),A))))).upper();
        StepSizeType hlip = static_cast<StepSizeType>(cast_exact(config.lipschitz_tolerance()/lipschitz));
        h=min(hlip,h);
    }

    IntervalDomainType T = to_time_bounds(t,t+h);
    auto wD=D+(config.domain_widening()-1u)*(D-D.midpoint());
    UpperBoxType B=D;

    Bool success=false;
    StepSizeType hprev=h*1.5_dy;
    while(!success) {
        // Initially, have no estimate of the bounding box, so use the initial domain
        // and widen the estimated reached set accordingly
        // Revert to using D for the initial bounds estimate B when reducing the step size
        // since otherwise B could be huge from the previous step
        auto wf=config.starting_widening()*f;
        B=this->_formula(wf,D,T,A,wD);
        for(Nat i=0; i<config.expansion_steps(); ++i) {
            // Check if a flow step from D over time range T stays within D
            UpperBoxType Br=this->_formula(f,D,T,A,B);
            if(not definitely(is_bounded(Br))) {
                success=false;
                break;
            } else if(refines(Br,B)) {
                B=Br;
                success=true;
                break;
            } else {
                wf=config.expanding_widening()*f;
                B=this->_formula(wf,D,T,A,B);
            }
        }
        if(!success) {
            StepSizeType hnew=hlf(hprev);
            hprev=h;
            h=hnew;
            if (h < config.minimum_step_size())
                ARIADNE_THROW(BoundingNotFoundException,"EulerBounder::_compute","The step size is lower than the minimum (" << config.minimum_step_size() << ") allowed, bounding could not be found.");

            T = to_time_bounds(t,t+h);
        }
    }

    for(Nat i=0; i<config.refinement_steps(); ++i) {
        B = this->_formula(f,D,T,A,B);
    }

    return std::make_pair(h,B);
}

UpperBoxType EulerBounder::_formula(ValidatedVectorMultivariateFunction const& f, BoxDomainType const& D, IntervalDomainType const& T, BoxDomainType const& A, UpperBoxType const& B) const {
    UpperIntervalType const& rT=reinterpret_cast<UpperIntervalType const&>(T);
    UpperBoxType const& rA=reinterpret_cast<UpperBoxType const&>(A);
    UpperIntervalType const rH=rT-T.lower();

    const bool is_autonomous = (f.argument_size() == D.dimension()+A.dimension());
    UpperBoxType dom = is_autonomous ? product(B,rA) : product(B,rT,rA);

    return D+rH*cast_vector(apply(f,dom));
}


} // namespace Ariadne;
