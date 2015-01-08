/***************************************************************************
 *            hybrid_enclosure.cc
 *
 *  Copyright  2009-10  Pieter Collins
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

#include "function/functional.h"
#include "config.h"

#include "numeric/numeric.h"
#include "algebra/vector.h"
#include "expression/operators.h"
#include "function/function.h"
#include "function/constraint.h"
#include "function/polynomial.h"
#include "function/taylor_function.h"
#include "geometry/box.h"
#include "geometry/grid_set.h"
#include "hybrid/hybrid_time.h"
#include "hybrid/discrete_event.h"
#include "hybrid/discrete_location.h"

#include "solvers/linear_programming.h"
#include "solvers/nonlinear_programming.h"
#include "solvers/constraint_solver.h"
#include "geometry/affine_set.h"

#include "output/graphics_interface.h"
#include "hybrid/hybrid_enclosure.h"
#include "hybrid/hybrid_set.h"


namespace Ariadne {


OutputStream& operator<<(OutputStream& os, const EnclosureVariableType& evt) {
    switch (evt) {
        case INITIAL: return os << "x";
        case TEMPORAL: return os << "t";
        //case CROSSING: return os << "t";
        //case STEP: return os << "h";
        case PARAMETER: return os << "a";
        case INPUT: return os << "u";
        case NOISE: return os << "v";
        case ERROR: return os << "e";
        case UNKNOWN: default: return os << "s";
    }
}

template<class T> StringType str(const T& t) { StringStream ss; ss<<t; return ss.str(); }

List<String> variable_names(const List<EnclosureVariableType>& vt) {
    std::map<EnclosureVariableType,Nat> counts;
    List<String> result;
    for(Nat i=0; i!=vt.size(); ++i) {
        result.append( str(vt[i]) + str(counts[vt[i]]++) );
    }
    return result;
}

struct Variables2d {
    RealVariable _x,_y;
    Variables2d(const RealVariable& x, const RealVariable& y) : _x(x), _y(y) { }
    RealVariable const& x_variable() const { return this->_x; };
    RealVariable const& y_variable() const { return this->_y; };
};


//-------------- HybridEnclosure -----------------------------------------//

HybridEnclosure::~HybridEnclosure() {
}

HybridEnclosure::HybridEnclosure()
    : _location(""), _events(), _space(), _set(), _variables()
{
}

HybridEnclosure::HybridEnclosure(const HybridBoundedConstraintSet& hybrid_set,
                                 const RealSpace& space,
                                 const IntervalFunctionModelFactoryInterface& factory)
    : _location(hybrid_set.location()), _events(), _space(space.variable_names()), _set(),
      _variables(space.dimension(),INITIAL)
{
    BoundedConstraintSet euclidean_set=hybrid_set.euclidean_set(this->_location,space);
    this->_set=Enclosure(euclidean_set,factory);
}

HybridEnclosure::HybridEnclosure(const DiscreteLocation& location, const RealSpace& space,
                                 const ExactBox& box, const IntervalFunctionModelFactoryInterface& factory)
    : _location(location), _events(), _space(space.variable_names()), _set(box,factory),
      _variables(box.dimension(),INITIAL)
{
}

HybridEnclosure::HybridEnclosure(const HybridBox& hbox, const IntervalFunctionModelFactoryInterface& factory)
    : _location(hbox.location()), _events(), _space(hbox.space().variable_names()), _set(hbox.continuous_set(),factory),
      _variables(hbox.continuous_set().dimension(),INITIAL)
{
}



HybridEnclosure::HybridEnclosure(const DiscreteLocation& location, const RealSpace& spc, const Enclosure& set)
    : _location(location), _events(), _space(spc.variable_names()), _set(set),
      _variables(catenate(List<EnclosureVariableType>(set.dimension(),INITIAL),List<EnclosureVariableType>(set.number_of_parameters()-set.dimension(),UNKNOWN)))
{
}


HybridEnclosure* HybridEnclosure::clone() const {
    return new HybridEnclosure(*this);
}

const RealSpace
HybridEnclosure::space() const
{
    return RealSpace(this->_space);
}

List<DiscreteEvent> const&
HybridEnclosure::previous_events() const
{
    return this->_events;
}

IntervalFunctionModelFactoryInterface const&
HybridEnclosure::function_factory() const
{
    return this->_set.function_factory();
}

ValidatedVectorFunctionModel const&
HybridEnclosure::space_function() const
{
    return this->_set.space_function();
}

ValidatedScalarFunctionModel const&
HybridEnclosure::time_function() const
{
    return this->_set.time_function();
}

ValidatedScalarFunctionModel const&
HybridEnclosure::dwell_time_function() const
{
    return this->_set.dwell_time_function();
}

Void HybridEnclosure::set_time_function(const ValidatedScalarFunctionModel& time_function)
{
    ARIADNE_NOT_IMPLEMENTED;
    ARIADNE_ASSERT_MSG(Ariadne::subset(this->parameter_domain(),time_function.domain()),
                       "Domain of "<<time_function<<" does not contain parameter domain "<<this->parameter_domain());
    ARIADNE_ASSERT_MSG(this->parameter_domain()==time_function.domain(),
                       "Domain of "<<time_function<<" does not equal parameter domain "<<this->parameter_domain());
    //this->_set._time_function=time_function;
}

UpperBox
HybridEnclosure::space_bounding_box() const
{
    ARIADNE_LOG(8,"space_codomain="<<this->space_function().codomain()<<" space_range="<<apply(this->space_function(),this->_set.reduced_domain())<<"\n");
    //return this->space_function()(this->_set.reduced_domain());
    return make_exact_box(this->_set.bounding_box());
}

UpperInterval
HybridEnclosure::time_range() const
{
    ARIADNE_LOG(8,"time_codomain="<<this->time_function().codomain()<<" time_range="<<apply(this->time_function(),this->_set.reduced_domain())<<"\n");
    //return this->time_function().codomain();
    return apply(this->time_function(),this->_set.reduced_domain());
}

UpperInterval
HybridEnclosure::dwell_time_range() const
{
    ARIADNE_LOG(8,"dwell_time_codomain="<<this->dwell_time_function().codomain()<<
                  " dwell_time_range="<<apply(this->dwell_time_function(),this->_set.reduced_domain())<<"\n");
    return apply(this->dwell_time_function(),this->_set.reduced_domain());
}

Nat
HybridEnclosure::number_of_constraints() const
{
    return this->_set.number_of_constraints();
}

Nat
HybridEnclosure::number_of_parameters() const
{
    return this->_set.number_of_parameters();
}

const ExactBox
HybridEnclosure::parameter_domain() const
{
    return this->_set.domain();
}

HybridUpperBox
HybridEnclosure::bounding_box() const
{
    return HybridUpperBox(this->_location,this->_space,make_exact_box(this->space_bounding_box()));
}


Void HybridEnclosure::new_parameter(ExactInterval ivl, EnclosureVariableType vt)
{
    this->_set.new_parameter(ivl);
    this->_variables.append(vt);
}

Void HybridEnclosure::new_variable(ExactInterval ivl, EnclosureVariableType vt)
{
    this->_set.new_variable(ivl);
    this->_variables.append(vt);
}

Void HybridEnclosure::new_invariant(DiscreteEvent event, ValidatedScalarFunction constraint_function) {
    this->_set.new_negative_state_constraint(constraint_function);
}


Void HybridEnclosure::new_activation(DiscreteEvent event, ValidatedScalarFunction constraint_function) {
    this->_set.new_positive_state_constraint(constraint_function);
}

Void HybridEnclosure::new_guard(DiscreteEvent event, ValidatedScalarFunction constraint_function) {
    this->_set.new_zero_state_constraint(constraint_function);
}

Void HybridEnclosure::new_parameter_constraint(DiscreteEvent event, ValidatedConstraint constraint) {
    this->_set.new_parameter_constraint(constraint);
}

Void HybridEnclosure::new_state_constraint(DiscreteEvent event, ValidatedConstraint constraint) {
    this->_set.new_state_constraint(constraint);
}


Void HybridEnclosure::new_constraint(DiscreteEvent event, ValidatedConstraint constraint) {
    this->new_state_constraint(event,constraint);
}


Void HybridEnclosure::apply_reset(DiscreteEvent event, DiscreteLocation target, RealSpace space, const ValidatedVectorFunction& map)
{
    ARIADNE_ASSERT_MSG(map.argument_size()==this->dimension(),"*this="<<*this<<", event="<<event<<", space="<<space<<", target="<<target<<", map="<<map);
    this->_events.append(event);
    this->_location=target;
    this->_space=space.variable_names();
    this->_set.apply_map(map);
}

Void HybridEnclosure::apply_fixed_evolve_step(const ValidatedVectorFunctionModel& phi, const ExactFloat& elps)
{
    this->_set.apply_fixed_evolve_step(phi,elps);
}

Void HybridEnclosure::apply_spacetime_evolve_step(const ValidatedVectorFunctionModel& phi, const ValidatedScalarFunctionModel& elps)
{
    this->_set.apply_spacetime_evolve_step(phi,elps);
}

Void HybridEnclosure::apply_spacetime_reach_step(const ValidatedVectorFunctionModel& phi, const ValidatedScalarFunctionModel& elps)
{
    this->_set.apply_spacetime_reach_step(phi,elps);
}


Void HybridEnclosure::apply_evolve_step(const ValidatedVectorFunctionModel& phi, const ValidatedScalarFunctionModel& elps)
{
    this->_set.apply_parameter_evolve_step(phi,elps);
}

Void HybridEnclosure::apply_finishing_evolve_step(const ValidatedVectorFunctionModel& phi, const ValidatedScalarFunctionModel& omega)
{
    this->_set.apply_finishing_parameter_evolve_step(phi,omega);
}


Void HybridEnclosure::apply_reach_step(const ValidatedVectorFunctionModel& phi, const ValidatedScalarFunctionModel& elps)
{
    this->_set.apply_parameter_reach_step(phi,elps);
}

Void HybridEnclosure::apply_full_reach_step(const ValidatedVectorFunctionModel& phi)
{
    this->_set.apply_full_reach_step(phi);
}



Void HybridEnclosure::bound_time(Real tmax) {
    if(this->time_range().upper()>ValidatedNumber(tmax).lower()) {
        this->_set.new_negative_parameter_constraint(this->time_function()-ValidatedNumber(tmax));
    }
}

Void HybridEnclosure::bound_time(ValidatedScalarFunction tmax) {
    this->_set.new_negative_parameter_constraint(this->time_function()-this->_set.function_factory().create(this->_set.domain(),tmax));
}

Void HybridEnclosure::set_time(Real time)
{
    this->_set.new_zero_parameter_constraint(this->time_function()-ValidatedNumber(time));
}

Void HybridEnclosure::set_time(ValidatedScalarFunction time)
{
    this->_set.new_zero_parameter_constraint(this->time_function()-this->function_factory().create(this->_set.domain(),time));
}


Void HybridEnclosure::set_maximum_time(DiscreteEvent event, Float final_time)
{
    this->_set.new_negative_parameter_constraint(this->time_function()-ValidatedNumber(final_time)); // Deprecated
}

Void HybridEnclosure::new_time_step_bound(DiscreteEvent event, ValidatedScalarFunction constraint) {
    ARIADNE_NOT_IMPLEMENTED; // Deprecated
}

Void HybridEnclosure::set_step_time(ExactFloat time)
{
    ARIADNE_NOT_IMPLEMENTED; // Deprecated
}



const DiscreteLocation& HybridEnclosure::location() const {
    return this->_location;
}

/*
ValidatedConstrainedImageSet HybridEnclosure::continuous_set() const {
    return ValidatedConstrainedImageSet(this->_set.domain(),this->space_function(),this->_constraints);
}
*/


const Enclosure&
HybridEnclosure::continuous_set() const {
    return this->_set;
}


Nat HybridEnclosure::dimension() const {
    return this->space_function().result_size();
}

Tribool HybridEnclosure::empty() const {
    return this->_set.empty();
}

Tribool HybridEnclosure::inside(const HybridBox& hbx) const {
    if(this->_location==hbx.location()) { return this->continuous_set().inside(hbx.continuous_set()); }
    else { return this->continuous_set().empty(); }
}

Tribool HybridEnclosure::separated(const HybridBox& hbx) const {
    if(this->_location==hbx.location()) { return this->continuous_set().separated(hbx.continuous_set()); }
    else { return true; }
}

Tribool HybridEnclosure::satisfies(EffectiveConstraint c) const
{
    return this->continuous_set().satisfies(c);
}


List<ValidatedConstraint> HybridEnclosure::constraints() const {
    return this->continuous_set().constraints();
}

Void
HybridEnclosure::restrict(const ExactBox& subdomain)
{
    this->_set.restrict(subdomain);
}

Void
HybridEnclosure::recondition()
{
    this->uniform_error_recondition();
    this->kuhn_recondition();
}

Void
HybridEnclosure::uniform_error_recondition()
{
    Nat old_number_of_parameters = this->number_of_parameters();
    this->_set.uniform_error_recondition();
    ExactIntervalVector new_variables = project(this->parameter_domain(),range(old_number_of_parameters,this->number_of_parameters()));
    this->_variables.concatenate(List<EnclosureVariableType>(new_variables.size(),ERROR));
    this->_check();
}


Void
HybridEnclosure::kuhn_recondition()
{
    this->_set.kuhn_recondition();
    this->_check();
}


List<HybridEnclosure>
HybridEnclosure::split() const
{
    List<ExactIntervalVector> subdomains = this->continuous_set().splitting_subdomains_zeroth_order();
    List<HybridEnclosure> result(subdomains.size(),*this);
    for(Nat i=0; i!=result.size(); ++i) {
        result[i].restrict(subdomains[i]);
    }
    return result;
}


Void check_subset(const ExactIntervalVector& dom1, const ExactIntervalVector& dom2, const char* msg)
{
    if(dom1.size()!=dom2.size()) {
        ARIADNE_FAIL_MSG(msg<<" size("<<dom1<<")!=size("<<dom2<<")");
    } else if(!subset(dom1,dom2)) {
        ARIADNE_FAIL_MSG(msg<<" !subset("<<dom1<<","<<dom2<<")");
    }
}

Void
HybridEnclosure::_check() const
{
    //this->_set.check();
    const ExactIntervalVector& reduced_domain = this->_set.reduced_domain();
    check_subset(reduced_domain,this->_set.domain(),"domain");
    check_subset(reduced_domain,this->space_function().domain(),"function domain");
    for(Nat i=0; i!=this->_set.constraints().size(); ++i) {
        check_subset(reduced_domain,this->_set.constraint(i).function().domain(),"constraint");
    }
    check_subset(reduced_domain,this->time_function().domain(),"time");
}

Void HybridEnclosure::draw(CanvasInterface& canvas, const Set<DiscreteLocation>& locations, const Variables2d& axes) const
{
    Projection2d proj=Ariadne::projection(this->space(),axes);
    this->continuous_set().draw(canvas,proj);
}

OutputStream& HybridEnclosure::write(OutputStream& os) const
{
    return os << "HybridEnclosure"
              << "( variables = " << variable_names(this->_variables)
              << ", events=" << this->_events
              << ", location=" << this->_location
              << ", space=" << this->space()
              << ", range=" << apply(this->space_function(),this->_set.domain())
              << ", domain=" << this->_set.domain()
              << ", subdomain=" << this->_set.reduced_domain()
              << ", empty=" << Ariadne::empty(this->_set.reduced_domain())
              << ", state=" << (this->_set.space_function())
              << ", constraints=" << (this->_set.constraints())
              << ", time="<< (this->_set.time_function())
              << ")";
}

OutputStream& HybridEnclosure::print(OutputStream& os) const
{
    return os << "HybridEnclosure"
              << "( events=" << this->_events
              << ", location=" << this->_location
              << ", range=" << apply(this->space_function(),this->_set.domain())
              << ", domain=" << this->_set.domain()
              << ", subdomain=" << this->_set.reduced_domain()
              << ", subdomain=" << this->_set.reduced_domain()
              << ", #constraints=" << this->_set.constraints().size()
              << ", time_range="<<this->time_range() << ")";
}

OutputStream& HybridEnclosure::repr(OutputStream& os) const
{
    return os << "HybridEnclosure"
              << "( " << this->_location
              << ", " << this->_set.domain()
              << ", " << this->_set.reduced_domain()
              << ", " << this->_set.space_function()
              << ", " << this->_set.constraints()
              << ", " << this->time_function() << ")";
}


Void HybridEnclosure::adjoin_outer_approximation_to(HybridGridTreeSet& hgts, Int depth) const {
    const Enclosure& set = this->continuous_set();
    GridTreeSet& paving = hgts[this->location()];
    set.adjoin_outer_approximation_to(paving,depth);
}

HybridGridTreeSet outer_approximation(const ListSet<HybridEnclosure>& hls, const HybridGrid& g, Int depth) {
    HybridGridTreeSet result(g);
    for(ListSet<HybridEnclosure>::ConstIterator iter=hls.begin(); iter!=hls.end(); ++iter) {
        result[iter->location()].adjoin_outer_approximation(iter->continuous_set(),depth);
    }
    return result;
}

Void
draw(FigureInterface& figure, const ListSet<HybridEnclosure>& hels) {
    for(ListSet<HybridEnclosure>::ConstIterator iter=hels.begin(); iter!=hels.end(); ++iter) {
        draw(figure,iter->continuous_set());
    }
}

ListSet<HybridEnclosure::ContinuousStateSetType>
ListSet<HybridEnclosure>::operator[](const DiscreteLocation& loc) const
{
    ListSet<HybridEnclosure::ContinuousStateSetType> result;
    for(ConstIterator iter=this->begin(); iter!=this->end(); ++iter) {
        if(iter->location()==loc) {
            result.adjoin(iter->continuous_set());
        }
    }
    return result;
}




} // namespace Ariadne