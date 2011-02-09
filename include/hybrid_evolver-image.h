/***************************************************************************
 *            hybrid_evolver-image.h
 *
 *  Copyright  2007-10  Alberto Casagrande, Pieter Collins, Luca Geretti
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

/*! \file hybrid_evolver-image.h
 *  \brief Evolver for hybrid systems using image sets.
 */

#ifndef ARIADNE_HYBRID_EVOLVER_IMAGE_H
#define ARIADNE_HYBRID_EVOLVER_IMAGE_H

#include <string>
#include <vector>
#include <list>
#include <iostream>

#include <boost/smart_ptr.hpp>

#include "tuple.h"

#include "hybrid_set.h"

#include "hybrid_automaton.h"
#include "evolver_interface.h"
#include "evolver_base.h"
#include "evolution_parameters.h"

#include "logging.h"

namespace Ariadne {

template<class Sys, class BS> class Evolver;

class VectorTaylorFunction;
class TaylorSet;
class HybridAutomaton;
template<class ES> class Orbit;

class EvolutionParameters;
class TaylorModel;
template<class MDL> class CalculusInterface;

class HybridTime;

class DiscreteEvent;

/*! \brief A class for computing the evolution of a hybrid system.
 *
 * The actual evolution steps are performed by the HybridEvolver class.
 */
class ImageSetHybridEvolver
    : public EvolverBase<HybridAutomaton,HybridTaylorSet>
    , public Loggable
{
    typedef ScalarFunction ScalarFunctionType;
    typedef VectorFunction VectorFunctionType;
    typedef Vector<Interval> BoxType;
    typedef VectorTaylorFunction FunctionModelType;
    typedef VectorTaylorFunction MapModelType;
    typedef VectorTaylorFunction FlowModelType;
    typedef ScalarTaylorFunction ConstraintModelType;
    typedef TaylorModel TimeModelType;
    typedef TaylorSet FlowSetModelType;
    typedef TaylorSet SetModelType;
    typedef TaylorSet TimedSetModelType;
  public:
    typedef ContinuousEvolutionParameters EvolutionParametersType;
    typedef HybridAutomaton::TimeType TimeType;
    typedef int IntegerType;
    typedef Float RealType;
    typedef std::vector<DiscreteEvent> EventListType;
    typedef HybridAutomaton SystemType;
    typedef TaylorSet ContinuousEnclosureType;
    typedef HybridBasicSet<TaylorSet> HybridEnclosureType;
    typedef HybridEnclosureType EnclosureType;
    typedef Orbit<EnclosureType> OrbitType;
    typedef ListSet<EnclosureType> EnclosureListType;
    typedef Float ContinuousTimeType;
  public:

    //! \brief Default constructor.
    ImageSetHybridEvolver();

    //! \brief Construct from parameters using a default integrator.
    ImageSetHybridEvolver(const EvolutionParametersType& parameters);

    /*! \brief Make a dynamically-allocated copy. */
    ImageSetHybridEvolver* clone() const { return new ImageSetHybridEvolver(*this); }

    //@{
    //! \name Parameters controlling the evolution.

    //! \brief A reference to the parameters controlling the evolution.
    EvolutionParametersType& parameters() { return *this->_parameters; }
	//! \brief A constant reference to the parameters controlling the evolution.
    const EvolutionParametersType& parameters() const { return *this->_parameters; }

    //@}

    //@{
    //! \name Evolution using abstract sets.

    //! \brief Compute an approximation to the orbit set using the given semantics.
    Orbit<EnclosureType> orbit(const SystemType& system, const EnclosureType& initial_set, const TimeType& time, Semantics semantics=UPPER_SEMANTICS) const;

    //! \brief Compute an approximation to the orbit set for upper semantics, with continuous evolution only.
    Orbit<EnclosureType> upper_orbit_continuous(const SystemType& system, const EnclosureType& initial_set, const TimeType& time) const;

    //! \brief Compute an approximation to the evolution set using the given semantics.
    EnclosureListType evolve(const SystemType& system, const EnclosureType& initial_set, const TimeType& time, Semantics semantics=UPPER_SEMANTICS) const {
        EnclosureListType final; EnclosureListType reachable; EnclosureListType intermediate;
        this->_evolution(final,reachable,intermediate,system,initial_set,time,semantics,false);
        return final; }

    //! \brief Compute an approximation to the evolution set under the given semantics.
    EnclosureListType reach(const SystemType& system, const EnclosureType& initial_set, const TimeType& time, Semantics semantics=UPPER_SEMANTICS) const {
        EnclosureListType final; EnclosureListType reachable; EnclosureListType intermediate;
        this->_evolution(final,reachable,intermediate,system,initial_set,time,semantics,true);
        reachable.adjoin(final);
        return reachable; }

    //! \brief Compute an approximation to the evolution set under the given semantics, returning the reached and final sets.
    std::pair<EnclosureListType,EnclosureListType> reach_evolve(const SystemType& system, const EnclosureType& initial_set, const TimeType& time, Semantics semantics=LOWER_SEMANTICS) const {
        EnclosureListType final; EnclosureListType reachable; EnclosureListType intermediate;
        this->_evolution(final,reachable,intermediate,system,initial_set,time,semantics,true);
        reachable.adjoin(final);
        return make_pair<EnclosureListType,EnclosureListType>(reachable,final); }

    //! \brief Compute an approximation to the evolution set under the given semantics, returning the reached and final sets, and the information
    //! on having disproved.
    tuple<EnclosureListType,EnclosureListType,DisproveData> lower_chain_reach_evolve_disprove(const SystemType& system, const EnclosureType& initial_set,
																					  const TimeType& time, const HybridBoxes& disprove_bounds,
																					  const bool& skip_if_disproved) const {
        EnclosureListType final; EnclosureListType reachable; EnclosureListType intermediate;
        DisproveData falsInfo = this->_lower_evolution_disprove(final,reachable,intermediate,system,initial_set,time,
														   true,disprove_bounds,skip_if_disproved);
        reachable.adjoin(final);
        return make_tuple<EnclosureListType,EnclosureListType,DisproveData>(reachable,final,falsInfo); }

  protected:
    virtual void _evolution(EnclosureListType& final, EnclosureListType& reachable, EnclosureListType& intermediate,
                            const SystemType& system, const EnclosureType& initial, const TimeType& time,
                            Semantics semantics, bool reach) const;

    virtual DisproveData _lower_evolution_disprove(EnclosureListType& final, EnclosureListType& reachable,
										   EnclosureListType& intermediate, const SystemType& system,
										   const EnclosureType& initial, const TimeType& time, bool reach,
										   const HybridBoxes& disprove_bounds, bool skip_if_disproved) const;

    virtual void _upper_evolution_continuous(EnclosureListType& final, EnclosureListType& reachable, EnclosureListType& intermediate,
                            const SystemType& system, const EnclosureType& initial, const TimeType& time, bool reach) const;

    typedef tuple<DiscreteState, EventListType, SetModelType, TimeModelType> HybridTimedSetType;
    virtual void _evolution_step(std::list< HybridTimedSetType >& working_sets,
                                  EnclosureListType& final, EnclosureListType& reachable, EnclosureListType& intermediate,
                                  const SystemType& system, const HybridTimedSetType& current_set, const TimeType& time,
                                  Semantics semantics, bool reach) const;

    virtual void _upper_evolution_continuous_step(std::list< HybridTimedSetType >& working_sets,
                                  				  EnclosureListType& final, EnclosureListType& reachable, EnclosureListType& intermediate,
                                  				  const SystemType& system, const HybridTimedSetType& current_set, const TimeType& time, 
												  bool reach) const;

    virtual DisproveData _lower_evolution_disprove_step(std::list< HybridTimedSetType >& working_sets,
												EnclosureListType& final, EnclosureListType& reachable, EnclosureListType& intermediate,
												const SystemType& system, const HybridTimedSetType& current_set, const TimeType& time,
												bool reach, const HybridBoxes& disprove_bounds) const;

  protected:
    TimeModelType crossing_time(VectorFunction guard, const FlowSetModelType& flow_set) const;

    Interval normal_derivative(VectorFunction guard, const FlowSetModelType& flow_set, const TimeModelType& crossing_time) const;

    void compute_initially_active_events(std::map<DiscreteEvent,tribool>&,
                                         const std::map<DiscreteEvent,VectorFunction>&,
                                         const ContinuousEnclosureType&) const;

    void compute_flow_model(FlowSetModelType&, BoxType&, Float&, VectorFunction, 
                            const SetModelType&, const TimeModelType&, Float) const;

    void compute_crossing_time_and_direction(TimeModelType&, Interval&,
                                             VectorFunction guard, const FlowSetModelType& flow_set) const;

    void compute_blocking_events(std::map<DiscreteEvent,TimeModelType>&, std::set<DiscreteEvent>&,
                                 const std::map<DiscreteEvent,VectorFunction>& guards,
                                 const FlowSetModelType& flow_set_model) const;

    void compute_blocking_time(std::set<DiscreteEvent>&, TimeModelType&,
                               const std::map<DiscreteEvent,TimeModelType>&) const;

    void compute_activation_events(std::map<DiscreteEvent,tuple<tribool,TimeModelType,tribool> >&,
                                  const std::map<DiscreteEvent,VectorFunction>& activations,
                                  const FlowSetModelType& flow_set_model) const;

    void compute_activation_times(std::map<DiscreteEvent,tuple<TimeModelType,TimeModelType> >&,
                                  const std::map<DiscreteEvent,VectorFunction>& activations,
                                  const FlowSetModelType& flow_set_model,
                                  const TimeModelType& blocking_time_model,
                                  const Semantics sematics) const;

  private:

    bool _check_bounds(const uint& numDivisions,
					   const TaylorSet& reachable_set,
					   const Box& disprove_bounds) const;

  protected:
    // Special events
    static const DiscreteEvent starting_event;
    static const DiscreteEvent finishing_event;
    static const DiscreteEvent blocking_event;

 private:
    boost::shared_ptr< EvolutionParametersType > _parameters;
    boost::shared_ptr< CalculusInterface<TaylorModel> > _toolbox;
};



} // namespace Ariadne

#endif // ARIADNE_HYBRID_EVOLVER_IMAGE_H