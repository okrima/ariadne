/***************************************************************************
 *            reachability_analyser.h
 *
 *  Copyright  2006-10  Alberto Casagrande, Pieter Collins, Luca Geretti
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
 
/*! \file reachability_analyser.h
 *  \brief Methods for computing abstract reachable sets.
 */

#ifndef ARIADNE_REACHABILITY_ANALYSER_H
#define ARIADNE_REACHABILITY_ANALYSER_H

#include <boost/smart_ptr.hpp>
#include <cstdarg>
#include <config.h>

#include "hybrid_set_interface.h"
#include "evolver_interface.h"
#include "discretiser_interface.h"
#include "reachability_analyser_interface.h"

#include "orbit.h"
#include "grid_set.h"
#include "hybrid_set.h"
#include "graphics.h"

#include "discretiser.h"
#include "hybrid_evolver.h"

#include "logging.h"
#include "parametric.h"

namespace Ariadne {
 
template<class ES> class Orbit;

class DiscreteState;

template<class BS> class HybridBasicSet;
typedef HybridBasicSet<Box> HybridBox;
typedef std::map<DiscreteState,Vector<Float> > HybridFloatVector;
typedef std::map<RealConstant,int,ConstantComparator<Real> > RealConstantIntMap;

class HybridGrid;
class HybridGridCell;
class HybridGridTreeSet;

template<class ES> class HybridListSet;
template<class ES> class HybridDiscretiser;

/** Provides a convenient structure for systems and their relative info for verification. */
struct SystemVerificationInfo
{
public:
	HybridAutomaton& _system;
	HybridImageSet& _initial_set;
	HybridBoxes& _domain;
	HybridBoxes _safe_region;
	std::vector<uint> _projection;

	HybridAutomaton& getSystem() { return _system; }
	HybridImageSet& getInitialSet() { return _initial_set; }
	HybridBoxes& getDomain() { return _domain; }
	HybridBoxes& getSafeRegion() { return _safe_region; }
	std::vector<uint>& getProjection() { return _projection; }

	SystemVerificationInfo(HybridAutomaton& system,
						   HybridImageSet& initial_set,
						   HybridBoxes& domain,
						   HybridBoxes& safe_region);

	SystemVerificationInfo(HybridAutomaton& system,
						   HybridImageSet& initial_set,
						   HybridBoxes& domain,
						   std::vector<uint>& projection);

	SystemVerificationInfo(HybridAutomaton& system,
						   HybridImageSet& initial_set,
						   HybridBoxes& domain,
						   HybridBoxes& safe_region,
						   std::vector<uint>& projection);

	virtual std::ostream& write(std::ostream&) const;

private:

	/** \brief Checks fields consistency */
	void _check_fields() const;

};

inline std::ostream& operator<<(std::ostream& os, const SystemVerificationInfo& verInfo) {
    return verInfo.write(os); }

/*! \brief A class for performing reachability analysis on a hybrid system-
 */
class HybridReachabilityAnalyser
    : public Loggable
{
  private:
    boost::shared_ptr< DiscreteEvolutionParameters > _parameters;
	boost::shared_ptr< DiscreteEvolutionStatistics > _statistics;
    boost::shared_ptr< HybridDiscretiser<HybridEvolver::ContinuousEnclosureType> > _discretiser;
  public:
    typedef DiscreteEvolutionParameters EvolutionParametersType;
	typedef DiscreteEvolutionStatistics	EvolutionStatisticsType;
    typedef HybridAutomaton SystemType;
    typedef SystemType::StateSpaceType StateSpaceType;
    typedef SystemType::TimeType TimeType;
    typedef HybridGridTreeSet SetApproximationType;
    typedef HybridEvolver::EnclosureType EnclosureType;
    typedef HybridEvolver::ContinuousEnclosureType ContinuousEnclosureType;
    typedef std::pair<SetApproximationType,bool> (HybridReachabilityAnalyser::*UpperChainReachFuncPtr)(const SystemType&, const HybridImageSet&) const;
  public:
    //@{
    //! \name Constructors and destructors
    /*! \brief Virtual destructor */
    virtual ~HybridReachabilityAnalyser();

    /*! \brief Construct from a method for evolving basic sets. */
    HybridReachabilityAnalyser(const HybridDiscretiser<HybridEvolver::ContinuousEnclosureType>& discretiser);

    /*! \brief Construct from evolution parameters and a method for evolving basic sets. */
    template<class HybridEnclosureType>
    HybridReachabilityAnalyser(const EvolutionParametersType& parameters,
                               const EvolverInterface<HybridAutomaton,HybridEnclosureType>& evolver);

    template<class HybridEnclosureType>
    HybridReachabilityAnalyser(const EvolverInterface<HybridAutomaton,HybridEnclosureType>& evolver);

    /*! \brief Make a dynamically-allocated copy. */
    virtual HybridReachabilityAnalyser* clone() const { return new HybridReachabilityAnalyser(*this); }
    //@}
  
    //@{ 
    //! \name Methods to set and get the parameters controlling the accuracy
    /*! \brief The parameters controlling the accuracy. */
    const EvolutionParametersType& parameters() const { return *this->_parameters; }
    /*! \brief A reference to the parameters controlling the accuracy. */
    EvolutionParametersType& parameters() { return *this->_parameters; }
    //@}

    //@{ 
    //! \name Methods to set and get the statistics related to the analyses
    /*! \brief The statistics stemming from the analyses. */
    const EvolutionStatisticsType& statistics() const { return *this->_statistics; }
    /*! \brief A reference to the statistics stemming from the analyses. */
    EvolutionStatisticsType& statistics() { return *this->_statistics; }
    //@}
  
    //@{
    //! \name Evaluation of systems on abstract sets
    /*! \brief Compute a lower-approximation to the set obtained by evolving \a system for \a time starting in \a initial_set. */
    virtual SetApproximationType lower_evolve(const SystemType& system,
                                              const HybridImageSet& initial_set, 
                                              const TimeType& time) const;
  
    /*! \brief Compute a lower-approximation to the reachable set of \a system starting in \a initial_set up to \a time (discrete part only). */
    virtual SetApproximationType
    lower_reach(const SystemType& system, 
                const HybridImageSet& initial_set, 
                const TimeType& time) const;
  
    /*! \brief Compute a lower-approximation to the reachable and evolved sets of \a system starting in \a initial_set up to \a time. */
    virtual std::pair<SetApproximationType,SetApproximationType>
    lower_reach_evolve(const SystemType& system, 
                       const HybridImageSet& initial_set, 
                       const TimeType& time) const;
  
    /*! \brief Compute an approximation to the set obtained by iterating \a time times \a system starting in \a initial_set. */
    virtual SetApproximationType upper_evolve(const SystemType& system,
                                              const HybridImageSet& initial_set, 
                                              const TimeType& time) const;
  
    /*! \brief Compute an approximation to the reachable set of \a system starting in \a initial_set iterating at most \a time times. */
    virtual SetApproximationType upper_reach(const SystemType& system,
                                             const HybridImageSet& initial_set, 
                                             const TimeType& timeType) const;
  
    /*! \brief Compute an approximation to the reachable and evolved sets of \a system starting in \a initial_set iterating at most \a time times. */
    virtual std::pair<SetApproximationType,SetApproximationType>
    upper_reach_evolve(const SystemType& system, 
                       const HybridImageSet& initial_set, 
                       const TimeType& time) const;
  
    /*! \brief Compute an outer-approximation to the chain-reachable set of \a system starting in \a initial_set. */
    virtual SetApproximationType chain_reach(const SystemType& system,
                                             const HybridImageSet& initial_set) const;

    /*! \brief Compute an outer-approximation to the chain-reachable set of \a system starting in \a initial_set and remaining in \a bounding_domain. \deprecated */
    virtual SetApproximationType chain_reach(const SystemType& system,
                                             const HybridImageSet& initial_set,
                                             const HybridBoxes& bounding_domain) const;

    /*! \brief Compute an outer-approximation to the chain-reachable set of \a system starting in \a initial_set, with
     * upper semantics; the method performs discretisation before transitions, then checks activations on the discretised cells.
     * \return The reach set and a flag notifying if the result is valid (i.e. it has not been restricted due to
     * the bounding domain and it is safe). */
    virtual std::pair<SetApproximationType,bool> upper_chain_reach(SystemType& system,
																   const HybridImageSet& initial_set) const;

    /*! \brief Compute an outer-approximation to the chain-reachable set of \a system starting in \a initial_set, with
         * lower semantics; the method performs periodical discretisations and checks the new reached region for inclusion
         * The resulting set is a subset of the outer-approximation of the whole evolution set.
         * \return The reach set and the falsification information. */
    virtual std::pair<SetApproximationType,DisproveData> lower_chain_reach(SystemType& system,
																	   const HybridImageSet& initial_set) const;
  
    /*! \brief Compute an outer-approximation to the viability kernel of \a system within \a bounding_set. */
    virtual SetApproximationType viable(const SystemType& system,
                                        const HybridImageSet& bounding_set) const;
  
    /*! \brief Attempt to verify that the reachable set of \a system starting in \a initial_set remains in \a safe_box. */
    virtual tribool verify(SystemType& system,
                           const HybridImageSet& initial_set);

	/*! \brief Attempt to verify that the reachable set of a system starting in an initial_set remains in a safe region inside a domain.
	 * \details This is done in an iterative way by tuning the evolution/analysis parameters. The \a verInfo contains all the information
	 * necessary for verification.
	 */
	tribool verify_iterative(SystemVerificationInfo& verInfo);

	/*! \brief Performs iterative verification where \a parameter is substituted into the system.
	 */
	tribool verify_iterative(SystemVerificationInfo& verInfo, const RealConstant& parameter);

	/*! \brief Performs iterative verification where the singleton \a value is substituted into the system for the given \a parameter.
	 */
	tribool verify_iterative(SystemVerificationInfo& verInfo, const RealConstant& parameter, const Float& value);

	/*! \brief Compute an underapproximation of the safety/unsafety intervals of \a parameter (defined as an interval) for the automaton
		\a system starting in \a initial_set, where the safe region is \a safe inside \a domain. 
        \details The procedure uses the bisection method. The parameter is assumed as having separable safe and unsafe intervals in its range.
        The tolerance in [0 1] is a percentage of the parameter interval width and is used to provide a termination condition for the
		bisection search.
        \return The intervals of safety and unsafety. */
	std::pair<Interval,Interval> parametric_verification_1d_bisection(SystemVerificationInfo& verInfo,
										 							  const RealConstant& parameter,
										 							  const Float& tolerance);

	/**
	 * \brief Performs a parametric verification on two parameters \a xParam, \a yParam,
	 * discretizing with \a numPointsPerAxis points for each axis.
	 * \details The procedure uses the bisection method. Saves the results in a file called "<systemName>-<xName>-<yName>" and
	 * generates a "<systemName>-<xName>-<yName>.png" plot, where <systemName> is the name of the system,
	 * <xName> is the name of xParam and <yName> is the name of yParam.
	 */
	Parametric2DBisectionResults parametric_verification_2d_bisection(SystemVerificationInfo& verInfo,
								 const RealConstant& xParam,
								 const RealConstant& yParam,
								 const Float& tolerance,
								 const unsigned& numPointsPerAxis);

	/**
	 * \brief Performs a parametric verification on a set of two parameters \a params, by using bisection.
	 */
	Parametric2DBisectionResults parametric_verification_2d_bisection(SystemVerificationInfo& verInfo,
												   const RealConstantSet& params,
												   const Float& tolerance,
												   const unsigned& numPointsPerAxis);

	/**
	 * \brief Performs a parametric verification on a set of parameters \a params, by partitioning the parameters space.
	 * \details The \a minPartitioningRatio variable determines the minimum granularity for partitioning any parameter (expressed as a
	 * percentage in respect to the range of a given parameter). The values in \a params are substituted into the system, the latter
	 * being restored to its initial conditions by the end of the method.
	 */
	ParametricPartitioningOutcomeList parametric_verification_partitioning(SystemVerificationInfo& verInfo,
													   const RealConstantSet& params,
													   const Float& minPartitioningRatio);

	/**
	 * \brief Performs dominance checking.
	 * \details Verifies if the \a dominating system dominates the \a dominated system.
	 */
	tribool dominance(SystemVerificationInfo& dominating,
					  SystemVerificationInfo& dominated);

	/**
	 * \brief Performs dominance checking with \a parameter substituted into the \a dominating system.
	 * \details Verifies if the \a dominating system dominates the \a dominated system.
	 */
	tribool dominance(SystemVerificationInfo& dominating,
					  SystemVerificationInfo& dominated,
					  const RealConstant& parameter);

	/**
	 * \brief Performs dominance checking with \a parameter substituted into the \a dominating system with a value of \a value.
	 * \details Verifies if the \a dominating system dominates the \a dominated system.
	 */
	tribool dominance(SystemVerificationInfo& dominating,
					  SystemVerificationInfo& dominated,
					  const RealConstant& parameter,
					  const Float& value);

	/*! \brief Compute an underapproximation of the dominating/non-dominating intervals of \a parameter for the dominance problem.
        \details The parameter is varied on the \a dominating system alone. The procedure uses the bisection method. The parameter is assumed as having separable dominating and non-dominating intervals in its range.
        The tolerance in [0 1] is a percentage of the parameter interval width and is used to provide a termination condition for the
		bisection search.
        \return The intervals of safety and unsafety. */
	std::pair<Interval,Interval> parametric_dominance_1d_bisection(SystemVerificationInfo& dominating,
			  	  	  	  	  	  	  	  	  	  	  	  	  	   SystemVerificationInfo& dominated,
										 							  const RealConstant& parameter,
										 							  const Float& tolerance);
	/**
	 * \brief Performs a parametric dominance checking on two parameters \a xParam, \a yParam,
	 * discretizing with \a numPointsPerAxis points for each axis.
	 * \details The procedure uses the bisection method. Saves the results in a file called "<dominatingName>&<dominatedName>-<xName>-<yName>" and
	 * generates a "<dominatingName>&<dominatedName>-<xName>-<yName>.png" plot, where <systemName> is the name of the system,
	 * <xName> is the name of xParam and <yName> is the name of yParam.
	 */
	Parametric2DBisectionResults parametric_dominance_2d_bisection(SystemVerificationInfo& dominating,
 	  	  	  	  	  	   	   	   	   	   SystemVerificationInfo& dominated,
								 const RealConstant& xParam,
								 const RealConstant& yParam,
								 const Float& tolerance,
								 const unsigned& numPointsPerAxis);

	/**
	 * \brief Performs a parametric dominance checking on a set of two parameters \a params, by using bisection.
	 */
	Parametric2DBisectionResults parametric_dominance_2d_bisection(SystemVerificationInfo& dominating,
 	  	  	  	  	  	   	   	   	   	   SystemVerificationInfo& dominated,
												   const RealConstantSet& params,
												   const Float& tolerance,
												   const unsigned& numPointsPerAxis);

	/**
	 * \brief Performs a parametric dominance checking on a set of parameters \a dominating_params of the \a dominating system, by partitioning
	 * the parameters space.
	 * \details The \a minPartitioningRatio variable determines the minimum granularity for partitioning any parameter (expressed as a
	 * percentage in respect to the range of a given parameter). The values in \a dominating_params are substituted into the \a dominating
	 * system alone, the latter being restored to its initial conditions by the end of the method.
	 */
	ParametricPartitioningOutcomeList parametric_dominance_partitioning(SystemVerificationInfo& dominating,
													 	   SystemVerificationInfo& dominated,
													 	   const RealConstantSet& dominating_params,
													 	   const Float& minPartitioningRatio);

    //@}

  public:

	// Determine if the upper/lower reached regions of verify_iterative (and of all the functions using it) must be saved into figures (false by default)
	bool plot_verify_results;
	// Determine if the reached region (resulting from the chain_reach_grid() procedure) is to be dumped into disk for those locations not involved in
	// the evolution, and reloaded when required again. Only partially implemented at the moment.
	bool chain_reach_dumping;
	// The reduction in the number of logical cores used in multithreading (down from the maximum concurrency of the machine) (zero by default)
	uint free_cores;
  
  public:

    typedef HybridTime T;
    typedef HybridAutomaton Sys;
    typedef HybridListSet<Box> BxLS;
    typedef HybridGrid Gr;
    typedef HybridGridCell GC;
    typedef HybridGridTreeSet GCLS;
    typedef HybridGridTreeSet GTS;
    typedef HybridOpenSetInterface OpSI;
    typedef HybridOvertSetInterface OvSI;
    typedef HybridCompactSetInterface CoSI;
    typedef std::map<RealConstant,int,ConstantComparator<Real> > RealConstantMap;

  private:

    /*! \brief Generates a list of hybrid enclosures from the \a initial_set, depending on the minimum cell size
     * given by the \a grid. */
    list<HybridBasicSet<TaylorSet> > _split_initial_set(const HybridImageSet initial_set,
										   const HybridGrid grid) const;

    /*! \brief Generates the target hybrid enclosures given an enclosure \a encl and a transition \a trans,
     * eventually splitting if the original enclosure has any width larger than the values in \a minTargetCellWidths.
     * Also, removes target enclosures lying outside the \a target_bounding, consequently invalidating the verification by setting
     * \a isValid to false. */
    void _splitAndCreateTargetEnclosures(bool& isValid,
										 std::list<EnclosureType>& initial_enclosures,
    									 const ContinuousEnclosureType& encl,
    									 const Vector<Float>& minTargetCellWidths,
    									 const Box& target_bounding,
    									 const DiscreteTransition& trans_it,
    									 const TaylorCalculus& tc) const;

    void _splitTargetEnclosures(bool& isValid,
							    std::list<EnclosureType>& initial_enclosures,
							    const DiscreteState& target_loc,
							    const ContinuousEnclosureType& target_encl,
							    const Vector<Float>& minTargetCellWidths,
							    const Box& target_bounding) const;

    /*! \brief Gets for each non-singleton constant the factor determining the number of chunks its interval should be split into.
     *
     * \details Splits until the deviation of the derivatives is reasonably low in respect to the deviation calculated at the midpoint. This
     * limit value is expressed as a percentage using \a tolerance.
     *
     * @param system The system to get the accessible constants from.
     * @param targetRatioPerc The derivative widths ratio percentage to reach before termination.
     *
     * @return A split factor for each non-singleton accessible constant of the \a system.
     */
    RealConstantIntMap _getSplitFactorsOfConstants(HybridAutomaton& system, const RealConstantSet& locked_constants,
												  const Float& targetRatioPerc) const;

    /*! \brief Gets the set of all the split intervals from the stored split factors.
     *  \details Orders the list elements by first picking the leftmost subintervals, followed by the rightmost and then
     *  all the remaining from right to left.
     */
    std::list<RealConstantSet> _getSplitConstantsSet() const;

    /*! \brief Gets the best constant among the \a working_constants of the \a system to split, in terms of
     * relative reduction of derivative widths compared to some \a referenceWidths.
     */
    RealConstant _getBestConstantToSplit(SystemType& system, const RealConstantSet& working_constants,
    							         const HybridFloatVector& referenceWidths) const;

    /*! \brief Helper function to get the maximum value of the derivative width ratio \f$ (w-w^r)/w_r \f$, where the \f$ w^r \f$ values
     * are stored in \a referenceWidths and the \f$ w \f$ values are obtained from the \a system.
     */
    Float _getMaxDerivativeWidthRatio(const HybridAutomaton& system, const HybridFloatVector& referenceWidths) const;

    /*! \brief Helper function to get the widths of the derivatives from the \a system */
    HybridFloatVector _getDerivativeWidths(const HybridAutomaton& system) const;

    // Helper functions for operators on lists of sets.
    GTS _upper_reach(const Sys& sys, const GTS& set, const T& time, const int accuracy) const;
    GTS _upper_evolve(const Sys& sys, const GTS& set, const T& time, const int accuracy) const;
    std::pair<GTS,GTS> _upper_reach_evolve(const Sys& sys, const GTS& set, const T& time, const int accuracy) const;
    std::pair<GTS,GTS> _upper_reach_evolve_continuous(const Sys& sys, const list<EnclosureType>& initial_enclosures, const T& time, const int accuracy) const;
    std::pair<SetApproximationType,bool> _upper_chain_reach(SystemType& system, const HybridImageSet& initial_set,UpperChainReachFuncPtr func) const;
    std::pair<SetApproximationType,bool> _upper_chain_reach_forward(const SystemType& system, const HybridImageSet& initial_set) const;

    /*! \brief Checks if a break is to be issued.
     * \details Read the \a flagToCheck for validity: if false, sets \a flagToUpdate to false and returns true ("do break"), otherwise
     * just returns false.
     */
    bool _upper_chain_reach_break(const bool& flagToCheck, bool& flagToUpdate) const;

    /*! \brief Pushes the enclosures from \a reachCells into \a destination.
     * \details Ignores enclosures that lie outside the domain.
     * \return True iff any enclosure has been ignored due to lying outside the domain.
     */
    bool _upper_chain_reach_pushReachCells(const HybridGridTreeSet& reachCells,
    									   const SystemType& system,
    									   std::list<EnclosureType>& destination) const;

    /*! \brief Pushes the enclosures from the \a source enclosure into the \a destination enclosure list, for all \a transitions.
     * \return True iff any enclosure has been ignored due to lying outside the domain.
     */
    bool _upper_chain_reach_pushTargetEnclosures(const std::list<DiscreteTransition>& transitions,
												 const ContinuousEnclosureType& source,
												 const HybridGrid& grid,
												 std::list<EnclosureType>& destination) const;

    /*! \brief Pushes the enclosures from the \a source enclosure into the \a destination enclosure list for a specific transition \a trans.
     * \details Splits the \a source until the enclosure is definitely active for \a trans or the minimum allowed target cell widths \a minTargetCellWidths has been reached.
     * \return True iff any enclosure has been ignored due to lying outside the domain.
     */
    bool _upper_chain_reach_pushTargetEnclosuresOfTransition(const DiscreteTransition& trans,
    														 const ContinuousEnclosureType& source,
    														 const Vector<Float>& minTargetCellWidths,
    														 std::list<EnclosureType>& destination) const;

    /*! \brief Pushes the target enclosure from the \a source enclosure into the \a destination enclosure list for a specific transition \a trans.
     * \return True iff the target enclosure lies outside the domain.
     */
    bool _upper_chain_reach_pushTransitioningEnclosure(const DiscreteTransition& trans,
    												   const ContinuousEnclosureType& source,
    												   const TaylorCalculus& tc,
    												   std::list<EnclosureType>& destination) const;

    /*! \brief Pushes the enclosures from the \a finalCells tree set into the \a destination enclosure list.
     * \return True iff any cell enclosure has been ignored due to lying outside the domain.
     */
    bool _upper_chain_reach_pushFinalCells(const HybridGridTreeSet& finalCells,
    									   std::list<EnclosureType>& destination) const;

    std::pair<SetApproximationType,DisproveData> _lower_chain_reach(const SystemType& system,
    																const HybridImageSet& initial_set) const;

	/*! \brief Attempt to verify that the reachable set of a system starting in an initial set remains in a safe region inside a domain,
		in an iterative way by tuning the evolution/analysis parameters. In addition, the constants in the \a locked_constants set
		are not allowed to be split */
	tribool _verify_iterative(SystemVerificationInfo& verInfo,
							 const RealConstantSet& locked_constants);

	/*! \brief Prove that the automaton \a system starting in \a initial_set definitely remains in the safe region. */
	bool _prove(SystemType& system, const HybridImageSet& initial_set);
	/*! \brief Disprove that the automaton \a system starting in \a initial_set definitely remains in the safe region. */
	bool _disprove(SystemType& system, const HybridImageSet& initial_set);

	/*! \brief Get the hybrid maximum absolute derivatives of \system given a previously computed chain reach statistics. 
		\details ASSUMPTION: the continuous variables are preserved in order and quantity between discrete states. */
	HybridFloatVector _getHybridMaximumAbsoluteDerivatives(const SystemType& system) const;

	/*! \brief Set the lock to grid time of \system given a previously computed chain reach statistics.
		\details The value is taken as the maximum over the times required by any variable on any location to cover a distance equal to
		the domain width of the location, moving at the maximum absolute derivative.
		ASSUMPTION: the continuous variables are preserved in order and quantity between discrete states. */
	void _setLockToGridTime(const SystemType& system) const;

	/*! \brief Get the hybrid maximum integration step size, under the assumption that given the maximum derivatives \a hmad,
		all variables in a step must cover a length greater than a length determined by the \a hgrid. */
	std::map<DiscreteState,Float> _getHybridMaximumStepSize(const HybridFloatVector& hmad, const HybridGrid& hgrid);

	/*! \brief Set the hybrid maximum integration step size, under the assumption that given the maximum derivatives \a hmad,
		all variables in a step must cover a length greater than a length determined by the \a hgrid. The value is equal for all
		locations and corresponds to the largest integration step among the locations. */
	void _setEqualizedHybridMaximumStepSize(const HybridFloatVector& hmad, const HybridGrid& hgrid);

	/*! \brief Set the maximum enclosure cell from the hybrid grid \a hgrid. */
	void _setMaximumEnclosureCell(const HybridGrid& hgrid);

	/*! \brief Get the hybrid grid given the maximum derivative \a hmad and the bounding domain parameter, where the grid is chosen differently for each location.
	 * \details The grid is chosen so that each cell is included into the domain corresponding to its location. */
	HybridGrid _getLooselyConstrainedHybridGrid(const HybridFloatVector& hmad) const;

	/*! \brief Get the hybrid grid given the maximum derivative \a hmad and the bounding domain parameter, where the grid is chosen differently for each location.
	 * \details The grid is chosen so that each cell is included into the domains for all locations. */
	HybridGrid _getHybridGrid(const HybridFloatVector& hmad) const;

	/*! \brief Get the hybrid grid given the maximum derivative \a hmad and the bounding domain parameter, where the grid is chosen equally for all locations.*/
	HybridGrid _getEqualizedHybridGrid(const HybridFloatVector& hmad) const;

	/*! \brief Set the initial evolution parameters and the grid given the automaton \a system, the bounding domain \a domain and the safe region \a safe.*/
	void _setInitialParameters(SystemType& system, 
							   const HybridBoxes& domain,
							   const HybridBoxes& safe,
							   const RealConstantSet& locked_constants);

	/*! \brief Tune the parameters for the next verification iteration, given previous statistics. */
	void _tuneIterativeStepParameters(SystemType& system);

	/* \brief Processes the \a positive_int and \a negative_int initial intervals based on the lower and upper results.
	 * \return A variable determining if we must proceed further with bisection refining, and another variable determining
	 * the bound where positive values are present.
	 */
	std::pair<bool,bool> _process_initial_bisection_results(Interval& positive_int,
														    Interval& negative_int,
														    const Interval& parameter_range,
														    const tribool& lower_result,
														    const tribool& upper_result) const;

	/*! \brief Shows the verification \a result into the standard output. */
	void _log_verification_result(const tribool& result) const;

	/* \brief Processes the \a result in order to update the \a positive_int interval, possibly updating \a negative_int too. */
	void _process_positive_bisection_result(const tribool& result,
											Interval& positive_int,
										    Interval& negative_int,
											const Float& current_value,
											const bool& safeOnBottom) const;

	/*! \brief Processes the \a result in order to update the \a negative_int interval, possibly updating \a positive_int too. */
	void _process_negative_bisection_result(const tribool& result,
											Interval& positive_int,
										    Interval& negative_int,
											const Float& current_value,
											const bool& safeOnBottom) const;

	/*! \brief Converts the positive/negative search intervals into positive/negative bounds.
	 * \details The result is obtained by knowing the range of the parameter \a parameter_range and the side where
	 * positive values hold, deduced from \a positiveOnBottom. */
	std::pair<Interval,Interval> _pos_neg_bounds_from_search_intervals(const Interval& positive_int,
												 	 	 	 	 	 	 	  const Interval& negative_int,
												 	 	 	 	 	 	 	  const Interval& parameter_range,
												 	 	 	 	 	 	 	  const bool& positiveOnBottom) const;

	/*! \brief Performs one verification sweep along the X axis if \a sweepOnX is true, the Y axis otherwise. */
	void _parametric_verification_2d_bisection_sweep(Parametric2DBisectionResults& results,
						  	  	  	    SystemVerificationInfo& verInfo,
						  	  	  	    RealConstant xParam,
						  	  	  	    RealConstant yParam,
						  	  	  	    const Float& tolerance,
						  	  	  	    const uint& numPointsPerAxis,
						  	  	  	    bool sweepOnX);

	/*! \brief Performs one dominance sweep along the X axis if \a sweepOnX is true, the Y axis otherwise. */
	void _parametric_dominance_2d_bisection_sweep(Parametric2DBisectionResults& results,
						  	  	  	    		  SystemVerificationInfo& dominating,
						  	  	  	    		  SystemVerificationInfo& dominated,
						  	  	  	    		  RealConstant xParam,
						  	  	  	    		  RealConstant yParam,
						  	  	  	    		  const Float& tolerance,
						  	  	  	    		  const uint& numPointsPerAxis,
						  	  	  	    		  bool sweepOnX);

	/*! \brief Set the parameters for the next dominance iteration, given a bundle of information around a system and a set of constants
	 * that must be ignore when choosing the splitting factors of the system. */
	void _setDominanceParameters(SystemVerificationInfo& systemBundle, const RealConstantSet& lockedConstants);

	/*! \brief Helper function to perform dominance in the more general case when some \a dominatingLockedConstants are enforced. */
	tribool _dominance(SystemVerificationInfo& dominating,
					   SystemVerificationInfo& dominated,
					   const RealConstantSet& dominatingLockedConstants);

	/*! \brief Performs the positive part of dominance checking. */
	bool _dominance_positive(SystemVerificationInfo& dominating,
							 SystemVerificationInfo& dominated,
							 const RealConstantSet& dominatingLockedConstants);

	/*! \brief Performs the negative part of dominance checking. */
	bool _dominance_negative(SystemVerificationInfo& dominating,
							 SystemVerificationInfo& dominated,
							 const RealConstantSet& dominatingLockedConstants);

	/*! \brief Handles splitting the parameters \a current_params, if this is allowed by the given \a tolerance and by a given
	 * verification \a outcome.
	 * \details Additionally saves the un-splittable results into \a output_list. */
	void _split_parameters_set(const tribool& outcome,
							   std::list<RealConstantSet>& working_list,
							   ParametricPartitioningOutcomeList& output_list,
							   RealConstantSet& current_params,
							   const RealConstantSet& working_params,
							   const Float& tolerance) const;
};

template<class HybridEnclosureType>
HybridReachabilityAnalyser::
HybridReachabilityAnalyser(const EvolverInterface<HybridAutomaton,HybridEnclosureType>& evolver)
    : _parameters(new EvolutionParametersType())
	, _statistics(new EvolutionStatisticsType())
    , _discretiser(new HybridDiscretiser<typename HybridEnclosureType::ContinuousStateSetType>(evolver))
{
	this->plot_verify_results = false;
	this->chain_reach_dumping = false;
	this->free_cores = 0;
}


template<class HybridEnclosureType>
HybridReachabilityAnalyser::
HybridReachabilityAnalyser(const EvolutionParametersType& parameters,
                           const EvolverInterface<HybridAutomaton,HybridEnclosureType>& evolver)
    : _parameters(new EvolutionParametersType(parameters))
	, _statistics(new EvolutionStatisticsType())
    , _discretiser(new HybridDiscretiser<typename HybridEnclosureType::ContinuousStateSetType>(evolver))
{
	this->plot_verify_results = false;
	this->chain_reach_dumping = false;
	this->free_cores = 0;
}

} // namespace Ariadne

#endif // ARIADNE_REACHABILITY_ANALYSER_H
