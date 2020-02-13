/***************************************************************************
 *            solvers/bounder.hpp
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

/*! \file solvers/bounder.hpp
 *  \brief Classes for finding the bounds of the flow for differential equations.
 */

#ifndef ARIADNE_BOUNDER_HPP
#define ARIADNE_BOUNDER_HPP

#include "../utility/typedefs.hpp"
#include "../utility/attribute.hpp"
#include "../algebra/algebra.hpp"
#include "../function/domain.hpp"
#include "../function/function_model.hpp"
#include "../output/logging.hpp"

#include "configuration_interface.hpp"
#include "integrator_interface.hpp"

namespace Ariadne {

class BoundingNotFoundException : public std::runtime_error {
  public:
    BoundingNotFoundException(const String& str) : std::runtime_error(str) { }
};


class BounderConfiguration : public ConfigurationInterface {
  public:
    typedef Dyadic RealType;
  public:
    virtual ~BounderConfiguration() = default;
    BounderConfiguration();
  public:
    RealType const& minimum_step_size() const { return _minimum_step_size; }
    BounderConfiguration& set_minimum_step_size(RealType value) {
        assert(0<=value && value<=_maximum_step_size); _minimum_step_size=value; return *this; }

    RealType const& maximum_step_size() const { return _maximum_step_size; }
    BounderConfiguration& set_maximum_step_size(RealType value) {
        assert(_minimum_step_size<=value); _maximum_step_size=value; return *this; }

    RealType const& lipschitz_tolerance() const { return _lipschitz_tolerance; }
    BounderConfiguration& set_lipschitz_tolerance(RealType value) {
        assert(0<=value); _lipschitz_tolerance=value; return *this; }
  private:
    RealType _minimum_step_size;
    RealType _maximum_step_size;
    RealType _lipschitz_tolerance;
  protected:
    virtual OutputStream& _write(OutputStream& os) const override;
};

//! \ingroup SolverModule EvaluationModule
//! \brief Interface for classes calculating the bounds of a flow.
class BounderInterface {
  public:
    //! \brief Compute a pair \f$(h,B)\f$ consisting of a bound \a B for the flow
    //! of \f$\dt{x}=f(x)\f$ starting in \f$dom\f$  for time step \f$h\leq h_{sug}\f$.
    //! <br>
    //! Arguments: \f$f\f$ is the \a vector_field, \f$dom\f$ is the \a state_domain and \f$h_{sug}\f$ is the suggested time step.
    virtual Pair<StepSizeType,UpperBoxType> compute(ValidatedVectorMultivariateFunction const& f, BoxDomainType const& D, StepSizeType const& hsug) const = 0;

    //! \brief Compute a pair \f$(h,B)\f$ consisting of a bound \a B for the flow
    //! of \f$\dt{x}=f(x,t)\f$ starting in \f$dom\f$  for time step \f$h\leq h_{sug}\f$.
    virtual Pair<StepSizeType,UpperBoxType> compute(ValidatedVectorMultivariateFunction const& differential_equation, BoxDomainType const& state_domain, StepSizeType const& initial_time, StepSizeType const& suggested_time_step) const = 0;

    //! \brief Compute a pair \f$(h,B)\f$ consisting of a bound \a B for the flow
    //! of \f$\dt{x}=f(x,a)\f$ starting in \f$D\f$ with \f$a\in A\f$ for time step \f$h\leq h_{sug}\f$.
    //! <br>
    //! Arguments: \f$f\f$ is the \a vector_field, \f$D\f$ is the \a state_domain and \f$h_{sug}\f$ is the suggested time step.
    virtual Pair<StepSizeType,UpperBoxType> compute(ValidatedVectorMultivariateFunction const& vector_field, BoxDomainType const& state_domain, BoxDomainType const& parameter_domain, StepSizeType const& suggested_time_step) const = 0;

    //! \brief Compute a pair \f$(h,B)\f$ consisting of a bound \a B for the flow
    //! of \f$\dt{x}=f(x,t,a)\f$ starting in \f$D\f$ with \f$a\in A\f$ for time step \f$h\leq h_{sug}\f$.
    //! Arguments: \f$f\f$ is the \a vector_field, \f$dom\f$ is the \a state_domain and \f$h_{sug}\f$ is the suggested time step.
    virtual Pair<StepSizeType,UpperBoxType> compute(ValidatedVectorMultivariateFunction const& differential_equation, BoxDomainType const& state_domain, StepSizeType const& initial_time, BoxDomainType const& parameter_domain, StepSizeType const& suggested_time_step) const = 0;

  public:
    virtual OutputStream& _write(OutputStream& os) const = 0;
    virtual BounderInterface* clone() const = 0;

    friend inline OutputStream& operator<<(OutputStream& os, BounderInterface const& bounder) {
        bounder._write(os); return os; }

    virtual ~BounderInterface() = default;
    virtual BounderConfiguration& configuration() = 0;
    virtual BounderConfiguration const& configuration() const = 0;
};


class BounderBaseConfiguration : public BounderConfiguration { };

class BounderBase : public BounderInterface, Loggable {
  protected:
    SharedPointer<BounderBaseConfiguration> _configuration;
  protected:
    BounderBase(BounderBaseConfiguration* config);
  public:
    virtual Pair<StepSizeType,UpperBoxType> compute(ValidatedVectorMultivariateFunction const& f, BoxDomainType const& D, StepSizeType const& hsug) const override;
    virtual Pair<StepSizeType,UpperBoxType> compute(ValidatedVectorMultivariateFunction const& f, BoxDomainType const& D, StepSizeType const& t, StepSizeType const& hsug) const override;
    virtual Pair<StepSizeType,UpperBoxType> compute(ValidatedVectorMultivariateFunction const& f, BoxDomainType const& D, BoxDomainType const& A, StepSizeType const& hsug) const override;
    virtual Pair<StepSizeType,UpperBoxType> compute(ValidatedVectorMultivariateFunction const& f, BoxDomainType const& D, StepSizeType const& t, BoxDomainType const& A, StepSizeType const& hsug) const override = 0;
};


class ExpandContractBounderConfiguration : public BounderBaseConfiguration {
  public:
    ExpandContractBounderConfiguration();
  public:
    RealType const& domain_widening() const { return _domain_widening; }
    ExpandContractBounderConfiguration& set_domain_widening(RealType value) { _domain_widening=value; return *this; }
    RealType const& starting_widening() const { return _starting_widening; }
    ExpandContractBounderConfiguration& set_starting_widening(RealType value) { _starting_widening=value; return *this; }
    RealType const& expanding_widening() const { return _expanding_widening; }
    ExpandContractBounderConfiguration& set_expanding_widening(RealType value) { _expanding_widening=value; return *this; }
    Nat const& expansion_steps() const { return _expansion_steps; }
    ExpandContractBounderConfiguration& set_expansion_steps(Nat value) { _expansion_steps=value; return *this; }
    Nat const& refinement_steps() const { return _refinement_steps; }
    ExpandContractBounderConfiguration& set_refinement_steps(Nat value) { _refinement_steps=value; return *this; }
  private:
    RealType _domain_widening;
    RealType _starting_widening;
    RealType _expanding_widening;
    Nat _expansion_steps;
    Nat _refinement_steps;
  protected:
    virtual OutputStream& _write(OutputStream& os) const override;
};

class ExpandContractBounder : public BounderBase {
  public:
    using BounderBase::compute;
    // Compute the bounds on the reach set of dx/dt = f(x,t,a) starting at time t, for a suggested step size of h.
    virtual Pair<StepSizeType,UpperBoxType> compute(ValidatedVectorMultivariateFunction const& f, BoxDomainType const& D, StepSizeType const& t, BoxDomainType const& A, StepSizeType const& hsug) const override;
  private:
    virtual OutputStream& _write(OutputStream& os) const override = 0;
  protected:
    ExpandContractBounder();
    ExpandContractBounder(ExpandContractBounderConfiguration* config);
  public:
    virtual ExpandContractBounderConfiguration& configuration() override {
        return static_cast<ExpandContractBounderConfiguration&>(*_configuration); }
    virtual ExpandContractBounderConfiguration const& configuration() const override {
        return static_cast<ExpandContractBounderConfiguration const&>(*_configuration); }
  private: public:
    // Estimate a bounding box for the flow of f starting in D at T.lower() to time T.upper() with parameters A, assuming flow remains in B,
    // but widening the flow vectors by VECTOR_WIDENING.
    virtual UpperBoxType _formula(ValidatedVectorMultivariateFunction const& f, BoxDomainType const& D, IntervalDomainType const& T, BoxDomainType const& A, UpperBoxType const& B) const = 0;
};

struct EulerBounderConfiguration : public ExpandContractBounderConfiguration { };

//! \brief A flow bounder which uses the estimate \f$\phi(D,T)\subset D+hf(B,T)\f$ to test refinement.
class EulerBounder final : public ExpandContractBounder {
  public:
    virtual OutputStream& _write(OutputStream& os) const override { return os << "EulerBounder(configuation="<<this->configuration()<<")"; }
    virtual EulerBounder* clone() const override { return new EulerBounder(*this); }
  private: public:
    virtual UpperBoxType _formula(ValidatedVectorMultivariateFunction const& f, BoxDomainType const& D, IntervalDomainType const& T, BoxDomainType const& A, UpperBoxType const& B) const override;
};




class BounderHandle {
  private:
    SharedPointer<BounderInterface> _impl;
  public:
    BounderHandle(BounderInterface const& bounder) : _impl(bounder.clone()) { }
    BounderHandle(BounderHandle const& other) : _impl(other._impl) { }
    BounderHandle& operator=(BounderHandle const& other) { _impl = other._impl; return *this; }

    operator BounderInterface const& () const { return *_impl; }

    BounderHandle* clone() const { return new BounderHandle(*this); }

    Pair<StepSizeType,UpperBoxType> compute(ValidatedVectorMultivariateFunction const& f, BoxDomainType const& D, StepSizeType const& hsug) const {
        return _impl->compute(f,D,hsug);
    }

    Pair<StepSizeType,UpperBoxType> compute(ValidatedVectorMultivariateFunction const& f, BoxDomainType const& D, BoxDomainType const& A, StepSizeType const& hsug) const {
        return _impl->compute(f,D,A,hsug);
    }

    Pair<StepSizeType,UpperBoxType> compute(ValidatedVectorMultivariateFunction const& f, BoxDomainType const& D, StepSizeType const& t, StepSizeType const& hsug) const {
        return _impl->compute(f,D,t,hsug);
    }

    Pair<StepSizeType,UpperBoxType> compute(ValidatedVectorMultivariateFunction const& f, BoxDomainType const& D, StepSizeType const& t, BoxDomainType const& A, StepSizeType const& hsug) const {
        return _impl->compute(f,D,t,A,hsug);
    }

    friend OutputStream& operator<<(OutputStream& os, BounderHandle const& b) { return b._impl->_write(os); }

};



} // namespace Ariadne

#endif /* ARIADNE_BOUNDER_HPP */
