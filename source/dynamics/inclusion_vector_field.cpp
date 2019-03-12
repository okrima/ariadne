/***************************************************************************
 *            inclusion_vector_field.cpp
 *
 *  Copyright  2008-18  Luca Geretti
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

#include "inclusion_vector_field.hpp"
#include "../symbolic/expression_set.hpp"
#include "../function/symbolic_function.hpp"

namespace Ariadne {

typedef Pair<Nat,EffectiveFormula> CoordinateFormulaPair;
typedef List<CoordinateFormulaPair> CoordinateFormulaPairs;

List<Nat> input_indices(SizeType num_variables, SizeType num_inputs) {
    List<Nat> result;
    for (Nat i : range(num_variables,num_variables+num_inputs)) { result.append(i); }
    return result;
}

BoxDomainType input_bounds_to_domain(RealVariablesBox const& inputs) {
    List<IntervalDomainType> result;

    auto vars = inputs.variables();
    for (auto v : vars) {
        result.push_back(cast_exact(IntervalDomainType(inputs[v].lower().get_d(),inputs[v].upper().get_d())));
    }
    return Vector<IntervalDomainType>(result);
}

Pair<CoordinateFormulaPair,ExactIntervalType> centered_coordinate_transformation(Nat const& i, ExactIntervalType const& bounds) {
    if (same(bounds.lower(),-bounds.upper())) return Pair<CoordinateFormulaPair,ExactIntervalType>({i,EffectiveFormula::coordinate(i)},bounds);
    else return Pair<CoordinateFormulaPair,ExactIntervalType>({i,EffectiveFormula::coordinate(i)+EffectiveFormula::constant(EffectiveNumber(bounds.midpoint()))},ExactIntervalType(cast_exact(bounds.lower()-bounds.midpoint()),cast_exact(bounds.upper()-bounds.midpoint())));
}


Pair<CoordinateFormulaPairs,BoxDomainType> centered_coordinates_transformation(Nat const& dim, BoxDomainType const& inputs) {
    CoordinateFormulaPairs assignments;
    BoxDomainType new_bounds(inputs.size());
    for (auto i : range(inputs.size())) {
        auto tr = centered_coordinate_transformation(i+dim,inputs[i]);
        assignments.push_back(tr.first);
        new_bounds[i] = tr.second;
    }
    return Pair<CoordinateFormulaPairs,BoxDomainType>(assignments,new_bounds);
}

Void incorporate_additive_inputs_coefficients(Vector<EffectiveFormula>& transformed_formulae, BoxDomainType& transformed_inputs) {

    SizeType n = transformed_formulae.size();
    SizeType m = transformed_inputs.size();
    for (SizeType j : range(m)) {
        for (SizeType i : range(n)) {
            EffectiveFormula der = simplify(derivative(transformed_formulae[i],n+j));
            if (not identical(der,EffectiveFormula::zero()) and not identical(der,EffectiveFormula::constant(1))) {
                transformed_inputs[j] = cast_exact(transformed_inputs[j] * der.val());
                transformed_formulae[i] = simplify(substitute(transformed_formulae[i],n+j,EffectiveFormula::zero())+EffectiveFormula::coordinate(n+j));
            }
        }
    }
}

EffectiveVectorFormulaFunction noise_independent_component(EffectiveVectorFormulaFunction const& function, SizeType num_inputs) {

    CoordinateFormulaPairs substitutions;
    for (auto i : range(function.result_size(),function.result_size()+num_inputs)) {
        substitutions.append({i,EffectiveFormula::zero()});
    }

    return EffectiveVectorFormulaFunction(function.argument_size(),simplify(substitute(function._formulae,substitutions)));
}

Vector<EffectiveVectorMultivariateFunction> input_derivatives(EffectiveVectorFormulaFunction const& function, SizeType num_inputs) {

    Vector<EffectiveVectorMultivariateFunction> result(num_inputs);

    SizeType n = function.result_size();

    for (auto j : range(num_inputs)) {
        Vector<EffectiveFormula> derivative_formulae(n);
        for (auto i : range(n)) {
            derivative_formulae[i] = simplify(derivative(function._formulae[i],n+j));
        }
        result[j] = EffectiveVectorFormulaFunction(function.argument_size(),derivative_formulae);
    }

    return result;
}


Void InclusionVectorField::_acquire_and_assign_properties() {

    List<Nat> input_indices = Ariadne::input_indices(this->dimension(),this->number_of_inputs());

    const EffectiveVectorFormulaFunction& ff = dynamic_cast<const EffectiveVectorFormulaFunction&>(_function.reference());

    _noise_independent_component = Ariadne::noise_independent_component(ff,this->number_of_inputs());
    _input_derivatives = Ariadne::input_derivatives(ff,this->number_of_inputs());

    _is_input_affine = is_affine_in(ff,input_indices);
    _is_input_additive = is_additive_in(ff,input_indices);
}

Void InclusionVectorField::_transform_and_assign(EffectiveVectorMultivariateFunction const& function, BoxDomainType const& inputs) {

    const EffectiveVectorFormulaFunction& ff = dynamic_cast<const EffectiveVectorFormulaFunction&>(function.reference());

    auto transformation = centered_coordinates_transformation(function.result_size(),inputs);
    CoordinateFormulaPairs centering_substitution = transformation.first;
    BoxDomainType transformed_inputs = transformation.second;

    Vector<EffectiveFormula> transformed_formulae = substitute(ff._formulae,centering_substitution);

    List<Nat> input_indices = Ariadne::input_indices(function.result_size(),inputs.size());
    if (is_additive_in(transformed_formulae,input_indices)) {
        incorporate_additive_inputs_coefficients(transformed_formulae,transformed_inputs);
    }

    _function = EffectiveVectorFormulaFunction(function.argument_size(),transformed_formulae);
    _inputs = transformed_inputs;
}

InclusionVectorField::InclusionVectorField(DottedRealAssignments const& dynamics, RealVariablesBox const& inputs)
{
    List<RealVariable> dyn_var_list = left_hand_sides(dynamics);
    List<RealVariable> inp_var_list;
    for (auto var : inputs.variables()) { inp_var_list.append(var); }

    Set<UntypedVariable> dyn_arg_set;
    for (auto expr : right_hand_sides(dynamics)) {
        dyn_arg_set.adjoin(expr.arguments());
    }

    Set<RealVariable> inp_var_set(inp_var_list);
    for (auto iv : inp_var_set) {
        if(not dyn_arg_set.contains(iv))
            ARIADNE_THROW(UnusedInputException,"InclusionVectorField(dynamics,inputs) with dynamics="<<dynamics<<" and inputs="<<inputs,"The input '" << iv << "' is not in the arguments of the dynamics");
    }

    Set<UntypedVariable> dyn_var_set(dyn_var_list);
    Set<UntypedVariable> inp_var_set_untyped(inp_var_set);
    for (auto av : dyn_arg_set) {
        if (not dyn_var_set.contains(av) and not inp_var_set_untyped.contains(av))
            ARIADNE_THROW(MissingInputException,"InclusionVectorField(dynamics,inputs) with dynamics="<<dynamics<<" and inputs="<<inputs,"The argument '" << av << "' of the dynamics is not among the inputs");
    }

    RealSpace var_spc(dyn_var_list);
    RealSpace inp_spc(inp_var_list);
    RealSpace spc = var_spc.adjoin(inp_spc);

    _variable_names = variable_names(left_hand_sides(dynamics));

    _transform_and_assign(make_function(spc,Vector<RealExpression>(right_hand_sides(dynamics))),input_bounds_to_domain(inputs));
    _acquire_and_assign_properties();
}

InclusionVectorField::InclusionVectorField(EffectiveVectorMultivariateFunction const& function, BoxDomainType const& inputs)
{
    if (function.argument_size() != function.result_size()+inputs.size())
        ARIADNE_THROW(FunctionArgumentsMismatchException,"InclusionVectorField(function,inputs) with function="<<function<<" and inputs="<<inputs,
                "Incompatible inputs size (" << inputs.size() << ") with the provided function (R^" << function.argument_size() << " -> R^" << function.result_size() << ")");

    if (dynamic_cast<const EffectiveVectorFormulaFunction*>(function.raw_pointer()) == nullptr)
        ARIADNE_THROW(NotFormulaFunctionException,"InclusionVectorField of EffectiveVectorMultivariateFunction","The function must be constructed from a Formula at the moment.");

    List<Identifier> variable_names;
    for (auto i : range(0,function.result_size()))
        variable_names.append(Identifier("x"+std::to_string(i)));
    _variable_names = variable_names;

    _transform_and_assign(function,inputs);
    _acquire_and_assign_properties();
}


RealSpace InclusionVectorField::state_space() const
{
    return real_space(this->_variable_names);
}

}