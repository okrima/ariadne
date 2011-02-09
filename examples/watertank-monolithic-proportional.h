/***************************************************************************
 *            watertank-monolithic-proportional.h
 *
 *  Copyright  2011  Luca Geretti
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

#ifndef WATERTANK_MONOLITHIC_PROPORTIONAL_H_
#define WATERTANK_MONOLITHIC_PROPORTIONAL_H_

#include <cstdarg>
#include "ariadne.h"

namespace Ariadne {

HybridAutomaton getWatertankMonolithicProportional()
{
	HybridAutomaton system("watertank-mono-pr");

    // Set the system parameters
	RealConstant a("a",0.02); // The constant defining the decrease rate of the tank level
	RealConstant tau("tau",1.25); // The characteristic time for the opening/closing of the valve
	RealConstant ref("ref",6.75); // A reference tank level
	RealConstant bfp("bfp",Interval(0.3,0.32863)); // The product beta*f(p)
	RealConstant Kp("Kp",0.6); // The gain of the proportional controller
	RealConstant delta("delta",Interval(-0.1,0.1)); // An indeterminacy in guards evaluation

	// System variables
	RealVariable x("x"); // water level
	RealVariable y("y"); // valve level
	List<RealVariable> varlist;
	varlist.append(x);
	varlist.append(y);

	// Accessible constants
	system.register_accessible_constant(delta);
	system.register_accessible_constant(a);
	system.register_accessible_constant(bfp);
	system.register_accessible_constant(tau);
	system.register_accessible_constant(Kp);

	// Constants
	ScalarFunction one=ScalarFunction::constant(2,1.0);
	ScalarFunction zero=ScalarFunction::constant(2,0.0);

	/// Build the Hybrid System

	/// Create four discrete states
	DiscreteState l1(1);      // Zero saturated
	DiscreteState l2(2);      // Stabilized
	DiscreteState l3(3);      // One saturated

	/// Create the discrete events
	DiscreteEvent e12(12);
	DiscreteEvent e21(21);
	DiscreteEvent e23(23);
	DiscreteEvent e32(32);

	/// Create the dynamics

	RealExpression x_d = -a*x+bfp*y;
	RealExpression y_towardszero_d = -y/tau;
	RealExpression y_controlled_d = (Kp*(ref-x-delta)-y)/tau;
	RealExpression y_towardsone_d = (1-y)/tau;

	// Dynamics at the different modes
	List<RealExpression> exprlist;
	exprlist.append(x_d);
	exprlist.append(y_towardszero_d);
	VectorFunction zerosaturated_d(exprlist, varlist);
	exprlist[1] = y_controlled_d;
	VectorFunction controlled_d(exprlist, varlist);
	exprlist[1] = y_towardsone_d;
	VectorFunction onesaturated_d(exprlist, varlist);

	/// Create the reset
	IdentityFunction reset_id(2);

	/// Create the guards.
	/// Guards are true when f(x) = Ax + b > 0
	/// x <= ref - Delta
	RealExpression x_lesser_ref_minus_delta = -x-delta+ref;
	ScalarFunction guard12(x_lesser_ref_minus_delta,varlist);
	//ScalarAffineFunction guard12(Vector<Float>(4, -1.0,0.0,0.0,-1.0),Rif);
	/// x >= ref - Delta
	RealExpression x_greater_ref_minus_delta = x+delta-ref;
	ScalarFunction guard21(x_greater_ref_minus_delta,varlist);
	//ScalarAffineFunction guard21(Vector<Float>(4, 1.0,0.0,0.0,1.0),-Rif);
	/// x <= ref - 1/Kp - Delta
	RealExpression x_lesser_ref_kp_minus_delta = -x+ref-1.0/Kp-delta;
	ScalarFunction guard23(x_lesser_ref_kp_minus_delta,varlist);
	//ScalarAffineFunction guard23(Vector<Float>(4, -1.0,0.0,0.0,-1.0),(Rif-1.0/Kp));
	/// x >= ref - 1/Kp - Delta
	RealExpression x_greater_ref_kp_minus_delta = x-ref+1.0/Kp+delta;
	ScalarFunction guard32(x_greater_ref_kp_minus_delta,varlist);
	//ScalarAffineFunction guard32(Vector<Float>(4, 1.0,0.0,0.0,1.0),(1.0/Kp - Rif));

	/// Create the invariants.
	/// Invariants are true when f(x) = Ax + b < 0
	/// forced transitions do not need an explicit invariant,
	/// hence we do not need invariants

	/// Build the automaton
	system.new_mode(l1,zerosaturated_d);
	system.new_mode(l2,controlled_d);
	system.new_mode(l3,onesaturated_d);

	system.new_forced_transition(e12,l1,l2,reset_id,guard12);
	system.new_forced_transition(e21,l2,l1,reset_id,guard21);
	system.new_forced_transition(e23,l2,l3,reset_id,guard23);
	system.new_forced_transition(e32,l3,l2,reset_id,guard32);

	return system;
}


}

#endif /* WATERTANK_MONOLITHIC_HYSTERESIS_H_ */