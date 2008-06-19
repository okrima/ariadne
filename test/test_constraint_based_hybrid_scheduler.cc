/***************************************************************************
 *            test_constraint_based_hybrid_scheduler.cc
 *
 *  Copyright  2007  Pieter Collins
 *  Email  Pieter.Collins@cwi.nl
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

#include <iostream>
#include <fstream>
#include <string>

#include "test_float.h"
#include "geometry/set_interface.h"
#include "geometry/set_reference.h"
#include "geometry/hybrid_set.h"
#include "geometry/rectangle.h"
#include "geometry/empty_set.h"
#include "geometry/polyhedral_set.h"
#include "geometry/linear_constraint.h"
#include "system/affine_map.h"
#include "system/affine_vector_field.h"
#include "system/constraint_based_hybrid_automaton.h"
#include "evaluation/evolution_parameters.h"
#include "evaluation/applicator.h"
#include "evaluation/lohner_integrator.h"
#include "evaluation/affine_integrator.h"
#include "evaluation/detector.h"
#include "evaluation/constraint_based_hybrid_evolver.h"
#include "evaluation/constraint_based_hybrid_scheduler.h"
#include "output/epsstream.h"
#include "output/logging.h"

#include "test.h"

static const int verbosity = 0;

using namespace Ariadne;
using namespace std;

static const DiscreteState mode1_id(1);
static const DiscreteState mode2_id(2);
static const DiscreteEvent event3_id(3);
static const DiscreteEvent event4_id(4);
static const DiscreteEvent event5_id(5);
static const DiscreteEvent event6_id(6);
  
template<class R>
EvolutionParameters<R> 
construct_parameters() 
{
  EvolutionParameters<R> parameters;
  parameters.set_maximum_step_size(0.125);
  parameters.set_lock_to_grid_time(0.5);
  parameters.set_maximum_basic_set_radius(0.25);
  parameters.set_grid_length(0.125);
  
  return parameters;
}


template<class R>
ConstraintBasedHybridScheduler<R> 
construct_scheduler() 
{
  typedef Interval<R> I;
  typedef Zonotope<I,I> BS;

  Applicator<R> applicator;
  C1LohnerIntegrator<R> lohner_integrator; 
  Detector<R> detector;
  return ConstraintBasedHybridScheduler<R>(applicator,lohner_integrator,detector);
}


template<class R>
ConstraintBasedHybridAutomaton<R> 
construct_automaton() 
{
    AffineVectorField<R> dynamic1(Matrix<R>("[-0.375,-0.5;0.25,0.25]"),Vector<R>("[2,0.5]"));
    AffineVectorField<R> dynamic2(Matrix<R>("[0,0; 0,0]"),Vector<R>("[0.75,1]"));

    AffineMap<R> reset312(Matrix<R>("[0,1;-1,0]"),Vector<R>("[-2,3]"));
    LinearConstraint<R> guard312(Vector<R>("[1,0]"),less,R(2));
    AffineMap<R> reset412(Matrix<R>("[0,1;-1,0]"),Vector<R>("[-2,1]"));
    LinearConstraint<R> guard412(Vector<R>("[0,-1]"),less,R(2));
    AffineMap<R> reset512(Matrix<R>("[1,0;0,1]"),Vector<R>("[0,-3]"));
    LinearConstraint<R> activation512(Vector<R>("[0,1]"),less,R(2));
    LinearConstraint<R> invariant61(Vector<R>("[0,1]"),less,R(2.5));
    
    ConstraintBasedHybridAutomaton<R> automaton("");

    automaton.new_mode(mode1_id,dynamic1);
    automaton.new_mode(mode2_id,dynamic2);
    automaton.new_forced_transition(event3_id,mode1_id,mode2_id,reset312,guard312);
    automaton.new_forced_transition(event4_id,mode1_id,mode2_id,reset412,guard412);
    automaton.new_unforced_transition(event5_id,mode1_id,mode2_id,reset512,activation512);
    automaton.new_invariant(event6_id,mode1_id,invariant61);

    cout << "automaton = " << flush;
    cout << automaton << endl << endl;

    return automaton;
}




template<class R>
class TestConstraintBasedHybridScheduler
{
  typedef Interval<R> I;
  typedef TimeModelHybridBasicSet< Zonotope<I,I> > timed_set_type;
 public:
  ConstraintBasedHybridAutomaton<R> automaton;
  ConstraintBasedHybridScheduler<R> scheduler;
  EvolutionParameters<R> parameters;
  Box<R> bounding_box;
    

  TestConstraintBasedHybridScheduler()
    : automaton(construct_automaton<R>()), 
      scheduler(construct_scheduler<R>()),
      parameters(construct_parameters<R>()), 
      bounding_box("[-3,3]x[-3,3]")
  { }

  epsfstream& write_invariants(epsfstream& eps) const {
    Box<R> guard3_rectangle("[2,3]x[-3,3]");
    Box<R> guard4_rectangle("[-3,3]x[-3,-2]");
    Box<R> activation5_rectangle("[-3,3]x[2,3]");
    Box<R> invariant6_rectangle("[-3,3]x[2.5,3]");
    eps << fill_colour(cyan) << activation5_rectangle;
    eps << fill_colour(black) << invariant6_rectangle;
    eps << fill_colour(magenta) << guard4_rectangle;
    eps << fill_colour(magenta) << guard3_rectangle;
    return eps;
  }

  void plot(const char* name, timed_set_type initial_set, const std::vector<timed_set_type> reach_set) const
  {
    epsfstream eps;
    std::string filename=std::string("test_constraint_based_hybrid_scheduler-")+name+".eps";
    eps.open(filename.c_str(),bounding_box);
    write_invariants(eps);
    assert(reach_set.size()>=1);
    eps << fill_colour(green);
    eps << reach_set.back().continuous_state_set();
    eps << fill_colour(white);
    for(size_type i=0; i!=reach_set.size()-1; ++i) {
      eps << reach_set[i].continuous_state_set();
    }
    eps << fill_colour(yellow);
    eps << initial_set.continuous_state_set();
    eps.close();
    return;
  }

  int test_reachability_step() {
    cout << "\n" << __FUNCTION__ << "\n";
    time_type step_size=0.25;
    timed_set_type initial_set(mode1_id,Zonotope<I,I>("{(0.5,-0.375),[0.125,0,0.015625;0,0.125,0.015625]}"));
    std::vector<timed_set_type> reach_set=scheduler.evolution_step(automaton,initial_set,step_size,step_size,upper_semantics,compute_reachable_set);
    cout << "\nreach_set=" << reach_set << endl;
    plot("reach",initial_set,reach_set);
    return 0;
  }
  
  int test_guard_step() {
    cout << "\n" << __FUNCTION__ << "\n";
    time_type step_size=0.25;
    time_type maximum_step_size=0.25;
    timed_set_type initial_set(mode1_id,Zonotope<I,I>("{(1.90,-0.4),[0.05,0,0.01;0,0.05,0.01]}"));
    std::vector<timed_set_type> reach_set=scheduler.evolution_step(automaton,initial_set,step_size,step_size,upper_semantics,compute_reachable_set);
    cout << "\nreach_set=" << reach_set << endl;
    plot("guard",initial_set,reach_set);
    return 0;
  }
  
  int test_activation_step() {
    cout << "\n" << __FUNCTION__ << "\n";
    time_type step_size=0.25;
    timed_set_type initial_set(mode1_id,Zonotope<I,I>(Box<R>("[-1.45,-1.35]x[2.05,2.15]")));
    std::vector<timed_set_type> reach_set=scheduler.evolution_step(automaton,initial_set,step_size,step_size,upper_semantics,compute_reachable_set);
    cout << "\nreach_set=" << reach_set << endl;
    plot("activation",initial_set,reach_set);
    return 0;
  }
  
  int test_invariant_step() {
    cout << "\n" << __FUNCTION__ << "\n";
    time_type step_size=0.25;
    TimeModelHybridBasicSet< Zonotope<I,I> > initial_set(mode1_id,Zonotope<I,I>(Box<R>("[1.35,1.45]x[2.35,2.45]")));
    std::vector<timed_set_type> reach_set=scheduler.evolution_step(automaton,initial_set,step_size,step_size,upper_semantics,compute_reachable_set);
    cout << "\nreach_set=" << reach_set << endl;
    plot("invariant",initial_set,reach_set);
    return 0;
  }
  
  int test_initial_guard_step() {
    cout << "\n" << __FUNCTION__ << "\n";
    time_type step_size=0.25;
    TimeModelHybridBasicSet< Zonotope<I,I> > initial_set(mode1_id,Zonotope<I,I>(Box<R>("[1.95,2.05]x[1.35,1.45]")));
    std::vector<timed_set_type> reach_set=scheduler.evolution_step(automaton,initial_set,step_size,step_size,upper_semantics,compute_reachable_set);
    cout << "\nreach_set=" << reach_set << endl;
    plot("initial_guard",initial_set,reach_set);
    return 0;
  }
  
  int test_corner_collision_step() {
    cout << "\n" << __FUNCTION__ << "\n";
    time_type step_size=0.25;
    TimeModelHybridBasicSet< Zonotope<I,I> > initial_set(mode1_id,Zonotope<I,I>(Box<R>("[1.65,1.75]x[2.35,2.45]")));
    std::vector<timed_set_type> subdivided_set=scheduler.evolution_step(automaton,initial_set,step_size,step_size,upper_semantics,compute_reachable_set);
    std::vector<timed_set_type> reach_set;
    cout << "\nsubdivided_set=" << subdivided_set << endl;
    for(uint i=0; i!=subdivided_set.size(); ++i) { 
      std::vector<timed_set_type> reach_subset=scheduler.evolution_step(automaton,subdivided_set[i],step_size,step_size,upper_semantics,compute_reachable_set);
      reach_set.insert(reach_set.end(),reach_subset.begin(),reach_subset.end());
    }
    cout << "reach_set=" << reach_set << endl;
    plot("corner_collision",initial_set,reach_set);
    return 0;
  }
  
  int test() {
    test_reachability_step();
    test_guard_step();
    test_activation_step();
    test_invariant_step();
    test_initial_guard_step();
    //test_corner_collision_step();
    return 0;
  }

};



int main() {
  set_hybrid_evolver_verbosity(::verbosity);

  TestConstraintBasedHybridScheduler<Flt>().test();
  cerr << "INCOMPLETE ";
  return 0;
}

