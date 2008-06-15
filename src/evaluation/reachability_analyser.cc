/***************************************************************************
 *            reachability_analyser.cc
 *
 *  Copyright  2006-8 Pieter Collins
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

#include "numeric/integer.h"
#include "numeric/rational.h"
#include "numeric/float.h"

#include "system/transition_system_interface.h"
#include "system/map.h"
#include "system/vector_field.h"

#include "evaluation/reachability_analyser.h"
#include "evaluation/reachability_analyser.code.h"

namespace Ariadne {
  namespace Evaluation {
    using namespace Numeric;
    using namespace Geometry;
    using namespace System;

#ifdef ENABLE_FLOAT64
    template class ReachabilityAnalyser< TransitionSystemInterface<Integer,GridApproximationScheme<Float64> >, GridApproximationScheme<Float64> >;
    template class ReachabilityAnalyser< TransitionSystemInterface<Rational,GridApproximationScheme<Float64> >, GridApproximationScheme<Float64> >;
    template class ReachabilityAnalyser< Map<Float64>, GridApproximationScheme<Float64> >;
    template class ReachabilityAnalyser< VectorField<Float64>, GridApproximationScheme<Float64> >;
#endif
  
#ifdef ENABLE_FLOATMP
    template class ReachabilityAnalyser< TransitionSystemInterface<Integer,GridApproximationScheme<FloatMP> >, GridApproximationScheme<FloatMP> >;
    template class ReachabilityAnalyser< TransitionSystemInterface<Rational,GridApproximationScheme<FloatMP> >, GridApproximationScheme<FloatMP> >;
    template class ReachabilityAnalyser< Map<FloatMP>, GridApproximationScheme<FloatMP> >;
    template class ReachabilityAnalyser< VectorField<FloatMP>, GridApproximationScheme<FloatMP> >;
#endif

  }
}