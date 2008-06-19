/***************************************************************************
 *            grid_cell_list_set.cc
 *
 *  Copyright  2006-7  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, Pieter.Collins@cwi.nl
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

#include "geometry/grid_cell_list_set.h"
#include "geometry/grid_cell_list_set.code.h"

#include "numeric/float.h"

namespace Ariadne {
  

    
    
#ifdef ENABLE_FLOAT64
    template class GridCellListSet<Float64>;
  //    template class GridCellListSetIterator<Float64>;
#endif

#ifdef ENABLE_FLOATMP
    template class GridCellListSet<FloatMP>;
  //    template class GridCellListSetIterator<FloatMP>;
#endif
    
  
} // namespace Ariadne
