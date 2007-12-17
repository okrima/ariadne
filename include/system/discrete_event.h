/***************************************************************************
 *            discrete_event.h
 *
 *  Copyright  2007  Pieter Collins
 *  Pieter.Collins@cwi.nl
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
 
#ifndef ARIADNE_DISCRETE_EVENT_H
#define ARIADNE_DISCRETE_EVENT_H

#include "base/types.h"

namespace Ariadne {
  namespace System {

    class DiscreteEvent;
    std::ostream& operator<<(std::ostream&, const DiscreteEvent&);

    /*! \brief A class representing the identifier of a discrete mode. */
    class DiscreteEvent {
     public:
      //! \brief Construct from an identifier. 
      explicit DiscreteEvent(const id_type& id) : _id(id) { }
      id_type id() const { return this->_id; }
      bool operator==(const DiscreteEvent& other) const { return this->_id==other._id; }
      bool operator!=(const DiscreteEvent& other) const { return this->_id!=other._id; }
      bool operator<(const DiscreteEvent& other) const { return this->_id<other._id; }
     private:
      friend std::ostream& operator<<(std::ostream&, const DiscreteEvent&);
     private:
      int _id;
    };

    inline std::ostream& operator<<(std::ostream& os, const DiscreteEvent& event) { return os << event._id; }
  }

}

#endif 