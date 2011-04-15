/***************************************************************************
 *            hybrid_space.h
 *
 *  Copyright 2008-11  Pieter Collins
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

/*! \file hybrid_space.h
 *  \brief Sets in hybrid spaces.
 */

#ifndef ARIADNE_HYBRID_SPACE_H
#define ARIADNE_HYBRID_SPACE_H

#include <map>

#include "container.h"
#include "stlio.h"
#include "space.h"
#include "discrete_location.h"
#include "hybrid_set_interface.h"


namespace Ariadne {

class LocationException : public std::runtime_error {
  public:
    LocationException(const std::string& what) : std::runtime_error(what) { }
};

class HybridGridTreeSet;

class HybridSpaceInterface
{
  public:
    virtual bool has_location(const DiscreteLocation& q) const  = 0;
    virtual RealSpace operator[](const DiscreteLocation& q) const = 0;
  public:
    virtual HybridSpaceInterface* clone() const = 0;
    virtual std::ostream& write(std::ostream& os) const = 0;
  public:
    friend std::ostream& operator<<(std::ostream& os, const HybridSpaceInterface& hsp) { return hsp.write(os); }
};

//! \ingroup HybridModule
//! \brief A hybrid space \f$\bigsqcup_{q\in Q} \R^{d_q}\f$ with discrete states \f$Q\f$.
class HybridSpace
{
  public:
    //! \brief The interface satisified by bounded sets in the space.
    typedef HybridBoundedSetInterface BoundedSetInterfaceType;
    //! \brief The interface satisified by overt sets in the space.
    typedef HybridOvertSetInterface OvertSetInterfaceType;
    //! \brief The interface satisified by over sets in the space.
    typedef HybridOpenSetInterface OpenSetInterfaceType;
    //! \brief The interface satisified by closed sets in the space.
    typedef HybridClosedSetInterface ClosedSetInterfaceType;
    //! \brief The interface satisified by compact sets in the space.
    typedef HybridCompactSetInterface CompactSetInterfaceType;
    //! \brief The interface satisified by regular sets in the space.
    typedef HybridRegularSetInterface RegularSetInterfaceType;
    //! \brief The interface satisified by located sets in the space.
    typedef HybridLocatedSetInterface LocatedSetInterfaceType;
    //! \brief The type of approximations to sets in the space.
    typedef HybridGridTreeSet SetApproximationType;
  public:
    HybridSpace(const HybridSpaceInterface& hspc) : _ptr(hspc.clone()) { }
    HybridSpace(const HybridSpaceInterface* hspc_ptr) : _ptr(hspc_ptr) { }

    bool has_location(const DiscreteLocation& q) const { return this->_ptr->has_location(q); }
    RealSpace operator[](const DiscreteLocation& q) const { return this->_ptr->operator[](q); }

    operator const HybridSpaceInterface& () const { return *_ptr; }

    friend std::ostream& operator<<(std::ostream& os, const HybridSpace& hsp) { return os << *hsp._ptr; }
  private:
    shared_ptr<const HybridSpaceInterface> _ptr;
};


//! \ingroup HybridModule
//! \brief A hybrid space \f$\bigsqcup_{q\in Q} \R^{d_q}\f$ with discrete states \f$Q\f$.
class MonolithicHybridSpace
    : public HybridSpaceInterface
{
  public:
    typedef Map<DiscreteLocation, RealSpace >::const_iterator const_iterator;

    MonolithicHybridSpace() : _locations() { }

    MonolithicHybridSpace* clone() const { return new MonolithicHybridSpace(*this); }

    void new_location(const DiscreteLocation& q, const RealSpace& spc) { this->_locations.insert(q,spc); }

    bool has_location(const DiscreteLocation& q) const { return _locations.has_key(q); }

    RealSpace operator[](const DiscreteLocation& q) const {
        const_iterator iter = this->_locations.find(q);
        if(iter==this->_locations.end()) {
            ARIADNE_THROW(LocationException,"HybridSpace[DiscreteLocation q]","Space has no location "<<q); }
        return iter->second; }

    Set<DiscreteLocation> locations() const { return this->_locations.keys(); }
    const_iterator begin() const { return this->_locations.begin(); }
    const_iterator end() const { return this->_locations.end(); }

    std::ostream& write(std::ostream& os) const { return os << "HybridSpace( " << this->_locations << " )"; }
  private:
    Map< DiscreteLocation, RealSpace > _locations;
};



}

#endif // ARIADNE_HYBRID_SPACE_H
