/***************************************************************************
 *            box.h
 *
 *  Mon 2 May 2005
 *  Copyright 2005  Alberto Casagrande, Pieter Collins
 *  Email casagrande@dimi.uniud.it, Pieter.Collins@cwi.nl
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
 
/*! \file box.h
 *  \brief Rectangles and cuboids.
 */

#ifndef ARIADNE_RECTANGLE_H
#define ARIADNE_RECTANGLE_H

#include <iosfwd>

#include "base/array.h"
#include "base/iterator.h"
#include "base/tribool.h"

#include "numeric/arithmetic.h"
#include "numeric/function.h"
#include "numeric/interval.h"

#include "linear_algebra/vector.h"

#include "combinatoric/declarations.h"

#include "geometry/exceptions.h"
#include "geometry/point.h"
#include "geometry/rectangle_expression.h"

namespace Ariadne {
  namespace Geometry {

    class basic_set_tag;
    template<class R> class PointList;
    template<class BS> class ListSet;
    template<class R> class BoxVerticesIterator;
    
    /*! \brief A box of arbitrary dimension.
     * 
     * The most important geomtric object in %Ariadne, the %Box class 
     * describes a box, cuboid or hypercuboid in Euclidean space.
     * Boxes are closed under most major geometric operations, and geometric
     * predicates and operations are simple to compute.
     * Boxes sets require little data to describe, and are used as cells in grid-based and lattice based sets,
     * making them ubiquitous in storage-representations.
     * Further, boxes are easily and exactly convertible to the other major polyhedral
     * set representations, namely, Rectangle, Zonotope,  Polytope and Polyhedron.
     * 
     * Boxes are decribed by the lower and upper bounds in each of their
     * dimensions, as accessed by the lower_bound(dimension_type) const and upper_bound(dimension_type) const methods. 
     * 
     * Boxes are by default ordered by the lexicographic order on lower-left corner.
     * If these are equal, then points are ordered by the lexicographic order 
     * on the upper-right corner.
     *
     * A box is described by the literal "[a0,b0]x[a1,b1]x...".
     * 
     * Boxes are always bounded, but may be empty. A zero-dimensional box is considered to be non-empty.
     *
     * \b Storage: A %Box of dimension \a d is specified by \a 2d real 
     * number giving the lower and upper bounds in each coordinate.
     * The lower bound in the \a i th coordinate is the \a 2i th element, and the 
     * upper bound is the \a 2i+1 st element.
     */
    template<class R>
    class Box 
      : public RectangleExpression< Box<R> >
    {
      typedef typename Numeric::traits<R>::interval_type I;
      typedef typename Numeric::traits<R>::arithmetic_type A;
     private:
      array<R> _data;
     public:
      /*! \brief A tag describing the type of set. */
      typedef basic_set_tag set_category;
      /*! \brief The type of denotable real number used for the corners. */
      typedef R real_type;
      /*! \brief The type used for the corners. */
      typedef R value_type;
      /*! \brief The type of denotable point contained by the box. */
      typedef Point<R> state_type;
      /*! \brief An iterator to the vertices of the box. */
      typedef BoxVerticesIterator<R> vertices_const_iterator;
     public:
      //@{
      //! \name Constructors
      /*! \brief Construct an empty box with dimension \a d. */
      explicit Box(size_type d=0);
      
      /*! \brief Construct a box from a range of interval values. */
      template<class ForwardIterator> explicit Box(ForwardIterator b, ForwardIterator e);
      
      /*! \brief Construct from a one-dimensional C-style array of the form <code>{ l1,u2, l2,u2, ... , ln,un }</code>. */
      template<class RR> explicit Box(const dimension_type& d, const RR* ptr);
      
      /*! \brief Construct from a two-dimensional C-style array of the form <code>{ {l1,u2}, {l2,u2}, ... , {ln,un} }</code>. */
      template<class RR> explicit Box(const dimension_type& d, const RR ary[][2]);
      
      /*! \brief Construct from an array of intervals. */
      explicit Box(const Base::array< Numeric::Interval<R> >& a);
      
      /*! \brief Construct from a std::vector of intervals. */
      explicit Box(const std::vector< Numeric::Interval<R> >& v);

      /*! \brief Construct a degenerate box from a single point. */
      explicit Box(const Point<R>& pt);
      
      /*! \brief Construct a box from an interval point. */
      explicit Box(const Point< Numeric::Interval<R> >& pt);;
      
      /*! \brief Construct from two corners. */
      explicit Box(const Point<R>& pt1, const Point<R>& pt2);
      
      /*! \brief Construct from a string literal. */
      explicit Box(const std::string& s);
      
      /*! \brief Construct from an interval vector. */
      explicit Box(const LinearAlgebra::Vector< Numeric::Interval<R> >& iv);
      
      /*! \brief Convert from a box expression. */
      template<class E> Box(const RectangleExpression<E>& r);
    
    
      /*! \brief Copy constructor. */
      Box(const Box<R>& r);
    
      /*! \brief Copy assignment operator. */
      Box<R>& operator=(const Box<R>& r);

      /*! \brief Assign from a box expression. */
      template<class E> Box<R>& operator=(const RectangleExpression<E>& r);
    
      //@}
      
      
      //@{
      //! \name Conversion operators
      /*! \brief Convert to an interval point. */
      operator Point< Numeric::Interval<R> >() const;
      //@}
      
      
      //{@
      //! \name Comparison operators
      /*! \brief The equality operator */
      bool operator==(const Box<R>& A) const;
      
      /*! \brief The inequality operator */
      bool operator!=(const Box<R>& A) const;
      //@}
      

      //@{
      //! \name Data access
      /*! \brief Returns a reference to the array of data. */
      array<R>& data();
     
      /*! \brief Returns a constant reference to the array of data. */
      const array<R>& data() const;
     
      /*! \brief The lower bound of the \a i th coordinate */
      const R& lower_bound(dimension_type i) const;
      
      /*! \brief A reference to the lower bound of the \a i th coordinate */
      R& lower_bound(dimension_type i);
      
      /*! \brief The upper bound of the \a i th coordinate */
      const R& upper_bound(dimension_type i) const;
      
      /*! \brief A reference to the upper bound of the \a i th coordinate */
      R& upper_bound(dimension_type i);
      
      /*! \brief Returns the projection onto the \a i th coordinate (unchecked). */
      Numeric::Interval<R>& operator[] (dimension_type i);
     
      /*! \brief The projection onto the \a i th coordinate (unchecked). */
      const Numeric::Interval<R>& operator[] (dimension_type i) const;
      
      /*! \brief The interval of values in the \a i th coordinate. */
      const Numeric::Interval<R>& interval(dimension_type i) const;
      
      /*! \brief The lower corner. */
      Point<R> lower_corner() const;
      
      /*! \brief The upper corner. */
      Point<R> upper_corner() const;
      
      /*! \brief The set of position vectors of the box. */
      LinearAlgebra::Vector< Numeric::Interval<R> > position_vectors() const;
      //@}
      
      
      //@{ 
      //! \name Modifying operations
      /*! \brief Makes the box empty. */
      void clear();
      
      /*! \brief Sets the \a i th interval. */
      void set_interval(dimension_type i, Numeric::Interval<R> x);
      
      /*! \brief Sets the lower bound of the \a i th coordinate to \a r. */
      void set_lower_bound(dimension_type i, const R& l);
      
      /*! \brief Sets the upper bound of the \a i th coordinate to \a u. */
      void set_upper_bound(dimension_type i, const R& u);

      /*! \brief Expand the %Box by \a delta in each direction. (Deprecated, use neighbourhood(...) instead) */
      Box<R>& expand_by(const R& delta);

      /*! \brief Return a copy of the %Box expanded by \a delta in each direction. (Deprecated, use neighbourhood(...) instead) */
      Box<R> expand(const R& delta) const;
      /*! \brief Return a copy of the %Box expanded by \a delta in each direction. */
      Box<R> neighbourhood(const R& delta) const;
      //@}
      
      
      //@{
      //! \name Box geometric operations
      /*! \brief The dimension of the Euclidean space the box lies in. */
      dimension_type dimension() const;
      
      /*! \brief True if the box is empty. A zero-dimensional box is considered empty. */
      tribool empty() const;
      
      /*! \brief The centre. */
      Point<R> centre() const;
      
      /*! \brief The radius in the sup norm. */
      R radius() const;
      
      /*! \brief An approximation to the volume. */
      R volume() const;
      
      /*! \brief Compute a quadrant of the Box determined by \a q.
       *  \a q is a binary word such that the ith bit of q is 0 if the lower half
       *  of the box in the ith coordinate is used, and 1 if the upper
       *  half is used.
       */
      Box<R> quadrant(const Combinatoric::BinaryWord& q) const;
      
      /*! \brief The vertices of the box. */
      PointList<R> vertices() const;

      /*! \brief An iterator to the first vertex of the box. */
      vertices_const_iterator vertices_begin() const;
      /*! \brief An iterator to the end vertex of the box. */
      vertices_const_iterator vertices_end() const;

      /*! \brief The number of vertices. */
      size_type number_of_vertices() const;
        
      /*! \brief The \a i th vertex. */
      Point<R> vertex(size_type i) const;
        
      /*! \brief Tests if \a point is included into a box. */
      tribool contains(const Point<R>& pt) const;
      
      /*! \brief Checks for boundedness. */
      tribool bounded() const;
      
      /*! \brief A box containing the given box; returns a copy. */
      Box bounding_box() const;
      //@}
      
#ifdef DOXYGEN
      //@{ 
      //! \name Binary geometric predicates
      /*! \brief SetInterface equality operator. */
      friend tribool equal(const Box<R>& A, const Box<R>& B) const;
      /*! \brief Tests disjointness with \a r. */
      friend tribool disjoint(const Box<R>& A, const Box<R>& B) const;
      /*! \brief Tests if the box is a subset of another box \a r. */
      friend tribool subset(const Box<R>& A, const Box<R>& B) const;
      //@}

      //@{ 
      //! \name Binary geometric operations
      /*! \brief The intersection of \a A and \a B. */
      friend Box<R> intersection(const Box<R>& A, const Box<R>& B); 
      /*! \brief The closure of the intersection of the interiors of \a A and \a B. */
      friend Box<R> regular_intersection(const Box<R>& A, const Box<R>& B); 

      /*! \brief The smallest box containing \a A and \a B. */
      friend Box<R> rectangular_hull(const Box<R>& A, const Box<R>& B); 
      /*! \brief The smallest box containing \a A and \a B. */
      friend Box<R> rectangular_hull(const Box<R>& A, const Point<R>& B); 
      /*! \brief The smallest box containing \a A and \a B. */
      friend Box<R> rectangular_hull(const Point<R>& A, const Box<R>& B); 
      /*! \brief The smallest box containing \a A and \a B. */
      friend Box<R> rectangular_hull(const Point<R>& A, const Point<R>& B); 

      /*! \brief The componentwise sum of rectangles \a A and \a B. */
      friend Box<R> minkowski_sum(const Box<R>& A, const Box<R>& B); 
      /*! \brief The componentwise difference of rectangles \a A and \a B. */
      friend Box<R> minkowski_difference(const Box<R>& A, const Box<R>& B); 
      
      /*! \brief The difference between two rectangles. */
      friend LinearAlgebra::Vector< Numeric::Interval<R> > operator-(const Box<R>& A, const Box& B);
      /*! \brief Adds a vector to a box. */
      friend Box<R> operator+(const Box<R>& r, const LinearAlgebra::Vector<R>& v);
      /*! \brief Adds an interval vector to a box. */
      friend Box<R> operator+(const Box<R>& r, const LinearAlgebra::Vector< Numeric::Interval<R> >& v);
      /*! \brief Subtracts a vector from a box. */
      friend Box<R> operator-(const Box<R>& r, const LinearAlgebra::Vector<R>& v);
      /*! \brief Subtracts an interval vector from a box. */
      friend Box<R> operator-(const Box<R>& r, const LinearAlgebra::Vector< Numeric::Interval<R> >& v);
      //@}
#endif
      
      //@{ 
      //! \name Input/output operations
      /*! \brief The name of the class. */
      static std::string name();
      /*! \brief Write to an output stream. */
      std::ostream& write(std::ostream& os) const;
      /*! \brief Read from an input stream. */
      std::istream& read(std::istream& is);
      //@}
    };
    
    
    
    template<class R>
    class BoxVerticesIterator 
      : public boost::iterator_facade<BoxVerticesIterator<R>,
                                      Point<R>,
                                      boost::forward_traversal_tag,
                                      Point<R> const&,
                                      Point<R> const*
                                     >
    {
     public:
      BoxVerticesIterator(const Box<R>& r, const bool end);
      bool equal(const BoxVerticesIterator<R>& other) const;
      const Point<R>& dereference() const;
      void increment();
     private:
      const Box<R>* _r; 
      long unsigned int _i; 
      bool _parity; 
      Point<R> _pt;
    };
    
    
    
    
    
    template<class R> tribool equal(const Box<R>& A, const Box<R>& B);
      
    template<class R> tribool disjoint(const Box<R>& A, const Box<R>& B);
      
    template<class R> tribool subset(const Box<R>& A, const Box<R>& B);
    
    template<class R> Box<R> closed_intersection(const Box<R>& A, const Box<R>& B);
      
    template<class R> Box<R> open_intersection(const Box<R>& A, const Box<R>& B);
      
    template<class R> Box<R> rectangular_hull(const Box<R>& A, const Box<R>& B);
      
    
    
    template<class R>
    Geometry::Box<R> 
    operator+(const Geometry::Box<R>& r, 
              const LinearAlgebra::Vector<R>& v);
    
    template<class R>
    Geometry::Box<R> 
    operator-(const Geometry::Box<R>& r, 
              const LinearAlgebra::Vector<R>& v);
    
    template<class R>
    std::istream&
    operator>>(std::istream& is, Box<R>& r);

    template<class R>
    std::ostream&
    operator<<(std::ostream& os, const Box<R>& r);

    
  }
}

#include "box.inline.h"

#endif /* ARIADNE_RECTANGLE_H */