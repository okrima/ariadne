/***************************************************************************
 *            parallelotope.h
 *
 *  6 January 2006
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, pieter.collins@cwi.nl
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
 *  along with this program; if not, write to bouthe Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */
 
/*! \file parallelotope.h
 *  \brief Parallelotopes.
 */

#ifndef _ARIADNE_PARALLELOTOPE_H
#define _ARIADNE_PARALLELOTOPE_H

#include <iosfwd>

#include "../declarations.h"

#include "../numeric/interval.h"

#include "../linear_algebra/vector.h"
#include "../linear_algebra/matrix.h"

#include "../geometry/point.h"
#include "../geometry/polyhedron.h"

namespace Ariadne {
  namespace Geometry {
 
    template<> 
    inline bool is_a<Parallelotope,Parallelotope>() { return true; }
    template<> 
    inline bool is_a<Parallelotope,Zonotope>() { return true; }
    template<>
    inline bool is_a<Parallelotope,Polyhedron>() { return true; }

    /* Forward declaration of friends. */
    template<typename R> std::ostream& operator<<(std::ostream&, const Parallelotope<R>&);
    template<typename R> std::istream& operator>>(std::istream&, Parallelotope<R>&);

    /*! \brief A parallelotope of arbitrary dimension.
     */
    template <typename R>
    class Parallelotope {
     private:
      typedef typename numerical_traits<R>::field_extension_type F;
     public:
      /*! \brief The type of denotable real number used for the corners. */
      typedef R real_type;
      /*! \brief The type of denotable point contained by the rectangle. */
      typedef Point<R> state_type;
      /*! \brief The type of Matrix giving principal directions. */
      typedef LinearAlgebra::Vector<R> Vector_type;
      /*! \brief The type of Matrix giving principal directions. */
      typedef LinearAlgebra::Matrix<R> Matrix_type;
     public:
      /*! \brief Default constructor constructs an empty parallelotope of dimension \a n. */
      explicit Parallelotope(size_type n = 0)
        : _centre(n),  _generators(n,n) { }
      
      /*! \brief Construct from centre and directions. */
      explicit Parallelotope(const Vector_type& c, const Matrix_type& m)
        : _centre(c), _generators(m)
      {
        if (m.size1()!=m.size2()) {
          throw std::domain_error(
              "The the Matrix of principal directions is not a square Matrix");
        }
        
        if (c.size()!=m.size1()) {
          throw std::domain_error("The centre and directions have different dimensions.");
        }
      }
      
      /*! \brief Construct from centre and directions. */
      explicit Parallelotope(const state_type& c, const Matrix_type& m)
        : _centre(c), _generators(m)
      {
        if (m.size1()!=m.size2()) {
          throw std::domain_error(
              "The the Matrix of principal directions is not a square Matrix");
        }
        
        if (c.dimension()!=m.size1()) {
          throw std::domain_error("The centre and directions have different dimensions.");
        }
      }
       
      /*! \brief Construct from a Rectangle. */
      explicit Parallelotope(const Rectangle<real_type>& r)
        : _centre(r.dimension()), _generators(r.dimension(),r.dimension())
      {
        for(size_type i=0; i!=dimension(); ++i) {
          _centre[i] = (r.lower_bound(i)+r.upper_bound(i))/2;
          _generators(i,i) = (r.upper_bound(i)-r.lower_bound(i))/2;
        }
      }
      
      /*! \brief Construct from a string literal. */
      explicit Parallelotope(const std::string& s)
        : _centre(), _generators()
      {
        std::stringstream ss(s);
        ss >> *this;
      }
      
      /*! \brief Copy constructor. */
      Parallelotope(const Parallelotope<R>& original)
        : _centre(original._centre),
          _generators(original._generators)
      { }
      
      /*! \brief Copy assignment operator. */
      Parallelotope<R>& operator=(const Parallelotope<R>& original) {
        if(this != &original) {
          this->_centre = original._centre;
          this->_generators = original._generators;
        }
        return *this;
      }
      
      /*! \brief The equality operator */
      bool operator==(const Parallelotope<R>& A) const;
      
      /*! \brief The inequality operator */
      bool operator!=(const Parallelotope<R>& A) const; 

      /*! \brief The dimension of the Euclidean space the parallelotope lies in. */
      size_type dimension() const {
        return (this->_centre).dimension();
      }
      
      /*! \brief True if the parallelotope is empty. */
      bool empty() const {
        // FIXME: This is probably ok since we're not talking about interior.
        return false;
      }
      
      /*! \brief True if the paralleltope has empty interior. */
      bool empty_interior() const {
        using namespace Ariadne::LinearAlgebra;

        return !independent_rows(this->_generators);     
      }
      
      /*! \brief The centre of the parallelotope. */
      state_type centre() const {
        return this->_centre;
      }
      
      /*! \brief The radius of the parallelotope. */
      real_type radius() const {
        return this->_generators.norm();
      }
      
      /*! \brief The \a n th of principle direction. */
      Vector_type generator(size_type n) const {
        return column(this->_generators,n);
      }
      
      /*! \brief The Matrix of principle directions. */
      Matrix_type generators() const {
        return this->_generators;
      }
     
      /*! \brief A rectangle containing the given parallelotope. */
      Rectangle<R> bounding_box() const;
      
      /*! \brief The \a i th vertex. */
      state_type vertex(const size_type& i) const;
      
      /*! \brief The vertices of the parallelotope. */
      std::vector< state_type > vertices() const;
      
      /*! \brief Convert to a zonotope. */
      operator Zonotope<R> () const;
      
      /*! \brief Convert to a polyhedron. */
      operator Polyhedron<R> () const;
       
      /*! \brief Tests if the parallelotope contains \a point. */
      bool contains(const state_type& point) const;
      
      /*! \brief Tests if the parallelotope contains a \a rectangle. */
      bool contains(const Rectangle<R>& rect) const; 
    
      /*! \brief Tests if the interior of the parallelotope contains \a point. */
      bool interior_contains(const state_type& point) const;

      /*! \brief Subdivide into two smaller pieces. */
      ListSet<R,Geometry::Parallelotope> divide() const;
      
      /*! \brief Subdivide into smaller pieces in each dimension. */
      ListSet<R,Geometry::Parallelotope> subdivide() const;
      
      /*! \brief Tests disjointness from a rectangle. */
      bool disjoint(const Rectangle<R>& r) const;
      
      friend std::istream& operator>> <> (std::istream& is, Parallelotope<R>& r);
     private:
      void compute_linear_inequalities(Matrix_type&, Vector_type&, Vector_type&) const;
      LinearAlgebra::Vector<F> coordinates(const state_type& s) const;
     private:
      /* Parallelotope's centre. */
      state_type _centre;
      /* Parallelotope's principal directions. */
      Matrix_type _generators;
      
    };
 
    /*! \brief Performs the Minkoswi sum of two parallelotopes */
    template<typename R> 
    Parallelotope<R> 
    minkowski_sum(const Parallelotope<R>& A, const Parallelotope<R>& B) {
        throw std::domain_error("minkowski_sum(const Parallelotope<R>& A, const Parallelotope<R>& B): not implemented.");

	return Parallelotope<R>(A.dimension());
    }

    /*! \brief Performs the Minkoswi sum of a zonotope and a basic set */
    template<typename R, template <typename> class BS> 
    inline Parallelotope<R> 
    minkowski_sum(const Parallelotope<R>& A, const BS<R>& B) {
      return minkowski_sum(A, Parallelotope<R>(B) );
    }

    /*! \brief Tests disjointness */
    template <typename R>
    inline
    bool 
    disjoint(const Parallelotope<R>& A, const Parallelotope<R>& B) 
    {
      return disjoint(Polyhedron<R>(A),Polyhedron<R>(B));
    }
    
    /*! \brief Tests disjointness */
    template <typename R>
    inline
    bool 
    disjoint(const Parallelotope<R>& A, const Rectangle<R>& B)
    {
      return A.disjoint(B);
    }
    
    /*! \brief Tests disjointness */
    template <typename R>
    inline
    bool 
    disjoint(const Rectangle<R>& A, const Parallelotope<R>& B) 
    {
      return disjoint(B,A);
    }
    
    /*! \brief Tests intersection of interiors */
    template <typename R>
    inline
    bool 
    interiors_intersect(const Parallelotope<R>& A,
                        const Parallelotope<R>& B) 
    {
      return interiors_intersect(Polyhedron<R>(A), Polyhedron<R>(B));
    }
    
    /*! \brief Tests intersection of interiors */
    template <typename R>
    inline
    bool 
    interiors_intersect(const Parallelotope<R>& A,
                        const Rectangle<R>& B) 
    {
      return interiors_intersect(Polyhedron<R>(A), Polyhedron<R>(B));
    }
    
    /*! \brief Tests intersection of interiors */
    template <typename R>
    inline
    bool
    interiors_intersect(const Rectangle<R>& A,
                        const Parallelotope<R>& B) 
    {
      return interiors_intersect(Polyhedron<R>(A), Polyhedron<R>(B));
    }
    
    
    /*! \brief Tests inclusion of \a A in the interior of \a B. */
    template <typename R>
    inline
    bool
    inner_subset(const Parallelotope<R>& A,
                 const Parallelotope<R>& B) 
    {
      return inner_subset(Polyhedron<R>(A), Polyhedron<R>(B));
    }

    /*! \brief Tests inclusion of \a A in the interior of \a B. */
    template <typename R>
    inline
    bool
    inner_subset(const Parallelotope<R>& A,
                 const Rectangle<R>& B) 
    {
      return inner_subset(Polyhedron<R>(A), Polyhedron<R>(B));
    }

    /*! \brief Tests inclusion of \a A in the interior of \a B. */
    template <typename R>
    inline
    bool
    inner_subset(const Rectangle<R>& A,
                 const Parallelotope<R>& B) 
    {
      return inner_subset(Polyhedron<R>(A), Polyhedron<R>(B));
    }

    /*! \brief Tests inclusion */
    template <typename R>
    inline
    bool 
    subset(const Parallelotope<R>& A, 
           const Parallelotope<R>& B) 
    {
      return subset(Polyhedron<R>(A), Polyhedron<R>(B));
    }
    
    /*! \brief Tests inclusion */
    template <typename R>
    inline
    bool 
    subset(const Parallelotope<R>& A, 
           const Rectangle<R>& B) 
    {
      return subset(Polyhedron<R>(A), Polyhedron<R>(B));
    }
    
    /*! \brief Tests inclusion */
    template <typename R>
    inline
    bool 
    subset(const Rectangle<R>& A, 
           const Parallelotope<R>& B) 
    {
      return subset(Polyhedron<R>(A), Polyhedron<R>(B));
    }
    

    /*! \brief Tests inclusion in an open cover.  */
    template <typename R>
    inline
    bool 
    subset_of_open_cover(const Parallelotope<R>& A,
                         const ListSet<R, Parallelotope >& cover) 
    {
      throw std::domain_error("subset_of_open_cover(Parallelotope, std::vector<Parallelotope>) not implemented");
    }

    /*! \brief Scale a parallelotope. */
    template<typename R>
    inline
    Geometry::Parallelotope<R> 
    scale(const Geometry::Parallelotope<R>& p, const R& scale_factor) {

      const Geometry::Point<R>& centre=p.centre();
      const LinearAlgebra::Matrix<R>& generators=p.generators();
      
      Geometry::Point<R> new_centre(p.dimension());

      for(size_type i=0; i!=p.dimension(); ++i) {
        new_centre[i]=scale_factor*centre[i];
      }

      return Geometry::Parallelotope<R>(new_centre, scale_factor*generators);
    }

    
    /*! \brief Tests inclusion of \a A om the interior of \a B. */
    template <typename R>
    inline
    bool 
    inner_subset(const Parallelotope<R>& A,
                 const ListSet<R,Parallelotope>& B) 
    {
      throw std::domain_error("subset_of_closed_cover(Parallelotope, std::vector<Parallelotope>) not implemented");
    }

  }
}

#endif /* _ARIADNE_PARALLELOTOPE_H */
