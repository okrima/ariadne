/***************************************************************************
 *            polyhedron.cc
 *
 *  Copyright  2006-8  Alberto Casagrande, Pieter Collins
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
 *  along with this program; if not, write to bouthe Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#include <iostream>
#include <string>
#include <sstream>

#include <list>
#include <set>
#include <vector>
#include <valarray>

#include "polyhedron.h"

#include "macros.h"
#include "exceptions.h"

#include "numeric.h"

#include "vector.h"
#include "matrix.h"

#include "function.h"

#include "box.h"
#include "polytope.h"



namespace Ariadne {


extern int global_verbosity;
int verbosity=global_verbosity;


void ddconv(std::vector< Vector<Float> >&, const std::vector< Vector<Float> >&);

Polyhedron polyhedron(const Box& bx);
Polyhedron polyhedron(const Polytope& p);
Polytope polytope(const Polyhedron& p);



Polyhedron::Polyhedron()
    : ConstraintSet(Box(0),ConstantFunction(Vector<Float>(0,1.0),0))
{
}


Polyhedron::Polyhedron(uint d)
    : ConstraintSet( Box(0), ConstantFunction(Vector<Float>(0,1.0),d) )
{
}


Polyhedron::Polyhedron(const Matrix<Float>& A, const Vector<Float>& b)
    : ConstraintSet( Box(A.row_size(),Interval(0,inf())), AffineFunction(-A,b) )
{
}




Polyhedron::Polyhedron(const Box& bx)
    : ConstraintSet( Ariadne::polyhedron(bx) )
{
}

Polyhedron::Polyhedron(const Polytope& p)
    : ConstraintSet( Ariadne::polyhedron(p) )
{
}



Matrix<Float>
Polyhedron::A() const
{
    return -dynamic_cast<const AffineFunction&>(this->function()).A();
}

Vector<Float>
Polyhedron::b() const
{
    return dynamic_cast<const AffineFunction&>(this->function()).b();
}



tribool
Polyhedron::empty() const
{
    ARIADNE_NOT_IMPLEMENTED;
}


tribool
Polyhedron::bounded() const
{
    ARIADNE_NOT_IMPLEMENTED;
}


Polyhedron
Polyhedron::halfspace(size_t i) const
{
    ARIADNE_ASSERT(i<this->number_of_constraints());
    uint d=this->dimension();
    Matrix<Float> A=project(this->A(),range(i,i+1),range(0,d));
    Vector<Float> b=project(this->b(),range(i,i+1));
    return Polyhedron(A,b);
}





Polyhedron 
intersection(const Polyhedron& plhd1, const Polyhedron& plhd2)
{
    ARIADNE_ASSERT(plhd1.dimension()==plhd2.dimension());
    uint d=plhd1.dimension();
    size_t nc1=plhd1.number_of_constraints();
    size_t nc2=plhd2.number_of_constraints();
    Matrix<Float> A(nc1+nc2,d);
    project(A,range(0,nc1),range(0,d)) = plhd1.A();
    project(A,range(nc1,nc1+nc2),range(0,d)) = plhd1.A();
    Vector<Float> b=join(plhd1.b(),plhd2.b());
    return Polyhedron(A,b);
}


Polyhedron
polyhedron(const Box& bx)
{
    Vector<Float> b; Matrix<Float> A;
    uint d=bx.size();
    for(uint i=0; i!=d; ++i) {
        A[i][i]=-1.0;
        b[i]=-bx[i].lower();
        A[i+d][i]=1.0;
        b[i+d]=bx[i].upper();
    }
    return Polyhedron(A,b);
}



Polyhedron
polyhedron(const Polytope& pltp)
{
    ARIADNE_NOT_IMPLEMENTED;
}



Polytope
polytope(const Polyhedron& pltp)
{
    ARIADNE_NOT_IMPLEMENTED;
}







std::ostream& 
Polyhedron::write(std::ostream& os) const
{
    //return os << "Polyhedron( A=" << this->A() << ", b=" << this->b() << " )";
    const Matrix<Float> A=this->A();
    const Vector<Float> b=this->b();
    os << "Polyhedron( constraints=";
    uint d=this->dimension();
    size_t nc=this->number_of_constraints();
    for(size_t i=0; i!=nc; ++i) {
        os << ( i==0 ? "[" : "," );
        for(size_t j=0; j!=d; ++j) {
            os << ( j==0 ? "(" : ",");
            os << A[i][j]; 
        }
        os << ":" << b[i] << ")";
    }
    os << "] )";
    return os;
}

std::istream& 
operator>>(std::istream& is, Polyhedron& p) 
{
    std::vector< std::vector<Float> > Alst;
    std::vector< Float > Blst;
  
    std::vector<Float> a;
    Float b;
  
    char c;
    is >> c;
    assert(c=='[');
  
    c=is.peek();
    while(c=='[') {
        // Read constraint ax<=b in form [a_1,a_2,...,a_n;b];
        read_sequence(is,a,'[',';',',');
        is >> b;
        is >> c;
        assert(c==']');
        Alst.push_back(a);
        Blst.push_back(b);
    }
  
    size_t m=Alst.size();
    size_t n=Alst[0].size();
    Matrix<Float> A(m,n);
    Vector<Float> B(m);
    for(uint i=0; i!=m; ++i) {
        for(size_t j=0; j!=n; ++j) {
            A[i][j]=Alst[i][j];
        }
        B[i]=Blst[i];
    }
  
    p=Polyhedron(A,B);
  
    return is;
}




} // namespace Ariadne

