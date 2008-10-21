/***************************************************************************
 *            test_graphics.cc
 *
 *  Copyright 2008  Pieter Collins
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
 
#include "function.h"
#include "graphics.h"
#include "box.h"
#include "zonotope.h"
#include "approximate_taylor_model.h"
#include "set.h"

using namespace Ariadne;


struct RadiusSquare : FunctionData<1,2,1> {
  template<class R, class A, class P>
  void compute(R& r, const A& x, const P& p) const {
    r[0]=sqr(x[0])+sqr(x[1])-sqr(p[0]);
  }
};
                   


int main(int argc, char **argv) 
{

    Box bx1(2); bx1[0]=Interval(-0.2,0.2); bx1[1]=Interval(-0.1,0.10);
    Box bx2(2); bx2[0]=Interval(0.1,0.3); bx2[1]=Interval(0.05,0.15);
    Box bx3(2); bx3[0]=Interval(0.2,0.4); bx3[1]=Interval(0.10,0.25);
    Box bx4(2); bx4[0]=Interval(0.25,0.5); bx4[1]=Interval(0.20,0.50);
    Box bx5(2); bx5[0]=Interval(0.4,0.8); bx5[1]=Interval(0.40,1.1);
    double z1cdata[]={0.15,0.6}; double z1gdata[]={0.05,0.0,0.05, 0.0,0.05,0.05};
    Vector<Float> z1c(2,z1cdata);
    Matrix<Float> z1g(2,3,z1gdata);
    Zonotope z1(z1c,z1g);
    Vector<Float> ts1c=z1c-Vector<Float>(2,Float(0.25));
    Matrix<Float> ts1g=z1g;
    AffineFunction afn1(ts1g,ts1c);
    ApproximateTaylorModel ts1(Box(3,Interval(-1,1)),afn1,1,0);

    Function<RadiusSquare> radius(Vector<Float>(1u,0.5));
    ConstraintSet cs1(Box(1u,Interval(-1,0)),radius);
    

    Graphic g;
    g.set_fill_colour(Colour(0.5,1.0,1.0));
    g.plot(bx1);
    g.plot(bx2);
    g.plot(bx3);
    g.plot(bx4);
    g.plot(bx5);
    g.set_fill_colour(Colour(0.0,0.5,0.5));
    g.plot(polytope(z1));
    g.plot(polytope(ts1));
    g.write("test_graphics-1");
    
    g.set_fill_colour(Colour(1.0,1.0,0.5));
    g.display();

    g.clear();

    g.set_fill_colour(Colour(1.0,0.5,1.0));
    g.plot(bx1);
    g.plot(bx2);
    g.plot(bx5);
    g.write("test_graphics-2");

}

