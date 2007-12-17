/***************************************************************************
 *            test_taylor_derivative.cc
 *
 *  Copyright  2007  Alberto Casagrande, Pieter Collins
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
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#include <cassert>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>

#include "test_float.h"
#include "numeric/rational.h"
#include "linear_algebra/vector.h"
#include "function/taylor_series.h"
#include "function/taylor_variable.h"
#include "function/taylor_derivative.h"

#include "test.h"

using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::Function;
using namespace std;

template<class X>
class TestTaylorDerivative {
 private:
  TaylorDerivative<X> x1,x2,x3;
 public:
  TestTaylorDerivative() {
    double a1[15]={ 2.0, 1.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    double a2[15]={ 3.0, 1.0, 0.0, 0.25, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    double a3[15]={ 2.0, 1.0, 0.0, 0.125, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    x1=TaylorDerivative<X>(1,2,4,a1);
    x2=TaylorDerivative<X>(1,2,4,a2);
    x3=TaylorDerivative<X>(1,1,4,a3);
    
    ARIADNE_TEST_CALL(test_degree());
    ARIADNE_TEST_CALL(test_add());
    ARIADNE_TEST_CALL(test_sub());
    ARIADNE_TEST_CALL(test_mul());
    ARIADNE_TEST_CALL(test_div());
    ARIADNE_TEST_CALL(test_compose());
  }

  void test_degree() {
    ARIADNE_TEST_ASSERT(x1.degree()==4);
  }

  void test_add() {
    cout << x1 << "+" << x2 << " = " << x1+x2 << std::endl;
    //assert((x1+x2)==TaylorDerivative<X>("[3,2,0,0]"));
  }

  void test_sub() {
    cout << x1 << "-" << x2 << " = " << x1-x2 << std::endl;
    //assert((x1-x2)==TaylorDerivative<X>("[-1,0,0,0]"));
  }

  void test_mul() {
    X c=2;
    cout << x1 << "*" << c << " = " << x1*c << std::endl;
    cout << c << "*" << x1 << " = " << c*x1 << std::endl;
    //assert((x1*x2)==TaylorDerivative<X>("[2,3,2,0]"));
  }

  void test_div() {
    X c=2;
    cout << x1 << "/" << c << " = " << x1/c << std::endl;
  }

  void test_compose() {
    //double ax[10] = { 3.0, 1.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    double ax[10] = { 3.0, 1.0, 0.0, 0.0, 0.125, 0.25, 0.0, 0.0, 0.0, 0.0 };
    double ay[4] = { 1.0, -1.0, 0.5, -0.25 };
    double aid[4] = { ax[0], 1.0, 0.0, 0.0 };
    TaylorDerivative<X> x(1,2,3,ax);
    TaylorVariable<X> y(1,3,ay);
    TaylorDerivative<X> id(1,1,3,aid);
    cout << "x=" << x << endl;
    cout << "y=" << y << endl;
    cout << "compose(y,x)=" << compose(y,x) << endl;
    cout << "compose(id,x)=" << compose(id,x) << endl;
    cout << "compose(id,x)-x=" << compose(id,x)-x << endl;
  }


};


int main() {
  TestTaylorDerivative<Rational> t1;
  cout << "INCOMPLETE " << flush;
  return 0;
}