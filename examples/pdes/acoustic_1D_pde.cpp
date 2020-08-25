#include "numeric/numeric.hpp"
#include "numeric/real.hpp"
#include "algebra/algebra.hpp"
#include "algebra/vector.hpp"
#include "utility/array.hpp"

#include "output/gnuplot.hpp"
#include "dynamics/1D_pde.hpp"

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <unistd.h>

using namespace Ariadne;
using namespace std;

int main() {

    SizeType Nx = 10;                // # point in space

    //pde_solution_1D solution;
    Parameter1D stringModel;

    stringModel.length = 1.0;  // Length
    stringModel.tension = 100000;
    stringModel.mass = 4;
    FloatMP x0 = 0.85*stringModel.length;  // Point of max amplitube - Triangular IC
    stringModel.wavelength = 2*stringModel.length;
    //Real c = string.frequency*string.wavelength;
    FloatMP amp = 0.8;
    stringModel.damping = 100;

    stringModel.c = sqrt((stringModel.tension/(stringModel.mass/stringModel.length)));

    stringModel.frequency = stringModel.c/stringModel.wavelength; // Frequency

    //Real T = string.wavelength/c;
    stringModel.CourantNumber = 0.8;
    //Real C2 = pow(string.CourantNumber, 2);

    FloatMP k = 2*pi_opp()/stringModel.wavelength;
    FloatMP omega = 2*pi_opp()*stringModel.frequency;

    // Set function initial position condition
    //std::function<FloatMP(FloatMP)> PosInitSteady = [&](FloatMP x){return sin(k*x);};  //Steady-state wave
   
    std::function<FloatMP(FloatMP)> PosInitTri = [&](FloatMP x){ //Triangular
        if (x.get_d() < x0.get_d())
            return amp*(x/x0);
        else
            return amp*(stringModel.length - x)/(stringModel.length - x0);
    };

    std::function<FloatMP(FloatMP, FloatMP)> source = [&](FloatMP x, FloatMP t){return (sin(k*x)*cos(omega*t));};

    auto data = pde_1Dsolver(PosInitTri, source, stringModel, Nx+1);
  

    //auto VelInitCond = [&](Real x){return 0;};            // Set function initial velocity position
   
    Gnuplot gp;
    GnuplotCanvas canvas = GnuplotCanvas();

    Image2D image;

    canvas.setTerminal(gp, image, _png, "string-SpaceEvol");
    canvas.setMultiplot(gp, true);
    canvas.setTitle(gp, "String Evolution with triangular IC");
    canvas.setXLabel(gp, "x - Space");
    canvas.setYLabel(gp, "Amplitude");
    canvas.setRange2D(image, 0, Nx+1, -1, 1);

    canvas.plotTensor2D(gp, image, data);


    return 0;

}